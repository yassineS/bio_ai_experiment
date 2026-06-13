package samtools

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// TestDecodeKey is a table-driven check of the key->action map, covering
// ordinary keys, control bytes, and the ANSI arrow/Home escape sequences. It
// needs no TTY: decodeKey is a pure function.
func TestDecodeKey(t *testing.T) {
	cases := []struct {
		name string
		seq  []byte
		want tvAction
	}{
		{"q quits", []byte{'q'}, tvActQuit},
		{"ctrl-c quits", []byte{keyCtrlC}, tvActQuit},
		{"lone esc quits", []byte{keyEsc}, tvActQuit},
		{"help", []byte{'?'}, tvActHelp},
		{"goto g", []byte{'g'}, tvActGoto},
		{"goto slash", []byte{'/'}, tvActGoto},
		{"h left", []byte{'h'}, tvActLeft},
		{"l right", []byte{'l'}, tvActRight},
		{"j row down", []byte{'j'}, tvActRowDown},
		{"k row up", []byte{'k'}, tvActRowUp},
		{"H small left", []byte{'H'}, tvActSmallLeft},
		{"L small right", []byte{'L'}, tvActSmallRight},
		{"space page right", []byte{' '}, tvActPageRight},
		{"backspace page left", []byte{keyBackspace}, tvActPageLeft},
		{"ctrl-h big left", []byte{keyCtrlH}, tvActBigLeft},
		{"ctrl-l big right", []byte{keyCtrlL}, tvActBigRight},
		{"zero start", []byte{'0'}, tvActStart},
		{"dot toggle", []byte{'.'}, tvActToggleDot},
		{"i toggle ins", []byte{'i'}, tvActToggleIns},
		{"r toggle name", []byte{'r'}, tvActToggleName},
		{"m color mapq", []byte{'m'}, tvActColorMapQ},
		{"b color baseq", []byte{'b'}, tvActColorBaseQ},
		{"n color nucl", []byte{'n'}, tvActColorNucl},
		{"N color none", []byte{'N'}, tvActColorNone},
		{"unknown key", []byte{'Z'}, tvActNone},
		{"empty", nil, tvActNone},
		// Escape sequences (CSI form).
		{"up arrow", []byte{keyEsc, '[', 'A'}, tvActRowUp},
		{"down arrow", []byte{keyEsc, '[', 'B'}, tvActRowDown},
		{"right arrow", []byte{keyEsc, '[', 'C'}, tvActRight},
		{"left arrow", []byte{keyEsc, '[', 'D'}, tvActLeft},
		{"home CSI", []byte{keyEsc, '[', 'H'}, tvActStart},
		{"home tilde", []byte{keyEsc, '[', '1', '~'}, tvActStart},
		// SS3 form (ESC O ...).
		{"up arrow ss3", []byte{keyEsc, 'O', 'A'}, tvActRowUp},
		{"home ss3", []byte{keyEsc, 'O', 'H'}, tvActStart},
		{"unknown escape", []byte{keyEsc, '[', 'Z'}, tvActNone},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := decodeKey(tc.seq); got != tc.want {
				t.Errorf("decodeKey(%v) = %d, want %d", tc.seq, got, tc.want)
			}
		})
	}
}

// TestApplyAction checks the action->state transitions: movement, paging,
// jumps, row scroll, toggles, and the control flags. applyAction is pure.
func TestApplyAction(t *testing.T) {
	base := tvViewState{Chrom: "chr1", LeftPos0: 100, Width: 40, Height: 24, RowShift: 5}

	t.Run("left/right move one column", func(t *testing.T) {
		if got := applyAction(base, tvActLeft); got.LeftPos0 != 99 {
			t.Errorf("left: LeftPos0 = %d, want 99", got.LeftPos0)
		}
		if got := applyAction(base, tvActRight); got.LeftPos0 != 101 {
			t.Errorf("right: LeftPos0 = %d, want 101", got.LeftPos0)
		}
	})

	t.Run("small/big/page moves", func(t *testing.T) {
		if got := applyAction(base, tvActSmallLeft); got.LeftPos0 != 80 {
			t.Errorf("smallLeft = %d, want 80", got.LeftPos0)
		}
		if got := applyAction(base, tvActSmallRight); got.LeftPos0 != 120 {
			t.Errorf("smallRight = %d, want 120", got.LeftPos0)
		}
		if got := applyAction(base, tvActBigRight); got.LeftPos0 != 1100 {
			t.Errorf("bigRight = %d, want 1100", got.LeftPos0)
		}
		if got := applyAction(base, tvActPageRight); got.LeftPos0 != 140 {
			t.Errorf("pageRight = %d, want 140 (width 40)", got.LeftPos0)
		}
		if got := applyAction(base, tvActPageLeft); got.LeftPos0 != 60 {
			t.Errorf("pageLeft = %d, want 60", got.LeftPos0)
		}
	})

	t.Run("clamp at zero", func(t *testing.T) {
		st := tvViewState{LeftPos0: 5, Width: 40}
		if got := applyAction(st, tvActBigLeft); got.LeftPos0 != 0 {
			t.Errorf("bigLeft from 5 = %d, want 0 (clamped)", got.LeftPos0)
		}
		st.RowShift = 0
		if got := applyAction(st, tvActRowUp); got.RowShift != 0 {
			t.Errorf("rowUp from 0 = %d, want 0 (clamped)", got.RowShift)
		}
	})

	t.Run("start jumps to zero", func(t *testing.T) {
		if got := applyAction(base, tvActStart); got.LeftPos0 != 0 {
			t.Errorf("start: LeftPos0 = %d, want 0", got.LeftPos0)
		}
	})

	t.Run("row scroll", func(t *testing.T) {
		if got := applyAction(base, tvActRowDown); got.RowShift != 6 {
			t.Errorf("rowDown: RowShift = %d, want 6", got.RowShift)
		}
		if got := applyAction(base, tvActRowUp); got.RowShift != 4 {
			t.Errorf("rowUp: RowShift = %d, want 4", got.RowShift)
		}
	})

	t.Run("quit sets flag", func(t *testing.T) {
		if got := applyAction(base, tvActQuit); !got.Quit {
			t.Error("quit: Quit flag not set")
		}
	})

	t.Run("goto and help set control flags", func(t *testing.T) {
		if got := applyAction(base, tvActGoto); !got.PendingGoto || got.ShowHelp || got.Quit {
			t.Errorf("goto: flags = %+v", got)
		}
		if got := applyAction(base, tvActHelp); !got.ShowHelp || got.PendingGoto {
			t.Errorf("help: flags = %+v", got)
		}
	})

	t.Run("toggles flip", func(t *testing.T) {
		st := applyAction(base, tvActToggleDot)
		if !st.IsDot {
			t.Error("dot toggle did not set IsDot")
		}
		if st2 := applyAction(st, tvActToggleDot); st2.IsDot {
			t.Error("second dot toggle did not clear IsDot")
		}
		if st := applyAction(base, tvActToggleIns); !st.HideInserts {
			t.Error("ins toggle did not set HideInserts")
		}
		if st := applyAction(base, tvActToggleName); !st.ShowName {
			t.Error("name toggle did not set ShowName")
		}
	})

	t.Run("color modes", func(t *testing.T) {
		for act, want := range map[tvAction]tvColorMode{
			tvActColorMapQ:  tvColorMapQ,
			tvActColorBaseQ: tvColorBaseQ,
			tvActColorNucl:  tvColorNucl,
			tvActColorNone:  tvColorNone,
		} {
			if got := applyAction(base, act); got.ColorMode != want {
				t.Errorf("color action %d: ColorMode = %v, want %v", act, got.ColorMode, want)
			}
		}
	})

	t.Run("control flags cleared on next non-control action", func(t *testing.T) {
		st := applyAction(base, tvActGoto) // sets PendingGoto
		st = applyAction(st, tvActRight)   // should clear it
		if st.PendingGoto {
			t.Error("PendingGoto not cleared by a subsequent movement action")
		}
	})
}

// TestParseGotoRegion checks the goto-prompt region parser.
func TestParseGotoRegion(t *testing.T) {
	cases := []struct {
		in        string
		wantChrom string
		wantPos0  int
		wantOK    bool
	}{
		{"chr1:1000", "chr1", 999, true},
		{"chr2:1", "chr2", 0, true},
		{"chrX:1,000,000", "chrX", 999999, true},
		{"500", "", 499, true},    // bare position on current contig
		{"chr3", "chr3", 0, true}, // bare contig -> start
		{"  chr1:50  ", "chr1", 49, true},
		{"", "", 0, false},                           // empty -> not ok
		{"chr1:0", "", 0, false},                     // pos < 1 invalid
		{"chr1:abc", "", 0, false},                   // malformed
		{"HLA-A*01:01:100", "HLA-A*01:01", 99, true}, // last colon splits
	}
	for _, tc := range cases {
		t.Run(tc.in, func(t *testing.T) {
			gotChrom, gotPos, gotOK := parseGotoRegion(tc.in)
			if gotOK != tc.wantOK || gotChrom != tc.wantChrom || gotPos != tc.wantPos0 {
				t.Errorf("parseGotoRegion(%q) = (%q,%d,%v), want (%q,%d,%v)",
					tc.in, gotChrom, gotPos, gotOK, tc.wantChrom, tc.wantPos0, tc.wantOK)
			}
		})
	}
}

// TestTvStatusLine spot-checks that the status line includes the position
// (1-based) and the active mode/toggles.
func TestTvStatusLine(t *testing.T) {
	st := tvViewState{Chrom: "chr1", LeftPos0: 99, Width: 60, ColorMode: tvColorNucl, IsDot: true, HideInserts: true}
	line := tvStatusLine(st)
	for _, want := range []string{"chr1:100", "w=60", "nucleotide", "dot=on", "ins=hidden"} {
		if !strings.Contains(line, want) {
			t.Errorf("status line %q missing %q", line, want)
		}
	}
}

// fakeTerminal is a scripted tvTerminal for exercising the interactive loop
// without a TTY: it replays a queue of key events, captures all writes, and
// returns a fixed size. ReadLine pops the next queued "line".
type fakeTerminal struct {
	keys   [][]byte
	lines  []string
	cols   int
	rows   int
	out    bytes.Buffer
	closed bool
}

func (f *fakeTerminal) ReadKey() ([]byte, error) {
	if len(f.keys) == 0 {
		return nil, io.EOF
	}
	k := f.keys[0]
	f.keys = f.keys[1:]
	return k, nil
}

func (f *fakeTerminal) ReadLine(prompt string) (string, error) {
	if len(f.lines) == 0 {
		return "", io.EOF
	}
	l := f.lines[0]
	f.lines = f.lines[1:]
	return l, nil
}

func (f *fakeTerminal) Size() (int, int) { return f.cols, f.rows }

func (f *fakeTerminal) Write(p []byte) error {
	f.out.Write(p)
	return nil
}

func (f *fakeTerminal) Close() error { f.closed = true; return nil }

// TestRunTviewInteractiveLoop drives the full control loop against a scripted
// fake terminal (no TTY): it sends a few movement keys, opens the help screen,
// performs a goto, and quits, asserting the loop renders frames and exits
// cleanly. It reuses the live-fixture builder so the renderer runs for real.
func TestRunTviewInteractiveLoop(t *testing.T) {
	bin := upstreamSamtoolsBinary(t)
	bamPath, refPath := buildTviewFixture(t, bin, tviewParitySAM, tviewParityRef)

	fake := &fakeTerminal{
		cols: 40,
		rows: 24,
		keys: [][]byte{
			{'l'},              // right one column
			{keyEsc, '[', 'C'}, // right arrow
			{'?'},              // help screen...
			{' '},              // ...any key returns
			{'g'},              // goto prompt
			{'q'},              // quit
		},
		lines: []string{"chr1:30"},
	}
	opts := TviewOptions{Input: bamPath, Reference: refPath, Position: "chr1:1"}
	if err := RunTviewInteractive(fake, opts); err != nil {
		t.Fatalf("RunTviewInteractive: %v", err)
	}
	out := fake.out.String()
	// The help screen must have been shown.
	if !strings.Contains(out, "key bindings") {
		t.Error("help screen text not found in output")
	}
	// After the goto to chr1:30 a status line for that position must appear.
	if !strings.Contains(out, "chr1:30") {
		t.Errorf("expected a frame at chr1:30 after goto; output:\n%s", out)
	}
	// The clear-and-home escape must have been emitted for each redraw.
	if !strings.Contains(out, ansiClearHome) {
		t.Error("expected ANSI clear/home escape in output")
	}
}

// TestRunTviewInteractiveEOF checks that the loop exits cleanly when input ends
// (no keys queued => immediate EOF after the initial draw).
func TestRunTviewInteractiveEOF(t *testing.T) {
	bin := upstreamSamtoolsBinary(t)
	bamPath, refPath := buildTviewFixture(t, bin, tviewParitySAM, tviewParityRef)

	fake := &fakeTerminal{cols: 40, rows: 24}
	opts := TviewOptions{Input: bamPath, Reference: refPath, Position: "chr1:1"}
	if err := RunTviewInteractive(fake, opts); err != nil {
		t.Fatalf("RunTviewInteractive on empty input: %v", err)
	}
	if fake.out.Len() == 0 {
		t.Error("expected at least the initial frame to be drawn")
	}
}

// TestRenderTviewFrameParity confirms the interactive frame's alignment grid
// matches the existing Tview text backend for the same window: the read rows
// (after the status header and CR/LF normalisation) are identical to `-d T`
// output. This reuses the existing tview fixtures.
func TestRenderTviewFrameParity(t *testing.T) {
	bin := upstreamSamtoolsBinary(t)
	bamPath, refPath := buildTviewFixture(t, bin, tviewParitySAM, tviewParityRef)

	st := tvViewState{Chrom: "chr1", LeftPos0: 0, Width: 40, Height: 100, ColorMode: tvColorMapQ}
	opts := TviewOptions{Input: bamPath, Reference: refPath}
	frame, err := renderTviewFrame(st, opts)
	if err != nil {
		t.Fatalf("renderTviewFrame: %v", err)
	}

	// Strip the leading clear/home and the status header line; convert CRLF to
	// LF; the remainder must equal the text backend's grid.
	body := strings.TrimPrefix(string(frame), ansiClearHome)
	nl := strings.IndexByte(body, '\n')
	if nl < 0 {
		t.Fatalf("frame had no status line")
	}
	grid := strings.ReplaceAll(body[nl+1:], "\r\n", "\n")

	var want bytes.Buffer
	if err := Tview(TviewOptions{Input: bamPath, Reference: refPath, Position: "chr1:1", Width: 40, Mode: TviewText}, &want); err != nil {
		t.Fatalf("Tview text: %v", err)
	}
	if grid != want.String() {
		t.Errorf("interactive frame grid != -d T output\n--- frame ---\n%s\n--- want ---\n%s", grid, want.String())
	}
}
