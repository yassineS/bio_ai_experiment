package samtools

// tview_interactive.go implements the `-d C` interactive alignment viewer as a
// pure-Go control loop on top of the existing tview frame renderer (tview.go).
// Upstream's `-d C` is an ncurses UI; this port avoids that dependency by
// driving a small raw-mode terminal directly (see tview_tty_linux.go) and
// reusing the SAME tvScreen grid the text/HTML backends already build.
//
// The design splits cleanly into two halves:
//
//   - This file: pure, TTY-free logic. Keystrokes are decoded into a tvAction
//     (decodeKey), and a tvAction applied to a tvViewState produces the next
//     state (applyAction). Both are ordinary functions with no I/O, so they are
//     fully unit-testable without a terminal. RunTviewInteractive ties the
//     state machine to a terminal abstraction (the tvTerminal interface) and
//     the frame renderer.
//
//   - tview_tty_linux.go (and the tview_tty_other.go stub): the tiny,
//     untestable OS layer — raw-mode termios set/restore and the TIOCGWINSZ
//     window-size ioctl. It is intentionally minimal.
//
// # Key bindings (ported from bam_tview_curses.c::curses_loop)
//
//   Left  / h           scroll one column left
//   Right / l           scroll one column right
//   Up    / k           scroll one read row up    (row_shift)
//   Down  / j           scroll one read row down  (row_shift)
//   H                   page 20 columns left
//   L                   page 20 columns right
//   space               page one screen-width right
//   backspace / DEL     page one screen-width left
//   Ctrl-H (^H)         jump 1000 columns left
//   Ctrl-L (^L)         jump 1000 columns right
//   0 / Home            jump to the start of the contig (position 0)
//   g / /               prompt for a chr:pos region and jump there
//   m / b / n / N       colour mode: mapping-q / base-q / nucleotide / none
//   . (dot)             toggle base-vs-dot rendering for matches
//   r                   toggle by-read-name colouring
//   i                   toggle insertion sub-columns
//   ?                   show the help screen (any key returns)
//   q / Esc / Ctrl-C    quit
//
// Note: upstream binds j/k (and Up/Down) to read-row scrolling and h/l (and
// Left/Right) to column scrolling; this port keeps that mapping. `0`/Home (jump
// to contig start) is an ergonomic addition. The colour-mode toggles (m/b/n/N)
// and `r` track state shown in the status line but do not change the rendered
// glyphs: the text grid this port draws is colourless, so the toggles are inert
// beyond the status report (the same way `-d T` never emits ANSI colour).

import (
	"fmt"
	"io"
	"strconv"
)

// tvColorMode selects which per-read attribute would drive colouring,
// mirroring upstream's TV_COLOR_* enum. The text grid is colourless, so the
// value is tracked and surfaced in the status line but otherwise inert.
type tvColorMode int

const (
	// tvColorMapQ colours by mapping quality (upstream default).
	tvColorMapQ tvColorMode = iota
	// tvColorBaseQ colours by base quality.
	tvColorBaseQ
	// tvColorNucl colours by nucleotide.
	tvColorNucl
	// tvColorNone disables colouring.
	tvColorNone
)

// String returns the short human-readable name of the colour mode for the
// status line.
func (m tvColorMode) String() string {
	switch m {
	case tvColorBaseQ:
		return "base-q"
	case tvColorNucl:
		return "nucleotide"
	case tvColorNone:
		return "none"
	default:
		return "mapping-q"
	}
}

// tvAction is the decoded result of a keystroke: the abstract command the
// interactive loop should apply, independent of the exact key pressed.
type tvAction int

// The interactive actions. Movement actions adjust the window; toggle actions
// flip a render flag; the control actions (quit/goto/help) are surfaced via the
// state's control flags for the loop to act on.
const (
	tvActNone tvAction = iota
	tvActQuit
	tvActLeft       // one column left
	tvActRight      // one column right
	tvActRowUp      // scroll reads up
	tvActRowDown    // scroll reads down
	tvActPageLeft   // page one screen-width left
	tvActPageRight  // page one screen-width right
	tvActSmallLeft  // 20 columns left
	tvActSmallRight // 20 columns right
	tvActBigLeft    // 1000 columns left
	tvActBigRight   // 1000 columns right
	tvActStart      // jump to contig start (pos 0)
	tvActGoto       // prompt for chr:pos and jump
	tvActHelp       // show help screen
	tvActToggleDot  // toggle base-vs-dot
	tvActToggleIns  // toggle insertion columns
	tvActToggleName // toggle by-read-name colouring
	tvActColorMapQ
	tvActColorBaseQ
	tvActColorNucl
	tvActColorNone
)

// ANSI / control byte constants used by the key decoder.
const (
	keyEsc       = 0x1b // ESC
	keyCtrlC     = 0x03
	keyCtrlH     = 0x08 // ^H
	keyCtrlL     = 0x0c
	keyBackspace = 0x7f // DEL
)

// decodeKey maps a key event to a tvAction. A key event is one byte for an
// ordinary key, or a multi-byte ANSI escape sequence (ESC '[' final, or ESC 'O'
// final) for arrow/Home keys. It is a pure function: the entire interactive key
// map lives here so it can be unit-tested without a terminal.
func decodeKey(seq []byte) tvAction {
	if len(seq) == 0 {
		return tvActNone
	}
	if seq[0] == keyEsc {
		if len(seq) == 1 {
			// A lone ESC quits, matching upstream's '\033' case.
			return tvActQuit
		}
		switch {
		case isArrow(seq, 'A'):
			return tvActRowUp
		case isArrow(seq, 'B'):
			return tvActRowDown
		case isArrow(seq, 'C'):
			return tvActRight
		case isArrow(seq, 'D'):
			return tvActLeft
		case isHome(seq):
			return tvActStart
		default:
			return tvActNone
		}
	}
	switch seq[0] {
	case 'q', keyCtrlC:
		return tvActQuit
	case '?':
		return tvActHelp
	case 'g', '/':
		return tvActGoto
	case 'h':
		return tvActLeft
	case 'l':
		return tvActRight
	case 'j':
		return tvActRowDown
	case 'k':
		return tvActRowUp
	case 'H':
		return tvActSmallLeft
	case 'L':
		return tvActSmallRight
	case ' ':
		return tvActPageRight
	case keyBackspace:
		return tvActPageLeft
	case keyCtrlH:
		return tvActBigLeft
	case keyCtrlL:
		return tvActBigRight
	case '0':
		return tvActStart
	case '.':
		return tvActToggleDot
	case 'i':
		return tvActToggleIns
	case 'r':
		return tvActToggleName
	case 'm':
		return tvActColorMapQ
	case 'b':
		return tvActColorBaseQ
	case 'n':
		return tvActColorNucl
	case 'N':
		return tvActColorNone
	default:
		return tvActNone
	}
}

// isArrow reports whether seq is an arrow escape ending in final ('A'..'D'),
// accepting both the CSI form (ESC '[' final) and the SS3 form (ESC 'O' final).
func isArrow(seq []byte, final byte) bool {
	return len(seq) == 3 && seq[0] == keyEsc &&
		(seq[1] == '[' || seq[1] == 'O') && seq[2] == final
}

// isHome reports whether seq is a Home-key escape (ESC [ H, ESC O H, or
// ESC [ 1 ~).
func isHome(seq []byte) bool {
	if len(seq) == 3 && seq[0] == keyEsc && (seq[1] == '[' || seq[1] == 'O') && seq[2] == 'H' {
		return true
	}
	if len(seq) == 4 && seq[0] == keyEsc && seq[1] == '[' && seq[2] == '1' && seq[3] == '~' {
		return true
	}
	return false
}

// tvViewState is the full mutable state of the interactive viewer: the window
// (contig, left position, width, height), the read-row scroll offset, the
// render toggles, and the control flags the loop reads back after applyAction.
type tvViewState struct {
	Chrom    string // displayed contig
	LeftPos0 int    // 0-based left edge of the window
	Width    int    // display width in columns (mcol)
	Height   int    // terminal height in rows (mrow)
	RowShift int    // first read row to show (scroll offset, >= 0)

	HideInserts bool        // -i toggle (true = insertion columns hidden)
	IsDot       bool        // '.' toggle (cosmetic; tracked in status line)
	ShowName    bool        // 'r' toggle (cosmetic; tracked in status line)
	ColorMode   tvColorMode // m/b/n/N colour mode (cosmetic; tracked in status)

	// Control flags consumed by the loop after applyAction. applyAction clears
	// them at the start of each dispatch so a returned state reflects only the
	// action just applied.
	Quit        bool // set by tvActQuit
	PendingGoto bool // set by tvActGoto: the loop should prompt for a region
	ShowHelp    bool // set by tvActHelp: the loop should show the help screen
}

// applyAction returns the next view state after applying act to st. It is a
// pure function (no I/O): scroll/page/jump actions adjust LeftPos0 (clamped at
// 0); row actions adjust RowShift (clamped at 0); toggles flip the matching
// flag; goto/help/quit set their control flag for the loop to act on.
func applyAction(st tvViewState, act tvAction) tvViewState {
	// Start each dispatch with the transient control flags cleared.
	st.Quit = false
	st.PendingGoto = false
	st.ShowHelp = false

	switch act {
	case tvActQuit:
		st.Quit = true
	case tvActLeft:
		st.LeftPos0--
	case tvActRight:
		st.LeftPos0++
	case tvActSmallLeft:
		st.LeftPos0 -= 20
	case tvActSmallRight:
		st.LeftPos0 += 20
	case tvActBigLeft:
		st.LeftPos0 -= 1000
	case tvActBigRight:
		st.LeftPos0 += 1000
	case tvActPageLeft:
		st.LeftPos0 -= st.Width
	case tvActPageRight:
		st.LeftPos0 += st.Width
	case tvActStart:
		st.LeftPos0 = 0
	case tvActRowUp:
		st.RowShift--
	case tvActRowDown:
		st.RowShift++
	case tvActGoto:
		st.PendingGoto = true
	case tvActHelp:
		st.ShowHelp = true
	case tvActToggleDot:
		st.IsDot = !st.IsDot
	case tvActToggleIns:
		st.HideInserts = !st.HideInserts
	case tvActToggleName:
		st.ShowName = !st.ShowName
	case tvActColorMapQ:
		st.ColorMode = tvColorMapQ
	case tvActColorBaseQ:
		st.ColorMode = tvColorBaseQ
	case tvActColorNucl:
		st.ColorMode = tvColorNucl
	case tvActColorNone:
		st.ColorMode = tvColorNone
	case tvActNone:
		// no-op
	}
	// Clamp upstream-style: pos and row_shift never go negative.
	if st.LeftPos0 < 0 {
		st.LeftPos0 = 0
	}
	if st.RowShift < 0 {
		st.RowShift = 0
	}
	return st
}

// tvStatusLine builds the one-line status/header shown above the alignment
// frame: contig, 1-based position, width, colour mode, and the active toggles,
// plus a hint for help. It is pure so it can be unit-tested.
func tvStatusLine(st tvViewState) string {
	dot := "off"
	if st.IsDot {
		dot = "on"
	}
	ins := "shown"
	if st.HideInserts {
		ins = "hidden"
	}
	name := "off"
	if st.ShowName {
		name = "on"
	}
	return fmt.Sprintf("%s:%d  w=%d  color=%s  dot=%s  ins=%s  byname=%s   [? help, q quit]",
		st.Chrom, st.LeftPos0+1, st.Width, st.ColorMode, dot, ins, name)
}

// tvHelpText is the help screen shown for '?'. It lists every key binding.
const tvHelpText = "samtools tview - interactive viewer key bindings\r\n\r\n" +
	"  Left  / h      scroll one column left\r\n" +
	"  Right / l      scroll one column right\r\n" +
	"  Up    / k      scroll read rows up\r\n" +
	"  Down  / j      scroll read rows down\r\n" +
	"  H              page 20 columns left\r\n" +
	"  L              page 20 columns right\r\n" +
	"  space          page one screen left-to-right\r\n" +
	"  backspace      page one screen right-to-left\r\n" +
	"  Ctrl-H         jump 1000 columns left\r\n" +
	"  Ctrl-L         jump 1000 columns right\r\n" +
	"  0 / Home       jump to start of contig\r\n" +
	"  g  /           go to a chr:pos region\r\n" +
	"  m b n N        colour by map-q / base-q / nucleotide / none\r\n" +
	"  .              toggle base-vs-dot\r\n" +
	"  i              toggle insertion columns\r\n" +
	"  r              toggle by-read-name colour\r\n" +
	"  ?              this help\r\n" +
	"  q / Esc        quit\r\n\r\n" +
	"Press any key to return.\r\n"

// tvTerminal abstracts the small set of terminal operations the interactive
// loop needs. The real implementation (linuxTerminal) wraps stdin/stdout in
// raw mode; tests supply a scripted fake so the loop logic is exercised without
// a TTY. Implementations set up raw mode in their constructor and restore the
// terminal in Close.
type tvTerminal interface {
	// ReadKey returns the next key event: a single byte for an ordinary key,
	// or a full ANSI escape sequence (starting with ESC) for arrow/Home keys.
	// It returns io.EOF when input ends.
	ReadKey() ([]byte, error)
	// ReadLine reads a line of input (for the goto prompt), returning the
	// entered text without the trailing newline. prompt is shown first.
	ReadLine(prompt string) (string, error)
	// Size returns the current terminal (width, height) in (columns, rows).
	Size() (int, int)
	// Write writes raw bytes (already containing any ANSI escapes) to the
	// terminal output.
	Write(p []byte) error
	// Close restores the terminal to its original (cooked) mode.
	Close() error
}

// ansiClearHome clears the screen (ESC [ 2 J) and homes the cursor (ESC [ H).
const ansiClearHome = "\x1b[2J\x1b[H"

// renderTviewFrame renders the current frame for st into a byte slice ready to
// write to the terminal: a leading clear-and-home, the status line, then the
// alignment grid (read rows scrolled by st.RowShift, limited to the terminal
// height), with CR/LF line endings (required in raw mode, where the terminal
// does not translate '\n'). It reuses the existing tview grid builder.
func renderTviewFrame(st tvViewState, opts TviewOptions) ([]byte, error) {
	screen, err := tviewRenderWindowAt(opts, st.Chrom, st.LeftPos0, st.Width)
	if err != nil {
		return nil, err
	}
	var buf []byte
	buf = append(buf, ansiClearHome...)
	buf = append(buf, tvStatusLine(st)...)
	buf = append(buf, '\r', '\n')

	rows := screen.rows
	// One line is reserved for the status header; show as many grid rows as
	// fit. Rows 0..tvMinAlnRow (ruler/reference/consensus) are always shown;
	// RowShift scrolls only the read rows (index > tvMinAlnRow).
	avail := st.Height - 1
	if avail < 3 {
		avail = 3
	}
	shown := 0
	for y := 0; y < len(rows) && shown < avail; y++ {
		if y > tvMinAlnRow && (y-(tvMinAlnRow+1)) < st.RowShift {
			continue
		}
		row := rows[y]
		for x := 0; x < screen.mcol; x++ {
			buf = append(buf, row[x].ch)
		}
		buf = append(buf, '\r', '\n')
		shown++
	}
	return buf, nil
}

// RunTviewInteractive runs the interactive `-d C` viewer loop against term,
// starting from opts. It seeds the view state from opts and the terminal size,
// draws the initial frame, then loops: read a key, decode it to an action,
// apply it, handle goto/help/quit, and redraw. It returns when the user quits
// or input ends. The caller (the CLI) verifies stdin/stdout are TTYs and builds
// term in raw mode before calling this; term.Close restores the terminal.
func RunTviewInteractive(term tvTerminal, opts TviewOptions) error {
	// Resolve the starting contig/position with the same logic the
	// non-interactive backend uses (so a bare `-d C` starts where `-d T`
	// would).
	chrom, left0, err := tviewResolveStartForOpts(opts)
	if err != nil {
		return err
	}

	w, h := term.Size()
	if opts.Width > 0 {
		w = opts.Width
	}
	st := tvViewState{
		Chrom:       chrom,
		LeftPos0:    left0,
		Width:       w,
		Height:      h,
		HideInserts: opts.HideInserts,
		ColorMode:   tvColorMapQ,
	}

	redraw := func() error {
		// Re-query size each frame so a resized terminal is honoured.
		cw, ch := term.Size()
		if opts.Width > 0 {
			cw = opts.Width
		}
		st.Width, st.Height = cw, ch
		frame, ferr := renderTviewFrame(st, opts)
		if ferr != nil {
			return ferr
		}
		return term.Write(frame)
	}

	if err := redraw(); err != nil {
		return err
	}

	for {
		seq, rerr := term.ReadKey()
		if rerr == io.EOF {
			return nil
		}
		if rerr != nil {
			return rerr
		}
		st = applyAction(st, decodeKey(seq))

		switch {
		case st.Quit:
			return nil
		case st.ShowHelp:
			if werr := term.Write([]byte(ansiClearHome + tvHelpText)); werr != nil {
				return werr
			}
			if _, herr := term.ReadKey(); herr != nil && herr != io.EOF {
				return herr
			}
		case st.PendingGoto:
			line, gerr := term.ReadLine("Go to (chr:pos): ")
			if gerr != nil && gerr != io.EOF {
				return gerr
			}
			if nc, np, ok := parseGotoRegion(line); ok {
				if nc != "" {
					st.Chrom = nc
				}
				st.LeftPos0 = np
				if st.LeftPos0 < 0 {
					st.LeftPos0 = 0
				}
			}
		}

		if err := redraw(); err != nil {
			return err
		}
	}
}

// parseGotoRegion parses the goto input ("chr:pos", "chr", or a bare position)
// into a (chrom, 0-based left position, ok) triple. A bare number is a position
// on the current contig (chrom == ""); a bare name is a contig jumped to at its
// start. Returns ok=false for empty input. It is pure so it can be unit-tested.
func parseGotoRegion(s string) (string, int, bool) {
	s = trimSpace(s)
	if s == "" {
		return "", 0, false
	}
	// chr:pos — split on the last ':' so contig names containing ':' parse.
	for i := len(s) - 1; i >= 0; i-- {
		if s[i] == ':' {
			chrom := s[:i]
			pos1, err := strconv.Atoi(removeCommas(s[i+1:]))
			if err != nil || pos1 < 1 {
				return "", 0, false
			}
			return chrom, pos1 - 1, true
		}
	}
	// Bare number => position on the current contig.
	if pos1, err := strconv.Atoi(removeCommas(s)); err == nil && pos1 >= 1 {
		return "", pos1 - 1, true
	}
	// Otherwise a bare contig name, jumped to at its start.
	return s, 0, true
}

// removeCommas strips ',' grouping separators from a numeric string (samtools
// accepts "1,000,000"-style positions).
func removeCommas(s string) string {
	out := make([]byte, 0, len(s))
	for i := 0; i < len(s); i++ {
		if s[i] != ',' {
			out = append(out, s[i])
		}
	}
	return string(out)
}

// trimSpace trims leading/trailing ASCII spaces, tabs, CR and LF.
func trimSpace(s string) string {
	i := 0
	for i < len(s) && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n') {
		i++
	}
	j := len(s)
	for j > i && (s[j-1] == ' ' || s[j-1] == '\t' || s[j-1] == '\r' || s[j-1] == '\n') {
		j--
	}
	return s[i:j]
}
