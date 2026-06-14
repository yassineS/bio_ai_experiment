package mosdepth

// Parity tests for mosdepth against the upstream Nim test suite at
// reference_code/mosdepth/tests/ and reference_code/mosdepth/functional-tests.sh.
//
// Methodology: each upstream `run <name> $exe ...` invocation in
// functional-tests.sh is mirrored as a `TestParity_<Name>` here. The Go
// port is driven through Run (or OpenAndRun) on the same in-tree BAM
// fixture, and the resulting per-base / regions / summary / dist files
// are diffed against the inline `assert_equal` heredoc string from the
// upstream script.
//
// One intentional difference from upstream:
//
//  1. Indexes. Like upstream, we now emit a `.csi` alongside each
//     bgzipped BED output (built with min_shift=14, depth=5 to match
//     htslib's tbx_index_build). Parity tests assert on the bgzipped
//     data files and the plain-text summary / distribution files; the
//     `.csi` is validated structurally and via round-trip query rather
//     than byte-diffed against upstream's (see TestParity_IndexFiles_Csi
//     and TestRunCsiReadable).
//
// Default-mode overlap-pair detection is now implemented: where the two
// mates of a read pair overlap the same reference bases, the shared bases
// are counted once, exactly as upstream mosdepth does (see addRecords /
// addOverlapCorrection in coverage.go). Every default-mode case below
// therefore asserts byte-for-byte against the values upstream's
// functional-tests.sh expects, with no t.Skip.

import (
	"compress/gzip"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// fixtureDir returns the absolute path to tools/mosdepth/testdata/parity.
func fixtureDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity"))
	if err != nil {
		t.Fatalf("Abs: %v", err)
	}
	return abs
}

// runParity invokes OpenAndRun against bamFile in fixtureDir with opts
// and returns the prefix used so the caller can read individual outputs
// back.
func runParity(t *testing.T, bamFile string, opts Options) string {
	t.Helper()
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "t")
	opts.Prefix = prefix
	bamPath := filepath.Join(fixtureDir(t), bamFile)
	if err := OpenAndRun(bamPath, opts); err != nil {
		t.Fatalf("OpenAndRun(%s): %v", bamPath, err)
	}
	return prefix
}

// readGzLines decompresses a gzip-wrapped (BGZF-compatible) file and
// returns its lines, with the trailing blank line trimmed.
func readGzLines(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip.NewReader: %v", err)
	}
	gr.Multistream(true)
	data, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("read gz: %v", err)
	}
	return strings.Split(strings.TrimRight(string(data), "\n"), "\n")
}

// readPlainLines reads a plain-text file and returns its lines with the
// trailing blank line trimmed.
func readPlainLines(t *testing.T, path string) []string {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return strings.Split(strings.TrimRight(string(data), "\n"), "\n")
}

// linesWithPrefix returns the subset of lines that start with one of the
// provided prefixes.
func linesWithPrefix(lines []string, prefixes ...string) []string {
	out := []string{}
	for _, ln := range lines {
		for _, p := range prefixes {
			if strings.HasPrefix(ln, p) {
				out = append(out, ln)
				break
			}
		}
	}
	return out
}

// equalLines reports whether two string slices are element-wise equal.
func equalLines(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// TestParity_OverlapM_DefaultPerBase mirrors `run overlapM $exe t tests/ovl.bam`.
// In default mode mosdepth de-duplicates the bases covered by both mates of a
// read pair, so the MT per-base output collapses to a single depth-1 run over
// [0,80). Upstream asserts exactly `MT 0 80 1 / MT 80 16569 0`
// (functional-tests.sh lines 28-29); we now match it byte-for-byte.
func TestParity_OverlapM_DefaultPerBase(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	got := linesWithPrefix(readGzLines(t, prefix+".per-base.bed.gz"), "MT\t")
	want := []string{
		"MT\t0\t80\t1",
		"MT\t80\t16569\t0",
	}
	if !equalLines(got, want) {
		t.Fatalf("MT default-mode per-base mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(got, "\n"))
	}
}

// TestParity_OverlapFastMode_MT mirrors
// `run overlapFastMode $exe t --fast-mode tests/ovl.bam`. The MT
// per-base output matches upstream byte-for-byte because --fast-mode
// skips the CIGAR walk; our default-mode and upstream's --fast-mode line
// up exactly. We restrict to chrom MT via --chrom so the test only has
// to sweep one reference (ovl.bam has ~80 references).
func TestParity_OverlapFastMode_MT(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		FastMode:    true,
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	got := linesWithPrefix(readGzLines(t, prefix+".per-base.bed.gz"), "MT\t")
	want := []string{
		"MT\t0\t6\t1",
		"MT\t6\t42\t2",
		"MT\t42\t80\t1",
		"MT\t80\t16569\t0",
	}
	if !equalLines(got, want) {
		t.Fatalf("MT per-base mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(got, "\n"))
	}
}

// TestParity_OverlapFastMode_Chr1Zero mirrors
// `assert_equal "$(zgrep -w ^1 t.per-base.bed.gz)" "1 0 249250621 0"` —
// chromosome 1 has no reads, so the per-base output is a single
// whole-chromosome zero-depth run.
func TestParity_OverlapFastMode_Chr1Zero(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{FastMode: true, ExcludeFlag: DefaultExcludeFlag})
	lines := readGzLines(t, prefix+".per-base.bed.gz")
	want := "1\t0\t249250621\t0"
	for _, ln := range lines {
		if ln == want {
			return
		}
	}
	t.Fatalf("missing chr1 zero-depth row.\nlines containing 1\\t:\n%s", strings.Join(linesWithPrefix(lines, "1\t"), "\n"))
}

// TestParity_OverlapM_SummaryMT mirrors
// `cat t.mosdepth.summary.txt | grep 'MT'` and `... grep 'total'` — the
// per-chromosome MT row and the total row. With overlap-pair de-duplication
// the 80 covered MT bases sit at depth 1, so the mean prints as 0.00
// (80/16569). Upstream asserts `MT 16569 80 0.00 0 1` and
// `total 16569 80 0.00 0 1` (functional-tests.sh lines 31-32).
func TestParity_OverlapM_SummaryMT(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	lines := readPlainLines(t, prefix+".mosdepth.summary.txt")
	var mt, total string
	for _, ln := range lines {
		if strings.HasPrefix(ln, "MT\t") {
			mt = ln
		}
		if strings.HasPrefix(ln, "total\t") {
			total = ln
		}
	}
	if mt != "MT\t16569\t80\t0.00\t0\t1" {
		t.Fatalf("summary MT row mismatch.\nwant: %q\ngot:  %q", "MT\t16569\t80\t0.00\t0\t1", mt)
	}
	if total != "total\t16569\t80\t0.00\t0\t1" {
		t.Fatalf("summary total row mismatch.\nwant: %q\ngot:  %q", "total\t16569\t80\t0.00\t0\t1", total)
	}
}

// TestParity_BigWindow mirrors
// `run big_window $exe t tests/ovl.bam --by 100000000`. The single 100Mb
// window spans all of MT; with overlap-pair de-duplication the 80 covered
// bases sit at depth 1, so the window mean is 80/16569 = 0.0048 → 0.00.
// Upstream asserts exactly `MT 0 16569 0.00` (functional-tests.sh line 62);
// we now match it byte-for-byte.
func TestParity_BigWindow(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		ByWindow:    100000000,
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	mt := linesWithPrefix(readGzLines(t, prefix+".regions.bed.gz"), "MT\t")
	want := []string{"MT\t0\t16569\t0.00"}
	if !equalLines(mt, want) {
		t.Fatalf("MT regions mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(mt, "\n"))
	}
}

// TestParity_LengthFilter_Excluded mirrors `--min-frag-len 81` and
// `--max-frag-len 79`. Both exclude all MT reads (TLEN==80), leaving the
// whole chromosome at depth 0.
func TestParity_LengthFilter_Excluded(t *testing.T) {
	cases := []struct {
		name string
		opts Options
	}{
		{"min-frag-len 81", Options{MinFragLen: 81, Chrom: "MT", ExcludeFlag: DefaultExcludeFlag}},
		{"max-frag-len 79", Options{MaxFragLen: 79, Chrom: "MT", ExcludeFlag: DefaultExcludeFlag}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			prefix := runParity(t, "ovl.bam", tc.opts)
			mt := linesWithPrefix(readGzLines(t, prefix+".per-base.bed.gz"), "MT\t")
			want := []string{"MT\t0\t16569\t0"}
			if !equalLines(mt, want) {
				t.Fatalf("want %v got %v", want, mt)
			}
		})
	}
}

// TestParity_LengthFilter_Included_FastMode mirrors
// `--min-frag-len 80 --max-frag-len 80` under --fast-mode. The TLEN==80
// reads are kept; output identical to the unfiltered fast-mode run.
func TestParity_LengthFilter_Included_FastMode(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		MinFragLen:  80,
		MaxFragLen:  80,
		FastMode:    true,
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	mt := linesWithPrefix(readGzLines(t, prefix+".per-base.bed.gz"), "MT\t")
	want := []string{
		"MT\t0\t6\t1",
		"MT\t6\t42\t2",
		"MT\t42\t80\t1",
		"MT\t80\t16569\t0",
	}
	if !equalLines(mt, want) {
		t.Fatalf("MT per-base mismatch under fast-mode + frag-len 80.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(mt, "\n"))
	}
}

// TestParity_ThresholdByBED mirrors
// `run threshold_test_by $exe --by tests/track.bed -T 0,1,2 -c MT t tests/ovl.bam`.
// With overlap-pair de-duplication the region MT:2..80 (width 78) is depth 1
// everywhere, so >=0 and >=1 each cover all 78 bases while >=2 covers none.
// Upstream asserts exactly `MT 2 80 aregion 78 78 0` (functional-tests.sh
// line 120); we now match it byte-for-byte.
func TestParity_ThresholdByBED(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		ByBED:       filepath.Join(fixtureDir(t), "track.bed"),
		Thresholds:  []int{0, 1, 2},
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	lines := readGzLines(t, prefix+".thresholds.bed.gz")
	if len(lines) != 2 {
		t.Fatalf("expected header + 1 row, got %d:\n%s", len(lines), strings.Join(lines, "\n"))
	}
	wantHdr := "#chrom\tstart\tend\tregion\t0X\t1X\t2X"
	if lines[0] != wantHdr {
		t.Errorf("header mismatch.\nwant: %q\ngot:  %q", wantHdr, lines[0])
	}
	wantRow := "MT\t2\t80\taregion\t78\t78\t0"
	if lines[1] != wantRow {
		t.Fatalf("row mismatch.\nwant: %q\ngot:  %q", wantRow, lines[1])
	}
}

// TestParity_TrackHeader mirrors
// `run track_header $exe --by tests/track.bed t tests/ovl.bam`. With
// overlap-pair de-duplication the region MT:2..80 is depth 1 over all 78
// bases, so its mean is exactly 1.00. Upstream asserts
// `MT 2 80 aregion 1.00` (functional-tests.sh line 137); we match it.
func TestParity_TrackHeader(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		ByBED:       filepath.Join(fixtureDir(t), "track.bed"),
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	got := linesWithPrefix(readGzLines(t, prefix+".regions.bed.gz"), "MT\t")
	want := []string{"MT\t2\t80\taregion\t1.00"}
	if !equalLines(got, want) {
		t.Fatalf("regions mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(got, "\n"))
	}
}

// TestParity_UnorderedBED mirrors
// `run unordered_bed $exe --by tests/unordered.bed t tests/ovl.bam` —
// regions.bed.gz must have exactly 2 lines.
func TestParity_UnorderedBED(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		ByBED:       filepath.Join(fixtureDir(t), "unordered.bed"),
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	lines := readGzLines(t, prefix+".regions.bed.gz")
	if len(lines) != 2 {
		t.Fatalf("expected 2 rows in regions.bed.gz, got %d:\n%s",
			len(lines), strings.Join(lines, "\n"))
	}
}

// TestParity_ReadGroupFilter_Match mirrors
// `run test_read_group $exe -n tt tests/ovl.bam -R GT04008021_119`.
// All ovl.bam reads have that RG; the dist file should contain MT and a
// total row.
func TestParity_ReadGroupFilter_Match(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		ReadGroups:  []string{"GT04008021_119"},
		NoPerBase:   true,
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	lines := readPlainLines(t, prefix+".mosdepth.global.dist.txt")
	mt := linesWithPrefix(lines, "MT\t")
	if len(mt) == 0 {
		t.Errorf("no MT rows in dist file")
	}
	hasTotal := false
	for _, ln := range lines {
		if strings.HasPrefix(ln, "total\t") {
			hasTotal = true
			break
		}
	}
	if !hasTotal {
		t.Errorf("missing total row in dist file:\n%s", strings.Join(lines, "\n"))
	}
}

// TestParity_ReadGroupFilter_Missing mirrors
// `run test_missing_read_group $exe -n tt tests/ovl.bam -R MISSING` —
// no read matches "MISSING".
func TestParity_ReadGroupFilter_Missing(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		ReadGroups:  []string{"MISSING"},
		NoPerBase:   true,
		ExcludeFlag: DefaultExcludeFlag,
		Chrom:       "MT",
	})
	lines := readPlainLines(t, prefix+".mosdepth.global.dist.txt")
	want := []string{"MT\t0\t1.00", "total\t0\t1.00"}
	if !equalLines(lines, want) {
		t.Fatalf("dist mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(lines, "\n"))
	}
}

// TestParity_MissingChrom_StrictFail mirrors
// `run missing_chrom $exe -c nonexistent --by 20000 t tests/ovl.bam`. Upstream
// writes `[mosdepth] chromosome nonexistent not found` to stderr and exits 1.
// We now return a ChromNotFoundError carrying the same message instead of
// silently emitting empty outputs.
func TestParity_MissingChrom_StrictFail(t *testing.T) {
	tmp := t.TempDir()
	bamPath := filepath.Join(fixtureDir(t), "ovl.bam")
	err := OpenAndRun(bamPath, Options{
		Prefix:      filepath.Join(tmp, "t"),
		Chrom:       "nonexistent",
		ByWindow:    20000,
		ExcludeFlag: DefaultExcludeFlag,
	})
	if err == nil {
		t.Fatalf("expected an error for --chrom nonexistent, got nil")
	}
	var cnf *ChromNotFoundError
	if !errors.As(err, &cnf) {
		t.Fatalf("expected *ChromNotFoundError, got %T: %v", err, err)
	}
	want := "[mosdepth] chromosome nonexistent not found"
	if err.Error() != want {
		t.Fatalf("error message mismatch.\nwant: %q\ngot:  %q", want, err.Error())
	}
}

// TestParity_BadFragLenBounds mirrors
// `run bad_frag_len_filter $exe t tests/ovl.bam --min-frag-len 10 --max-frag-len 9`.
// Upstream writes `[mosdepth] error --max-frag-len was lower than
// --min-frag-len.` to stderr and exits 2. We surface ErrBadFragLenBounds with
// the same message before any output is produced.
func TestParity_BadFragLenBounds(t *testing.T) {
	tmp := t.TempDir()
	bamPath := filepath.Join(fixtureDir(t), "ovl.bam")
	err := OpenAndRun(bamPath, Options{
		Prefix:      filepath.Join(tmp, "t"),
		MinFragLen:  10,
		MaxFragLen:  9,
		ExcludeFlag: DefaultExcludeFlag,
	})
	if !errors.Is(err, ErrBadFragLenBounds) {
		t.Fatalf("expected ErrBadFragLenBounds, got %v", err)
	}
	want := "[mosdepth] error --max-frag-len was lower than --min-frag-len."
	if err.Error() != want {
		t.Fatalf("error message mismatch.\nwant: %q\ngot:  %q", want, err.Error())
	}
	// No output files should have been created.
	if _, statErr := os.Stat(filepath.Join(tmp, "t.mosdepth.summary.txt")); statErr == nil {
		t.Errorf("summary file should not exist after a bad-frag-len error")
	}
}

// TestParity_NoPerBase — `-n` / `--no-per-base` suppresses the per-base
// output file.
func TestParity_NoPerBase(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		NoPerBase:   true,
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
	})
	if _, err := os.Stat(prefix + ".per-base.bed.gz"); !os.IsNotExist(err) {
		t.Errorf("per-base.bed.gz should not exist (err=%v)", err)
	}
	if _, err := os.Stat(prefix + ".mosdepth.summary.txt"); err != nil {
		t.Errorf("summary.txt should exist: %v", err)
	}
}

// TestParity_ChromRestrict — `-c MT` only emits MT rows.
func TestParity_ChromRestrict(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		Chrom:       "MT",
		ExcludeFlag: DefaultExcludeFlag,
		FastMode:    true,
	})
	pb := readGzLines(t, prefix+".per-base.bed.gz")
	for _, ln := range pb {
		if !strings.HasPrefix(ln, "MT\t") {
			t.Errorf("non-MT row leaked into per-base.bed.gz: %q", ln)
		}
	}
	sum := readPlainLines(t, prefix+".mosdepth.summary.txt")
	if len(sum) != 3 {
		t.Errorf("summary rows: got %d (%v), want 3", len(sum), sum)
	}
}

// TestParity_MAPQFilter exercises `-Q 60`. The MAPQ floor of 60 drops the
// MAPQ-21 mate (the 32S42M read) and keeps only the MAPQ-60 read (74M over
// [6,80)). With just one read of the pair surviving, there is no mate to
// de-duplicate, so the surviving read's coverage stands alone: depth 0 over
// [0,6), depth 1 over [6,80). This now matches upstream byte-for-byte (the
// overlap detector only fires when both mates pass the filters).
func TestParity_MAPQFilter(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		Chrom:       "MT",
		MinMAPQ:     60,
		ExcludeFlag: DefaultExcludeFlag,
	})
	got := linesWithPrefix(readGzLines(t, prefix+".per-base.bed.gz"), "MT\t")
	want := []string{
		"MT\t0\t6\t0",
		"MT\t6\t80\t1",
		"MT\t80\t16569\t0",
	}
	if !equalLines(got, want) {
		t.Fatalf("MAPQ60 MT per-base mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(got, "\n"))
	}
}

// TestParity_FlagExclude mirrors `-F 4` (exclude only unmapped reads). Both MT
// reads pass, so overlap-pair de-duplication collapses their shared bases to a
// single depth-1 run over [0,80) — identical to the default overlapM result.
func TestParity_FlagExclude(t *testing.T) {
	prefix := runParity(t, "ovl.bam", Options{
		Chrom:       "MT",
		ExcludeFlag: 4,
	})
	got := linesWithPrefix(readGzLines(t, prefix+".per-base.bed.gz"), "MT\t")
	want := []string{
		"MT\t0\t80\t1",
		"MT\t80\t16569\t0",
	}
	if !equalLines(got, want) {
		t.Fatalf("flag-exclude MT per-base mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(got, "\n"))
	}
}

// TestParity_FragmentMode mirrors
// `run fragment_mode $exe t --fragment-mode tests/full-fragment-pairs.bam`.
// Implemented: full byte-for-byte validation against the upstream binary
// lives in TestUpstream_FragmentMode_Parity. This case pins the exact
// fragment-coverage runs so a regression fails even when the upstream binary
// is unavailable.
func TestParity_FragmentMode(t *testing.T) {
	prefix := runParity(t, "full-fragment-pairs.bam", Options{
		FragmentMode: true,
		ExcludeFlag:  DefaultExcludeFlag,
	})
	got := readGzLines(t, prefix+".per-base.bed.gz")
	want := []string{
		"chr22:20000000-23000000\t0\t17318\t0",
		"chr22:20000000-23000000\t17318\t17320\t1",
		"chr22:20000000-23000000\t17320\t17420\t2",
		"chr22:20000000-23000000\t17420\t17756\t1",
		"chr22:20000000-23000000\t17756\t52130\t0",
		"chr22:20000000-23000000\t52130\t52135\t1",
		"chr22:20000000-23000000\t52135\t52235\t2",
		"chr22:20000000-23000000\t52235\t52546\t1",
		"chr22:20000000-23000000\t52546\t3000001\t0",
	}
	if !equalLines(got, want) {
		t.Fatalf("fragment-mode per-base mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(got, "\n"))
	}
}

// TestParity_Quantized mirrors `-q 0:1:4`. Implemented: byte-for-byte
// validation against the upstream binary lives in
// TestUpstream_Quantize_Parity. This case pins the MT segments so a
// regression fails without the upstream binary.
func TestParity_Quantized(t *testing.T) {
	quants, err := ParseQuantize("0:1:4")
	if err != nil {
		t.Fatalf("ParseQuantize: %v", err)
	}
	prefix := runParity(t, "ovl.bam", Options{
		FastMode:    true,
		Chrom:       "MT",
		Quantize:    quants,
		ExcludeFlag: DefaultExcludeFlag,
	})
	got := linesWithPrefix(readGzLines(t, prefix+".quantized.bed.gz"), "MT\t")
	want := []string{
		"MT\t0\t80\t1:4",
		"MT\t80\t16569\t0:1",
	}
	if !equalLines(got, want) {
		t.Fatalf("quantize MT mismatch.\nwant:\n%s\ngot:\n%s",
			strings.Join(want, "\n"), strings.Join(got, "\n"))
	}
}

// TestParity_D4RoundTrip mirrors `--d4` against a real fixture BAM. Upstream
// mosdepth is Nim, so there is no C binary to diff against; instead we assert
// the per-base depths decoded from the D4 file match the per-base BED depths
// produced for the same input (round-trip parity).
func TestParity_D4RoundTrip(t *testing.T) {
	bamFile := "ovl.bam"
	// ovl.bam declares the full GRCh37 reference but its reads only cover
	// MT. A dense D4 track over the whole genome would be multi-gigabyte, so
	// we scope this round-trip to MT via --chrom. This is purely a test
	// constraint; the writer handles all declared chromosomes.
	bedPrefix := runParity(t, bamFile, Options{ExcludeFlag: DefaultExcludeFlag, Chrom: "MT"})
	d4Prefix := runParity(t, bamFile, Options{ExcludeFlag: DefaultExcludeFlag, D4Output: true, Chrom: "MT"})

	bedLines := readGzLines(t, bedPrefix+".per-base.bed.gz")
	// Collapse the BED into per-chrom dense arrays keyed by the end coord of
	// the last run (the chromosome length).
	lengths := map[string]int{}
	for _, ln := range bedLines {
		f := strings.Split(ln, "\t")
		if len(f) < 4 {
			continue
		}
		end, err := strconv.Atoi(f[2])
		if err != nil {
			continue
		}
		if end > lengths[f[0]] {
			lengths[f[0]] = end
		}
	}

	r, err := openD4Reader(d4Prefix + ".per-base.d4")
	if err != nil {
		t.Fatalf("openD4Reader: %v", err)
	}
	for chrom, length := range lengths {
		want := make([]int32, length)
		for _, ln := range bedLines {
			f := strings.Split(ln, "\t")
			if len(f) < 4 || f[0] != chrom {
				continue
			}
			start, _ := strconv.Atoi(f[1])
			end, _ := strconv.Atoi(f[2])
			depth, _ := strconv.Atoi(f[3])
			for p := start; p < end && p < length; p++ {
				want[p] = int32(depth)
			}
		}
		got, err := r.chromDepths(chrom)
		if err != nil {
			t.Fatalf("chromDepths(%q): %v", chrom, err)
		}
		// The D4 array spans the full reference length; compare the prefix
		// the BED covers (BED omits trailing zero runs past the last event
		// only when the run is non-zero — emitRuns always flushes to refLen,
		// so lengths should match, but guard anyway).
		if len(got) < length {
			t.Fatalf("%s: D4 length %d < BED length %d", chrom, len(got), length)
		}
		for i := 0; i < length; i++ {
			if got[i] != want[i] {
				t.Fatalf("%s pos %d: D4 depth %d, BED depth %d", chrom, i, got[i], want[i])
			}
		}
	}
}

// TestParity_EmptyTids mirrors
// `run empty_tids $exe t -n --thresholds 1,5 --by tests/empty-tids.bed tests/empty-tids.bam`.
func TestParity_EmptyTids(t *testing.T) {
	prefix := runParity(t, "empty-tids.bam", Options{
		ByBED:       filepath.Join(fixtureDir(t), "empty-tids.bed"),
		NoPerBase:   true,
		Thresholds:  []int{1, 5},
		ExcludeFlag: DefaultExcludeFlag,
	})
	if _, err := os.Stat(prefix + ".thresholds.bed.gz"); err != nil {
		t.Errorf("thresholds.bed.gz missing: %v", err)
	}
	if _, err := os.Stat(prefix + ".mosdepth.summary.txt"); err != nil {
		t.Errorf("summary.txt missing: %v", err)
	}
}

// TestParity_IndexFiles_Csi confirms that, like upstream mosdepth, we now
// emit a .csi (not a .tbi) alongside each bgzipped BED output, and that the
// .csi is a structurally valid index of its companion .bed.gz.
func TestParity_IndexFiles_Csi(t *testing.T) {
	// Restrict to the small MT contig so the window sweep stays cheap; the
	// CSI's structure and round-trip behaviour are contig-independent.
	prefix := runParity(t, "ovl.bam", Options{
		Chrom:       "MT",
		ByWindow:    100,
		FastMode:    true,
		ExcludeFlag: DefaultExcludeFlag,
	})
	dataPath := prefix + ".regions.bed.gz"
	if _, err := os.Stat(dataPath + ".csi"); err != nil {
		t.Fatalf(".csi missing: %v", err)
	}
	if _, err := os.Stat(dataPath + ".tbi"); err == nil {
		t.Errorf("unexpected .tbi emitted; upstream emits .csi only")
	}
	csi, err := tabix.ReadCSIFile(dataPath + ".csi")
	if err != nil {
		t.Fatalf("ReadCSIFile: %v", err)
	}
	if len(csi.Refs) == 0 {
		t.Fatalf("csi has no reference entries")
	}
}
