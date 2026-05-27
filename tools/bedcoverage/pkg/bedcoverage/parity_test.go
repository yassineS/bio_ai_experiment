package bedcoverage

// Parity tests against the upstream `bedtools coverage` test suite.
//
// Cases are mirrored from reference_code/bedtools/test/coverage/test-coverage.sh.
// Inputs are vendored under tools/bedcoverage/testdata/parity/ (literal copies
// of the upstream fixtures). Expected outputs are the inline heredoc strings
// from the upstream script, copied verbatim.
//
// Tests for options we do not implement (BAM input, `-split`, `-sorted` fast
// path, mean-as-float32 precision) are wrapped in t.Skip with a one-line
// rationale.

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bamtobed"
)

func readParityFixture(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runParity(t *testing.T, aFile, bFile string, opts Options) []byte {
	t.Helper()
	a := readParityFixture(t, aFile)
	b := readParityFixture(t, bFile)
	var buf bytes.Buffer
	if _, err := Coverage(bytes.NewReader(a), bytes.NewReader(b), &buf, opts); err != nil {
		t.Fatalf("Coverage: %v", err)
	}
	return buf.Bytes()
}

// coverage.t1 — BAM as A input via `-abam`. Upstream emits one BED12
// row per BAM alignment with the standard coverage suffix
// (count/bp/len/frac). t10/t11 already exercise BAM as B via
// pkg/bamtobed.FromBAM; closing this test needs the inverse (BAM as A
// in BED12 form, since A flows verbatim to the output). Roughly 60-90
// LOC: a new pkg/bamtobed.FromBAMBED12 + wiring into the parity test.
// Deferred — tracked in docs/PARITY_ROADMAP.md.
func TestParity_Coverage_T1_BAMInput(t *testing.T) {
	t.Skip("unimplemented: BAM-as-A via -abam (needs pkg/bamtobed.FromBAMBED12); tracked in docs/PARITY_ROADMAP.md")
}

// coverage.t2 — defaults: A, B BED; per-A count + bp + len + frac.
func TestParity_Coverage_T2_Default(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{})
	want := readParityFixture(t, "t2_default.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t3 — -counts.
func TestParity_Coverage_T3_Counts(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{Mode: ModeCounts})
	want := readParityFixture(t, "t3_counts.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t4 — -hist (per-A histogram + "all" footer).
func TestParity_Coverage_T4_Hist(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{Mode: ModeHist})
	want := readParityFixture(t, "t4_hist.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t5 — -d (per-base depth).
func TestParity_Coverage_T5_Depth(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{Mode: ModeDepth})
	want := readParityFixture(t, "t5_d.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t6 — -mean. Upstream prints with float32 precision
// ("1.3200001", "5.5599999"); our port uses float64 and prints the
// mathematically-equal "1.32" / "5.56". This test asserts that
// intentional non-divergence by reading the upstream expected text,
// parsing each row's mean column with strconv.ParseFloat, and
// asserting our output parses to the same value within float32 epsilon.
// Anything else (chrom, start, end, count, length, strand) must match
// exactly.
func TestParity_Coverage_T6_Mean(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{Mode: ModeMean})
	// The upstream expected output is the inline `echo` block from
	// reference_code/bedtools/test/coverage/test-coverage.sh:
	upstream := []string{
		"chr1\t20\t70\t6\t25\t+\t2.0000000",
		"chr1\t50\t100\t1\t25\t-\t2.2000000",
		"chr1\t200\t250\t3\t25\t+\t1.3200001",
		"chr2\t80\t130\t5\t25\t-\t3.0799999",
		"chr2\t150\t200\t4\t25\t+\t5.5599999",
		"chr2\t180\t230\t2\t25\t-\t3.4600000",
	}
	gotLines := strings.Split(strings.TrimRight(string(got), "\n"), "\n")
	if len(gotLines) != len(upstream) {
		t.Fatalf("row count mismatch: got %d, want %d\n%s", len(gotLines), len(upstream), got)
	}
	const eps = 1e-6
	for i, want := range upstream {
		gotF := strings.Split(gotLines[i], "\t")
		wantF := strings.Split(want, "\t")
		if len(gotF) != len(wantF) {
			t.Errorf("row %d field count mismatch: got %d, want %d (%q vs %q)",
				i, len(gotF), len(wantF), gotLines[i], want)
			continue
		}
		// First six columns must match byte-for-byte.
		for j := 0; j < 6; j++ {
			if gotF[j] != wantF[j] {
				t.Errorf("row %d col %d: got %q, want %q", i, j, gotF[j], wantF[j])
			}
		}
		// Final mean column: parse both as float64, compare.
		gotV, err := strconv.ParseFloat(gotF[6], 64)
		if err != nil {
			t.Errorf("row %d: parse got mean %q: %v", i, gotF[6], err)
			continue
		}
		wantV, err := strconv.ParseFloat(wantF[6], 64)
		if err != nil {
			t.Errorf("row %d: parse want mean %q: %v", i, wantF[6], err)
			continue
		}
		if diff := gotV - wantV; diff > eps || diff < -eps {
			t.Errorf("row %d mean: got %v, want %v (diff %v)", i, gotV, wantV, diff)
		}
	}
}

// coverage.t7 — -s (same strand).
func TestParity_Coverage_T7_SameStrand(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{SameStrand: true})
	want := readParityFixture(t, "t7_s.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t8 — -S (opposite strand).
func TestParity_Coverage_T8_OppositeStrand(t *testing.T) {
	got := runParity(t, "a.bed", "b.bed", Options{OppositeStrand: true})
	want := readParityFixture(t, "t8_S.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// coverage.t10 — A=BED, B=BAM with `-split`. The BAM is converted to BED
// blocks (one per CIGAR M-run, N breaks blocks) via bamtobed.FromBAMSplit
// before being fed to Coverage.
func TestParity_Coverage_T10_BAMSplit(t *testing.T) {
	a := readParityFixture(t, "c.bed")
	bRaw := readParityFixture(t, "three_blocks_match.bam")
	bBed := bamtobed.FromBAMSplit(bytes.NewReader(bRaw))
	var buf bytes.Buffer
	if _, err := Coverage(bytes.NewReader(a), bBed, &buf, Options{}); err != nil {
		t.Fatalf("Coverage: %v", err)
	}
	want := readParityFixture(t, "t10_bam_split.expected")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// coverage.t11 — A=BED, B=BAM, no `-split`. Full reference footprint of
// each alignment is counted via bamtobed.FromBAM.
func TestParity_Coverage_T11_BAMNoSplit(t *testing.T) {
	a := readParityFixture(t, "c.bed")
	bRaw := readParityFixture(t, "three_blocks_match.bam")
	bBed := bamtobed.FromBAM(bytes.NewReader(bRaw))
	var buf bytes.Buffer
	if _, err := Coverage(bytes.NewReader(a), bBed, &buf, Options{}); err != nil {
		t.Fatalf("Coverage: %v", err)
	}
	want := readParityFixture(t, "t11_bam_nosplit.expected")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}
