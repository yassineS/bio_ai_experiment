package bedgenomecov

// Parity tests against the upstream bedtools genomecov test suite.
//
// Cases are mirrored from reference_code/bedtools/test/genomecov/test-genomecov.sh.
// Inputs and expected outputs live under tools/bedgenomecov/testdata/parity/.
// BAM/SAM/CRAM input is supported (RunBAM / the `-ibam` flag): the upstream
// SAM fixtures and the htsutil-built BAMs (y.bam/empty.bam/merged.bam) plus the
// empty CRAM and its reference are vendored under testdata/parity/aln, and every
// alignment case asserts byte-for-byte against the upstream `bedtools genomecov`
// golden generated from bedtools v2.31.1. The BED-input tests (t11/t12/t13) and
// a BED12 `-split` case are exercised against goldens too.

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func readGenomecovParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runGenomecovParity(t *testing.T, bedFile, genomeFile string, opts Options) []byte {
	t.Helper()
	bed := readGenomecovParity(t, bedFile)
	gen := readGenomecovParity(t, genomeFile)
	g, err := ReadGenome(bytes.NewReader(gen))
	if err != nil {
		t.Fatalf("ReadGenome %s: %v", genomeFile, err)
	}
	var out bytes.Buffer
	if err := Run(bytes.NewReader(bed), g, &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.Bytes()
}

// runAlnParity reads a vendored SAM/BAM/CRAM fixture from testdata/parity/aln,
// runs it through RunBAM (the `-ibam` path, genome taken from the @SQ header),
// and asserts byte-for-byte equality against the vendored upstream golden.
// reference, when non-empty, names a CRAM decode FASTA under the same dir.
func runAlnParity(t *testing.T, fixture, golden, reference string, opts Options) {
	t.Helper()
	dir := filepath.Join("..", "..", "testdata", "parity", "aln")
	in, err := os.ReadFile(filepath.Join(dir, fixture))
	if err != nil {
		t.Fatalf("read fixture %s: %v", fixture, err)
	}
	if reference != "" {
		opts.CRAMReference = filepath.Join(dir, reference)
	}
	var out bytes.Buffer
	if err := RunBAM(bytes.NewReader(in), &out, opts); err != nil {
		t.Fatalf("RunBAM(%s): %v", fixture, err)
	}
	want, err := os.ReadFile(filepath.Join(dir, golden))
	if err != nil {
		t.Fatalf("read golden %s: %v", golden, err)
	}
	if !bytes.Equal(out.Bytes(), want) {
		t.Fatalf("mismatch for %s.\nwant:\n%s\ngot:\n%s", fixture, want, out.Bytes())
	}
}

// genomecov.t1..t10 — BAM/SAM input via RunBAM (`-ibam`), genome from the @SQ
// header. The upstream SAM fixtures are vendored under testdata/parity/aln (the
// htsutil-built BAMs are vendored too, for y.bam/merged.bam); every case
// asserts byte-for-byte against the upstream `bedtools genomecov` golden.
// t6/t7 each exposed and fixed a real parity bug on our side (the CIGAR-D
// split and the -dz 0-based offset).
func TestParity_Genomecov_T1to10_BAMInputs(t *testing.T) {
	t.Run("t1_three_blocks_bg", func(t *testing.T) {
		runAlnParity(t, "three_blocks.sam", "t1_three_blocks_bg.expected.txt", "",
			Options{Mode: ModeBedGraph})
	})
	t.Run("t2_three_blocks_bg_split", func(t *testing.T) {
		runAlnParity(t, "three_blocks.sam", "t2_three_blocks_bg_split.expected.txt", "",
			Options{Mode: ModeBedGraph, Split: true})
	})
	t.Run("t3_three_blocks_bga_split", func(t *testing.T) {
		runAlnParity(t, "three_blocks.sam", "t3_three_blocks_bga_split.expected.txt", "",
			Options{Mode: ModeBedGraphAll, Split: true})
	})
	t.Run("t4_merged_bga", func(t *testing.T) {
		runAlnParity(t, "merged.bam", "t4_merged_bga.expected.txt", "",
			Options{Mode: ModeBedGraphAll})
	})
	t.Run("t5_merged_bga_split", func(t *testing.T) {
		runAlnParity(t, "merged.bam", "t5_merged_bga_split.expected.txt", "",
			Options{Mode: ModeBedGraphAll, Split: true})
	})
	t.Run("t6_three_blocks_dz_split", func(t *testing.T) {
		runAlnParity(t, "three_blocks.sam", "t6_three_blocks_dz_split.expected.txt", "",
			Options{Mode: ModePerBaseNonZero, Split: true})
	})
	t.Run("t7_sam_w_del_bg", func(t *testing.T) {
		runAlnParity(t, "sam-w-del.sam", "t7_sam_w_del_bg.expected.txt", "",
			Options{Mode: ModeBedGraph})
	})
	t.Run("t8_y_hist", func(t *testing.T) {
		runAlnParity(t, "y.bam", "t8_y_hist.expected.txt", "",
			Options{Mode: ModeHistogram})
	})
	t.Run("t9_y_bg", func(t *testing.T) {
		runAlnParity(t, "y.bam", "t9_y_bg.expected.txt", "",
			Options{Mode: ModeBedGraph})
	})
	t.Run("t10_y_bga", func(t *testing.T) {
		runAlnParity(t, "y.bam", "t10_y_bga.expected.txt", "",
			Options{Mode: ModeBedGraphAll})
	})
}

// genomecov.t11 — histogram (default) over y.bed, including chroms in the
// genome that have no coverage.
func TestParity_Genomecov_T11_Histogram(t *testing.T) {
	got := runGenomecovParity(t, "y.bed", "genome.txt", Options{Mode: ModeHistogram, Scale: 1.0})
	want := readGenomecovParity(t, "t11_hist.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t12 — `-bg`: non-zero runs of constant depth.
func TestParity_Genomecov_T12_BedGraph(t *testing.T) {
	got := runGenomecovParity(t, "y.bed", "genome.txt", Options{Mode: ModeBedGraph, Scale: 1.0})
	want := readGenomecovParity(t, "t12_bg.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t13 — `-bga`: include zero-depth runs too.
func TestParity_Genomecov_T13_BedGraphAll(t *testing.T) {
	got := runGenomecovParity(t, "y.bed", "genome.txt", Options{Mode: ModeBedGraphAll, Scale: 1.0})
	want := readGenomecovParity(t, "t13_bga.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov -split on BED12 input: a 3-block record contributes coverage
// over each block ([0,10)/[20,30)/[40,50)) rather than the whole [0,50)
// span. Upstream's -split tests use BAM, but -split is equally valid for
// BED12 input; the expected output is generated from bedtools v2.31.1.
func TestParity_Genomecov_SplitBED12(t *testing.T) {
	got := runGenomecovParity(t, "split_blocks.bed", "split.genome",
		Options{Mode: ModeBedGraphAll, Scale: 1.0, Split: true})
	want := readGenomecovParity(t, "split_bga.expected.bed")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t14 — `-pc` paired-end coverage and t15 `-fs` fragment size are
// BAM-only. Both are implemented (RunBAM via -pc / -fs) and covered by the
// SAM-input unit tests TestRunBAM_PairedCoverage / TestRunBAM_FragmentSize
// plus a live cross-check against bedtools v2.31.1; the upstream BAM fixtures
// (built by htsutil) are not vendored.
func TestParity_Genomecov_T14_PairedEnd(t *testing.T) {
	runAlnParity(t, "pair-chip.sam", "t14_pairchip_pc.expected.txt", "",
		Options{Mode: ModeBedGraph, PairedCoverage: true})
}
func TestParity_Genomecov_T15_FragmentSize(t *testing.T) {
	runAlnParity(t, "chip.sam", "t15_chip_fs.expected.txt", "",
		Options{Mode: ModeBedGraph, FragmentSize: 100})
}
func TestParity_Genomecov_T16_EmptyBAM(t *testing.T) {
	runAlnParity(t, "empty.bam", "t16_empty_bam.expected.txt", "",
		Options{Mode: ModeHistogram})
}

// genomecov.t17 — an empty CRAM. CRAM input is now wired (RunBAM auto-detects
// the format via alnio.NewReaderWithReference); the empty.cram fixture decodes
// against the vendored test_ref.fa reference, matching upstream's
// `CRAM_REFERENCE=test_ref.fa bedtools genomecov -ibam empty.cram`.
func TestParity_Genomecov_T17_EmptyCRAM(t *testing.T) {
	runAlnParity(t, "empty.cram", "t17_empty_cram.expected.txt", "test_ref.fa",
		Options{Mode: ModeHistogram})
}

// genomecov.t18 — a deep SAM (1,000,000 reads stacked at c1:1, 100M each).
// The 33 MB SAM is generated at test time (mirroring upstream's mk-deep.py)
// rather than vendored; only the first per-base line is asserted, matching
// upstream's `bedtools genomecov -d -ibam deep.sam | head -1`.
func TestParity_Genomecov_T18_DeepSAM(t *testing.T) {
	var sb strings.Builder
	sb.WriteString("@HD\tVN:1.0\tSO:coordinate\n")
	sb.WriteString("@SQ\tSN:c1\tAS:genome.txt\tLN:100\n")
	for i := 0; i < 1000000; i++ {
		fmt.Fprintf(&sb, "r%d\t0\tc1\t1\t100\t100M\t*\t0\t0\t*\t*\n", i)
	}
	var out bytes.Buffer
	if err := RunBAM(strings.NewReader(sb.String()), &out, Options{Mode: ModePerBase}); err != nil {
		t.Fatalf("RunBAM(deep): %v", err)
	}
	firstLine := out.String()
	if idx := strings.IndexByte(firstLine, '\n'); idx >= 0 {
		firstLine = firstLine[:idx+1]
	}
	want, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", "aln", "t18_deep_d_head1.expected.txt"))
	if err != nil {
		t.Fatalf("read golden: %v", err)
	}
	if firstLine != string(want) {
		t.Fatalf("mismatch.\nwant: %q\ngot:  %q", want, firstLine)
	}
}
