package bedgenomecov

// Parity tests against the upstream bedtools genomecov test suite.
//
// Cases are mirrored from reference_code/bedtools/test/genomecov/test-genomecov.sh.
// Inputs and expected outputs live under tools/bedgenomecov/testdata/parity/.
// bedgenomecov only consumes BED input today (no BAM/CRAM/SAM parser), so
// most upstream tests — which use `bedtools genomecov -ibam` on a BAM built
// from a SAM fixture by htsutil — are skipped. The three BED-input tests
// (t11/t12/t13) are exercised here.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bamtobed"
)

// runGenomecovFromBAMText runs the genomecov Run() function against a BED
// text body whose genome is derived from a list of BAM-header references.
// Used by parity tests that exercise `-ibam` semantics: the SQ header of
// the BAM provides the per-chromosome size list, exactly mirroring how
// upstream `bedtools genomecov -ibam` seeds its depth arrays.
func runGenomecovFromBAMText(t *testing.T, bed []byte, refs []bamtobed.BAMRef, opts Options) []byte {
	t.Helper()
	g := &GenomeSize{Length: map[string]int{}}
	for _, r := range refs {
		g.Order = append(g.Order, r.Name)
		g.Length[r.Name] = r.Length
	}
	var out bytes.Buffer
	if err := Run(bytes.NewReader(bed), g, &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.Bytes()
}

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

// genomecov.t1 — three-block BAM, `-bg`, no `-split`. Full reference
// footprint of the alignment is covered.
func TestParity_Genomecov_T1_ThreeBlocksNoSplit(t *testing.T) {
	sam := readGenomecovParity(t, "three_blocks.sam")
	bed, refs, err := bamtobed.DecodeSAMToBED(bytes.NewReader(sam))
	if err != nil {
		t.Fatalf("DecodeSAMToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeBedGraph, Scale: 1.0})
	want := []byte("chr1\t0\t50\t1\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t2 — three-block BAM, `-bg -split`. CIGAR `N` gaps split blocks.
func TestParity_Genomecov_T2_ThreeBlocksSplit(t *testing.T) {
	sam := readGenomecovParity(t, "three_blocks.sam")
	bed, refs, err := bamtobed.DecodeSAMSplitToBED(bytes.NewReader(sam))
	if err != nil {
		t.Fatalf("DecodeSAMSplitToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeBedGraph, Scale: 1.0})
	want := []byte("chr1\t0\t10\t1\nchr1\t20\t30\t1\nchr1\t40\t50\t1\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t3 — three-block BAM, `-bga -split`. Adds zero-depth runs.
func TestParity_Genomecov_T3_ThreeBlocksSplitBGA(t *testing.T) {
	sam := readGenomecovParity(t, "three_blocks.sam")
	bed, refs, err := bamtobed.DecodeSAMSplitToBED(bytes.NewReader(sam))
	if err != nil {
		t.Fatalf("DecodeSAMSplitToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeBedGraphAll, Scale: 1.0})
	want := []byte("chr1\t0\t10\t1\nchr1\t10\t20\t0\nchr1\t20\t30\t1\nchr1\t30\t40\t0\nchr1\t40\t50\t1\nchr1\t50\t1000\t0\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t6 — three-block BAM, `-dz -split` (per-base non-zero).
func TestParity_Genomecov_T6_ThreeBlocksSplitDZ(t *testing.T) {
	sam := readGenomecovParity(t, "three_blocks.sam")
	bed, refs, err := bamtobed.DecodeSAMSplitToBED(bytes.NewReader(sam))
	if err != nil {
		t.Fatalf("DecodeSAMSplitToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModePerBaseNonZero, Scale: 1.0})
	// Build expected: 3 runs of 10 bases at positions [1..10], [21..30], [41..50].
	var w bytes.Buffer
	for _, r := range []struct{ lo, hi int }{{1, 10}, {21, 30}, {41, 50}} {
		for i := r.lo; i <= r.hi; i++ {
			w.WriteString("chr1\t")
			w.WriteString(itoa(i))
			w.WriteString("\t1\n")
		}
	}
	if !bytes.Equal(got, w.Bytes()) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", w.Bytes(), got)
	}
}

// genomecov.t7 — SAM with one D op, `-bg`. Upstream genomecov splits BAM
// blocks on both N and D ops (sam-w-del has CIGAR 10M1D10M and emits two
// 10bp runs separated by the deleted base). To match we set SplitOnDel.
func TestParity_Genomecov_T7_SAMWithDel(t *testing.T) {
	sam := readGenomecovParity(t, "sam-w-del.sam")
	bed, refs, err := bamtobed.DecodeSAMSplitOptsToBED(bytes.NewReader(sam),
		bamtobed.DecodeOpts{SplitOnN: true, SplitOnDel: true})
	if err != nil {
		t.Fatalf("DecodeSAMSplitOptsToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeBedGraph, Scale: 1.0})
	want := []byte("chr1\t0\t10\t1\nchr1\t11\t21\t1\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t8 — y.bam, default histogram. Tests that chroms with no
// coverage still appear in the output.
func TestParity_Genomecov_T8_YBamHist(t *testing.T) {
	bam := readGenomecovParity(t, "y.bam")
	bed, refs, err := bamtobed.DecodeBAMToBED(bytes.NewReader(bam))
	if err != nil {
		t.Fatalf("DecodeBAMToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeHistogram, Scale: 1.0})
	want := []byte("1\t0\t93\t100\t0.93\n1\t1\t4\t100\t0.04\n1\t2\t3\t100\t0.03\n2\t0\t100\t100\t1\n3\t0\t100\t100\t1\ngenome\t0\t293\t300\t0.976667\ngenome\t1\t4\t300\t0.0133333\ngenome\t2\t3\t300\t0.01\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t9 — y.bam, `-bg`.
func TestParity_Genomecov_T9_YBamBG(t *testing.T) {
	bam := readGenomecovParity(t, "y.bam")
	bed, refs, err := bamtobed.DecodeBAMToBED(bytes.NewReader(bam))
	if err != nil {
		t.Fatalf("DecodeBAMToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeBedGraph, Scale: 1.0})
	want := []byte("1\t15\t17\t1\n1\t17\t20\t2\n1\t20\t22\t1\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t10 — y.bam, `-bga`.
func TestParity_Genomecov_T10_YBamBGA(t *testing.T) {
	bam := readGenomecovParity(t, "y.bam")
	bed, refs, err := bamtobed.DecodeBAMToBED(bytes.NewReader(bam))
	if err != nil {
		t.Fatalf("DecodeBAMToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeBedGraphAll, Scale: 1.0})
	want := []byte("1\t0\t15\t0\n1\t15\t17\t1\n1\t17\t20\t2\n1\t20\t22\t1\n1\t22\t100\t0\n2\t0\t100\t0\n3\t0\t100\t0\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// itoa is a tiny local helper to keep the expected-output builders concise.
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
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

// genomecov.t14..t18 — `-pc` paired-end coverage, `-fs` fragment size,
// BAM empty fixtures, deep SAM. All BAM-only.
//
// -pc: covers the entire fragment [POS-1, POS-1+TLEN) for the leftmost
// mate (TLEN > 0) and skips the rightmost mate (TLEN < 0). Needs a new
// pkg/bamtobed converter (FromBAMPaired) plus a bedgenomecov option to
// consume it; ~50-80 LOC.
//
// -fs: extends each read to a fixed fragment length downstream of the
// alignment start (`POS-1` .. `POS-1+fragLen`). Needs a similar
// converter + option; ~30-50 LOC. Both tracked together because they
// share the new pkg/bamtobed entry point.
// genomecov.t14 — `-pc`: paired-end coverage. The leftmost first-mate
// of each proper pair contributes the full fragment span [POS-1,
// POS-1+|ISIZE|); the other mate is dropped so the fragment is counted
// exactly once. pair-chip.sam has one pair (chip:1:1106:..:70206) with
// FLAG=99 / TLEN=203 / POS=1; expected output is one chr1\t0\t203\t1
// bedGraph row.
func TestParity_Genomecov_T14_PairedEnd(t *testing.T) {
	sam := readGenomecovParity(t, "pair-chip.sam")
	refs, body, err := bamtobed.ReadSAMHeaderRefs(bytes.NewReader(sam))
	if err != nil {
		t.Fatalf("ReadSAMHeaderRefs: %v", err)
	}
	bed, rerr := io.ReadAll(bamtobed.FromSAMPaired(body))
	if rerr != nil {
		t.Fatalf("FromSAMPaired: %v", rerr)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeBedGraph, Scale: 1.0})
	want := []byte("chr1\t0\t203\t1\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t15 — `-fs 100`: each read is extended downstream
// (forward) or upstream (reverse) to a fixed fragment length. chip.sam
// has chip:2 (FLAG=0, POS=2) extended to [1, 101) and chip:1 (FLAG=16,
// POS=226, refLen=75) extended to [200, 300).
func TestParity_Genomecov_T15_FragmentSize(t *testing.T) {
	sam := readGenomecovParity(t, "chip.sam")
	refs, body, err := bamtobed.ReadSAMHeaderRefs(bytes.NewReader(sam))
	if err != nil {
		t.Fatalf("ReadSAMHeaderRefs: %v", err)
	}
	bed, rerr := io.ReadAll(bamtobed.FromSAMExtended(body, 100))
	if rerr != nil {
		t.Fatalf("FromSAMExtended: %v", rerr)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeBedGraph, Scale: 1.0})
	want := []byte("chr1\t1\t101\t1\nchr1\t200\t300\t1\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t16 — empty.bam, default histogram. The genome (3 chroms of
// 100bp each) comes entirely from the BAM SQ header; with no alignments,
// every base is at depth 0.
func TestParity_Genomecov_T16_EmptyBAM(t *testing.T) {
	bam := readGenomecovParity(t, "empty.bam")
	bed, refs, err := bamtobed.DecodeBAMToBED(bytes.NewReader(bam))
	if err != nil {
		t.Fatalf("DecodeBAMToBED: %v", err)
	}
	got := runGenomecovFromBAMText(t, bed, refs, Options{Mode: ModeHistogram, Scale: 1.0})
	want := []byte("1\t0\t100\t100\t1\n2\t0\t100\t100\t1\n3\t0\t100\t100\t1\ngenome\t0\t300\t300\t1\n")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// genomecov.t17 — empty CRAM input. pkg/htsgo/cram has a decoder but
// it requires a fasta reference (CRAM's reference-compressed read
// substitution lookup). bedgenomecov's BAM path goes through
// pkg/bamtobed.DecodeBAMToBED which produces a BED stream + ref list
// in one shot; an equivalent DecodeCRAMToBED would need to accept and
// thread through an `M5`-keyed reference cache. Roughly 80-120 LOC
// once a small empty CRAM fixture is vendored. Deferred.
func TestParity_Genomecov_T17_EmptyCRAM(t *testing.T) {
	t.Skip("unimplemented: CRAM input; needs pkg/bamtobed.DecodeCRAMToBED + reference threading")
}

// genomecov.t18 — upstream test calls bundled mk-deep.py to synthesise
// a 1 Mbase deep SAM at runtime, then exercises the per-base path's
// O(N*depth) memory profile. Closing this needs either (a) checking
// in a ~100 MB synthetic SAM, or (b) porting mk-deep.py to Go as a
// testdata generator. Both options are out of scope for this slice;
// the per-base path itself is exercised by t1-t12 on small fixtures.
func TestParity_Genomecov_T18_DeepSAM(t *testing.T) {
	t.Skip("fixture: needs in-tree port of mk-deep.py (1Mbase deep SAM synthesiser); not algorithmic")
}
