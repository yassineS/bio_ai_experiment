package bcftools

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// withTinyMpileupBlocks shrinks the block-streaming grid so a small
// synthetic input actually splits into multiple blocks (production uses
// 4 Mbp / 50 kbp, far larger than any unit-test reference). The left and
// right flanks are set to the same width (the production right flank equals
// the left flank). It restores the production values on cleanup.
func withTinyMpileupBlocks(t *testing.T, width, flank int) {
	t.Helper()
	withTinyMpileupBlocksLR(t, width, flank, flank)
}

// withTinyMpileupBlocksLR is withTinyMpileupBlocks with independently
// chosen left and right flanks, so a test can shrink one flank to zero to
// prove the OTHER flank is load-bearing (e.g. a right-edge indel whose
// downstream mate only co-resides because of the right flank).
func withTinyMpileupBlocksLR(t *testing.T, width, leftFlank, rightFlank int) {
	t.Helper()
	oldW, oldF, oldR := mpileupBlockWidth, mpileupBlockFlank, mpileupBlockRightFlank
	mpileupBlockWidth = width
	mpileupBlockFlank = leftFlank
	mpileupBlockRightFlank = rightFlank
	t.Cleanup(func() {
		mpileupBlockWidth = oldW
		mpileupBlockFlank = oldF
		mpileupBlockRightFlank = oldR
	})
}

// TestMpileupBlockedMatchesWholeBuffer is the core regression guard for
// task #51: the block-streaming indexed path (which re-reads each block's
// reads and emits only the block's columns) must produce byte-identical
// VCF records to the whole-buffer linear path, INCLUDING when a read pair
// straddles a block boundary. The left flank guarantees the straddling
// mate co-resides so the overlap-merge, BAQ and snapshot state at the
// block start matches a whole-contig run.
func TestMpileupBlockedMatchesWholeBuffer(t *testing.T) {
	// Tiny grid: 100bp blocks, 60bp flank — large enough that a 30bp read
	// plus its mate fit inside one flank, small enough that a 600bp
	// reference splits into several blocks.
	withTinyMpileupBlocks(t, 100, 60)

	dir := t.TempDir()
	const refLen = 600
	famPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(famPath, []byte(">chr1\n"+strings.Repeat("A", refLen)+"\n"), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	if err := os.WriteFile(famPath+".fai",
		[]byte(fmt.Sprintf("chr1\t%d\t6\t%d\t%d\n", refLen, refLen, refLen+1)), 0o644); err != nil {
		t.Fatalf("write fai: %v", err)
	}

	// Build coordinate-sorted reads. Several proper read pairs are placed so
	// that at least one pair straddles each 100bp block edge (100, 200, ...):
	// a forward mate ending just past the edge and its reverse mate starting
	// just before it, overlapping. The overlapping bases exercise
	// tweakOverlapQuality across the boundary. A C at column ~95/195/... makes
	// a variant-ish site near each edge.
	seq30 := func(varAt int) string {
		b := []byte(strings.Repeat("A", 30))
		if varAt >= 0 && varAt < 30 {
			b[varAt] = 'C'
		}
		return string(b)
	}
	qual30 := strings.Repeat("?", 30)
	lines := []string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:600",
		"@RG\tID:rg1\tSM:sample1",
	}
	// For each block edge E in {100,200,300,400,500}, lay an overlapping
	// proper pair: mate1 (fwd) at E-15 spanning [E-15, E+15), mate2 (rev) at
	// E-5 spanning [E-5, E+25). They overlap on [E-5, E+15). FLAG 99 = paired
	// +proper+mate-reverse+read1; FLAG 147 = paired+proper+reverse+read2.
	recID := 0
	mk := func(pos int, flag int, mpos, tlen, varCol int) string {
		recID++
		return fmt.Sprintf("p%d\t%d\tchr1\t%d\t60\t30M\t=\t%d\t%d\t%s\t%s\tRG:Z:rg1",
			recID, flag, pos, mpos, tlen, seq30(varCol), qual30)
	}
	type rec struct {
		pos  int
		line string
	}
	var recs []rec
	for _, e := range []int{100, 200, 300, 400, 500} {
		m1pos := e - 15 // 1-based
		m2pos := e - 5
		// var column inside mate1 at the edge (col index in the 30bp read)
		recs = append(recs, rec{m1pos, mk(m1pos, 99, m2pos, 45, 15)})  // C at edge-ish
		recs = append(recs, rec{m2pos, mk(m2pos, 147, m1pos, -45, 5)}) // C overlapping
	}
	// A standalone read far from any edge for baseline coverage.
	recs = append(recs, rec{50, fmt.Sprintf("s1\t0\tchr1\t50\t60\t30M\t*\t0\t0\t%s\t%s\tRG:Z:rg1", seq30(-1), qual30)})
	recs = append(recs, rec{350, fmt.Sprintf("s2\t0\tchr1\t350\t60\t30M\t*\t0\t0\t%s\t%s\tRG:Z:rg1", seq30(-1), qual30)})

	// Sort by pos (coordinate order) so the BAM is valid and the index works.
	for i := 1; i < len(recs); i++ {
		for j := i; j > 0 && recs[j-1].pos > recs[j].pos; j-- {
			recs[j-1], recs[j] = recs[j], recs[j-1]
		}
	}
	for _, r := range recs {
		lines = append(lines, r.line)
	}
	lines = append(lines, "")
	samText := strings.Join(lines, "\n")

	bamPath := writeIndexedTestBAM(t, dir, samText)
	region := "chr1:1-600"

	// Indexed path → block-streaming driver (tiny grid forces multiple
	// blocks with a straddling pair on every edge).
	indexed := mpileupDataRecords(t, MpileupOptions{
		Inputs:   []string{bamPath},
		FastaRef: famPath,
		Regions:  []string{region},
	})

	// Linear whole-buffer path: copy the BAM without an index sidecar so
	// openBAMRegionReader returns nil and the linear scan runs.
	linDir := t.TempDir()
	linBAM := filepath.Join(linDir, "lin.bam")
	if err := copyFileForTest(bamPath, linBAM); err != nil {
		t.Fatalf("copy bam: %v", err)
	}
	linear := mpileupDataRecords(t, MpileupOptions{
		Inputs:   []string{linBAM},
		FastaRef: famPath,
		Regions:  []string{region},
	})

	if len(indexed) == 0 {
		t.Fatalf("block-streaming path produced no records")
	}
	if len(indexed) != len(linear) {
		t.Fatalf("record count differs: blocked=%d wholeBuffer=%d", len(indexed), len(linear))
	}
	for i := range indexed {
		if indexed[i] != linear[i] {
			t.Fatalf("record %d differs across block boundary:\n blocked    : %s\n wholeBuffer: %s",
				i, indexed[i], linear[i])
		}
	}
}

// TestMpileupBlockedIndelAcrossRightEdge is a CI-portable guard that the
// block driver runs the INDEL pass byte-identically to the whole-buffer path
// when a real, firing indel sits just LEFT of a block boundary and its reads
// extend DOWNSTREAM past that boundary into the indel-realignment window. It
// reuses the vendored indel-AD.2 fixture (a homopolymer insertion at 11:75
// with paired reads spanning well past the call), placing several tiny block
// grids so the boundary lands at/around the indel. With the right flank the
// downstream-extending reads co-reside and the block emits the same INDEL
// record (QS/PL/AD/IDV) as the linear scan. This locks the indel pass to the
// block driver; the right-flank's effect on the harder downstream-second-mate
// case is covered by TestMpileupBlockedRightEdgeIndelRealData (real GIAB
// data, skipped when the fixture is absent).
func TestMpileupBlockedIndelAcrossRightEdge(t *testing.T) {
	ref := mpileupFixture(t, "indel-AD.2.fa")
	mpileupFixture(t, "indel-AD.2.fa.fai")
	bam := mpileupFixture(t, "indel-AD.2.bam")
	mpileupFixture(t, "indel-AD.2.bam.bai")

	// Indexed copy (with .bai → block driver) and a no-index copy (→ linear).
	dir := t.TempDir()
	idxBAM := filepath.Join(dir, "in.bam")
	if err := copyFileForTest(bam, idxBAM); err != nil {
		t.Fatalf("copy bam: %v", err)
	}
	if err := copyFileForTest(bam+".bai", idxBAM+".bai"); err != nil {
		t.Fatalf("copy bai: %v", err)
	}
	linBAM := filepath.Join(t.TempDir(), "lin.bam")
	if err := copyFileForTest(bam, linBAM); err != nil {
		t.Fatalf("copy lin bam: %v", err)
	}

	region := "11:1-150"
	linear := mpileupDataRecords(t, MpileupOptions{
		Inputs: []string{linBAM}, FastaRef: ref, Regions: []string{region}, Annotate: "AD",
	})
	if len(linear) == 0 {
		t.Fatalf("linear path produced no records")
	}
	// Sanity: the fixture must actually fire an INDEL, else the test is vacuous.
	sawIndel := false
	for _, l := range linear {
		if strings.Contains(l, "\tINDEL") || strings.HasPrefix(l, "11\t75\t") && strings.Contains(l, "INDEL") {
			sawIndel = true
			break
		}
	}
	if !sawIndel {
		t.Fatalf("indel-AD.2 fixture produced no INDEL record; test would be vacuous")
	}

	// Several block widths place the boundary AT, just BEFORE, and just AFTER
	// the indel column (75), each with the production-shape left+right flank.
	// All must reproduce the linear records byte-for-byte.
	for _, w := range []int{74, 75, 76, 78, 80, 100} {
		withTinyMpileupBlocksLR(t, w, 60, 60)
		got := mpileupDataRecords(t, MpileupOptions{
			Inputs: []string{idxBAM}, FastaRef: ref, Regions: []string{region}, Annotate: "AD",
		})
		if len(got) != len(linear) {
			t.Fatalf("block width %d: record count blocked=%d wholeBuffer=%d", w, len(got), len(linear))
		}
		for i := range got {
			if got[i] != linear[i] {
				t.Fatalf("block width %d: record %d differs across right-edge indel:\n blocked    : %s\n wholeBuffer: %s",
					w, i, got[i], linear[i])
			}
		}
	}
}

// giabRealDataFixture locates the real GIAB BAM + hs37d5 reference used by
// the host/bioval validation environment, returning ("","") (so the caller
// skips) when they are not present. It mirrors the submodule guard on the
// live parity tests: the GIAB inputs are multi-GB and not vendored, so this
// regression runs only where the operator has staged them.
func giabRealDataFixture() (bam, ref string) {
	bamCandidates := []string{
		os.Getenv("BIOAI_GIAB_BAM"),
		"/work/_giab/giab_b37.bam",
	}
	refCandidates := []string{
		os.Getenv("BIOAI_GIAB_REF"),
		"/work/pipeline/.fixtures/giab/hs37d5.fa.gz",
	}
	pick := func(cs []string, needSibling string) string {
		for _, c := range cs {
			if c == "" {
				continue
			}
			if _, err := os.Stat(c); err != nil {
				continue
			}
			if needSibling != "" {
				if _, err := os.Stat(c + needSibling); err != nil {
					continue
				}
			}
			return c
		}
		return ""
	}
	bam = pick(bamCandidates, ".bai")
	ref = pick(refCandidates, "")
	if bam == "" || ref == "" {
		return "", ""
	}
	return bam, ref
}

// TestMpileupBlockedRightEdgeIndelRealData reproduces the exact byte-exact
// failure the right-flank fix closes (task #51, attempt-2): the GIAB indel
// at 20:30142484 whose per-allele QS / indel PL shifted when a 200 kbp block
// boundary landed on/just after it, because the block dropped the indel's
// downstream second mates (which the smart-overlap merge needs to converge to
// the whole-contig value). With the symmetric right flank the block driver's
// records match the whole-buffer (linear) path byte-for-byte across every
// grid alignment around the indel.
//
// It is skipped when the multi-GB GIAB BAM / hs37d5 reference are not staged
// (set BIOAI_GIAB_BAM / BIOAI_GIAB_REF, or stage them at the bioval paths) —
// the same "environmental, not a defect" guard the live upstream-parity tests
// use for the un-vendored submodule sources.
func TestMpileupBlockedRightEdgeIndelRealData(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping real-data GIAB regression in -short mode")
	}
	bam, ref := giabRealDataFixture()
	if bam == "" {
		t.Skip("GIAB BAM/reference not staged; set BIOAI_GIAB_BAM and BIOAI_GIAB_REF to enable the real-data right-edge indel regression")
	}

	const indel = 30142484 // 1-based; the homopolymer indel that shifted
	// Linear whole-buffer reference: a tight window around the indel with no
	// block grid effect (the whole-buffer driver reads it all at once).
	linRegion := fmt.Sprintf("20:%d-%d", indel-2000, indel+2000)
	linear := mpileupDataRecords(t, MpileupOptions{
		Inputs: []string{bam}, FastaRef: ref, Regions: []string{linRegion},
	})
	if len(linear) == 0 {
		t.Fatalf("linear path produced no records for %s", linRegion)
	}

	// Grid-shift the production 200 kbp block boundary onto / around the
	// indel: a window whose start S makes the boundary S-1+200000 land at
	// EDGE. For each EDGE near the indel, the block path must match the
	// whole-buffer path over the SAME window byte-for-byte.
	for _, edge := range []int{indel - 2, indel - 1, indel, indel + 1, indel + 2, indel + 4} {
		s := edge - mpileupBlockWidth // window start (1-based)
		e := edge + mpileupBlockWidth
		reg := fmt.Sprintf("20:%d-%d", s, e)
		blockOut := mpileupDataRecords(t, MpileupOptions{
			Inputs: []string{bam}, FastaRef: ref, Regions: []string{reg},
		})
		// Linear comparison over the SAME window (no index sidecar would be
		// needed, but the whole-buffer path is selected by NOT taking the
		// block grid; here we compare the block run against the same-window
		// whole-buffer run by disabling the block driver via a huge width).
		oldW := mpileupBlockWidth
		mpileupBlockWidth = 1 << 30 // one block covers the whole window → whole-buffer behaviour
		whole := mpileupDataRecords(t, MpileupOptions{
			Inputs: []string{bam}, FastaRef: ref, Regions: []string{reg},
		})
		mpileupBlockWidth = oldW

		if len(blockOut) != len(whole) {
			t.Fatalf("edge %d (%s): record count block=%d whole=%d", edge, reg, len(blockOut), len(whole))
		}
		for i := range blockOut {
			if blockOut[i] != whole[i] {
				t.Fatalf("edge %d (%s): record %d differs (right-edge indel not byte-exact):\n block: %s\n whole: %s",
					edge, reg, i, blockOut[i], whole[i])
			}
		}
	}

	// Region-end-on-indel clamp (the attempt-3 corner). When the requested
	// region's RIGHT boundary Y lands EXACTLY on the indel, the block
	// driver's right flank must NOT read past Y: the last emitted column at Y
	// must use exactly the reads a region-bounded query [X, Y] sees. Before
	// the clamp the block driver read mpileupBlockRightFlank bases of
	// downstream context PAST Y and the indel's per-allele QS / PL diverged
	// from the region-bounded value.
	//
	// The reference is the genuinely region-bounded result: one giant block
	// (so blockEnd == wEnd0 == Y) with the right flank forced to zero, so it
	// reads nothing past Y regardless of the clamp. The production-grid block
	// run (its real right flank, clamped to the region end) must reproduce it
	// byte-for-byte. Without the clamp the production run reads past Y on its
	// final block and this assertion fails.
	for _, endY := range []int{indel, indel + 1, indel + 2} {
		reg := fmt.Sprintf("20:%d-%d", indel-2000, endY)
		blockOut := mpileupDataRecords(t, MpileupOptions{
			Inputs: []string{bam}, FastaRef: ref, Regions: []string{reg},
		})
		oldW, oldR := mpileupBlockWidth, mpileupBlockRightFlank
		mpileupBlockWidth = 1 << 30 // one block over the whole region
		mpileupBlockRightFlank = 0  // read nothing past the region end Y
		regionBounded := mpileupDataRecords(t, MpileupOptions{
			Inputs: []string{bam}, FastaRef: ref, Regions: []string{reg},
		})
		mpileupBlockWidth, mpileupBlockRightFlank = oldW, oldR
		if len(blockOut) != len(regionBounded) {
			t.Fatalf("region-end Y=%d (%s): record count block=%d regionBounded=%d", endY, reg, len(blockOut), len(regionBounded))
		}
		for i := range blockOut {
			if blockOut[i] != regionBounded[i] {
				t.Fatalf("region-end Y=%d (%s): record %d differs (right flank read past region end):\n block        : %s\n regionBounded: %s",
					endY, reg, i, blockOut[i], regionBounded[i])
			}
		}
	}
}

// TestMpileupBlockedAnnotateMatchesWholeBuffer repeats the boundary check
// with -a AD,DP,SP and -a AD,ADF,ADR so the per-allele FORMAT counts (which
// depend on the overlap-merge and depth state) also stay byte-identical
// across block edges.
func TestMpileupBlockedAnnotateMatchesWholeBuffer(t *testing.T) {
	withTinyMpileupBlocks(t, 100, 60)

	dir := t.TempDir()
	const refLen = 400
	famPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(famPath, []byte(">chr1\n"+strings.Repeat("A", refLen)+"\n"), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	if err := os.WriteFile(famPath+".fai",
		[]byte(fmt.Sprintf("chr1\t%d\t6\t%d\t%d\n", refLen, refLen, refLen+1)), 0o644); err != nil {
		t.Fatalf("write fai: %v", err)
	}

	qual30 := strings.Repeat("?", 30)
	seqVar := func() string { b := []byte(strings.Repeat("A", 30)); b[10] = 'C'; return string(b) }
	lines := []string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:400",
		"@RG\tID:rg1\tSM:sample1",
	}
	recID := 0
	type rec struct {
		pos  int
		line string
	}
	var recs []rec
	for _, e := range []int{100, 200, 300} {
		recID++
		m1 := e - 15
		m2 := e - 5
		recs = append(recs, rec{m1, fmt.Sprintf("q%d\t99\tchr1\t%d\t60\t30M\t=\t%d\t45\t%s\t%s\tRG:Z:rg1", recID, m1, m2, seqVar(), qual30)})
		recs = append(recs, rec{m2, fmt.Sprintf("q%d\t147\tchr1\t%d\t60\t30M\t=\t%d\t-45\t%s\t%s\tRG:Z:rg1", recID, m2, m1, seqVar(), qual30)})
	}
	for i := 1; i < len(recs); i++ {
		for j := i; j > 0 && recs[j-1].pos > recs[j].pos; j-- {
			recs[j-1], recs[j] = recs[j], recs[j-1]
		}
	}
	for _, r := range recs {
		lines = append(lines, r.line)
	}
	lines = append(lines, "")
	samText := strings.Join(lines, "\n")

	bamPath := writeIndexedTestBAM(t, dir, samText)
	linDir := t.TempDir()
	linBAM := filepath.Join(linDir, "lin.bam")
	if err := copyFileForTest(bamPath, linBAM); err != nil {
		t.Fatalf("copy bam: %v", err)
	}

	for _, annot := range []string{"AD,DP,SP", "AD,ADF,ADR"} {
		indexed := mpileupDataRecords(t, MpileupOptions{
			Inputs: []string{bamPath}, FastaRef: famPath,
			Regions: []string{"chr1:1-400"}, Annotate: annot,
		})
		linear := mpileupDataRecords(t, MpileupOptions{
			Inputs: []string{linBAM}, FastaRef: famPath,
			Regions: []string{"chr1:1-400"}, Annotate: annot,
		})
		if len(indexed) != len(linear) || len(indexed) == 0 {
			t.Fatalf("-a %s: count blocked=%d wholeBuffer=%d", annot, len(indexed), len(linear))
		}
		for i := range indexed {
			if indexed[i] != linear[i] {
				t.Fatalf("-a %s record %d differs:\n blocked    : %s\n wholeBuffer: %s",
					annot, i, indexed[i], linear[i])
			}
		}
	}
}
