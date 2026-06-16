package cram

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// eqxRefFASTA is a tiny reference whose first bases let a read with an
// internal X (mismatch) and = (match) ops exercise the writer's per-base
// match/substitution encoding. chr1's prefix is ACGTACGT..., so a read of
// ACGAACGT against it has exactly one mismatch at position 4 (T→A).
const eqxRefFASTA = ">chr1\n" +
	"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
	">chr2\n" +
	"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n"

// eqxSAM carries reads whose CIGARs use the =/X ops an aligner run with
// --eqx emits in place of M: r1 is a perfect match (all =), r2 has one X
// substitution, r3 alternates = and X. The seq/pos/qual must survive a
// round-trip through our CRAM writer and decode identically under
// samtools (where, as with samtools' own reference-based CRAM, =/X
// collapse to M).
const eqxSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:64\n" +
	"@SQ\tSN:chr2\tLN:64\n" +
	"r1\t0\tchr1\t1\t60\t8=\t*\t0\t0\tACGTACGT\tIIIIIIII\n" +
	"r2\t0\tchr1\t1\t60\t3=1X4=\t*\t0\t0\tACGAACGT\tIIIIIIII\n" +
	"r3\t0\tchr1\t5\t60\t2=1X1=1X3=\t*\t0\t0\tCGTTCATA\tIIIIIIII\n" +
	"r4\t16\tchr2\t1\t60\t4=2X2=\t*\t0\t0\tTTTTCCGG\tIIIIIIII\n"

// TestEqXWriteParity proves our CRAM writer accepts =/X CIGAR ops (the
// --eqx aligner output) and that the records it writes decode under
// upstream `samtools view` to exactly the input SAM fields (seq, pos,
// quals, with =/X reconstructing as M just as samtools' own CRAM does).
// It cross-checks against a samtools-written CRAM of the same alignment:
// the decoded SAM must match even though the container bytes differ.
func TestEqXWriteParity(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	dir := t.TempDir()

	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(eqxRefFASTA), 0o644); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	if out, err := exec.Command(samtools, "faidx", refPath).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx: %v\n%s", err, out)
	}

	samPath := filepath.Join(dir, "eqx.sam")
	if err := os.WriteFile(samPath, []byte(eqxSAM), 0o644); err != nil {
		t.Fatalf("write sam: %v", err)
	}

	// SAM -> BAM with samtools so we read genuine upstream-encoded BAM.
	bamPath := filepath.Join(dir, "eqx.bam")
	cmd := exec.Command(samtools, "view", "-b", "-o", bamPath, samPath)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools view -b: %v\n%s", err, out)
	}

	// Read the BAM with our reader and write it to CRAM with our writer.
	header, records := readBAMRecords(t, bamPath)
	// Confirm the BAM really carries =/X ops, otherwise the test is vacuous.
	assertHasEqXOps(t, records)

	ourCRAM := filepath.Join(dir, "ours.cram")
	if err := writeCRAMFile(t, ourCRAM, header, records); err != nil {
		t.Fatalf("our CRAM write: %v", err)
	}

	// Decode our CRAM with samtools; compare to the input SAM fields.
	// =/X collapse to M in CRAM (both our writer and samtools' do this, as
	// htslib's encoder stores per-base features and the decoder rebuilds a
	// plain M run), so the input CIGAR is normalised =/X→M before
	// comparison; every other core field must match verbatim.
	gotUp := samtoolsViewRecords(t, samtools, ourCRAM)
	wantInput := normaliseEqXToM(readSAMBodyFields(t, samtools, samPath))
	assertSameCoreFields(t, "our-cram-vs-input", gotUp, wantInput)

	// Cross-check: samtools' own CRAM of the same alignment decodes the
	// same core fields (container bytes may differ; the SAM must not).
	smCRAM := filepath.Join(dir, "samtools.cram")
	cmd = exec.Command(samtools, "view", "-C", "-T", refPath, "-o", smCRAM, bamPath)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools view -C: %v\n%s", err, out)
	}
	gotSM := samtoolsViewRecords(t, samtools, smCRAM)
	assertSameCoreFields(t, "our-cram-vs-samtools-cram", gotUp, gotSM)

	// And our own reader decodes our CRAM identically to samtools.
	ourDecode := ourViewRecords(t, ourCRAM)
	assertSameCoreFields(t, "our-reader-vs-samtools", ourDecode, gotUp)
}

// readBAMRecords reads a BAM file with our BAM reader and returns its
// header and every record.
func readBAMRecords(t *testing.T, path string) (*sam.Header, []*sam.Record) {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open BAM %s: %v", path, err)
	}
	defer f.Close()
	br, err := sam.NewBAMReader(f)
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	var recs []*sam.Record
	for {
		rec, rerr := br.Read()
		if rec != nil {
			recs = append(recs, rec)
		}
		if rerr != nil {
			break
		}
	}
	return br.Header(), recs
}

// writeCRAMFile writes records to a CRAM file at path.
func writeCRAMFile(t *testing.T, path string, h *sam.Header, records []*sam.Record) error {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return WriteCRAM(f, h, records)
}

// assertHasEqXOps fails unless at least one record uses an = op and at
// least one uses an X op, so the parity assertion genuinely exercises the
// =/X write path.
func assertHasEqXOps(t *testing.T, records []*sam.Record) {
	t.Helper()
	sawEqual, sawMismatch := false, false
	for _, rec := range records {
		for _, op := range rec.Cigar {
			switch op.Op() {
			case sam.CigarEqual:
				sawEqual = true
			case sam.CigarMismatch:
				sawMismatch = true
			}
		}
	}
	if !sawEqual || !sawMismatch {
		t.Fatalf("BAM fixture lacks =/X ops (equal=%v mismatch=%v) — the =/X path is untested", sawEqual, sawMismatch)
	}
}

// readSAMBodyFields runs `samtools view` over a SAM/BAM file and returns
// the record lines (no header).
func readSAMBodyFields(t *testing.T, samtools, path string) []string {
	t.Helper()
	return samtoolsViewRecords(t, samtools, path)
}

// assertSameCoreFields compares the first eleven mandatory SAM fields
// (QNAME..QUAL) of two record-line sets. Auxiliary tags are ignored: a
// reference-free CRAM (ours) and a reference-based CRAM (samtools') differ
// only in regenerated MD/NM tags, never in the core alignment, so the
// core-field equality is the meaningful parity check.
func assertSameCoreFields(t *testing.T, label string, got, want []string) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s: got %d records, want %d", label, len(got), len(want))
	}
	for i := range want {
		g := coreFields(got[i])
		w := coreFields(want[i])
		if g != w {
			t.Fatalf("%s: record %d core-field mismatch:\n got=%q\nwant=%q", label, i, g, w)
		}
	}
}

// normaliseEqXToM rewrites the CIGAR column (column 6) of each record line
// so every run of =, X and M ops merges into a single M, matching how CRAM
// stores =/X (as per-base features rebuilt to M on decode). It is applied
// to the input SAM before comparing it to a CRAM-decoded line.
func normaliseEqXToM(lines []string) []string {
	out := make([]string, len(lines))
	for i, line := range lines {
		parts := strings.Split(line, "\t")
		if len(parts) > 5 {
			parts[5] = collapseEqXCigar(parts[5])
		}
		out[i] = strings.Join(parts, "\t")
	}
	return out
}

// collapseEqXCigar parses a CIGAR and merges adjacent =/X/M ops into M,
// leaving every other op unchanged. "*" passes through.
func collapseEqXCigar(s string) string {
	c, err := sam.ParseCigar(s)
	if err != nil || len(c) == 0 {
		return s
	}
	var merged sam.Cigar
	var run uint32
	flush := func() {
		if run > 0 {
			merged = append(merged, sam.CigarOp(run<<4|sam.CigarMatch))
			run = 0
		}
	}
	for _, op := range c {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			run += op.Length()
		default:
			flush()
			merged = append(merged, op)
		}
	}
	flush()
	return merged.String()
}

// coreFields returns the first eleven tab-separated SAM columns
// (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL),
// dropping the auxiliary tags.
func coreFields(line string) string {
	parts := strings.Split(line, "\t")
	if len(parts) > 11 {
		parts = parts[:11]
	}
	return strings.Join(parts, "\t")
}
