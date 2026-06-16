package cram

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// lossyNamesRefFASTA is a two-reference FASTA, enough to force a
// cross-reference pair (which CRAM stores DETACHED) alongside same-ref
// pairs (stored with a downstream-mate link, whose duplicate name a
// lossy-names CRAM drops).
const lossyNamesRefFASTA = ">chr1\n" +
	"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
	">chr2\n" +
	"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n"

// lossyNamesSAM has three kinds of record whose names a lossy-names CRAM
// treats differently:
//   - readA: a proper same-reference pair. Both names are reconstructable
//     from the mate link, so the encoder DROPS the name; the decoder
//     synthesises "<prefix>:<n>".
//   - readB: a cross-reference pair (mates on chr1 and chr2), stored
//     DETACHED, so the real name is kept — in the mate block, after MF and
//     before NS, the layout this change adds support for reading.
//   - readC: an unpaired singleton, also DETACHED, real name kept.
const lossyNamesSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:64\n" +
	"@SQ\tSN:chr2\tLN:64\n" +
	"readA\t99\tchr1\t1\t60\t8M\t=\t30\t37\tACGTACGT\tIIIIIIII\n" +
	"readB\t97\tchr1\t5\t60\t8M\tchr2\t10\t0\tACGTACGT\tIIIIIIII\n" +
	"readA\t147\tchr1\t30\t60\t8M\t=\t1\t-37\tACGTACGT\tIIIIIIII\n" +
	"readB\t145\tchr2\t10\t60\t8M\tchr1\t5\t0\tTTTTGGGG\tIIIIIIII\n" +
	"readC\t0\tchr2\t20\t60\t8M\t*\t0\t0\tGGGGCCCC\tIIIIIIII\n"

// TestLossyNamesReadParity proves our reader decodes a read-names-NOT-
// preserved (lossy_names) CRAM — a layout that previously errored with
// "detached mate with read names not preserved is not yet supported" —
// byte-for-byte against `samtools view` of the same CRAM. It checks both
// that the real names of DETACHED records are read from the mate block and
// that the synthesised "<prefix>:<n>" names of the dropped non-detached
// pair match samtools' own synthesis.
func TestLossyNamesReadParity(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	dir := t.TempDir()

	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(lossyNamesRefFASTA), 0o644); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	if out, err := exec.Command(samtools, "faidx", refPath).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx: %v\n%s", err, out)
	}

	samPath := filepath.Join(dir, "pairs.sam")
	if err := os.WriteFile(samPath, []byte(lossyNamesSAM), 0o644); err != nil {
		t.Fatalf("write sam: %v", err)
	}

	// Encode to CRAM with read names dropped (lossy_names=1).
	cramPath := filepath.Join(dir, "lossy.cram")
	cmd := exec.Command(samtools, "view", "-C", "-T", refPath,
		"--output-fmt-option", "lossy_names=1", "-o", cramPath, samPath)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools view -C lossy_names=1: %v\n%s", err, out)
	}

	// The fixture must genuinely be read-names-not-preserved and contain a
	// detached record, otherwise the test would pass vacuously.
	assertLossyNamesDetached(t, cramPath)

	// Decode with our reader (reference attached so MD/NM regenerate just
	// as `samtools view -T` does) and with samtools; the full record lines
	// — synthesised names, kept detached names, SEQ and tags — must match.
	got := ourViewRecordsWithRef(t, cramPath, refPath)
	want := samtoolsViewRecords(t, samtools, cramPath)
	if len(got) != len(want) {
		t.Fatalf("decoded %d records, samtools decoded %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("record %d mismatch:\n got=%q\nwant=%q", i, got[i], want[i])
		}
	}

	// Spot-check that at least one synthesised "<prefix>:" name and the kept
	// detached "readB"/"readC" names actually appear, proving both paths ran.
	joined := strings.Join(got, "\n")
	if !strings.Contains(joined, "lossy.cram:") {
		t.Fatalf("no synthesised <prefix>:<n> name decoded; got:\n%s", joined)
	}
	if !strings.Contains(joined, "readB") || !strings.Contains(joined, "readC") {
		t.Fatalf("kept detached names readB/readC missing; got:\n%s", joined)
	}
}

// ourViewRecordsWithRef decodes a reference-backed CRAM with our reader,
// attaching the given FASTA so reference-derived bases and MD/NM tags are
// reconstructed exactly as `samtools view -T ref` produces them. It
// returns the record body lines (no header).
func ourViewRecordsWithRef(t *testing.T, cramPath, refPath string) []string {
	t.Helper()
	rr, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords %s: %v", cramPath, err)
	}
	defer rr.Close()
	if err := rr.SetReferenceFASTA(refPath); err != nil {
		t.Fatalf("SetReferenceFASTA %s: %v", refPath, err)
	}
	var buf bytes.Buffer
	if err := rr.WriteSAM(&buf); err != nil {
		t.Fatalf("WriteSAM %s: %v", cramPath, err)
	}
	var recs []string
	for _, line := range splitNonEmptyLines(buf.String()) {
		if strings.HasPrefix(line, "@") {
			continue
		}
		recs = append(recs, line)
	}
	return recs
}

// assertLossyNamesDetached opens cramPath and confirms its compression
// header marks read names as NOT included — the precondition for the
// mate-block name path this change adds. (That a detached record exists in
// this exact fixture is what made the old code error; the parity loop in
// the caller would fail loudly were the detached names not decoded.)
func assertLossyNamesDetached(t *testing.T, cramPath string) {
	t.Helper()
	rd, err := Open(cramPath)
	if err != nil {
		t.Fatalf("Open %s: %v", cramPath, err)
	}
	defer rd.Close()
	conts, err := rd.Containers()
	if err != nil {
		t.Fatalf("Containers: %v", err)
	}
	sawNotPreserved := false
	for _, c := range conts {
		dc, derr := ParseDataContainer(c)
		if derr != nil {
			continue
		}
		if dc.Compression != nil && !dc.Compression.Preservation.ReadNamesIncluded {
			sawNotPreserved = true
		}
	}
	if !sawNotPreserved {
		t.Fatalf("fixture preserves read names — the lossy-names path is untested")
	}
}
