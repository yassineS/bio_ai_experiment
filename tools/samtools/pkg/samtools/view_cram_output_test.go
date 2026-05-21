package samtools

import (
	"bytes"
	"os"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
)

// TestViewCRAMOutputRoundTrip is the C9 correctness oracle: it runs the
// `view` code path with CRAM output over the upstream CRAM fixture, then
// runs `view` again over that emitted CRAM, and asserts the two record
// streams agree. SAM -> CRAM -> SAM through the CLI must be lossless for
// the round-trip-significant fields (the documented unmapped-read
// MAPQ/CIGAR caveat applies, so the per-record key is QNAME/FLAG/
// RNAME/POS, which CRAM stores for every record).
func TestViewCRAMOutputRoundTrip(t *testing.T) {
	path := openCRAMFixture(t)
	want := cramFixtureRecordKeys(t, path)

	// First pass: view the fixture, emitting CRAM.
	in, err := os.Open(path)
	if err != nil {
		t.Fatalf("open CRAM fixture: %v", err)
	}
	defer in.Close()
	var cramOut bytes.Buffer
	n, err := View(in, &cramOut, ViewOptions{OutputCRAM: true})
	if err != nil {
		t.Fatalf("View -> CRAM: %v", err)
	}
	if n != len(want) {
		t.Fatalf("View -> CRAM emitted %d records, want %d", n, len(want))
	}

	// The bytes must be a real CRAM stream the CRAM reader accepts.
	if _, err := cram.NewRecordReader(bytes.NewReader(cramOut.Bytes())); err != nil {
		t.Fatalf("emitted output is not valid CRAM: %v", err)
	}

	// Second pass: view the emitted CRAM, emitting text SAM.
	var samOut bytes.Buffer
	n2, err := View(bytes.NewReader(cramOut.Bytes()), &samOut, ViewOptions{})
	if err != nil {
		t.Fatalf("View(CRAM) -> SAM: %v", err)
	}
	if n2 != len(want) {
		t.Fatalf("second View emitted %d records, want %d", n2, len(want))
	}

	got := parseRecordKeys(samOut.String())
	if len(got) != len(want) {
		t.Fatalf("re-parsed %d records, want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("record %d: got %q, want %q", i, got[i], want[i])
		}
	}
}

// TestViewCRAMOutputFromSAM round-trips a synthetic SAM stream through
// CRAM output and back, exercising the writer adapter without the
// reference fixture so the C9 path is covered even when the submodule is
// absent.
func TestViewCRAMOutputFromSAM(t *testing.T) {
	samText := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:100000\n" +
		"r1\t0\tchr1\t100\t30\t6M\t*\t0\t0\tACGTAC\tIIIIII\n" +
		"r2\t16\tchr1\t200\t40\t4M\t*\t0\t0\tTTGG\tIIII\n" +
		"r3\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n"

	var cramOut bytes.Buffer
	if _, err := View(bytes.NewReader([]byte(samText)), &cramOut, ViewOptions{OutputCRAM: true}); err != nil {
		t.Fatalf("View(SAM) -> CRAM: %v", err)
	}
	if _, err := cram.NewRecordReader(bytes.NewReader(cramOut.Bytes())); err != nil {
		t.Fatalf("emitted output is not valid CRAM: %v", err)
	}

	var samOut bytes.Buffer
	if _, err := View(bytes.NewReader(cramOut.Bytes()), &samOut, ViewOptions{}); err != nil {
		t.Fatalf("View(CRAM) -> SAM: %v", err)
	}
	got := parseRecordKeys(samOut.String())
	want := []string{
		"r1\t0\tchr1\t100",
		"r2\t16\tchr1\t200",
		"r3\t4\t*\t0",
	}
	if len(got) != len(want) {
		t.Fatalf("got %d records, want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("record %d: got %q, want %q", i, got[i], want[i])
		}
	}
}

// TestViewCRAMOutputQualityBinning exercises the `--output-fmt-option
// qbin=8` path: viewing a SAM stream to CRAM with the binning option set
// must yield CRAM whose decoded qualities are the 8-level-binned input,
// while the default (no option) leaves them verbatim.
func TestViewCRAMOutputQualityBinning(t *testing.T) {
	// r1 carries raw qualities chosen to straddle several Illumina-8 bins.
	// SAM QUAL char '!' is Phred 0; the offset is 33.
	samText := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:100000\n" +
		"r1\t0\tchr1\t100\t30\t6M\t*\t0\t0\tACGTAC\t#)-9CI\n"
	// '#'=2 ')'=8 '-'=12 '9'=24 'C'=34 'I'=40 (after subtracting 33).
	rawQual := []byte{2, 8, 12, 24, 34, 40}
	wantBinned := cram.BinningIllumina8.BinQuality(rawQual)

	decodeQual := func(t *testing.T, opts ViewOptions) []byte {
		t.Helper()
		var cramOut bytes.Buffer
		if _, err := View(bytes.NewReader([]byte(samText)), &cramOut, opts); err != nil {
			t.Fatalf("View -> CRAM: %v", err)
		}
		rr, err := cram.NewRecordReader(bytes.NewReader(cramOut.Bytes()))
		if err != nil {
			t.Fatalf("NewRecordReader: %v", err)
		}
		recs, err := rr.ReadAll()
		if err != nil {
			t.Fatalf("ReadAll: %v", err)
		}
		if len(recs) != 1 {
			t.Fatalf("decoded %d records, want 1", len(recs))
		}
		return recs[0].Qual
	}

	// With qbin=8 the decoded qualities must be the binned values.
	binned := decodeQual(t, ViewOptions{OutputCRAM: true, CRAMQualityBinning: "8"})
	if !bytes.Equal(binned, wantBinned) {
		t.Errorf("qbin=8 decoded QUAL = %v, want binned %v", binned, wantBinned)
	}

	// Without the option the qualities must round-trip verbatim.
	plain := decodeQual(t, ViewOptions{OutputCRAM: true})
	if !bytes.Equal(plain, rawQual) {
		t.Errorf("default decoded QUAL = %v, want verbatim %v", plain, rawQual)
	}
}

// TestViewCRAMOutputBadBinningOption confirms an unknown binning value is
// surfaced as an error rather than silently ignored.
func TestViewCRAMOutputBadBinningOption(t *testing.T) {
	samText := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n"
	var out bytes.Buffer
	_, err := View(bytes.NewReader([]byte(samText)), &out, ViewOptions{
		OutputCRAM: true, CRAMQualityBinning: "bogus",
	})
	if err == nil {
		t.Fatal("View accepted an unknown CRAMQualityBinning value")
	}
}

// TestIndexFileCRAMWritesCRAI runs the `index` file path over a CRAM
// file and asserts a .crai is written that cram.ReadCRAI parses back into
// sane entries: one per slice, offsets within the file, and an overlap
// query that returns the expected slice.
func TestIndexFileCRAMWritesCRAI(t *testing.T) {
	path := openCRAMFixture(t)
	cramBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read CRAM fixture: %v", err)
	}

	// Copy the fixture into the temp dir so the .crai lands beside it.
	dir := t.TempDir()
	cramCopy := dir + "/in.cram"
	if err := os.WriteFile(cramCopy, cramBytes, 0o644); err != nil {
		t.Fatalf("copy CRAM fixture: %v", err)
	}

	if err := IndexFile(cramCopy, "", IndexOptions{}); err != nil {
		t.Fatalf("IndexFile(CRAM): %v", err)
	}
	craiPath := cramCopy + ".crai"
	if _, err := os.Stat(craiPath); err != nil {
		t.Fatalf("expected .crai at %s: %v", craiPath, err)
	}

	idx, err := cram.OpenCRAI(craiPath)
	if err != nil {
		t.Fatalf("OpenCRAI: %v", err)
	}
	if len(idx.Entries) == 0 {
		t.Fatal(".crai has no entries for a non-empty CRAM")
	}
	fileLen := int64(len(cramBytes))
	for i, e := range idx.Entries {
		if e.ContainerOffset < 0 || e.ContainerOffset >= fileLen {
			t.Errorf("entry %d ContainerOffset %d out of bounds [0,%d)", i, e.ContainerOffset, fileLen)
		}
		if e.SliceOffset <= 0 {
			t.Errorf("entry %d SliceOffset %d, want > 0", i, e.SliceOffset)
		}
		if e.SliceSize <= 0 || e.ContainerOffset+e.SliceOffset+e.SliceSize > fileLen {
			t.Errorf("entry %d slice [%d+%d+%d] runs past file end %d",
				i, e.ContainerOffset, e.SliceOffset, e.SliceSize, fileLen)
		}
	}

	// An overlap query against the first mapped slice must return it.
	first := idx.Entries[0]
	if first.RefID >= 0 {
		hits := idx.Query(first.RefID,
			int64(first.AlignmentStart)-1,
			int64(first.AlignmentStart)+int64(first.AlignmentSpan))
		found := false
		for _, h := range hits {
			if h == first {
				found = true
			}
		}
		if !found {
			t.Errorf("overlap query on the first slice did not return it (got %d hits)", len(hits))
		}
	}
}

// TestIndexFileCRAMExplicitOutput confirms the -o output path override is
// honoured for CRAM input.
func TestIndexFileCRAMExplicitOutput(t *testing.T) {
	path := openCRAMFixture(t)
	cramBytes, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read CRAM fixture: %v", err)
	}
	dir := t.TempDir()
	cramCopy := dir + "/in.cram"
	if err := os.WriteFile(cramCopy, cramBytes, 0o644); err != nil {
		t.Fatalf("copy CRAM fixture: %v", err)
	}
	out := dir + "/custom.crai"
	if err := IndexFile(cramCopy, out, IndexOptions{}); err != nil {
		t.Fatalf("IndexFile(CRAM, -o): %v", err)
	}
	if _, err := cram.OpenCRAI(out); err != nil {
		t.Fatalf("OpenCRAI(%s): %v", out, err)
	}
}

// parseRecordKeys re-parses View's text-SAM output and returns the
// QNAME/FLAG/RNAME/POS key of each record line, the same key shape
// cramFixtureRecordKeys uses.
func parseRecordKeys(samText string) []string {
	var keys []string
	for _, line := range splitNonEmptyLines(samText) {
		fields := splitTabFields(line)
		if len(fields) < 4 {
			continue
		}
		keys = append(keys, fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3])
	}
	return keys
}

// splitNonEmptyLines splits s on newlines, dropping empty lines.
func splitNonEmptyLines(s string) []string {
	var out []string
	start := 0
	for i := 0; i <= len(s); i++ {
		if i == len(s) || s[i] == '\n' {
			if i > start {
				out = append(out, s[start:i])
			}
			start = i + 1
		}
	}
	return out
}

// splitTabFields splits a line into tab-separated fields.
func splitTabFields(s string) []string {
	var out []string
	start := 0
	for i := 0; i <= len(s); i++ {
		if i == len(s) || s[i] == '\t' {
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	return out
}
