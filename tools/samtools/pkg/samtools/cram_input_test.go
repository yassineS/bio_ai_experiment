package samtools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// cramFixturePath is the upstream samtools reference-free v3.0 CRAM test
// file, vendored byte-for-byte into tools/samtools/testdata/parity/ (a copy
// of reference_code/samtools/test/dat/test_input_1_a.cram). Vendoring the
// fixture bytes keeps the parity suite self-contained — it does not depend on
// the samtools submodule being initialised. Its decoded SAM is the sibling
// test_input_1_a.sam.
var cramFixturePath = filepath.Join("..", "..", "testdata", "parity", "test_input_1_a.cram")

// openCRAMFixture returns the path to the vendored CRAM fixture. The fixture
// is committed under testdata, so a missing file is a hard failure rather
// than a skip.
func openCRAMFixture(t *testing.T) string {
	t.Helper()
	if _, err := os.Stat(cramFixturePath); err != nil {
		t.Fatalf("vendored CRAM fixture missing at %s: %v", cramFixturePath, err)
	}
	return cramFixturePath
}

// cramFixtureRecordKeys returns a stable per-record key (QNAME + FLAG +
// RNAME + POS) for every record cram.OpenRecords yields from the fixture.
// It is the oracle the samtools CRAM-input tests assert against.
func cramFixtureRecordKeys(t *testing.T, path string) []string {
	t.Helper()
	rr, err := cram.OpenRecords(path)
	if err != nil {
		t.Fatalf("cram.OpenRecords: %v", err)
	}
	defer rr.Close()
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("cram ReadAll: %v", err)
	}
	keys := make([]string, len(recs))
	for i, r := range recs {
		keys[i] = recordKey(r)
	}
	return keys
}

// recordKey is a stable identity string for a record: QNAME, FLAG, RNAME,
// POS joined by tabs. RNAME and POS are formatted the way SAM text does
// (an empty RNAME is "*"; POS is the 1-based coordinate) so the key
// matches View's re-parsed text output.
func recordKey(r *sam.Record) string {
	rname := r.RName
	if rname == "" {
		rname = "*"
	}
	return strings.Join([]string{
		r.QName,
		strconv.Itoa(int(r.Flag)),
		rname,
		strconv.Itoa(int(r.Pos)),
	}, "\t")
}

// TestViewAcceptsCRAM confirms the samtools view code path transparently
// decodes a CRAM input and emits exactly the records the CRAM reader
// yields for the same file.
func TestViewAcceptsCRAM(t *testing.T) {
	path := openCRAMFixture(t)
	want := cramFixtureRecordKeys(t, path)

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open CRAM fixture: %v", err)
	}
	defer f.Close()

	var out bytes.Buffer
	n, err := View(f, &out, ViewOptions{})
	if err != nil {
		t.Fatalf("View on CRAM: %v", err)
	}
	if n != len(want) {
		t.Fatalf("View emitted %d records, want %d", n, len(want))
	}

	// Re-parse the emitted SAM text and key each record. View emits text
	// SAM with no header, so each non-blank line is one record.
	var got []string
	for _, line := range strings.Split(strings.TrimRight(out.String(), "\n"), "\n") {
		if line == "" {
			continue
		}
		fields := strings.Split(line, "\t")
		// QNAME, FLAG, RNAME, POS.
		got = append(got, fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3])
	}
	if len(got) != len(want) {
		t.Fatalf("re-parsed %d records, want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("record %d: got %q, want %q", i, got[i], want[i])
		}
	}
}

// TestParity_View_CRAM_VsUpstream is the byte-for-byte CRAM-input parity
// test: it runs the Go port's CRAM-decoding View path and the live upstream
// `samtools view` on the same vendored reference-free CRAM file and asserts
// the emitted record bodies are identical. The fixture decodes without an
// external reference, so no `-T` reference is needed; the record stream
// (including the full aux column, which the Go CRAM decoder emits in the same
// order upstream does) must match exactly. Header @PG injection differs
// between the two (the Go port deliberately omits @PG), so the comparison is
// over the record bodies that `view` without `-h` produces.
func TestParity_View_CRAM_VsUpstream(t *testing.T) {
	path := openCRAMFixture(t)

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open CRAM fixture: %v", err)
	}
	defer f.Close()

	var out bytes.Buffer
	if _, err := View(f, &out, ViewOptions{}); err != nil {
		t.Fatalf("View on CRAM: %v", err)
	}

	want := upstreamSamtoolsRun(t, "view", path)
	if !bytes.Equal(out.Bytes(), want) {
		t.Errorf("CRAM view mismatch vs upstream samtools view.\nwant:\n%s\ngot:\n%s", want, out.String())
	}
}

// TestViewCountCRAM confirms `view -c` counts CRAM records.
func TestViewCountCRAM(t *testing.T) {
	path := openCRAMFixture(t)
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open CRAM fixture: %v", err)
	}
	defer f.Close()

	var out bytes.Buffer
	n, err := View(f, &out, ViewOptions{Count: true})
	if err != nil {
		t.Fatalf("View -c on CRAM: %v", err)
	}
	if n != 15 {
		t.Errorf("View -c counted %d records, want 15", n)
	}
}

// TestViewHeaderOnlyCRAM confirms `view -H` emits the CRAM embedded header.
func TestViewHeaderOnlyCRAM(t *testing.T) {
	path := openCRAMFixture(t)
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open CRAM fixture: %v", err)
	}
	defer f.Close()

	var out bytes.Buffer
	if _, err := View(f, &out, ViewOptions{HeaderOnly: true}); err != nil {
		t.Fatalf("View -H on CRAM: %v", err)
	}
	if !strings.Contains(out.String(), "@SQ\tSN:ref1") {
		t.Errorf("view -H output missing @SQ ref1 line:\n%s", out.String())
	}
}

// TestDepthAcceptsCRAM confirms samtools depth accepts a CRAM input
// without error and produces depth rows.
func TestDepthAcceptsCRAM(t *testing.T) {
	path := openCRAMFixture(t)
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open CRAM fixture: %v", err)
	}
	defer f.Close()

	var out bytes.Buffer
	if err := Depth([]io.Reader{f}, &out, DepthOptions{}); err != nil {
		t.Fatalf("Depth on CRAM: %v", err)
	}
	if out.Len() == 0 {
		t.Error("Depth on CRAM produced no output rows")
	}
}

// TestFastqAcceptsCRAM confirms samtools fastq accepts a CRAM input
// without error and emits FASTQ records.
func TestFastqAcceptsCRAM(t *testing.T) {
	path := openCRAMFixture(t)
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open CRAM fixture: %v", err)
	}
	defer f.Close()

	dir := t.TempDir()
	outPath := filepath.Join(dir, "out.fastq")
	counts, err := Fastq(f, FastqOptions{OutputPath: outPath})
	if err != nil {
		t.Fatalf("Fastq on CRAM: %v", err)
	}
	if counts.Output == 0 {
		t.Error("Fastq on CRAM emitted no records")
	}
	data, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("read fastq output: %v", err)
	}
	if !bytes.HasPrefix(data, []byte("@")) {
		t.Errorf("fastq output does not start with '@':\n%s", data)
	}
}

// TestMpileupAcceptsCRAM confirms samtools mpileup accepts a CRAM input
// without error.
func TestMpileupAcceptsCRAM(t *testing.T) {
	path := openCRAMFixture(t)
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open CRAM fixture: %v", err)
	}
	defer f.Close()

	var out bytes.Buffer
	if err := Mpileup([]io.Reader{f}, &out, MpileupOptions{}, nil, nil); err != nil {
		t.Fatalf("Mpileup on CRAM: %v", err)
	}
	// Mapped reads in the fixture cover several positions; pileup must
	// emit at least one row.
	if out.Len() == 0 {
		t.Error("Mpileup on CRAM produced no output")
	}
}
