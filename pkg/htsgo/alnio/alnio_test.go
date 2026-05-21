package alnio

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// samtoolsTestDir is the upstream samtools regression-test tree, vendored
// as a git submodule under reference_code/. Tests that depend on a fixture
// there t.Skip when the submodule is not initialised.
const samtoolsTestDir = "../../../reference_code/samtools/test"

// minimalSAM is a tiny well-formed SAM stream: a one-line header and one
// alignment record.
const minimalSAM = "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\n" +
	"r1\t0\tref\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n"

// TestNewReaderSAM routes a plain-text SAM stream through NewReader.
func TestNewReaderSAM(t *testing.T) {
	rd, err := NewReader(bytes.NewReader([]byte(minimalSAM)))
	if err != nil {
		t.Fatalf("NewReader(SAM): %v", err)
	}
	if got := len(rd.Header().Refs); got != 1 {
		t.Fatalf("header has %d refs, want 1", got)
	}
	rec, err := rd.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if rec.QName != "r1" {
		t.Errorf("QName = %q, want r1", rec.QName)
	}
	if _, err := rd.Read(); err != io.EOF {
		t.Errorf("second Read = %v, want io.EOF", err)
	}
}

// TestNewReaderGzipSAM confirms a plain-gzip-compressed SAM stream — which
// sam.NewReader cannot decode itself — is transparently decompressed.
func TestNewReaderGzipSAM(t *testing.T) {
	var buf bytes.Buffer
	gw := gzip.NewWriter(&buf)
	if _, err := gw.Write([]byte(minimalSAM)); err != nil {
		t.Fatalf("gzip write: %v", err)
	}
	if err := gw.Close(); err != nil {
		t.Fatalf("gzip close: %v", err)
	}
	rd, err := NewReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewReader(gzip SAM): %v", err)
	}
	rec, err := rd.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if rec.QName != "r1" {
		t.Errorf("QName = %q, want r1", rec.QName)
	}
}

// TestNewReaderCRAM routes a CRAM stream through NewReader and asserts the
// records match what cram.OpenRecords yields for the same file.
func TestNewReaderCRAM(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.cram")
	if _, err := os.Stat(path); err != nil {
		t.Skip("samtools submodule not initialised — fixture unavailable")
	}

	// Expected records via the CRAM reader directly.
	want, err := cram.OpenRecords(path)
	if err != nil {
		t.Fatalf("cram.OpenRecords: %v", err)
	}
	wantRecs, err := want.ReadAll()
	if err != nil {
		t.Fatalf("cram ReadAll: %v", err)
	}
	want.Close()

	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open fixture: %v", err)
	}
	defer f.Close()
	rd, err := NewReader(f)
	if err != nil {
		t.Fatalf("NewReader(CRAM): %v", err)
	}
	gotRecs, err := readAll(rd)
	if err != nil {
		t.Fatalf("alnio ReadAll: %v", err)
	}
	if len(gotRecs) != len(wantRecs) {
		t.Fatalf("alnio decoded %d records, cram decoded %d", len(gotRecs), len(wantRecs))
	}
	for i := range gotRecs {
		if got, w := recordText(t, gotRecs[i]), recordText(t, wantRecs[i]); got != w {
			t.Errorf("record %d mismatch:\n got: %s\nwant: %s", i, got, w)
		}
	}
}

// recordText renders one record as a SAM text line for comparison.
func recordText(t *testing.T, rec *sam.Record) string {
	t.Helper()
	var buf bytes.Buffer
	sw := sam.NewSAMWriter(&buf)
	if err := sw.Write(rec); err != nil {
		t.Fatalf("SAMWriter.Write: %v", err)
	}
	if err := sw.Close(); err != nil {
		t.Fatalf("SAMWriter.Close: %v", err)
	}
	return buf.String()
}

// TestOpenReaderCRAM exercises the path-based, closeable OpenReader entry
// point on the CRAM fixture.
func TestOpenReaderCRAM(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.cram")
	if _, err := os.Stat(path); err != nil {
		t.Skip("samtools submodule not initialised — fixture unavailable")
	}
	rd, err := OpenReader(path, "")
	if err != nil {
		t.Fatalf("OpenReader(CRAM): %v", err)
	}
	defer rd.Close()
	recs, err := readAll(rd)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) != 15 {
		t.Errorf("decoded %d records, want 15", len(recs))
	}
}

// TestOpenReaderSAM exercises OpenReader on a text SAM file.
func TestOpenReaderSAM(t *testing.T) {
	dir := t.TempDir()
	p := filepath.Join(dir, "x.sam")
	if err := os.WriteFile(p, []byte(minimalSAM), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	rd, err := OpenReader(p, "")
	if err != nil {
		t.Fatalf("OpenReader(SAM): %v", err)
	}
	defer rd.Close()
	recs, err := readAll(rd)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) != 1 {
		t.Errorf("decoded %d records, want 1", len(recs))
	}
}

// TestOpenReaderMissing reports a clear error for a non-existent path.
func TestOpenReaderMissing(t *testing.T) {
	if _, err := OpenReader(filepath.Join(t.TempDir(), "nope.bam"), ""); err == nil {
		t.Error("OpenReader on a missing file should error")
	}
}

// readAll drains a sam.Reader into a slice of records.
func readAll(rd sam.Reader) ([]*sam.Record, error) {
	var out []*sam.Record
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return out, err
		}
		out = append(out, rec)
	}
}
