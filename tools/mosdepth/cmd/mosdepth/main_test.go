package main

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
)

// makeBAM is a slim duplicate of the same helper inside the pkg/mosdepth
// tests; kept local so this test package doesn't depend on test-only
// exports.
func makeBAM(t *testing.T, refs []sam.Reference, recs []*sam.Record) []byte {
	t.Helper()
	hdr := &sam.Header{Refs: refs}
	for _, r := range refs {
		hdr.Lines = append(hdr.Lines, sam.HeaderLine{
			Tag: "SQ",
			Fields: []sam.HeaderField{
				{Tag: "SN", Value: r.Name},
				{Tag: "LN", Value: itoa(int(r.Length))},
			},
		})
	}
	var buf bytes.Buffer
	bw := sam.NewBAMWriter(&buf)
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, rec := range recs {
		if err := bw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes()
}

func itoa(v int) string {
	if v == 0 {
		return "0"
	}
	neg := v < 0
	if neg {
		v = -v
	}
	d := []byte{}
	for v > 0 {
		d = append([]byte{byte('0' + v%10)}, d...)
		v /= 10
	}
	if neg {
		return "-" + string(d)
	}
	return string(d)
}

func TestRun_HelpAndVersion(t *testing.T) {
	if rc := run([]string{"-h"}); rc != 0 {
		t.Errorf("help rc: %d", rc)
	}
	if rc := run([]string{"--version"}); rc != 0 {
		t.Errorf("version rc: %d", rc)
	}
}

func TestRun_BadPositional(t *testing.T) {
	if rc := run([]string{}); rc != 2 {
		t.Errorf("expected rc=2 for missing args, got %d", rc)
	}
}

func TestRun_BadFlag(t *testing.T) {
	if rc := run([]string{"--bogus-flag"}); rc != 2 {
		t.Errorf("expected rc=2 for bad flag, got %d", rc)
	}
}

func TestRun_Successful(t *testing.T) {
	dir := t.TempDir()
	cigar, err := sam.ParseCigar("5M")
	if err != nil {
		t.Fatal(err)
	}
	bam := makeBAM(t, []sam.Reference{{Name: "chr1", Length: 20}}, []*sam.Record{
		{QName: "r1", RName: "chr1", Pos: 1, Cigar: cigar, MapQ: 60, Seq: "AAAAA"},
	})
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bam, 0644); err != nil {
		t.Fatal(err)
	}
	prefix := filepath.Join(dir, "out")
	rc := run([]string{prefix, bamPath})
	if rc != 0 {
		t.Errorf("run rc: %d", rc)
	}
	// per-base file produced.
	if _, err := os.Stat(prefix + ".per-base.bed.gz"); err != nil {
		t.Errorf("per-base missing: %v", err)
	}
}

func TestRun_D4Rejected(t *testing.T) {
	dir := t.TempDir()
	cigar, _ := sam.ParseCigar("5M")
	bam := makeBAM(t, []sam.Reference{{Name: "chr1", Length: 20}}, []*sam.Record{
		{QName: "r1", RName: "chr1", Pos: 1, Cigar: cigar, MapQ: 60, Seq: "AAAAA"},
	})
	bamPath := filepath.Join(dir, "in.bam")
	os.WriteFile(bamPath, bam, 0644)
	prefix := filepath.Join(dir, "out")
	if rc := run([]string{"-d", prefix, bamPath}); rc != 1 {
		t.Errorf("D4 rc: got %d, want 1", rc)
	}
}

func TestRun_BadThresholds(t *testing.T) {
	if rc := run([]string{"-T", "abc", "/tmp/x", "/tmp/x.bam"}); rc != 1 {
		t.Errorf("bad thresholds rc: %d", rc)
	}
}

func TestSplitCSV(t *testing.T) {
	got := splitCSV(",a,,b,c,")
	want := []string{"a", "b", "c"}
	if strings.Join(got, ",") != strings.Join(want, ",") {
		t.Errorf("splitCSV: got %v, want %v", got, want)
	}
}

func TestUint8Clamp(t *testing.T) {
	if uint8Clamp(-1) != 0 || uint8Clamp(255) != 255 || uint8Clamp(300) != 255 || uint8Clamp(42) != 42 {
		t.Errorf("uint8Clamp behaves unexpectedly")
	}
}
