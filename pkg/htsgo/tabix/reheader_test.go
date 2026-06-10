package tabix

import (
	"bytes"
	"path/filepath"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// decodeBGZF decompresses bgzipped bytes into their plain text.
func decodeBGZF(t *testing.T, gz []byte) string {
	t.Helper()
	br, err := bgzip.NewReader(bytes.NewReader(gz))
	if err != nil {
		t.Fatalf("bgzf reader: %v", err)
	}
	var out bytes.Buffer
	if _, err := out.ReadFrom(br); err != nil {
		t.Fatalf("bgzf decode: %v", err)
	}
	return out.String()
}

func TestReheaderReplacesHeader(t *testing.T) {
	const body = "chr1\t100\t.\tA\tT\t.\t.\t.\nchr1\t200\t.\tC\tG\t.\t.\t.\n"
	in := writeBGZF(t, "in.vcf.gz",
		"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+body)
	const newHdr = "##fileformat=VCFv4.2\n##NEW=yes\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	hdr := writeFile(t, "hdr.txt", newHdr)

	var buf bytes.Buffer
	if err := Reheader(in, hdr, '#', &buf); err != nil {
		t.Fatalf("Reheader: %v", err)
	}
	got := decodeBGZF(t, buf.Bytes())
	want := newHdr + body
	if got != want {
		t.Errorf("reheader mismatch:\n got %q\nwant %q", got, want)
	}
}

func TestReheaderAppendsTrailingNewline(t *testing.T) {
	const body = "chr1\t100\t.\tA\tT\t.\t.\t.\n"
	in := writeBGZF(t, "in.vcf.gz", "#old\n"+body)
	// Header file with NO trailing newline must not merge into the first
	// data line.
	hdr := writeFile(t, "hdr.txt", "#new")

	var buf bytes.Buffer
	if err := Reheader(in, hdr, '#', &buf); err != nil {
		t.Fatalf("Reheader: %v", err)
	}
	got := decodeBGZF(t, buf.Bytes())
	want := "#new\n" + body
	if got != want {
		t.Errorf("trailing-newline handling: got %q want %q", got, want)
	}
}

func TestReheaderNoLeadingHeader(t *testing.T) {
	// File without any meta-char header: the whole content is body and is
	// preserved after the inserted header.
	const body = "chr1\t100\t.\tA\tT\t.\t.\t.\n"
	in := writeBGZF(t, "in.gz", body)
	hdr := writeFile(t, "hdr.txt", "#inserted\n")

	var buf bytes.Buffer
	if err := Reheader(in, hdr, '#', &buf); err != nil {
		t.Fatalf("Reheader: %v", err)
	}
	got := decodeBGZF(t, buf.Bytes())
	want := "#inserted\n" + body
	if got != want {
		t.Errorf("no-header case: got %q want %q", got, want)
	}
}

func TestReheaderEmptyHeaderFile(t *testing.T) {
	const body = "chr1\t100\t.\tA\tT\t.\t.\t.\n"
	in := writeBGZF(t, "in.gz", "#old\n"+body)
	hdr := writeFile(t, "hdr.txt", "")

	var buf bytes.Buffer
	if err := Reheader(in, hdr, '#', &buf); err != nil {
		t.Fatalf("Reheader: %v", err)
	}
	got := decodeBGZF(t, buf.Bytes())
	if got != body {
		t.Errorf("empty header: got %q want %q", got, body)
	}
}

func TestReheaderMissingHeaderFile(t *testing.T) {
	in := writeBGZF(t, "in.gz", "#old\nchr1\t1\t.\tA\tT\t.\t.\t.\n")
	var buf bytes.Buffer
	if err := Reheader(in, filepath.Join(t.TempDir(), "nope.txt"), '#', &buf); err == nil {
		t.Fatalf("expected error for missing header file")
	}
}

func TestReheaderMissingDataFile(t *testing.T) {
	hdr := writeFile(t, "hdr.txt", "#new\n")
	var buf bytes.Buffer
	if err := Reheader(filepath.Join(t.TempDir(), "nope.gz"), hdr, '#', &buf); err == nil {
		t.Fatalf("expected error for missing data file")
	}
}

func TestHeaderLen(t *testing.T) {
	cases := []struct {
		data string
		meta byte
		want int
	}{
		{"#a\n#b\nrow\n", '#', 6},
		{"row\n", '#', 0},
		{"#only", '#', 5}, // header line with no newline runs to EOF
		{"@h\nbody\n", '@', 3},
		{"", '#', 0},
	}
	for _, c := range cases {
		if got := headerLen([]byte(c.data), c.meta); got != c.want {
			t.Errorf("headerLen(%q)=%d want %d", c.data, got, c.want)
		}
	}
}
