package hfile

import (
	"bytes"
	"io"
	"os"
	"testing"
)

// TestNetworkBackends exercises the HTTP(S), S3 and GCS backends against real,
// publicly-readable objects. It is gated on HFILE_NETTEST=1 so it never runs
// in CI (which has no outbound network and must stay hermetic); run it
// manually with `HFILE_NETTEST=1 go test ./pkg/htsgo/hfile/ -run Network -v`.
func TestNetworkBackends(t *testing.T) {
	if os.Getenv("HFILE_NETTEST") == "" {
		t.Skip("set HFILE_NETTEST=1 to run live-network backend tests")
	}

	cases := []struct {
		name string
		url  string
	}{
		{"https", "https://raw.githubusercontent.com/samtools/htslib/develop/README.md"},
		{"s3", "s3://1000genomes/CHANGELOG"},
		{"gcs", "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			h, err := Open(c.url)
			if err != nil {
				t.Fatalf("Open(%s): %v", c.url, err)
			}
			defer h.Close()

			size, err := h.Size()
			if err != nil {
				t.Fatalf("Size(%s): %v", c.url, err)
			}
			if size <= 0 {
				t.Fatalf("Size(%s) = %d, want > 0", c.url, size)
			}
			t.Logf("%s: Size = %d bytes", c.name, size)

			// Ranged ReadAt: read a window from the middle.
			off := size / 2
			n := 64
			if int64(n) > size-off {
				n = int(size - off)
			}
			buf := make([]byte, n)
			got, err := h.ReadAt(buf, off)
			if err != nil && err != io.EOF {
				t.Fatalf("ReadAt(%s, off=%d): %v", c.url, off, err)
			}
			if got == 0 {
				t.Fatalf("ReadAt(%s) read 0 bytes", c.url)
			}
			t.Logf("%s: ReadAt(off=%d) -> %d bytes: %q", c.name, off, got, sample(buf[:got]))

			// Sequential Read of the first chunk via the Reader interface.
			head := make([]byte, 32)
			hn, err := io.ReadFull(h, head)
			if err != nil && err != io.EOF && err != io.ErrUnexpectedEOF {
				t.Fatalf("Read head(%s): %v", c.url, err)
			}
			t.Logf("%s: head[:%d] = %q", c.name, hn, sample(head[:hn]))
		})
	}
}

// sample trims a byte slice to a short, printable preview.
func sample(b []byte) string {
	const max = 48
	if len(b) > max {
		b = b[:max]
	}
	return string(bytes.Map(func(r rune) rune {
		if r < 32 || r > 126 {
			return '.'
		}
		return r
	}, b))
}
