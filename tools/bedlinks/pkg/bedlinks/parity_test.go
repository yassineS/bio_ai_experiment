package bedlinks

// Parity tests for `bedtools links`.
//
// Upstream ships no `links/` test subdir under `reference_code/bedtools/test/`,
// so expected outputs are derived directly from the upstream source
// (`reference_code/bedtools/src/linksBed/linksBed.cpp` + `linksMain.cpp`) by
// mechanically playing forward CreateLinks() and WriteURL() on the fixtures
// below.
//
// Fixtures live in testdata/parity/. Cases:
//
//   - t1: BED6 fixture, defaults (base/org/db).
//   - t2: BED6 fixture, custom -base/-org/-db (matches the upstream help
//         example: -base http://mymirror.example.edu -org mouse -db mm9).
//   - t3: BED3 fixture, defaults — exercises the bedType=3 branch which
//         emits no name/score/strand <td> blocks.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runParity(t *testing.T, inputFile string, opts Options) []byte {
	t.Helper()
	in := readParity(t, inputFile)
	var buf bytes.Buffer
	if _, err := Run(bytes.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return buf.Bytes()
}

// t1: BED6, defaults.
func TestParity_Links_T1_BED6_Defaults(t *testing.T) {
	got := runParity(t, "links.bed", Options{})
	want := readParity(t, "t1.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// t2: BED6, custom mirror (`-base ... -org mouse -db mm9`).
func TestParity_Links_T2_CustomMirror(t *testing.T) {
	got := runParity(t, "links.bed", Options{
		Base: "http://mymirror.example.edu",
		Org:  "mouse",
		DB:   "mm9",
	})
	want := readParity(t, "t2.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// t3: BED3, defaults — no name/score/strand <td> blocks per row.
func TestParity_Links_T3_BED3_Defaults(t *testing.T) {
	got := runParity(t, "links_bed3.bed", Options{})
	want := readParity(t, "t3.expected")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
