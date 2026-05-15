package bedlinks

import (
	"bytes"
	"errors"
	"strings"
	"testing"
)

func TestRun_BED3_Defaults(t *testing.T) {
	in := "chr1\t100\t200\nchr2\t1000\t1100\n"
	var buf bytes.Buffer
	n, err := Run(strings.NewReader(in), &buf, Options{})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 2 {
		t.Errorf("rows = %d, want 2", n)
	}
	out := buf.String()

	// Header invariants.
	if !strings.Contains(out, "<html>") || !strings.Contains(out, "</html>") {
		t.Errorf("missing html scaffolding: %s", out)
	}
	if !strings.Contains(out, "<title>stdin</title>") {
		t.Errorf("expected <title>stdin</title>, got: %s", out)
	}
	if !strings.Contains(out, "BED Entries from: stdin") {
		t.Errorf("missing 'BED Entries from: stdin' header: %s", out)
	}

	// URL uses 0-based start; visible text uses 1-based.
	wantURL := "<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg18&position=chr1:100-200>chr1:101-200</a>"
	if !strings.Contains(out, wantURL) {
		t.Errorf("missing chr1 row link.\nwant substring: %s\ngot:\n%s", wantURL, out)
	}
	wantURL2 := "<a href=http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db=hg18&position=chr2:1000-1100>chr2:1001-1100</a>"
	if !strings.Contains(out, wantURL2) {
		t.Errorf("missing chr2 row link: %s", out)
	}

	// BED3 rows should have NO name/score/strand <td> blocks.
	if strings.Contains(out, "<td>\n\n\t</td>") {
		t.Errorf("BED3 row should not emit extra <td>: %s", out)
	}
}

func TestRun_BED6_EmitsNameScoreStrand(t *testing.T) {
	in := "chr1\t100\t200\tgeneA\t500\t+\nchr1\t300\t400\tgeneB\t250\t-\n"
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(in), &buf, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	out := buf.String()
	// Each BED6 row emits four <td> blocks (link + name + score + strand).
	for _, want := range []string{"\t<td>\ngeneA\n\t</td>\n", "\t<td>\n500\n\t</td>\n", "\t<td>\n+\n\t</td>\n",
		"\t<td>\ngeneB\n\t</td>\n", "\t<td>\n250\n\t</td>\n", "\t<td>\n-\n\t</td>\n"} {
		if !strings.Contains(out, want) {
			t.Errorf("missing %q in:\n%s", want, out)
		}
	}
}

func TestRun_BED5_EmitsNameAndScoreOnly(t *testing.T) {
	in := "chr1\t100\t200\tgeneA\t500\n"
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(in), &buf, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	out := buf.String()
	if !strings.Contains(out, "\t<td>\ngeneA\n\t</td>\n") {
		t.Errorf("missing name <td>: %s", out)
	}
	if !strings.Contains(out, "\t<td>\n500\n\t</td>\n") {
		t.Errorf("missing score <td>: %s", out)
	}
	// No strand <td> on BED5.
	if strings.Contains(out, "<td>\n+\n") || strings.Contains(out, "<td>\n-\n") {
		t.Errorf("BED5 row should not emit strand <td>: %s", out)
	}
}

func TestRun_BED4_EmitsNameOnly(t *testing.T) {
	in := "chr1\t100\t200\tgeneA\n"
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(in), &buf, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	out := buf.String()
	if !strings.Contains(out, "\t<td>\ngeneA\n\t</td>\n") {
		t.Errorf("missing name <td>: %s", out)
	}
}

func TestRun_CustomBaseOrgDB(t *testing.T) {
	in := "chrM\t0\t100\n"
	var buf bytes.Buffer
	opts := Options{
		Base: "http://mymirror.example.edu",
		Org:  "mouse",
		DB:   "mm9",
	}
	if _, err := Run(strings.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	out := buf.String()
	want := "<a href=http://mymirror.example.edu/cgi-bin/hgTracks?org=mouse&db=mm9&position=chrM:0-100>chrM:1-100</a>"
	if !strings.Contains(out, want) {
		t.Errorf("missing custom URL.\nwant: %s\ngot:\n%s", want, out)
	}
}

func TestRun_CustomBedFileInTitle(t *testing.T) {
	in := "chr1\t1\t2\n"
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(in), &buf, Options{BedFile: "regions.bed"}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !strings.Contains(buf.String(), "<title>regions.bed</title>") {
		t.Errorf("expected custom title; got:\n%s", buf.String())
	}
}

func TestRun_EmptyInput_StillEmitsScaffold(t *testing.T) {
	var buf bytes.Buffer
	n, err := Run(strings.NewReader(""), &buf, Options{})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 0 {
		t.Errorf("rows = %d, want 0", n)
	}
	out := buf.String()
	for _, want := range []string{"<html>", "</html>", "<table border=", "</table>"} {
		if !strings.Contains(out, want) {
			t.Errorf("missing %q in empty-input output:\n%s", want, out)
		}
	}
	// No <tr> data rows.
	if strings.Contains(out, "<tr>") {
		t.Errorf("empty input should not emit any <tr>: %s", out)
	}
}

func TestRun_PropagatesBedParseError(t *testing.T) {
	in := "chr1\tnotanumber\t200\n"
	var buf bytes.Buffer
	_, err := Run(strings.NewReader(in), &buf, Options{})
	if err == nil {
		t.Fatalf("expected parse error from bed.Reader")
	}
}

func TestInferBedType_BED3to12(t *testing.T) {
	// Use minimal records that walk every branch of inferBedType /
	// writeRow's switch (3/4/5/6/9/12). bed.Reader populates optional
	// fields based on column count, so the columns below are crafted to
	// stop at each bedType boundary.
	//
	// BED7/8: same <td> layout as BED6 in upstream (the switch falls
	// through to default = link-only). BED11 is a degenerate case where
	// BlockCount is set but BlockStarts is not; same handling. We
	// exercise BED9 and BED12 explicitly to cover the case-6-9-12 arm.
	cases := []struct {
		desc      string
		in        string
		wantTDs   int // total <td> blocks (link + extras)
		wantWords []string
	}{
		{"BED3", "chr1\t0\t1\n", 1, nil},
		{"BED4", "chr1\t0\t1\tg\n", 2, []string{"\t<td>\ng\n\t</td>\n"}},
		{"BED5", "chr1\t0\t1\tg\t10\n", 3, []string{"\t<td>\ng\n\t</td>\n", "\t<td>\n10\n\t</td>\n"}},
		{"BED6", "chr1\t0\t1\tg\t10\t+\n", 4, []string{"\t<td>\n+\n\t</td>\n"}},
		// BED9: 9 columns — Name/Score/Strand + ThickStart/ThickEnd/RGB.
		{"BED9", "chr1\t0\t1\tg\t10\t+\t0\t1\t255,0,0\n", 4, []string{"\t<td>\n+\n\t</td>\n"}},
		// BED12: 12 columns — also exercises the BlockStarts branch in inferBedType.
		{"BED12", "chr1\t0\t1\tg\t10\t+\t0\t1\t255,0,0\t1\t1\t0\n", 4, []string{"\t<td>\n+\n\t</td>\n"}},
	}
	for _, tc := range cases {
		var buf bytes.Buffer
		if _, err := Run(strings.NewReader(tc.in), &buf, Options{}); err != nil {
			t.Fatalf("%s: Run: %v", tc.desc, err)
		}
		got := strings.Count(buf.String(), "<td>")
		if got != tc.wantTDs {
			t.Errorf("%s: <td> count = %d, want %d (output:\n%s)", tc.desc, got, tc.wantTDs, buf.String())
		}
		for _, w := range tc.wantWords {
			if !strings.Contains(buf.String(), w) {
				t.Errorf("%s: missing %q in:\n%s", tc.desc, w, buf.String())
			}
		}
	}
}

// errWriter fails on every Write.
type errWriter struct{}

func (errWriter) Write(_ []byte) (int, error) { return 0, errors.New("forced") }

func TestRun_WriteErrorPropagates(t *testing.T) {
	// bufio.Writer buffers; flooding it with many records guarantees a
	// flush, which surfaces the underlying writer error.
	in := strings.Repeat("chr1\t100\t200\n", 5000)
	_, err := Run(strings.NewReader(in), errWriter{}, Options{})
	if err == nil {
		t.Fatalf("expected write error to propagate")
	}
}
