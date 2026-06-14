package prinseq

import (
	"bytes"
	"image"
	"image/color"
	"image/png"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// gdFixturePath points at the real prinseq-lite .gd file vendored in
// testdata (produced by the upstream example1.fastq run). It exercises
// every sub-table the renderer knows about.
const gdFixturePath = "../../testdata/parity/graphdata_example1.gd"

func loadFixture(t *testing.T) *GDData {
	t.Helper()
	f, err := os.Open(gdFixturePath)
	if err != nil {
		t.Fatalf("opening fixture: %v", err)
	}
	defer f.Close()
	d, err := ParseGD(f)
	if err != nil {
		t.Fatalf("parsing fixture: %v", err)
	}
	return d
}

// TestParseGDExtractsSeries asserts the data-extraction layer pulls
// the exact series out of the .gd. These are the values the plots
// draw; byte-identity of the PNG is out of scope (different renderer
// than upstream Cairo/GD), so the data layer is what we pin.
func TestParseGDExtractsSeries(t *testing.T) {
	d := loadFixture(t)

	if d.NumSeqs != 12 {
		t.Fatalf("numseqs = %d, want 12", d.NumSeqs)
	}
	if d.NumBases != 1150 {
		t.Fatalf("numbases = %d, want 1150", d.NumBases)
	}
	if d.MaxLength != 200 {
		t.Fatalf("maxlength = %d, want 200", d.MaxLength)
	}
	if d.BinVal != 2 {
		t.Fatalf("binval = %d, want 2", d.BinVal)
	}
	if d.Format1 != "fastq" {
		t.Fatalf("format1 = %q, want fastq", d.Format1)
	}
	if d.Tail != 1 {
		t.Fatalf("tail = %d, want 1", d.Tail)
	}

	// counts.length: {"50":3,"200":1,"100":8}
	wantLen := map[int]int{50: 3, 200: 1, 100: 8}
	if !eqIntMap(d.Counts["length"], wantLen) {
		t.Fatalf("counts.length = %v, want %v", d.Counts["length"], wantLen)
	}
	// counts.gc: {"50":5,"53":1,"0":1,"49":1,"56":2,"57":2}
	wantGC := map[int]int{50: 5, 53: 1, 0: 1, 49: 1, 56: 2, 57: 2}
	if !eqIntMap(d.Counts["gc"], wantGC) {
		t.Fatalf("counts.gc = %v, want %v", d.Counts["gc"], wantGC)
	}
	// counts.ns: {"10":2}
	if !eqIntMap(d.Counts["ns"], map[int]int{10: 2}) {
		t.Fatalf("counts.ns = %v, want {10:2}", d.Counts["ns"])
	}
	// counts.tail5 / tail3
	if !eqIntMap(d.Counts["tail5"], map[int]int{50: 1}) {
		t.Fatalf("counts.tail5 = %v, want {50:1}", d.Counts["tail5"])
	}
	if !eqIntMap(d.Counts["tail3"], map[int]int{50: 1}) {
		t.Fatalf("counts.tail3 = %v, want {50:1}", d.Counts["tail3"])
	}

	// stats.length: mean 95.83 std 37.96 min 50 max 200 mode 100.
	sl := d.Stats["length"]
	if sl.Min != 50 || sl.Max != 200 || sl.Mode != 100 || sl.Modeval != 8 {
		t.Fatalf("stats.length unexpected: %+v", sl)
	}
	if !approx(sl.Mean, 95.83) || !approx(sl.Std, 37.96) {
		t.Fatalf("stats.length mean/std = %v/%v, want 95.83/37.96", sl.Mean, sl.Std)
	}

	// qualsmean histogram.
	wantQM := map[int]int{27: 1, 11: 1, 32: 1, 33: 2, 18: 1, 29: 1, 17: 4, 20: 1}
	if !eqIntMap(d.QualsMean, wantQM) {
		t.Fatalf("qualsmean = %v, want %v", d.QualsMean, wantQM)
	}

	// compldust / complentropy histograms.
	wantDust := map[int]int{4: 2, 3: 2, 23: 4, 100: 1, 2: 1, 5: 2}
	if !eqIntMap(d.ComplDust, wantDust) {
		t.Fatalf("compldust = %v, want %v", d.ComplDust, wantDust)
	}
	wantEnt := map[int]int{81: 1, 35: 4, 77: 1, 0: 1, 73: 3, 76: 2}
	if !eqIntMap(d.ComplEntropy, wantEnt) {
		t.Fatalf("complentropy = %v, want %v", d.ComplEntropy, wantEnt)
	}

	// dinucodds: 10 keys, values parsed as floats.
	if len(d.DinucOdds) != 10 {
		t.Fatalf("dinucodds has %d keys, want 10", len(d.DinucOdds))
	}
	if !approx(d.DinucOdds["CG"], 1.916620315) {
		t.Fatalf("dinucodds[CG] = %v, want ~1.9166", d.DinucOdds["CG"])
	}
	if !approx(d.DinucOdds["AATT"], 0.395437960) {
		t.Fatalf("dinucodds[AATT] = %v, want ~0.3954", d.DinucOdds["AATT"])
	}

	// dubscounts: {"1":{"0":1},"3":{"0":1}}; dubslength: {"100":{"0":4}}
	if d.DubsCounts[1][0] != 1 || d.DubsCounts[3][0] != 1 {
		t.Fatalf("dubscounts = %v", d.DubsCounts)
	}
	if d.DubsLength[100][0] != 4 {
		t.Fatalf("dubslength = %v", d.DubsLength)
	}

	// quals (relative) and qualsbin tables have 100 positions each.
	if len(d.Quals) != 100 {
		t.Fatalf("quals has %d positions, want 100", len(d.Quals))
	}
	if len(d.QualsBin) != 100 {
		t.Fatalf("qualsbin has %d positions, want 100", len(d.QualsBin))
	}
	// Spot-check one quals position (pos 0): min 0 p25 13 median 16 p75 33 max 37.
	q0 := d.Quals[0]
	if q0.Min != 0 || q0.P25 != 13 || q0.Median != 16 || q0.P75 != 33 || q0.Max != 37 {
		t.Fatalf("quals[0] = %+v", q0)
	}
}

// TestConvertOdToBinMatrix pins the binning helper against upstream's
// algorithm for the length histogram (binval 2, xmax 200).
func TestConvertOdToBinMatrix(t *testing.T) {
	d := loadFixture(t)
	matrix, xmax, ymax := convertOdToBinMatrix(d.Counts["length"], 1, false)
	// bin=getBinVal(200)=2, xmax=bin*100=200.
	if xmax != 200 {
		t.Fatalf("xmax = %d, want 200", xmax)
	}
	// Bins of width 2 over i=1..200. Length values 50 (count 3),
	// 100 (count 8), 200 (count 1). Bin index = ceil position.
	// i=50 falls in the 25th bin (i=49,50), i=100 in the 50th
	// (i=99,100), i=200 in the 100th (i=199,200).
	if len(matrix) != 100 {
		t.Fatalf("matrix len = %d, want 100", len(matrix))
	}
	total := 0
	for _, v := range matrix {
		total += v
	}
	if total != 12 {
		t.Fatalf("binned total = %d, want 12 (all reads)", total)
	}
	// ymax must be the max bin rounded up to a multiple of 4. The
	// busiest bin holds the 8 length-100 reads -> 8.
	if ymax != 8 {
		t.Fatalf("ymax = %d, want 8", ymax)
	}
}

// TestConvertToBoxValues pins the boxplot conversion: 100 ordered
// columns with the five-number summary preserved.
func TestConvertToBoxValues(t *testing.T) {
	d := loadFixture(t)
	matrix, ymax := convertToBoxValues(d.Quals, 4)
	if len(matrix) != 100 {
		t.Fatalf("box matrix len = %d, want 100", len(matrix))
	}
	// Ordered by position ascending.
	for i := 1; i < len(matrix); i++ {
		if matrix[i].pos < matrix[i-1].pos {
			t.Fatalf("box positions not sorted at %d", i)
		}
	}
	// First column corresponds to quals[0].
	if matrix[0].min != 0 || matrix[0].max != 37 || matrix[0].median != 16 {
		t.Fatalf("box[0] = %+v", matrix[0])
	}
	// ymax is rounded up to a multiple of 4; the max observed quality
	// across positions is 39 -> 40.
	if ymax != 40 {
		t.Fatalf("ymax = %d, want 40", ymax)
	}
}

// TestConvertToBarValues pins the bar (qualsmean) conversion.
func TestConvertToBarValues(t *testing.T) {
	d := loadFixture(t)
	matrix, xmax, ymax := convertToBarValues(d.QualsMean, 5, 1)
	// max key is 33; niceval 5 rounds xmax up to 35.
	if xmax != 35 {
		t.Fatalf("xmax = %d, want 35", xmax)
	}
	// Dense slice from start=1..35 => 35 entries.
	if len(matrix) != 35 {
		t.Fatalf("matrix len = %d, want 35", len(matrix))
	}
	// matrix[idx] where idx = value-1: value 17 has count 4 -> idx 16.
	if matrix[16] != 4 {
		t.Fatalf("matrix[16] (value 17) = %d, want 4", matrix[16])
	}
	// ymax rounded up to multiple of 4 from the max count (4) -> 4.
	if ymax != 4 {
		t.Fatalf("ymax = %d, want 4", ymax)
	}
}

// TestConvertOdToStackBinMatrix pins the stacked duplicate conversion.
func TestConvertOdToStackBinMatrix(t *testing.T) {
	d := loadFixture(t)
	matrix, xmax, ymax := convertOdToStackBinMatrix(d.DubsCounts, 5, 1, 100, false)
	if xmax != 100 {
		t.Fatalf("xmax = %d, want 100", xmax)
	}
	if len(matrix) != 5 {
		t.Fatalf("stack layers = %d, want 5", len(matrix))
	}
	// Layer 0 (exact dupl.) holds all the duplicate counts; total is 2
	// (precounts 1 and 3 each carry one exact-dup group).
	total := 0
	for _, v := range matrix[0] {
		total += v
	}
	if total != 2 {
		t.Fatalf("stack layer0 total = %d, want 2", total)
	}
	if ymax < 1 {
		t.Fatalf("ymax = %d, want >=1", ymax)
	}
}

// TestRenderGDGraphsSet asserts the renderer emits exactly the set of
// graphs (and filenames) upstream prinseq-graphs.pl produces for this
// .gd, in the same emission order.
func TestRenderGDGraphsSet(t *testing.T) {
	d := loadFixture(t)
	graphs := RenderGDGraphs(d)

	wantOrder := []string{
		"_ld.png", "_td5.png", "_td3.png", "_ns.png", "_gc.png",
		"_cd.png", "_ce.png", "_pm.png", "_pv.png", "_or.png",
		"_qd.png", "_qd2.png", "_qd3.png", "_df.png", "_dl.png", "_dm.png",
	}
	if len(graphs) != len(wantOrder) {
		var got []string
		for _, g := range graphs {
			got = append(got, g.Suffix)
		}
		t.Fatalf("got %d graphs %v, want %d %v", len(graphs), got, len(wantOrder), wantOrder)
	}
	for i, g := range graphs {
		if g.Suffix != wantOrder[i] {
			t.Fatalf("graph[%d] suffix = %q, want %q", i, g.Suffix, wantOrder[i])
		}
		if g.canvas == nil || g.canvas.img == nil {
			t.Fatalf("graph[%d] (%s) has nil canvas", i, g.Suffix)
		}
	}
}

// TestWriteGDGraphsValidPNG renders the full set to disk and decodes
// every PNG back, asserting valid dimensions and non-blank content.
func TestWriteGDGraphsValidPNG(t *testing.T) {
	d := loadFixture(t)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "ex1")

	written, err := WriteGDGraphs(d, prefix, true)
	if err != nil {
		t.Fatalf("WriteGDGraphs: %v", err)
	}
	// 16 PNGs + 1 HTML.
	if len(written) != 17 {
		t.Fatalf("wrote %d files, want 17", len(written))
	}

	pngCount := 0
	for _, path := range written {
		if strings.HasSuffix(path, ".html") {
			b, err := os.ReadFile(path)
			if err != nil {
				t.Fatalf("reading html %s: %v", path, err)
			}
			if !bytes.Contains(b, []byte("<img")) {
				t.Fatalf("html %s missing <img> tags", path)
			}
			continue
		}
		pngCount++
		f, err := os.Open(path)
		if err != nil {
			t.Fatalf("opening %s: %v", path, err)
		}
		img, format, err := image.Decode(f)
		f.Close()
		if err != nil {
			t.Fatalf("decoding %s: %v", path, err)
		}
		if format != "png" {
			t.Fatalf("%s decoded as %q, want png", path, format)
		}
		b := img.Bounds()
		if b.Dx() < 50 || b.Dy() < 50 {
			t.Fatalf("%s too small: %dx%d", path, b.Dx(), b.Dy())
		}
		if blankImage(img) {
			t.Fatalf("%s is blank (all one colour)", path)
		}
	}
	if pngCount != 16 {
		t.Fatalf("decoded %d PNGs, want 16", pngCount)
	}
}

// TestRenderEncodeDecodeRoundTrip independently confirms each canvas
// re-encodes/decodes through image/png and is non-blank.
func TestRenderEncodeDecodeRoundTrip(t *testing.T) {
	d := loadFixture(t)
	for _, g := range RenderGDGraphs(d) {
		var buf bytes.Buffer
		if err := png.Encode(&buf, g.canvas.img); err != nil {
			t.Fatalf("encoding %s: %v", g.Suffix, err)
		}
		img, err := png.Decode(&buf)
		if err != nil {
			t.Fatalf("decoding %s: %v", g.Suffix, err)
		}
		if blankImage(img) {
			t.Fatalf("graph %s round-trips to a blank image", g.Suffix)
		}
	}
}

// TestDinucOddsRowOrder pins the sorted-key flattening used for PCA.
func TestDinucOddsRowOrder(t *testing.T) {
	d := loadFixture(t)
	row := dinucOddsRow(d.DinucOdds)
	if len(row) != 10 {
		t.Fatalf("row len = %d, want 10", len(row))
	}
	keys := make([]string, 0, len(d.DinucOdds))
	for k := range d.DinucOdds {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for i, k := range keys {
		if !approx(row[i], d.DinucOdds[k]) {
			t.Fatalf("row[%d] = %v, want %v (key %s)", i, row[i], d.DinucOdds[k], k)
		}
	}
}

// TestParseGDErrors asserts malformed/empty input is rejected.
func TestParseGDErrors(t *testing.T) {
	if _, err := ParseGD(strings.NewReader("#only comments\n#nothing else\n")); err == nil {
		t.Fatal("expected error on comment-only input")
	}
	if _, err := ParseGD(strings.NewReader("not json at all")); err == nil {
		t.Fatal("expected error on non-JSON input")
	}
	// A minimal valid object parses.
	d, err := ParseGD(strings.NewReader(`{"numseqs":3,"numbases":9,"counts":{"length":{"3":3}}}`))
	if err != nil {
		t.Fatalf("unexpected error on minimal object: %v", err)
	}
	if d.NumSeqs != 3 || d.Counts["length"][3] != 3 {
		t.Fatalf("minimal parse wrong: %+v", d)
	}
}

func eqIntMap(a, b map[int]int) bool {
	if len(a) != len(b) {
		return false
	}
	for k, v := range b {
		if a[k] != v {
			return false
		}
	}
	return true
}

func approx(a, b float64) bool {
	d := a - b
	if d < 0 {
		d = -d
	}
	return d < 1e-6
}

// blankImage reports whether every pixel of img is the same colour
// (i.e. nothing was drawn beyond the background).
func blankImage(img image.Image) bool {
	b := img.Bounds()
	var first color.Color
	for y := b.Min.Y; y < b.Max.Y; y++ {
		for x := b.Min.X; x < b.Max.X; x++ {
			c := img.At(x, y)
			if first == nil {
				first = c
				continue
			}
			r1, g1, b1, a1 := first.RGBA()
			r2, g2, b2, a2 := c.RGBA()
			if r1 != r2 || g1 != g2 || b1 != b2 || a1 != a2 {
				return false
			}
		}
	}
	return true
}
