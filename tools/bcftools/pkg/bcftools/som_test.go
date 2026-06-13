package bcftools

import (
	"bytes"
	"io"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// --- helpers --------------------------------------------------------------

const somFixtureVCF = `##fileformat=VCFv4.2
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction MQ0">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Base-quality bias">
##INFO=<ID=MQB,Number=1,Type=Float,Description="MQ bias">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Read-pos bias">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation bias">
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	G	200	PASS	MQ=60;MQ0F=0;BQB=1.0;MQB=1.0;RPB=0.9;SGB=-0.6
1	200	.	C	T	180	PASS	MQ=58;MQ0F=0.01;BQB=0.95;MQB=0.98;RPB=0.8;SGB=-0.5
1	300	.	G	A	30	LowQual	MQ=20;MQ0F=0.4;BQB=0.2;MQB=0.1;RPB=0.1;SGB=-0.1
1	400	.	T	C	190	PASS	MQ=59;MQ0F=0;BQB=0.99;MQB=0.97;RPB=0.85;SGB=-0.55
1	500	.	A	T	25	LowQual	MQ=18;MQ0F=0.5;BQB=0.15;MQB=0.12;RPB=0.05;SGB=-0.08
1	600	.	G	C	210	PASS	MQ=60;MQ0F=0;BQB=1.0;MQB=1.0;RPB=0.92;SGB=-0.62
`

// somWriteFixture writes the VCF text to a temp file and returns its path.
func somWriteFixture(t *testing.T, text string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(p, []byte(text), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	return p
}

// --- SOM math unit tests --------------------------------------------------

// TestSomIdxToNdim checks the flat-index → grid-coordinate conversion on a
// hand-checkable 3x3 (nbin=3, ndim=2) map.
func TestSomIdxToNdim(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	s := newSom(2, 3, 1, 10, 1.0, rng)
	out := make([]int, 2)
	cases := []struct {
		idx      int
		row, col int
	}{
		{0, 0, 0},
		{1, 0, 1},
		{2, 0, 2},
		{3, 1, 0},
		{4, 1, 1},
		{8, 2, 2},
	}
	for _, c := range cases {
		s.idxToNdim(c.idx, out)
		if out[0] != c.row || out[1] != c.col {
			t.Errorf("idx %d: got (%d,%d), want (%d,%d)", c.idx, out[0], out[1], c.row, c.col)
		}
	}
}

// TestSomFindBMU checks the best-matching-unit search on a tiny hand-set map.
func TestSomFindBMU(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	s := newSom(2, 2, 2, 10, 1.0, rng) // 4 nodes, 2-D vectors
	// Hand-set the weights so node 2 is closest to (0.9, 0.9).
	s.w = []float64{
		0.0, 0.0, // node 0
		0.5, 0.5, // node 1
		0.9, 0.9, // node 2 (target)
		0.2, 0.8, // node 3
	}
	idx, dist := s.findBMU([]float64{0.9, 0.9})
	if idx != 2 {
		t.Fatalf("BMU idx = %d, want 2", idx)
	}
	if dist != 0 {
		t.Fatalf("BMU dist = %v, want 0", dist)
	}
	// A vector right at node 0.
	idx0, dist0 := s.findBMU([]float64{0.0, 0.0})
	if idx0 != 0 || dist0 != 0 {
		t.Fatalf("BMU(0,0) = (%d,%v), want (0,0)", idx0, dist0)
	}
}

// TestSomTrainSiteMovesBMUTowardInput verifies that presenting a vector
// pulls the BMU's weights toward that vector (the core SOM property).
func TestSomTrainSiteMovesBMUTowardInput(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	s := newSom(2, 2, 2, 1, 1.0, rng)
	vec := []float64{1.0, 1.0}
	before, _ := s.findBMU(vec)
	distBefore := vecDist(s.w[before*2:before*2+2], vec)
	s.trainSite(vec, true)
	distAfter := vecDist(s.w[before*2:before*2+2], vec)
	if distAfter >= distBefore {
		t.Fatalf("BMU did not move toward input: before=%v after=%v", distBefore, distAfter)
	}
	// A good site must have accumulated some count somewhere.
	total := 0.0
	for _, c := range s.c {
		total += c
	}
	if total <= 0 {
		t.Fatalf("good site accumulated no count (total=%v)", total)
	}
}

// TestSomTrainSiteBadDoesNotCount confirms a bad site shapes the map but
// does not accumulate counts (upstream's update_counts==0 path).
func TestSomTrainSiteBadDoesNotCount(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	s := newSom(2, 3, 2, 5, 1.0, rng)
	s.trainSite([]float64{0.3, 0.7}, false)
	for i, c := range s.c {
		if c != 0 {
			t.Fatalf("node %d got count %v from a bad site, want 0", i, c)
		}
	}
}

// TestSomNormCounts checks the max→1 scaling and the all-zero no-op.
func TestSomNormCounts(t *testing.T) {
	s := &som{size: 3, c: []float64{1, 2, 4}}
	s.normCounts()
	want := []float64{0.25, 0.5, 1.0}
	for i := range want {
		if math.Abs(s.c[i]-want[i]) > 1e-12 {
			t.Errorf("c[%d] = %v, want %v", i, s.c[i], want[i])
		}
	}
	// All-zero counts are left untouched (no divide-by-zero).
	z := &som{size: 2, c: []float64{0, 0}}
	z.normCounts()
	if z.c[0] != 0 || z.c[1] != 0 {
		t.Fatalf("all-zero normCounts mutated: %v", z.c)
	}
}

// TestSomGetScoreThreshold confirms getScore only considers nodes whose
// count clears the threshold, and returns +Inf when none do.
func TestSomGetScoreThreshold(t *testing.T) {
	s := &som{
		size: 2,
		kdim: 2,
		w:    []float64{0.0, 0.0, 1.0, 1.0},
		c:    []float64{0.5, 1.0},
	}
	// Threshold 0.9: only node 1 qualifies. Distance from (1,1) to node1=0.
	got := s.getScore([]float64{1.0, 1.0}, 0.9)
	if got != 0 {
		t.Fatalf("score = %v, want 0 (node1 exact)", got)
	}
	// Threshold above all counts: no node qualifies → +Inf.
	got = s.getScore([]float64{1.0, 1.0}, 2.0)
	if !math.IsInf(got, 1) {
		t.Fatalf("score = %v, want +Inf", got)
	}
}

func vecDist(a, b []float64) float64 {
	d := 0.0
	for i := range a {
		x := a[i] - b[i]
		d += x * x
	}
	return math.Sqrt(d)
}

// --- annotation extraction ------------------------------------------------

func TestAnnotValue(t *testing.T) {
	v := &vcf.Variant{
		Qual: 42.5,
		Info: map[string]string{"MQ": "60", "DP4": "1,2,3,4", "BAD": "x"},
	}
	cases := []struct {
		name string
		want float64
	}{
		{"QUAL", 42.5},
		{"MQ", 60},
		{"DP4", 1}, // first field of multi-value
		{"MISSING", 0},
		{"BAD", 0}, // unparsable → 0
	}
	for _, c := range cases {
		if got := annotValue(v, c.name); got != c.want {
			t.Errorf("annotValue(%q) = %v, want %v", c.name, got, c.want)
		}
	}
	// Missing QUAL (-1) reads as 0.
	if got := annotValue(&vcf.Variant{Qual: -1}, "QUAL"); got != 0 {
		t.Errorf("missing QUAL = %v, want 0", got)
	}
}

func TestAnnotExtractorNormalize(t *testing.T) {
	e := newAnnotExtractor([]string{"A", "B"})
	e.observe([]float64{0, 10})
	e.observe([]float64{100, 10}) // column B is degenerate (min==max)
	// Column A normalises 0→0, 50→0.5, 100→1; column B is degenerate→0.5.
	out := e.normalize([]float64{50, 999})
	if math.Abs(out[0]-0.5) > 1e-12 {
		t.Errorf("col A norm = %v, want 0.5", out[0])
	}
	if out[1] != 0.5 {
		t.Errorf("degenerate col B norm = %v, want 0.5", out[1])
	}
	// Out-of-range inputs clamp to [0,1].
	out = e.normalize([]float64{-50, 0})
	if out[0] != 0 {
		t.Errorf("below-min clamp = %v, want 0", out[0])
	}
	out = e.normalize([]float64{500, 0})
	if out[0] != 1 {
		t.Errorf("above-max clamp = %v, want 1", out[0])
	}
}

func TestSiteClass(t *testing.T) {
	good := &vcf.Variant{Filter: []string{"PASS"}}
	if siteClass(good, 2, 1) != 2 {
		t.Error("PASS should be good class")
	}
	dot := &vcf.Variant{Filter: []string{"."}}
	if siteClass(dot, 2, 1) != 2 {
		t.Error(". should be good class")
	}
	none := &vcf.Variant{Filter: nil}
	if siteClass(none, 2, 1) != 2 {
		t.Error("empty FILTER should be good class")
	}
	bad := &vcf.Variant{Filter: []string{"LowQual"}}
	if siteClass(bad, 2, 1) != 1 {
		t.Error("LowQual should be bad class")
	}
}

// --- on-disk map round-trip ----------------------------------------------

// TestMapRoundTripInMemory writes a model to a buffer and reads it back,
// asserting identical contents.
func TestMapRoundTripInMemory(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	m := newSom(2, 4, 3, 20, 1.0, rng)
	for i := range m.c {
		m.c[i] = rng.Float64()
	}
	orig := &somModel{
		annots: []string{"QUAL", "MQ", "BQB"},
		min:    []float64{0, 1, 2},
		max:    []float64{100, 60, 1},
		bmuTh:  0.9,
		maps:   []*som{m},
	}
	var buf bytes.Buffer
	if err := writeMaps(&buf, orig); err != nil {
		t.Fatalf("writeMaps: %v", err)
	}
	got, err := readMaps(&buf)
	if err != nil {
		t.Fatalf("readMaps: %v", err)
	}
	assertModelsEqual(t, orig, got)
}

// TestMapRoundTripOnDisk exercises the file path (writeMapFile/readMapFile).
func TestMapRoundTripOnDisk(t *testing.T) {
	rng := rand.New(rand.NewSource(11))
	m := newSom(2, 3, 2, 10, 1.0, rng)
	orig := &somModel{
		annots: []string{"QUAL", "MQ"},
		min:    []float64{0, 1},
		max:    []float64{99, 60},
		bmuTh:  0.8,
		maps:   []*som{m},
	}
	prefix := filepath.Join(t.TempDir(), "model")
	if err := writeMapFile(prefix, orig); err != nil {
		t.Fatalf("writeMapFile: %v", err)
	}
	got, err := readMapFile(prefix)
	if err != nil {
		t.Fatalf("readMapFile: %v", err)
	}
	assertModelsEqual(t, orig, got)
}

func TestReadMapsRejectsBadMagic(t *testing.T) {
	var buf bytes.Buffer
	buf.WriteString("NOPEXX")
	if _, err := readMaps(&buf); err == nil {
		t.Fatal("readMaps accepted bad magic")
	}
}

// TestReadMapsTruncated confirms a short file is reported rather than read
// past EOF — the symmetric guard to the upstream write bug.
func TestReadMapsTruncated(t *testing.T) {
	rng := rand.New(rand.NewSource(5))
	m := newSom(2, 3, 2, 10, 1.0, rng)
	orig := &somModel{annots: []string{"QUAL", "MQ"}, min: []float64{0, 0}, max: []float64{1, 1}, bmuTh: 0.9, maps: []*som{m}}
	var buf bytes.Buffer
	if err := writeMaps(&buf, orig); err != nil {
		t.Fatalf("writeMaps: %v", err)
	}
	full := buf.Bytes()
	for _, n := range []int{0, 3, 10, 20, len(full) - 1} {
		if n < 0 {
			continue
		}
		if _, err := readMaps(bytes.NewReader(full[:n])); err == nil {
			t.Errorf("readMaps accepted truncation at %d bytes", n)
		}
	}
}

// failWriter fails after letting limit bytes through, exercising the
// write-error paths that the upstream som_write_map fails to check.
type failWriter struct {
	limit int
	n     int
}

func (w *failWriter) Write(p []byte) (int, error) {
	if w.n+len(p) > w.limit {
		allowed := w.limit - w.n
		if allowed < 0 {
			allowed = 0
		}
		w.n = w.limit
		return allowed, io.ErrShortWrite
	}
	w.n += len(p)
	return len(p), nil
}

func TestWriteMapsReportsWriteError(t *testing.T) {
	rng := rand.New(rand.NewSource(8))
	m := newSom(2, 3, 2, 10, 1.0, rng)
	model := &somModel{annots: []string{"QUAL", "MQ"}, min: []float64{0, 0}, max: []float64{1, 1}, bmuTh: 0.9, maps: []*som{m}}
	for _, limit := range []int{2, 7, 15, 40} {
		if err := writeMaps(&failWriter{limit: limit}, model); err == nil {
			t.Errorf("writeMaps did not report a write failure at limit %d", limit)
		}
	}
}

func assertModelsEqual(t *testing.T, a, b *somModel) {
	t.Helper()
	if a.bmuTh != b.bmuTh {
		t.Errorf("bmuTh: %v vs %v", a.bmuTh, b.bmuTh)
	}
	if len(a.annots) != len(b.annots) {
		t.Fatalf("annot count: %d vs %d", len(a.annots), len(b.annots))
	}
	for i := range a.annots {
		if a.annots[i] != b.annots[i] || a.min[i] != b.min[i] || a.max[i] != b.max[i] {
			t.Errorf("annot %d differs: (%s,%v,%v) vs (%s,%v,%v)", i,
				a.annots[i], a.min[i], a.max[i], b.annots[i], b.min[i], b.max[i])
		}
	}
	if len(a.maps) != len(b.maps) {
		t.Fatalf("map count: %d vs %d", len(a.maps), len(b.maps))
	}
	for mi := range a.maps {
		am, bm := a.maps[mi], b.maps[mi]
		if am.ndim != bm.ndim || am.nbin != bm.nbin || am.kdim != bm.kdim || am.size != bm.size {
			t.Fatalf("map %d dims differ", mi)
		}
		for i := range am.w {
			if am.w[i] != bm.w[i] {
				t.Fatalf("map %d w[%d]: %v vs %v", mi, i, am.w[i], bm.w[i])
			}
		}
		for i := range am.c {
			if am.c[i] != bm.c[i] {
				t.Fatalf("map %d c[%d]: %v vs %v", mi, i, am.c[i], bm.c[i])
			}
		}
	}
}

// --- end-to-end train→classify ------------------------------------------

// TestTrainClassifyRoundtrip trains on the fixture and classifies it,
// asserting the map round-trips, scores are deterministic and in range,
// and good (PASS) sites score higher than bad (LowQual) sites.
func TestTrainClassifyRoundtrip(t *testing.T) {
	in := somWriteFixture(t, somFixtureVCF)
	prefix := filepath.Join(t.TempDir(), "som")
	annots := []string{"QUAL", "MQ", "MQ0F", "BQB", "MQB", "RPB", "SGB"}

	trainOpts := SomOptions{
		Action:         SomActionTrain,
		Prefix:         prefix,
		TrainingAnnots: annots,
		Size:           8,
		RandomSeed:     1,
	}
	res, err := SomTrain(in, trainOpts)
	if err != nil {
		t.Fatalf("SomTrain: %v", err)
	}
	if res.NSites != 6 {
		t.Errorf("NSites = %d, want 6", res.NSites)
	}
	if res.KDim != len(annots) {
		t.Errorf("KDim = %d, want %d", res.KDim, len(annots))
	}
	if res.MapSize != 64 {
		t.Errorf("MapSize = %d, want 64", res.MapSize)
	}

	// The .som file must exist and re-load.
	if _, err := os.Stat(prefix + ".som"); err != nil {
		t.Fatalf("map file not written: %v", err)
	}

	classifyOpts := SomOptions{
		Action:         SomActionClassify,
		Prefix:         prefix,
		TrainingAnnots: annots,
	}
	var buf bytes.Buffer
	scores, err := SomClassify(in, &buf, classifyOpts)
	if err != nil {
		t.Fatalf("SomClassify: %v", err)
	}
	if len(scores) != 6 {
		t.Fatalf("len(scores) = %d, want 6", len(scores))
	}
	for _, s := range scores {
		if math.IsNaN(s.Score) || s.Score < -1.0 || s.Score > 1.0 {
			t.Errorf("score for %s:%d = %v out of expected range", s.Chrom, s.Pos, s.Score)
		}
	}

	// Determinism: a second classify pass gives identical scores.
	scores2, err := SomClassify(in, nil, classifyOpts)
	if err != nil {
		t.Fatalf("SomClassify (2nd): %v", err)
	}
	for i := range scores {
		if scores[i].Score != scores2[i].Score {
			t.Errorf("classify non-deterministic at %d: %v vs %v", i, scores[i].Score, scores2[i].Score)
		}
	}

	// Good (PASS) sites should on average score higher than bad ones.
	// PASS sites are positions 100,200,400,600; LowQual are 300,500.
	goodSum, goodN, badSum, badN := 0.0, 0, 0.0, 0
	for _, s := range scores {
		switch s.Pos {
		case 300, 500:
			badSum += s.Score
			badN++
		default:
			goodSum += s.Score
			goodN++
		}
	}
	goodAvg := goodSum / float64(goodN)
	badAvg := badSum / float64(badN)
	if !(goodAvg > badAvg) {
		t.Errorf("expected good sites to score higher: good=%v bad=%v", goodAvg, badAvg)
	}
}

// TestTrainDeterministicWeights confirms two trains with the same seed
// produce byte-identical map files (the property upstream lacks).
func TestTrainDeterministicWeights(t *testing.T) {
	in := somWriteFixture(t, somFixtureVCF)
	p1 := filepath.Join(t.TempDir(), "a")
	p2 := filepath.Join(t.TempDir(), "b")
	opts := SomOptions{Action: SomActionTrain, TrainingAnnots: []string{"QUAL", "MQ"}, Size: 5, RandomSeed: 99}
	o1 := opts
	o1.Prefix = p1
	o2 := opts
	o2.Prefix = p2
	if _, err := SomTrain(in, o1); err != nil {
		t.Fatalf("train 1: %v", err)
	}
	if _, err := SomTrain(in, o2); err != nil {
		t.Fatalf("train 2: %v", err)
	}
	b1, _ := os.ReadFile(p1 + ".som")
	b2, _ := os.ReadFile(p2 + ".som")
	if !bytes.Equal(b1, b2) {
		t.Fatal("same-seed trains produced different map files")
	}
}

// TestSomDispatch checks the Som() top-level routing and default annots.
func TestSomDispatch(t *testing.T) {
	in := somWriteFixture(t, somFixtureVCF)
	prefix := filepath.Join(t.TempDir(), "d")
	// Train via Som() with default annotations.
	if err := Som(in, nil, SomOptions{Action: SomActionTrain, Prefix: prefix}); err != nil {
		t.Fatalf("Som train: %v", err)
	}
	var buf bytes.Buffer
	if err := Som(in, &buf, SomOptions{Action: SomActionClassify, Prefix: prefix}); err != nil {
		t.Fatalf("Som classify: %v", err)
	}
	if buf.Len() == 0 {
		t.Fatal("classify produced no output")
	}
	// No action is an error.
	if err := Som(in, nil, SomOptions{Prefix: prefix}); err == nil {
		t.Fatal("Som with no action should error")
	}
}

func TestSomTrainErrors(t *testing.T) {
	in := somWriteFixture(t, somFixtureVCF)
	// Missing prefix.
	if _, err := SomTrain(in, SomOptions{Action: SomActionTrain}); err == nil {
		t.Error("train without prefix should error")
	}
	// Missing input file.
	if _, err := SomTrain(filepath.Join(t.TempDir(), "nope.vcf"), SomOptions{Prefix: "x"}); err == nil {
		t.Error("train on missing file should error")
	}
	// Classify without a map file.
	if _, err := SomClassify(in, nil, SomOptions{Prefix: filepath.Join(t.TempDir(), "absent")}); err == nil {
		t.Error("classify without a map should error")
	}
	// Empty VCF.
	empty := somWriteFixture(t, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	if _, err := SomTrain(empty, SomOptions{Action: SomActionTrain, Prefix: filepath.Join(t.TempDir(), "e")}); err == nil {
		t.Error("train on empty VCF should error")
	}
}
