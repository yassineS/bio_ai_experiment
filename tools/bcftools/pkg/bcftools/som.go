// bcftools som — Self-Organizing Map (Kohonen map) variant classifier (port).
//
// Upstream reference: reference_code/bcftools/vcfsom.c. The tool trains a
// 2-D self-organizing map over per-site INFO annotation vectors and then
// scores new sites by their distance to the trained map. It has two modes,
// `--train` and `--classify`.
//
// # What this port changes vs upstream
//
//  1. Upstream's som_write_map (vcfsom.c:170) is broken: it checks
//     fwrite("SOMv1",5,1,fp)!=5, but fwrite returns the *element* count
//     (1 here), so the comparison is always true and `--train` calls
//     error() and exits 255 after truncating the `.som` file to 5 bytes.
//     As a result no usable map is ever written and `--classify` always
//     fails. This port writes the full map and validates byte counts
//     correctly, so the train→classify pipeline actually works.
//     See docs/UPSTREAM_BUGS.md#bcftools-som-write-map.
//
//  2. Upstream reads a pre-extracted annots.tab.gz table. This port reads
//     a VCF/BCF directly and extracts the per-site annotation vector from
//     the INFO fields named by -t/--training-annots, which is the pipeline
//     a Go re-implementation should expose. Each annotation column is
//     min/max normalised onto [0,1] before training/scoring so that no
//     single annotation dominates the Euclidean distance.
//
//  3. The on-disk map format is our own, clean, versioned binary format
//     (magic "SOMGO1"; little-endian) because upstream's is unusable. It
//     is documented in writeMaps / readMaps below.
//
//  4. Weight initialisation uses Go's math/rand with a caller-supplied
//     seed instead of glibc's random() (TYPE_3 additive-feedback PRNG).
//     Reproducing glibc's PRNG byte-for-byte only matters for matching a
//     tool that crashes; the Go RNG is deterministic for a given seed,
//     which is what the round-trip and determinism tests rely on.
//
// # The SOM math (ported verbatim from vcfsom.c)
//
//   - findBMU: the best-matching unit is the node whose weight vector has
//     the smallest squared Euclidean distance to the input vector.
//   - trainSite: on each presented vector the learning rate and radius
//     decay as exp(-t/nt); every node within the radius of the BMU (in
//     map-grid coordinates) is nudged toward the input vector by
//     influence = exp(-d²·0.5/radius)·learning_rate, where d is the grid
//     distance from the BMU. Good sites also accumulate that influence
//     into the per-node count c[].
//   - getScore: at classify time a site's raw score is the Euclidean
//     distance (sqrt of the squared distance) to the nearest node whose
//     normalised count c[] is at least the BMU threshold. The reported
//     score is 1 - rawScore/sqrt(kdim), so higher is "more like the
//     training (good) sites", exactly as upstream's do_classify prints.
package bcftools

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// SomAction selects the operating mode of Som.
type SomAction int

const (
	// SomActionNone is the zero value (no action selected).
	SomActionNone SomAction = iota
	// SomActionTrain trains a map and writes <prefix>.som.
	SomActionTrain
	// SomActionClassify loads <prefix>.som and scores each site.
	SomActionClassify
)

// DefaultTrainingAnnots is the default INFO annotation set used when
// -t/--training-annots is not given. It mirrors the kind of site-quality
// annotations bcftools/GATK pipelines feed a SOM/VQSR-style classifier.
var DefaultTrainingAnnots = []string{"QUAL", "MQ", "MQ0F", "BQB", "MQB", "RPB", "SGB"}

// Upstream default parameter values.
const (
	// SomDefaultSize is the default map edge length (-s/--size), so a
	// default 2-D map has SomDefaultSize² nodes.
	SomDefaultSize = 20
	// SomDefaultLearn is the default learning rate (-l/--learning-rate).
	SomDefaultLearn = 1.0
	// SomDefaultBmuThreshold is the default best-matching-unit count
	// threshold (-b/--bmu-threshold).
	SomDefaultBmuThreshold = 0.9
	// SomDefaultRandomSeed is the default RNG seed (-r/--random-seed).
	SomDefaultRandomSeed = 1
	// SomDefaultDimension is the default map dimensionality
	// (-d/--som-dimension). Upstream requires >= 2.
	SomDefaultDimension = 2
)

// SomOptions controls Som / SomTrain / SomClassify.
type SomOptions struct {
	// Action selects train or classify mode.
	Action SomAction
	// Prefix is the output/input file prefix; the map lives at
	// <Prefix>.som.
	Prefix string
	// TrainingAnnots names the INFO tags forming the per-site vector.
	// Empty means DefaultTrainingAnnots.
	TrainingAnnots []string
	// Size is the map edge length (nbin). 0 means SomDefaultSize.
	Size int
	// NDim is the map dimensionality. 0 means SomDefaultDimension.
	NDim int
	// NTrain is the effective number of training iterations used in the
	// learning-rate decay. 0 means "number of presented sites".
	NTrain int
	// Learn is the learning rate. 0 means SomDefaultLearn.
	Learn float64
	// BmuThreshold is the count threshold for classification nodes.
	// 0 means SomDefaultBmuThreshold.
	BmuThreshold float64
	// RandomSeed seeds the weight initialisation. 0 keeps it deterministic
	// at SomDefaultRandomSeed so round-trips are reproducible.
	RandomSeed int64
	// GoodClass / BadClass are the integer class labels distinguishing
	// good (training-target) sites from bad sites, matching upstream's
	// -good/-bad defaults (2 and 1). The port derives a per-site class
	// from a FILTER==PASS heuristic; see siteClass.
	GoodClass int
	BadClass  int
}

// SomResult is returned by SomTrain. It reports what was learned so the
// CLI and tests can assert on it without re-reading the map file.
type SomResult struct {
	// NSites is the number of sites read from the training VCF.
	NSites int
	// KDim is the dimension of the annotation vector (number of annots).
	KDim int
	// Annots is the resolved annotation list.
	Annots []string
	// MapSize is the number of nodes (size = nbin^ndim).
	MapSize int
}

// som is the trained Kohonen map for one cross-validation fold.
type som struct {
	ndim int // map dimensionality (2 = 2-D grid)
	nbin int // edge length of the grid
	size int // number of nodes = nbin^ndim
	kdim int // dimension of the input vectors

	nt    float64 // learning-cycle count used in the decay exp(-t/nt)
	t     int     // current learning cycle
	learn float64 // learning rate

	w []float64 // weights, size*kdim, row-major per node
	c []float64 // per-node counts (learning influence accumulated)

	div []float64 // index→grid helpers: div[i] = nbin^(ndim-i-1)
}

// newSom allocates and randomly initialises a map. Mirrors som_init in
// vcfsom.c, but seeds weights from Go's math/rand (see file header).
func newSom(ndim, nbin, kdim, ntrain int, learn float64, rng *rand.Rand) *som {
	s := &som{
		ndim:  ndim,
		nbin:  nbin,
		kdim:  kdim,
		nt:    float64(ntrain),
		learn: learn,
	}
	s.size = 1
	for i := 0; i < ndim; i++ {
		s.size *= nbin
	}
	s.w = make([]float64, s.size*kdim)
	s.c = make([]float64, s.size)
	for i := range s.w {
		// Weights start uniformly in [0,1), the same range the
		// per-column-normalised input vectors live in. (Upstream seeds
		// with the raw random() integer, which is then pulled down toward
		// the input range during training; starting in-range is cleaner
		// and makes the count threshold meaningful sooner.)
		s.w[i] = rng.Float64()
	}
	s.div = make([]float64, ndim)
	for i := 0; i < ndim; i++ {
		s.div[i] = math.Pow(float64(nbin), float64(ndim-i-1))
	}
	return s
}

// idxToNdim converts a flat node index to its grid coordinates. Verbatim
// port of som_idx_to_ndim.
func (s *som) idxToNdim(idx int, out []int) {
	out[0] = int(float64(idx) / s.div[0])
	sub := 0.0
	for i := 1; i < s.ndim; i++ {
		sub += float64(out[i-1]) * s.div[i-1]
		out[i] = int((float64(idx) - sub) / s.div[i])
	}
}

// findBMU returns the index of the best-matching unit for vec and the
// squared distance to it. Verbatim port of som_find_bmu.
func (s *som) findBMU(vec []float64) (int, float64) {
	minDist := math.Inf(1)
	minIdx := 0
	for i := 0; i < s.size; i++ {
		base := i * s.kdim
		dist := 0.0
		for k := 0; k < s.kdim; k++ {
			d := vec[k] - s.w[base+k]
			dist += d * d
		}
		if dist < minDist {
			minDist = dist
			minIdx = i
		}
	}
	return minIdx, minDist
}

// getScore returns the Euclidean distance from vec to the nearest node
// whose count is at least bmuTh. Verbatim port of som_get_score. When no
// node clears the threshold it returns +Inf (caller maps that to score 0).
func (s *som) getScore(vec []float64, bmuTh float64) float64 {
	minDist := math.Inf(1)
	for i := 0; i < s.size; i++ {
		if s.c[i] < bmuTh {
			continue
		}
		base := i * s.kdim
		dist := 0.0
		for k := 0; k < s.kdim; k++ {
			d := vec[k] - s.w[base+k]
			dist += d * d
		}
		if dist < minDist {
			minDist = dist
		}
	}
	return math.Sqrt(minDist)
}

// trainSite presents one vector to the map and updates the weights and,
// when updateCounts is set (good sites), the per-node counts. Verbatim
// port of som_train_site.
func (s *som) trainSite(vec []float64, updateCounts bool) {
	s.t++
	dt := math.Exp(-float64(s.t) / s.nt)
	learningRate := s.learn * dt
	radius := float64(s.nbin) * dt
	radius *= radius

	aIdx := make([]int, s.ndim)
	bIdx := make([]int, s.ndim)

	minIdx, _ := s.findBMU(vec)
	s.idxToNdim(minIdx, aIdx)

	for i := 0; i < s.size; i++ {
		s.idxToNdim(i, bIdx)
		dist := 0.0
		for j := 0; j < s.ndim; j++ {
			d := float64(aIdx[j] - bIdx[j])
			dist += d * d
		}
		if dist <= radius {
			influence := math.Exp(-dist*dist*0.5/radius) * learningRate
			base := i * s.kdim
			for k := 0; k < s.kdim; k++ {
				s.w[base+k] += influence * (vec[k] - s.w[base+k])
			}
			if updateCounts {
				s.c[i] += influence
			}
		}
	}
}

// normCounts scales the per-node counts so the maximum is 1. Verbatim port
// of som_norm_counts. A map that never accumulated any count (no good
// sites, or all-zero influence) is left untouched.
func (s *som) normCounts() {
	max := 0.0
	for _, v := range s.c {
		if v > max {
			max = v
		}
	}
	if max <= 0 {
		return
	}
	for i := range s.c {
		s.c[i] /= max
	}
}

// annotExtractor pulls the per-site annotation vector out of a VCF record
// and tracks the per-column min/max so the vectors can be normalised.
type annotExtractor struct {
	annots []string
	min    []float64
	max    []float64
}

func newAnnotExtractor(annots []string) *annotExtractor {
	e := &annotExtractor{
		annots: annots,
		min:    make([]float64, len(annots)),
		max:    make([]float64, len(annots)),
	}
	for i := range annots {
		e.min[i] = math.Inf(1)
		e.max[i] = math.Inf(-1)
	}
	return e
}

// vector extracts the raw (un-normalised) annotation vector for v. Missing
// or unparsable annotations are treated as 0, matching upstream's
// atof("")==0 behaviour for absent table cells.
func (e *annotExtractor) vector(v *vcf.Variant) []float64 {
	vec := make([]float64, len(e.annots))
	for i, name := range e.annots {
		vec[i] = annotValue(v, name)
	}
	return vec
}

// observe updates the per-column min/max from a raw vector.
func (e *annotExtractor) observe(vec []float64) {
	for i, x := range vec {
		if x < e.min[i] {
			e.min[i] = x
		}
		if x > e.max[i] {
			e.max[i] = x
		}
	}
}

// normalize maps a raw vector onto [0,1] per column using the observed
// min/max. A degenerate column (min==max, or never observed) maps to 0.5
// so it contributes no spurious distance.
func (e *annotExtractor) normalize(vec []float64) []float64 {
	out := make([]float64, len(vec))
	for i, x := range vec {
		lo, hi := e.min[i], e.max[i]
		switch {
		case math.IsInf(lo, 0) || math.IsInf(hi, 0) || hi <= lo:
			out[i] = 0.5
		default:
			v := (x - lo) / (hi - lo)
			if v < 0 {
				v = 0
			} else if v > 1 {
				v = 1
			}
			out[i] = v
		}
	}
	return out
}

// annotValue resolves a single annotation by name. QUAL is taken from the
// record's QUAL column; everything else is looked up as an INFO float.
func annotValue(v *vcf.Variant, name string) float64 {
	if name == "QUAL" {
		if v.Qual < 0 {
			return 0
		}
		return v.Qual
	}
	s, ok := v.Info[name]
	if !ok || s == "" {
		return 0
	}
	// Multi-value INFO (e.g. "1,2,3"): use the first numeric field.
	for i := 0; i < len(s); i++ {
		if s[i] == ',' {
			s = s[:i]
			break
		}
	}
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0
	}
	return f
}

// siteClass returns the upstream-style class label for a record. The port
// derives "good" from FILTER==PASS (or ".") and "bad" otherwise, which is
// the natural training signal when the SOM reads a VCF rather than a
// pre-labelled annots table.
func siteClass(v *vcf.Variant, goodClass, badClass int) int {
	if len(v.Filter) == 0 {
		return goodClass
	}
	for _, f := range v.Filter {
		if f == "PASS" || f == "." || f == "" {
			return goodClass
		}
	}
	return badClass
}

// readVariants opens path and returns all variants.
func readVariants(path string) ([]*vcf.Variant, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer in.Close()
	r := vcf.NewReader(in)
	if _, err := r.ReadHeader(); err != nil {
		return nil, err
	}
	var out []*vcf.Variant
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		out = append(out, v)
	}
	return out, nil
}

// resolveOpts fills in the upstream default values for any unset field.
func resolveOpts(opts SomOptions) SomOptions {
	if len(opts.TrainingAnnots) == 0 {
		opts.TrainingAnnots = DefaultTrainingAnnots
	}
	if opts.Size <= 0 {
		opts.Size = SomDefaultSize
	}
	if opts.NDim <= 0 {
		opts.NDim = SomDefaultDimension
	}
	if opts.Learn <= 0 {
		opts.Learn = SomDefaultLearn
	}
	if opts.BmuThreshold <= 0 {
		opts.BmuThreshold = SomDefaultBmuThreshold
	}
	if opts.RandomSeed == 0 {
		opts.RandomSeed = SomDefaultRandomSeed
	}
	if opts.GoodClass == 0 && opts.BadClass == 0 {
		opts.GoodClass = 2
		opts.BadClass = 1
	}
	return opts
}

// somModel bundles the trained map(s) with the annotation normalisation
// state so classification reproduces the training-time vector scaling.
type somModel struct {
	annots []string
	min    []float64
	max    []float64
	bmuTh  float64
	maps   []*som
}

// SomTrain reads the training VCF at path, extracts the annotation vectors,
// trains the map, writes <prefix>.som, and returns a summary. This is the
// `--train` entry point.
func SomTrain(path string, opts SomOptions) (SomResult, error) {
	opts = resolveOpts(opts)
	variants, err := readVariants(path)
	if err != nil {
		return SomResult{}, fmt.Errorf("bcftools som: read %s: %w", path, err)
	}
	if len(variants) == 0 {
		return SomResult{}, fmt.Errorf("bcftools som: no sites in %s", path)
	}

	ext := newAnnotExtractor(opts.TrainingAnnots)
	raw := make([][]float64, len(variants))
	classes := make([]int, len(variants))
	for i, v := range variants {
		raw[i] = ext.vector(v)
		ext.observe(raw[i])
		classes[i] = siteClass(v, opts.GoodClass, opts.BadClass)
	}

	ntrain := opts.NTrain
	if ntrain <= 0 {
		ntrain = len(variants)
	}
	rng := rand.New(rand.NewSource(opts.RandomSeed))
	m := newSom(opts.NDim, opts.Size, len(opts.TrainingAnnots), ntrain, opts.Learn, rng)

	// Present sites in file order. Bad sites still shape the map (matching
	// upstream's train_bad default) but only good sites accumulate counts.
	for i := range variants {
		norm := ext.normalize(raw[i])
		isGood := classes[i] == opts.GoodClass
		m.trainSite(norm, isGood)
	}
	m.normCounts()

	model := &somModel{
		annots: opts.TrainingAnnots,
		min:    ext.min,
		max:    ext.max,
		bmuTh:  opts.BmuThreshold,
		maps:   []*som{m},
	}
	if err := writeMapFile(opts.Prefix, model); err != nil {
		return SomResult{}, err
	}

	return SomResult{
		NSites:  len(variants),
		KDim:    len(opts.TrainingAnnots),
		Annots:  append([]string(nil), opts.TrainingAnnots...),
		MapSize: m.size,
	}, nil
}

// SomScore is one classified site.
type SomScore struct {
	Chrom string
	Pos   int
	Score float64 // 1 - rawDist/sqrt(kdim); higher means "more good-like"
}

// SomClassify loads <prefix>.som and scores each site of the VCF at path,
// writing one "<score>" line per site to out (matching upstream's
// do_classify). It also returns the scores so callers/tests can inspect
// them. This is the `--classify` entry point.
func SomClassify(path string, out io.Writer, opts SomOptions) ([]SomScore, error) {
	opts = resolveOpts(opts)
	model, err := readMapFile(opts.Prefix)
	if err != nil {
		return nil, err
	}
	if len(model.maps) == 0 {
		return nil, fmt.Errorf("bcftools som: %s.som has no maps", opts.Prefix)
	}
	variants, err := readVariants(path)
	if err != nil {
		return nil, fmt.Errorf("bcftools som: read %s: %w", path, err)
	}

	ext := &annotExtractor{annots: model.annots, min: model.min, max: model.max}
	maxScore := math.Sqrt(float64(model.maps[0].kdim))

	scores := make([]SomScore, 0, len(variants))
	for _, v := range variants {
		norm := ext.normalize(ext.vector(v))
		raw := avgScore(model.maps, norm, model.bmuTh)
		var score float64
		if math.IsInf(raw, 1) || maxScore == 0 {
			score = 0
		} else {
			score = 1.0 - raw/maxScore
		}
		scores = append(scores, SomScore{Chrom: v.Chrom, Pos: v.Pos, Score: score})
		if out != nil {
			if _, err := fmt.Fprintf(out, "%e\n", score); err != nil {
				return nil, err
			}
		}
	}
	return scores, nil
}

// avgScore averages the per-map raw scores, skipping maps that return no
// score (no node cleared the threshold). Mirrors upstream's MERGE_AVG
// default; with a single map it is just that map's score.
func avgScore(maps []*som, vec []float64, bmuTh float64) float64 {
	sum := 0.0
	n := 0
	for _, m := range maps {
		s := m.getScore(vec, bmuTh)
		if math.IsInf(s, 1) {
			continue
		}
		sum += s
		n++
	}
	if n == 0 {
		return math.Inf(1)
	}
	return sum / float64(n)
}

// Som is the top-level dispatch used by the CLI: it routes to SomTrain or
// SomClassify based on opts.Action.
func Som(path string, out io.Writer, opts SomOptions) error {
	switch opts.Action {
	case SomActionTrain:
		_, err := SomTrain(path, opts)
		return err
	case SomActionClassify:
		_, err := SomClassify(path, out, opts)
		return err
	default:
		return fmt.Errorf("bcftools som: no action selected (use --train or --classify)")
	}
}

// --- On-disk map format ---------------------------------------------------
//
// The file <prefix>.som is a clean, versioned binary format. All integers
// and floats are little-endian. Layout:
//
//	magic     : 6 bytes  "SOMGO1"
//	bmuTh     : float64  best-matching-unit count threshold
//	nAnnot    : int32    number of annotation columns (== kdim)
//	repeat nAnnot times:
//	    nameLen : int32   length of the annotation name in bytes
//	    name    : bytes   the annotation name (UTF-8)
//	    min     : float64 training-time per-column minimum
//	    max     : float64 training-time per-column maximum
//	nMap      : int32    number of maps (cross-validation folds; 1 here)
//	repeat nMap times:
//	    ndim   : int32    map dimensionality
//	    nbin   : int32    grid edge length
//	    kdim   : int32    input-vector dimension (== nAnnot)
//	    size   : int32    node count (== nbin^ndim)
//	    w      : size*kdim float64  node weights (row-major per node)
//	    c      : size      float64  per-node normalised counts
//
// Unlike upstream's format this stores the annotation names and the
// per-column normalisation bounds so classification reproduces the
// training-time vector scaling without re-reading the training VCF.

const somMagic = "SOMGO1"

func writeMapFile(prefix string, model *somModel) error {
	if prefix == "" {
		return fmt.Errorf("bcftools som: --prefix is required to write the map")
	}
	f, err := os.Create(prefix + ".som")
	if err != nil {
		return fmt.Errorf("bcftools som: create %s.som: %w", prefix, err)
	}
	defer f.Close()
	if err := writeMaps(f, model); err != nil {
		return fmt.Errorf("bcftools som: write %s.som: %w", prefix, err)
	}
	return f.Close()
}

// writeMaps serialises model to w using the format documented above. It
// validates every write so a short or failed write is reported rather than
// silently truncating the file (the upstream bug being fixed).
func writeMaps(w io.Writer, model *somModel) error {
	if err := writeBytes(w, []byte(somMagic)); err != nil {
		return err
	}
	if err := writeFloat(w, model.bmuTh); err != nil {
		return err
	}
	if err := writeInt32(w, int32(len(model.annots))); err != nil {
		return err
	}
	for i, name := range model.annots {
		if err := writeInt32(w, int32(len(name))); err != nil {
			return err
		}
		if err := writeBytes(w, []byte(name)); err != nil {
			return err
		}
		if err := writeFloat(w, model.min[i]); err != nil {
			return err
		}
		if err := writeFloat(w, model.max[i]); err != nil {
			return err
		}
	}
	if err := writeInt32(w, int32(len(model.maps))); err != nil {
		return err
	}
	for _, m := range model.maps {
		if err := writeInt32(w, int32(m.ndim)); err != nil {
			return err
		}
		if err := writeInt32(w, int32(m.nbin)); err != nil {
			return err
		}
		if err := writeInt32(w, int32(m.kdim)); err != nil {
			return err
		}
		if err := writeInt32(w, int32(m.size)); err != nil {
			return err
		}
		if err := writeFloats(w, m.w); err != nil {
			return err
		}
		if err := writeFloats(w, m.c); err != nil {
			return err
		}
	}
	return nil
}

func readMapFile(prefix string) (*somModel, error) {
	if prefix == "" {
		return nil, fmt.Errorf("bcftools som: --prefix is required to read the map")
	}
	f, err := os.Open(prefix + ".som")
	if err != nil {
		return nil, fmt.Errorf("bcftools som: open %s.som: %w", prefix, err)
	}
	defer f.Close()
	model, err := readMaps(f)
	if err != nil {
		return nil, fmt.Errorf("bcftools som: parse %s.som: %w", prefix, err)
	}
	return model, nil
}

// readMaps deserialises a model written by writeMaps. It validates the
// magic and every length so a corrupt file is reported rather than read
// past EOF.
func readMaps(r io.Reader) (*somModel, error) {
	magic := make([]byte, len(somMagic))
	if err := readFull(r, magic); err != nil {
		return nil, err
	}
	if string(magic) != somMagic {
		return nil, fmt.Errorf("bad magic %q (want %q)", magic, somMagic)
	}
	bmuTh, err := readFloat(r)
	if err != nil {
		return nil, err
	}
	nAnnot, err := readInt32(r)
	if err != nil {
		return nil, err
	}
	if nAnnot < 0 || nAnnot > 1<<20 {
		return nil, fmt.Errorf("implausible annotation count %d", nAnnot)
	}
	model := &somModel{
		bmuTh:  bmuTh,
		annots: make([]string, nAnnot),
		min:    make([]float64, nAnnot),
		max:    make([]float64, nAnnot),
	}
	for i := 0; i < int(nAnnot); i++ {
		nameLen, err := readInt32(r)
		if err != nil {
			return nil, err
		}
		if nameLen < 0 || nameLen > 1<<16 {
			return nil, fmt.Errorf("implausible annotation name length %d", nameLen)
		}
		name := make([]byte, nameLen)
		if err := readFull(r, name); err != nil {
			return nil, err
		}
		model.annots[i] = string(name)
		if model.min[i], err = readFloat(r); err != nil {
			return nil, err
		}
		if model.max[i], err = readFloat(r); err != nil {
			return nil, err
		}
	}
	nMap, err := readInt32(r)
	if err != nil {
		return nil, err
	}
	if nMap < 0 || nMap > 1<<16 {
		return nil, fmt.Errorf("implausible map count %d", nMap)
	}
	model.maps = make([]*som, nMap)
	for i := 0; i < int(nMap); i++ {
		m := &som{}
		var v int32
		if v, err = readInt32(r); err != nil {
			return nil, err
		}
		m.ndim = int(v)
		if v, err = readInt32(r); err != nil {
			return nil, err
		}
		m.nbin = int(v)
		if v, err = readInt32(r); err != nil {
			return nil, err
		}
		m.kdim = int(v)
		if v, err = readInt32(r); err != nil {
			return nil, err
		}
		m.size = int(v)
		if m.ndim < 1 || m.nbin < 1 || m.kdim < 1 || m.size < 1 {
			return nil, fmt.Errorf("implausible map dimensions ndim=%d nbin=%d kdim=%d size=%d", m.ndim, m.nbin, m.kdim, m.size)
		}
		m.w = make([]float64, m.size*m.kdim)
		if err := readFloats(r, m.w); err != nil {
			return nil, err
		}
		m.c = make([]float64, m.size)
		if err := readFloats(r, m.c); err != nil {
			return nil, err
		}
		m.div = make([]float64, m.ndim)
		for j := 0; j < m.ndim; j++ {
			m.div[j] = math.Pow(float64(m.nbin), float64(m.ndim-j-1))
		}
		model.maps[i] = m
	}
	return model, nil
}

// --- little-endian binary helpers ----------------------------------------

func writeBytes(w io.Writer, b []byte) error {
	n, err := w.Write(b)
	if err != nil {
		return err
	}
	if n != len(b) {
		return io.ErrShortWrite
	}
	return nil
}

func writeInt32(w io.Writer, v int32) error {
	var buf [4]byte
	binary.LittleEndian.PutUint32(buf[:], uint32(v))
	return writeBytes(w, buf[:])
}

func writeFloat(w io.Writer, v float64) error {
	var buf [8]byte
	binary.LittleEndian.PutUint64(buf[:], math.Float64bits(v))
	return writeBytes(w, buf[:])
}

func writeFloats(w io.Writer, vs []float64) error {
	buf := make([]byte, 8*len(vs))
	for i, v := range vs {
		binary.LittleEndian.PutUint64(buf[i*8:], math.Float64bits(v))
	}
	return writeBytes(w, buf)
}

func readFull(r io.Reader, b []byte) error {
	_, err := io.ReadFull(r, b)
	return err
}

func readInt32(r io.Reader) (int32, error) {
	var buf [4]byte
	if err := readFull(r, buf[:]); err != nil {
		return 0, err
	}
	return int32(binary.LittleEndian.Uint32(buf[:])), nil
}

func readFloat(r io.Reader) (float64, error) {
	var buf [8]byte
	if err := readFull(r, buf[:]); err != nil {
		return 0, err
	}
	return math.Float64frombits(binary.LittleEndian.Uint64(buf[:])), nil
}

func readFloats(r io.Reader, vs []float64) error {
	buf := make([]byte, 8*len(vs))
	if err := readFull(r, buf); err != nil {
		return err
	}
	for i := range vs {
		vs[i] = math.Float64frombits(binary.LittleEndian.Uint64(buf[i*8:]))
	}
	return nil
}
