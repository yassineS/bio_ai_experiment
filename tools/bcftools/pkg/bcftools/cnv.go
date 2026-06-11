// bcftools cnv — copy-number variation calling.
//
// This is a faithful Go port of upstream's vcfcnv.c. The caller is an
// HMM over a per-contig sweep with hidden copy-number states; emission
// probabilities are derived from each site's B-allele frequency (BAF)
// and log-R ratio (LRR) distributions. The generic HMM engine
// (Viterbi, forward-backward, Baum-Welch) is the shared hmm.go port,
// reused unchanged from bcftools roh.
//
// State set (upstream's N_STATES == 4):
//
//	CN0 — complete loss   (homozygous deletion)
//	CN1 — single-copy loss
//	CN2 — normal diploid
//	CN3 — single-copy gain
//
// In the single-sample case the HMM has 4 states. With a paired
// control sample (-c) the HMM has 4*4 == 16 states, the outer product
// of the query and control copy-number, with a transition matrix that
// favours both samples sharing the same state (the --same-prob prior).
//
// The emission model for each site combines a BAF term — a mixture of
// Gaussian peaks at the genotype-cluster B-allele frequencies, weighted
// by the population genotype frequencies fRR/fRA/fAA — and an LRR term,
// a Gaussian centred on the per-state expected log-R ratio. Their
// relative contributions are the --BAF-weight and --LRR-weight knobs.
//
// Output is a tab-delimited region summary mirroring upstream's
// summary.tab "RG" rows: one record per maximal run of constant copy
// number along a contig, with a phred-scaled quality from the
// forward-backward posterior.
//
// Upstream reference: reference_code/bcftools/vcfcnv.c and the shared
// reference_code/bcftools/HMM.c.

package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Copy-number state indices, matching vcfcnv.c's CN0..CN3 constants.
const (
	cnvCN0     = 0
	cnvCN1     = 1
	cnvCN2     = 2
	cnvCN3     = 3
	cnvNStates = 4
)

// cnvMissingBAF is upstream's arbitrary negative sentinel for a missing
// BAF measurement (parse_lrr_baf sets baf = -0.1).
const cnvMissingBAF = -0.1

// CNVOptions controls the behaviour of CNVFile / CNV. Field defaults
// follow upstream vcfcnv.c when zero values would be inappropriate;
// callers should set the documented defaults explicitly.
type CNVOptions struct {
	// QuerySample is upstream's -s/--query-sample. When empty and the
	// input is single-sample, that sample is used; a multi-sample
	// input without -s is an error.
	QuerySample string
	// ControlSample is upstream's -c/--control-sample. When set, the
	// HMM runs in paired mode (16 states).
	ControlSample string
	// OutputDir is upstream's -o/--output-dir. Recorded for parity;
	// this port streams the region summary to the caller's io.Writer.
	OutputDir string
	// PlotThreshold is upstream's -p/--plot-threshold. Recorded; this
	// port emits no plots.
	PlotThreshold float64

	// AberrantQuery, AberrantControl are upstream's -a/--aberrant
	// FLOAT[,FLOAT]: the fraction of aberrant cells in the query and
	// control mixtures. They shift the CN3 BAF peak means. Default 1.0.
	AberrantQuery   float64
	AberrantControl float64
	// BAFDev / LRRDev are upstream's -d/--BAF-dev and -k/--LRR-dev:
	// the expected per-sample BAF and LRR standard deviations. Upstream
	// defaults are 0.04 and 0.20. The control values default to the
	// query values.
	BAFDev        float64
	BAFDevControl float64
	LRRDev        float64
	LRRDevControl float64
	// LRRWeight is upstream's -l/--LRR-weight (lrr_bias). When zero the
	// LRR term is dropped entirely. Default 0.2.
	LRRWeight float64
	// BAFWeight is upstream's -b/--BAF-weight (baf_bias). Default 1.0.
	BAFWeight float64
	// LRRSmoothWin is upstream's -L/--LRR-smooth-win: the window of the
	// LRR moving-average smoother. Default 10.
	LRRSmoothWin int
	// ErrProb is upstream's -e/--err-prob: a uniform error floor added
	// to every non-CN0 emission. Default 1e-4.
	ErrProb float64
	// XYProb is upstream's -x/--xy-prob: the P(x|y) transition
	// probability to a different state. Default 1e-9.
	XYProb float64
	// SameProb is upstream's -P/--same-prob: the prior that the query
	// and control samples share the same state (paired mode only).
	// Default 0.5.
	SameProb float64
	// Optimize is upstream's -O/--optimize: when > 0 and < 1, the
	// per-contig fraction of aberrant cells in CN3 is estimated by
	// iterated forward-backward down to this floor.
	Optimize float64
	// BaumWelch is upstream's -W/--baum-welch convergence threshold.
	// When non-zero, the transition matrix is re-estimated per contig.
	BaumWelch float64

	// AFFile is upstream's --AF-file (CHR<TAB>POS<TAB>REF,ALT<TAB>AF). It
	// acts as a targets filter (sites whose position is absent are
	// skipped) and supplies the per-site non-reference allele frequency
	// used to recompute the genotype frequencies fRR/fRA/fAA under
	// Hardy-Weinberg. When empty the port uses the fixed default genotype
	// frequencies.
	AFFile string

	// Regions / Targets / RegionsFile / TargetsFile are post-filters.
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string
}

func (o *CNVOptions) bafWeight() float64 {
	if o.BAFWeight == 0 {
		return 1.0
	}
	return o.BAFWeight
}

func (o *CNVOptions) lrrWeight() float64 {
	if o.LRRWeight < 0 {
		return 0
	}
	if o.LRRWeight == 0 {
		return 0.2
	}
	return o.LRRWeight
}

func (o *CNVOptions) errProb() float64 {
	if o.ErrProb == 0 {
		return 1e-4
	}
	return o.ErrProb
}

func (o *CNVOptions) xyProb() float64 {
	if o.XYProb == 0 {
		return 1e-9
	}
	return o.XYProb
}

func (o *CNVOptions) sameProb() float64 {
	if o.SameProb == 0 {
		return 0.5
	}
	return o.SameProb
}

func (o *CNVOptions) lrrSmoothWin() int {
	if o.LRRSmoothWin == 0 {
		return 10
	}
	return o.LRRSmoothWin
}

// queryBafDev2 returns the squared default BAF deviation for the query
// sample (upstream's baf_dev2_dflt).
func (o *CNVOptions) queryBafDev2() float64 {
	d := o.BAFDev
	if d <= 0 {
		d = 0.04
	}
	return d * d
}

func (o *CNVOptions) controlBafDev2() float64 {
	d := o.BAFDevControl
	if d <= 0 {
		d = o.BAFDev
	}
	if d <= 0 {
		d = 0.04
	}
	return d * d
}

func (o *CNVOptions) queryLrrDev2() float64 {
	d := o.LRRDev
	if d <= 0 {
		d = 0.20
	}
	return d * d
}

func (o *CNVOptions) controlLrrDev2() float64 {
	d := o.LRRDevControl
	if d <= 0 {
		d = o.LRRDev
	}
	if d <= 0 {
		d = 0.20
	}
	return d * d
}

func (o *CNVOptions) queryCellFrac() float64 {
	if o.AberrantQuery <= 0 {
		return 1.0
	}
	return o.AberrantQuery
}

func (o *CNVOptions) controlCellFrac() float64 {
	if o.AberrantControl <= 0 {
		return 1.0
	}
	return o.AberrantControl
}

// CNVRow is one HMM-called copy-number region, mirroring an upstream
// summary.tab "RG" row. Start and End are 1-based inclusive contig
// coordinates. CNCall is "CN0".."CN3" for the query sample;
// ControlCNCall is set only in paired mode.
type CNVRow struct {
	Sample        string
	Chrom         string
	Start         int
	End           int
	NSites        int
	NHets         int
	Quality       float64
	CNCall        string
	ControlCNCall string
}

// CNVFile streams the VCF/BCF at path and writes the HMM copy-number
// region summary. It returns the number of region rows written.
func CNVFile(path string, w io.Writer, opts CNVOptions) (int, error) {
	reader, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("open %q: %w", path, err)
	}
	defer reader.Close()
	return CNV(reader, w, opts)
}

// cnvContig accumulates one contig's per-site BAF/LRR observations for
// the query and (optionally) control samples.
type cnvContig struct {
	chrom    string
	pos      []uint32
	queryBAF []float64
	queryLRR []float64
	ctrlBAF  []float64
	ctrlLRR  []float64
	// nonrefAF holds the per-site non-reference allele frequency from
	// --AF-file (vcfcnv.c's nonref_afs). It is nil when no AF-file is in
	// use, in which case the fixed default genotype frequencies apply.
	nonrefAF []float64
}

// CNV reads VCF records from r, runs the copy-number HMM per contig,
// and writes the region summary TSV to w. It returns the number of
// region rows written.
func CNV(r io.Reader, w io.Writer, opts CNVOptions) (int, error) {
	// --AF-file (vcfcnv.c). It acts as a targets filter (sites whose
	// CHROM:POS is absent from the file are skipped) AND supplies the
	// per-site non-reference allele frequency that drives the genotype
	// frequencies fRR/fRA/fAA. When the site's POS is listed but its
	// REF/ALT alleles do not match, upstream falls back to
	// nonref_af_dflt = 0.1 (vcfcnv.c:1194,:1257).
	var afFile *cnvAFFile
	if opts.AFFile != "" {
		af, ferr := loadCNVAFFile(opts.AFFile)
		if ferr != nil {
			return 0, ferr
		}
		afFile = af
	}

	vr := vcf.NewReader(r)
	hdr, err := vr.ReadHeader()
	if err != nil {
		return 0, fmt.Errorf("read header: %w", err)
	}

	query, control, err := cnvResolveSamples(hdr.Samples, opts)
	if err != nil {
		return 0, err
	}
	paired := control != ""

	regionSpecs := append([]string(nil), opts.Regions...)
	regionSpecs = append(regionSpecs, opts.Targets...)

	contigs := []*cnvContig{}
	var cur *cnvContig
	byChrom := map[string]*cnvContig{}

	for {
		v, err := vr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("read variant: %w", err)
		}
		if len(regionSpecs) > 0 && !regionMatches(v, regionSpecs) {
			continue
		}

		// --AF-file targets filter: skip sites whose position is not
		// listed (vcfcnv.c uses the AF-file as the SnpSift targets index).
		var siteAF float64
		if afFile != nil {
			af, listed := afFile.lookup(v)
			if !listed {
				continue
			}
			siteAF = af
		}

		qb, qBAFok := cnvSampleBAF(v, query)
		ql, qLRRok := cnvSampleLRR(v, query)
		if !qBAFok || (opts.lrrWeight() > 0 && !qLRRok) {
			// Upstream's parse_lrr_baf: a missing LRR (when LRR is
			// used) also voids the BAF; the site becomes a no-call.
			qb = cnvMissingBAF
			ql = 0
			qBAFok = false
		}
		var cb, cl float64
		var cBAFok bool
		if paired {
			var cLRRok bool
			cb, cBAFok = cnvSampleBAF(v, control)
			cl, cLRRok = cnvSampleLRR(v, control)
			if !cBAFok || (opts.lrrWeight() > 0 && !cLRRok) {
				cb = cnvMissingBAF
				cl = 0
				cBAFok = false
			}
			// Upstream skips the record only when neither sample has
			// a usable BAF.
			if !qBAFok && !cBAFok {
				continue
			}
		} else if !qBAFok {
			continue
		}

		if cur == nil || cur.chrom != v.Chrom {
			c, ok := byChrom[v.Chrom]
			if !ok {
				c = &cnvContig{chrom: v.Chrom}
				byChrom[v.Chrom] = c
				contigs = append(contigs, c)
			}
			cur = c
		}
		cur.pos = append(cur.pos, uint32(v.Pos-1))
		cur.queryBAF = append(cur.queryBAF, qb)
		cur.queryLRR = append(cur.queryLRR, ql)
		if paired {
			cur.ctrlBAF = append(cur.ctrlBAF, cb)
			cur.ctrlLRR = append(cur.ctrlLRR, cl)
		}
		if afFile != nil {
			cur.nonrefAF = append(cur.nonrefAF, siteAF)
		}
	}

	bw := newCNVTSVWriter(w, paired, query, control)
	if err := bw.WriteHeader(); err != nil {
		return 0, err
	}
	written := 0
	for _, c := range contigs {
		if len(c.pos) == 0 {
			continue
		}
		rows, err := cnvCallContig(c, query, control, opts)
		if err != nil {
			return written, err
		}
		for _, row := range rows {
			if err := bw.WriteRow(row); err != nil {
				return written, err
			}
			written++
		}
	}
	return written, nil
}

// cnvAFRecord is one CHR<TAB>POS<TAB>REF,ALT[,ALT2...]<TAB>AF entry from
// an --AF-file. ref is the first (reference) allele and alt is the
// comma-joined list of the remaining (alternate) alleles, so a match is
// tested against the record's full allele vector exactly as upstream's
// read_AF does. af is the alternate-allele frequency.
type cnvAFRecord struct {
	ref string
	alt string
	af  float64
}

// cnvAFFile indexes an --AF-file by "chrom\x00pos". A position present
// with no allele-match still counts as "listed" (so the targets filter
// keeps it), but uses the default non-reference AF.
type cnvAFFile struct {
	byPos map[string][]cnvAFRecord
}

// cnvNonrefAFDflt is upstream's nonref_af_dflt (vcfcnv.c:1257): the AF
// used when a site's position is listed but its alleles do not match.
const cnvNonrefAFDflt = 0.1

// loadCNVAFFile parses an --AF-file into a cnvAFFile. The format mirrors
// roh's: CHR<TAB>POS<TAB>REF,ALT<TAB>AF, one record per line, '#'
// comments skipped.
func loadCNVAFFile(path string) (*cnvAFFile, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("bcftools cnv: --AF-file: %w", err)
	}
	defer in.Close()
	out := &cnvAFFile{byPos: map[string][]cnvAFRecord{}}
	sc := bufio.NewScanner(in)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 4 {
			continue
		}
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			continue
		}
		// Column 3 is REF,ALT[,ALT2...]; split off the REF and keep the
		// remaining ALTs comma-joined so multiallelic entries can match
		// the record's full allele vector.
		alleles := strings.SplitN(fields[2], ",", 2)
		if len(alleles) != 2 {
			return nil, fmt.Errorf("bcftools cnv: --AF-file: expected two comma-separated alleles (REF,ALT) in column 3: %q", line)
		}
		af, err := strconv.ParseFloat(fields[3], 64)
		if err != nil {
			// A missing "." AF is allowed; treat it like a no-match (the
			// default is applied at lookup time).
			af = -1
		}
		key := cnvAFKey(fields[0], pos)
		out.byPos[key] = append(out.byPos[key], cnvAFRecord{ref: alleles[0], alt: alleles[1], af: af})
	}
	return out, sc.Err()
}

// lookup mirrors vcfcnv.c's read_AF + targets index: the second return
// value reports whether the record's position is listed at all (the
// targets filter), and the first is the non-reference AF to use — the
// matched value when the record's full allele vector (REF + all ALTs)
// equals a file entry, else the default cnvNonrefAFDflt.
func (a *cnvAFFile) lookup(v *vcf.Variant) (float64, bool) {
	recs, ok := a.byPos[cnvAFKey(v.Chrom, v.Pos)]
	if !ok {
		return 0, false
	}
	alt := strings.Join(v.Alt, ",")
	for _, r := range recs {
		if r.ref == v.Ref && r.alt == alt && r.af >= 0 {
			return r.af, true
		}
	}
	return cnvNonrefAFDflt, true
}

// cnvAFKey builds the chrom+position index key for cnvAFFile.
func cnvAFKey(chrom string, pos int) string {
	return chrom + "\x00" + strconv.Itoa(pos)
}

// cnvResolveSamples picks the query and control sample names, applying
// upstream's rule that a multi-sample input requires -s.
func cnvResolveSamples(all []string, opts CNVOptions) (query, control string, err error) {
	has := func(name string) bool {
		for _, s := range all {
			if s == name {
				return true
			}
		}
		return false
	}
	query = opts.QuerySample
	if query == "" {
		if len(all) == 0 {
			return "", "", fmt.Errorf("input has no samples")
		}
		if len(all) > 1 {
			return "", "", fmt.Errorf("multi-sample VCF, missing the -s/--query-sample option")
		}
		query = all[0]
	} else if !has(query) {
		return "", "", fmt.Errorf("the query sample %q was not found", query)
	}
	control = opts.ControlSample
	if control != "" && !has(control) {
		return "", "", fmt.Errorf("the control sample %q was not found", control)
	}
	return query, control, nil
}

// cnvSampleBAF returns the FORMAT/BAF for the named sample, falling
// back to a BAF synthesised from FORMAT/AD = REF,ALT.
func cnvSampleBAF(v *vcf.Variant, name string) (float64, bool) {
	s := findSample(v, name)
	if s == nil {
		return 0, false
	}
	if baf, ok := parseFloatField(s.Data, "BAF"); ok {
		return baf, true
	}
	if ad, ok := s.Data["AD"]; ok {
		if baf, bok := bafFromAD(ad); bok {
			return baf, true
		}
	}
	return 0, false
}

// cnvSampleLRR returns the FORMAT/LRR for the named sample.
func cnvSampleLRR(v *vcf.Variant, name string) (float64, bool) {
	s := findSample(v, name)
	if s == nil {
		return 0, false
	}
	return parseFloatField(s.Data, "LRR")
}

// findSample returns the per-sample FORMAT slot for the named sample,
// or nil when the sample is not present in the record.
func findSample(v *vcf.Variant, name string) *vcf.Sample {
	for i := range v.Samples {
		if v.Samples[i].Name == name {
			return &v.Samples[i]
		}
	}
	return nil
}

// bafFromAD synthesises BAF = ALT / (REF + ALT) from a FORMAT/AD string
// of the form "ref,alt". Used as a fallback when FORMAT/BAF is not
// emitted by the upstream pipeline.
func bafFromAD(raw string) (float64, bool) {
	if raw == "" || raw == "." {
		return 0, false
	}
	parts := strings.Split(raw, ",")
	if len(parts) < 2 {
		return 0, false
	}
	refN, err1 := strconv.ParseFloat(strings.TrimSpace(parts[0]), 64)
	altN, err2 := strconv.ParseFloat(strings.TrimSpace(parts[1]), 64)
	if err1 != nil || err2 != nil {
		return 0, false
	}
	if (refN + altN) <= 0 {
		return 0, false
	}
	return altN / (refN + altN), true
}

// parseFloatField extracts a FORMAT field by tag, returning (value, ok).
// Empty / "." values count as missing.
func parseFloatField(data map[string]string, tag string) (float64, bool) {
	s, ok := data[tag]
	if !ok || s == "" || s == "." {
		return 0, false
	}
	if i := strings.IndexByte(s, ','); i >= 0 {
		s = s[:i]
	}
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, false
	}
	return v, true
}

// cnvTSVWriter emits the HMM region summary.
type cnvTSVWriter struct {
	w       io.Writer
	paired  bool
	query   string
	control string
}

func newCNVTSVWriter(w io.Writer, paired bool, query, control string) *cnvTSVWriter {
	return &cnvTSVWriter{w: w, paired: paired, query: query, control: control}
}

// WriteHeader emits the column row, mirroring upstream's summary.tab
// "RG" header line.
func (w *cnvTSVWriter) WriteHeader() error {
	var s string
	if w.paired {
		s = fmt.Sprintf("# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:%s\t[6]Copy number:%s\t[7]Quality\t[8]nSites\t[9]nHETs\n",
			w.query, w.control)
	} else {
		s = fmt.Sprintf("# RG, Regions\t[2]Chromosome\t[3]Start\t[4]End\t[5]Copy number:%s\t[6]Quality\t[7]nSites\t[8]nHETs\n",
			w.query)
	}
	_, err := io.WriteString(w.w, s)
	return err
}

// WriteRow emits one "RG" region row.
func (w *cnvTSVWriter) WriteRow(r CNVRow) error {
	var err error
	if w.paired {
		_, err = fmt.Fprintf(w.w, "RG\t%s\t%d\t%d\t%s\t%s\t%.1f\t%d\t%d\n",
			r.Chrom, r.Start, r.End, r.CNCall, r.ControlCNCall, r.Quality, r.NSites, r.NHets)
	} else {
		_, err = fmt.Fprintf(w.w, "RG\t%s\t%d\t%d\t%s\t%.1f\t%d\t%d\n",
			r.Chrom, r.Start, r.End, r.CNCall, r.Quality, r.NSites, r.NHets)
	}
	return err
}

// --- HMM machinery (port of vcfcnv.c) -------------------------------

// cnvGaussParam is a Gaussian BAF peak: mean, squared deviation, and a
// normalisation factor from the truncated [0,1] support. It mirrors
// vcfcnv.c's gauss_param_t.
type cnvGaussParam struct {
	mean float64
	dev2 float64
	norm float64
}

// cnvNormProb evaluates a single truncated-Gaussian BAF peak, mirroring
// vcfcnv.c's norm_prob.
func cnvNormProb(baf float64, p cnvGaussParam) float64 {
	return math.Exp(-(baf-p.mean)*(baf-p.mean)*0.5/p.dev2) / p.norm / math.Sqrt(2*math.Pi*p.dev2)
}

// cnvNormCDF returns the probability mass of an N(mean,dev) Gaussian
// inside the unit interval [0,1], mirroring vcfcnv.c's norm_cdf.
func cnvNormCDF(mean, dev float64) float64 {
	top := 1 - 0.5*math.Erfc((1-mean)/(dev*math.Sqrt2))
	bot := 1 - 0.5*math.Erfc((0-mean)/(dev*math.Sqrt2))
	return top - bot
}

// Gaussian-peak slots within a sample's parameter table, mirroring the
// GAUSS_* macros in vcfcnv.c.
const (
	cnvPkCN1R   = 0
	cnvPkCN1A   = 1
	cnvPkCN2RR  = 2
	cnvPkCN2RA  = 3
	cnvPkCN2AA  = 4
	cnvPkCN3RRR = 5
	cnvPkCN3RRA = 6
	cnvPkCN3RAA = 7
	cnvPkCN3AAA = 8
	cnvNGauss   = 18
)

// cnvSample holds one sample's per-contig observations and the derived
// emission parameters, mirroring vcfcnv.c's sample_t.
type cnvSample struct {
	baf       []float64
	lrr       []float64
	bafDev2   float64 // current (possibly --optimize-updated) BAF dev^2
	bafDev2D  float64 // default BAF dev^2
	lrrDev2   float64
	cellFrac  float64 // current CN3 aberrant-cell fraction
	cellFracD float64 // default CN3 aberrant-cell fraction
	gauss     [cnvNGauss]cnvGaussParam
	pobs      [cnvNStates]float64
}

// setGaussParams fills a sample's BAF peak table, mirroring vcfcnv.c's
// set_gauss_params: peak means depend on the aberrant-cell fraction
// for CN3, and the truncated-Gaussian norm depends on the BAF dev.
func (s *cnvSample) setGaussParams() {
	for i := range s.gauss {
		s.gauss[i].dev2 = s.bafDev2
	}
	dev := math.Sqrt(s.bafDev2)

	s.gauss[cnvPkCN1R].mean = 0
	s.gauss[cnvPkCN1A].mean = 1
	s.gauss[cnvPkCN1R].norm = cnvNormCDF(s.gauss[cnvPkCN1R].mean, dev)
	s.gauss[cnvPkCN1A].norm = cnvNormCDF(s.gauss[cnvPkCN1A].mean, dev)

	s.gauss[cnvPkCN2RR].mean = 0
	s.gauss[cnvPkCN2RA].mean = 0.5
	s.gauss[cnvPkCN2AA].mean = 1
	s.gauss[cnvPkCN2RR].norm = cnvNormCDF(s.gauss[cnvPkCN2RR].mean, dev)
	s.gauss[cnvPkCN2RA].norm = cnvNormCDF(s.gauss[cnvPkCN2RA].mean, dev)
	s.gauss[cnvPkCN2AA].norm = cnvNormCDF(s.gauss[cnvPkCN2AA].mean, dev)

	s.gauss[cnvPkCN3RRR].mean = 0
	s.gauss[cnvPkCN3RRA].mean = 1.0 / (2 + s.cellFrac)
	s.gauss[cnvPkCN3RAA].mean = (1.0 + s.cellFrac) / (2 + s.cellFrac)
	s.gauss[cnvPkCN3AAA].mean = 1
	s.gauss[cnvPkCN3RRR].norm = cnvNormCDF(s.gauss[cnvPkCN3RRR].mean, dev)
	s.gauss[cnvPkCN3RRA].norm = cnvNormCDF(s.gauss[cnvPkCN3RRA].mean, dev)
	s.gauss[cnvPkCN3RAA].norm = cnvNormCDF(s.gauss[cnvPkCN3RAA].mean, dev)
	s.gauss[cnvPkCN3AAA].norm = cnvNormCDF(s.gauss[cnvPkCN3AAA].mean, dev)
}

// setObservedProb fills s.pobs with the per-state observation
// probability for site isite, mirroring vcfcnv.c's set_observed_prob.
// fRR/fRA/fAA are the population genotype frequencies; bafBias/lrrBias
// are the --BAF-weight / --LRR-weight knobs and errProb is --err-prob.
func (s *cnvSample) setObservedProb(isite int, fRR, fRA, fAA, bafBias, lrrBias, errProb float64) {
	baf := s.baf[isite]
	lrr := 0.0
	if lrrBias > 0 {
		lrr = s.lrr[isite]
	}

	if baf < 0 {
		// no call: either a technical issue or the call could not be
		// made because it is CN0.
		s.pobs[cnvCN0] = 0.5
		rest := (1.0 - s.pobs[cnvCN0]) / (cnvNStates - 1)
		for i := 1; i < cnvNStates; i++ {
			s.pobs[i] = rest
		}
		return
	}

	cn1Baf := cnvNormProb(baf, s.gauss[cnvPkCN1R])*(fRR+fRA*0.5) +
		cnvNormProb(baf, s.gauss[cnvPkCN1A])*(fAA+fRA*0.5)
	cn2Baf := cnvNormProb(baf, s.gauss[cnvPkCN2RR])*fRR +
		cnvNormProb(baf, s.gauss[cnvPkCN2RA])*fRA +
		cnvNormProb(baf, s.gauss[cnvPkCN2AA])*fAA
	cn3Baf := cnvNormProb(baf, s.gauss[cnvPkCN3RRR])*fRR +
		cnvNormProb(baf, s.gauss[cnvPkCN3RRA])*fRA*0.5 +
		cnvNormProb(baf, s.gauss[cnvPkCN3RAA])*fRA*0.5 +
		cnvNormProb(baf, s.gauss[cnvPkCN3AAA])*fAA

	norm := cn1Baf + cn2Baf + cn3Baf
	cn1Baf /= norm
	cn2Baf /= norm
	cn3Baf /= norm

	cn1Lrr := math.Exp(-(lrr + 0.45) * (lrr + 0.45) / s.lrrDev2)
	cn2Lrr := math.Exp(-(lrr - 0.00) * (lrr - 0.00) / s.lrrDev2)
	cn3Lrr := math.Exp(-(lrr - 0.30) * (lrr - 0.30) / s.lrrDev2)

	s.pobs[cnvCN0] = 0
	s.pobs[cnvCN1] = errProb + (1-bafBias+bafBias*cn1Baf)*(1-lrrBias+lrrBias*cn1Lrr)
	s.pobs[cnvCN2] = errProb + (1-bafBias+bafBias*cn2Baf)*(1-lrrBias+lrrBias*cn2Lrr)
	s.pobs[cnvCN3] = errProb + (1-bafBias+bafBias*cn3Baf)*(1-lrrBias+lrrBias*cn3Lrr)
}

// cnvHmm2cnState splits a paired-mode HMM state index into the query
// and control copy-number indices, mirroring vcfcnv.c's hmm2cn_state.
func cnvHmm2cnState(i int) (a, b int) {
	a = i / cnvNStates
	b = i - a*cnvNStates
	return a, b
}

// cnvInitTprob builds the HMM transition matrix, mirroring vcfcnv.c's
// init_tprob_matrix. ndim is 4 (single sample) or 16 (paired).
func cnvInitTprob(ndim int, ijProb, sameProb float64) ([]float64, error) {
	mat := make([]float64, ndim*ndim)
	if ndim == cnvNStates {
		pii := 1 - ijProb*(cnvNStates-1)
		if pii < ijProb {
			return nil, fmt.Errorf("-x set too high, P(x|x) < P(x|y): %e vs %e", pii, ijProb)
		}
		for j := 0; j < ndim; j++ {
			for i := 0; i < ndim; i++ {
				if i == j {
					mat[i*ndim+j] = pii
				} else {
					mat[i*ndim+j] = ijProb
				}
			}
		}
		return mat, nil
	}
	// Paired mode: interpret ij_prob as ii_prob so the behaviour is
	// closer to single-sample calling.
	pii := 1 - ijProb*(cnvNStates-1)
	ij := (1 - pii) / float64(ndim-1)
	for j := 0; j < ndim; j++ {
		ja, jb := cnvHmm2cnState(j)
		sum := 0.0
		for i := 0; i < ndim; i++ {
			ia, ib := cnvHmm2cnState(i)
			pa := ij
			if ja == ia {
				pa = pii
			}
			pb := ij
			if jb == ib {
				pb = pii
			}
			switch {
			case ia == ib && ja == jb:
				mat[i*ndim+j] = pa*pb - pa*pb*sameProb + math.Sqrt(pa*pb)*sameProb
			case ia == ib:
				mat[i*ndim+j] = pa * pb
			default:
				mat[i*ndim+j] = pa * pb * (1 - sameProb)
			}
			sum += mat[i*ndim+j]
		}
		for i := 0; i < ndim; i++ {
			mat[i*ndim+j] /= sum
		}
	}
	return mat, nil
}

// cnvInitIProbs builds the HMM initial-state probability vector,
// mirroring vcfcnv.c's init_iprobs.
func cnvInitIProbs(ndim int, sameProb float64) []float64 {
	probs := make([]float64, ndim)
	if ndim == cnvNStates {
		for i := 0; i < ndim; i++ {
			if i == cnvCN2 {
				probs[i] = 0.5
			} else {
				probs[i] = 0.5 / 3
			}
		}
		return probs
	}
	norm := 0.0
	for i := 0; i < ndim; i++ {
		ia, ib := cnvHmm2cnState(i)
		pa := 0.5 / 3
		if ia == cnvCN2 {
			pa = 0.5
		}
		pb := 0.5 / 3
		if ib == cnvCN2 {
			pb = 0.5
		}
		probs[i] = pa * pb
		if ia != ib {
			probs[i] *= 1 - sameProb
		}
		norm += probs[i]
	}
	for i := 0; i < ndim; i++ {
		probs[i] /= norm
	}
	return probs
}

// cnvSmoothData applies a centred moving-average smoother of window
// win to dat in place, mirroring vcfcnv.c's smooth_data.
func cnvSmoothData(dat []float64, win int) {
	if win <= 1 {
		return
	}
	ndat := len(dat)
	if ndat == 0 {
		return
	}
	k1 := win / 2
	k2 := win - k1
	buf := make([]float64, 0, win)
	sum := 0.0
	for i := 0; i < k2 && i < ndat; i++ {
		sum += dat[i]
		buf = append(buf, dat[i])
	}
	head := 0
	for i := 0; i < ndat; i++ {
		dat[i] = sum / float64(len(buf)-head)
		if i >= k1 {
			sum -= buf[head]
			head++
		}
		if i+k2 < ndat {
			sum += dat[i+k2]
			buf = append(buf, dat[i+k2])
		}
	}
}

// cnvCopyNumberState maps an HMM state index to its CN0..CN3 string for
// either the query (ismpl==0) or control (ismpl==1) sample.
func cnvCopyNumberState(paired bool, istate, ismpl int) string {
	code := []string{"CN0", "CN1", "CN2", "CN3", "CN4"}
	if !paired {
		return code[istate]
	}
	a, b := cnvHmm2cnState(istate)
	if ismpl == 1 {
		return code[b]
	}
	return code[a]
}

// bafLikelyHet reports whether a BAF value lies in the heterozygous
// band, mirroring vcfcnv.c's BAF_LIKELY_HET macro.
func bafLikelyHet(v float64) bool { return v > 0.25 && v < 0.75 }

// cnvEmissionProbs computes the [nsites*nstates] emission-probability
// array for one contig, mirroring vcfcnv.c's set_emission_probs. It
// reuses the per-sample observation probabilities and forms their
// outer product in paired mode. nonrefAF, when non-nil, supplies the
// per-site non-reference allele frequency from --AF-file; the genotype
// frequencies fRR/fRA/fAA are then recomputed per site under
// Hardy-Weinberg (vcfcnv.c:735-739) instead of using the fixed defaults.
func cnvEmissionProbs(query, control *cnvSample, nsites, nstates int, paired bool, opts CNVOptions, nonrefAF []float64) []float64 {
	// Fixed default genotype frequencies, used when no AF-file is given.
	fRR, fRA, fAA := 0.76, 0.14, 0.098
	bafBias := opts.bafWeight()
	lrrBias := opts.lrrWeight()
	errProb := opts.errProb()

	eprob := make([]float64, nsites*nstates)
	for i := 0; i < nsites; i++ {
		rr, ra, aa := fRR, fRA, fAA
		if nonrefAF != nil {
			// Hardy-Weinberg genotype frequencies from the site AF.
			af := nonrefAF[i]
			rr = (1 - af) * (1 - af)
			ra = 2 * af * (1 - af)
			aa = af * af
		}
		query.setObservedProb(i, rr, ra, aa, bafBias, lrrBias, errProb)
		dst := eprob[i*nstates : (i+1)*nstates]
		if paired {
			control.setObservedProb(i, rr, ra, aa, bafBias, lrrBias, errProb)
			for a := 0; a < cnvNStates; a++ {
				for b := 0; b < cnvNStates; b++ {
					dst[a*cnvNStates+b] = query.pobs[a] * control.pobs[b]
				}
			}
		} else {
			for a := 0; a < cnvNStates; a++ {
				dst[a] = query.pobs[a]
			}
		}
	}
	return eprob
}

// cnvUpdateSample re-estimates a sample's CN3 aberrant-cell fraction and
// BAF deviation from the current forward-backward posteriors, mirroring
// vcfcnv.c's update_sample_args. It returns true when the estimate has
// converged. ismpl is 0 for the query sample, 1 for the control.
func cnvUpdateSample(s *cnvSample, fwd []float64, nstates, nsites, ismpl int, paired bool, optimizeFrac float64) bool {
	probs := make([]float64, 0, nsites)
	for i := 0; i < nsites; i++ {
		baf := s.baf[i]
		if baf > 4.0/5 || baf < 1.0/5 {
			continue
		}
		row := fwd[i*nstates : (i+1)*nstates]
		probCN3 := 0.0
		switch {
		case !paired:
			probCN3 = row[cnvCN3]
		case ismpl == 0:
			for j := 0; j < cnvNStates; j++ {
				probCN3 += row[cnvCN3*cnvNStates+j]
			}
		default:
			for j := 0; j < cnvNStates; j++ {
				probCN3 += row[cnvCN3+j*cnvNStates]
			}
		}
		probs = append(probs, probCN3)
	}
	cnvSmoothData(probs, 50)

	meanCN3, normCN3 := 0.0, 0.0
	bafAADev2, normBafAADev2 := 0.0, 0.0
	k := 0
	for i := 0; i < nsites; i++ {
		baf := s.baf[i]
		if baf > 4.0/5 {
			bafAADev2 += (1.0 - baf) * (1.0 - baf)
			normBafAADev2++
			continue
		}
		if baf > 0.5 {
			baf = 1 - baf
		}
		if baf < 1.0/5 {
			continue
		}
		probCN3 := probs[k]
		k++
		meanCN3 += probCN3 * baf
		normCN3 += probCN3
	}
	if normCN3 == 0 {
		s.cellFrac = 1.0
		return true
	}
	meanCN3 /= normCN3

	bafDev2 := 0.0
	k = 0
	for i := 0; i < nsites; i++ {
		baf := s.baf[i]
		if baf > 4.0/5 {
			continue
		}
		if baf > 0.5 {
			baf = 1 - baf
		}
		if baf < 1.0/5 {
			continue
		}
		probCN3 := probs[k]
		k++
		bafDev2 += probCN3 * (baf - meanCN3) * (baf - meanCN3)
	}
	bafDev2 /= normCN3
	if normBafAADev2 > 0 {
		bafAADev2 /= normBafAADev2
	}
	if bafDev2 < bafAADev2 {
		bafDev2 = bafAADev2
	}
	maxMeanCN3 := 0.5 - math.Sqrt(bafDev2)*1.644854
	newFrac := 1.0/meanCN3 - 2
	if meanCN3 > maxMeanCN3 || newFrac < optimizeFrac {
		s.cellFrac = 1.0
		return true
	}
	if newFrac > 1 {
		newFrac = 1
	}
	converged := math.Abs(newFrac-s.cellFrac) < 1e-1
	if bafDev2 > 3*s.bafDev2D {
		bafDev2 = 3 * s.bafDev2D
	} else if bafDev2 < 0.5*s.bafDev2D {
		bafDev2 = 0.5 * s.bafDev2D
	}
	s.cellFrac = newFrac
	s.bafDev2 = bafDev2
	return converged
}

// cnvAvgIIProb returns the average diagonal (self-transition) entry of
// an n*n transition matrix, mirroring vcfcnv.c's avg_ii_prob.
func cnvAvgIIProb(n int, mat []float64) float64 {
	avg := 0.0
	for i := 0; i < n; i++ {
		avg += mat[i*n+i]
	}
	return avg / float64(n)
}

// cnvCallContig runs the HMM over one contig and returns its region
// summary rows, mirroring vcfcnv.c's cnv_flush_viterbi.
func cnvCallContig(c *cnvContig, query, control string, opts CNVOptions) ([]CNVRow, error) {
	paired := control != ""
	nsites := len(c.pos)
	if nsites == 0 {
		return nil, nil
	}
	nstates := cnvNStates
	if paired {
		nstates = cnvNStates * cnvNStates
	}

	qs := &cnvSample{
		baf:       c.queryBAF,
		lrr:       append([]float64(nil), c.queryLRR...),
		bafDev2D:  opts.queryBafDev2(),
		lrrDev2:   opts.queryLrrDev2(),
		cellFracD: opts.queryCellFrac(),
	}
	var cs *cnvSample
	if paired {
		cs = &cnvSample{
			baf:       c.ctrlBAF,
			lrr:       append([]float64(nil), c.ctrlLRR...),
			bafDev2D:  opts.controlBafDev2(),
			lrrDev2:   opts.controlLrrDev2(),
			cellFracD: opts.controlCellFrac(),
		}
	} else {
		cs = &cnvSample{}
	}

	// Smooth LRR to reduce noise (only when LRR contributes).
	if opts.lrrWeight() > 0 {
		cnvSmoothData(qs.lrr, opts.lrrSmoothWin())
		if paired {
			cnvSmoothData(cs.lrr, opts.lrrSmoothWin())
		}
	}

	// Reset per-contig estimates to their defaults.
	qs.cellFrac = qs.cellFracD
	qs.bafDev2 = qs.bafDev2D
	qs.setGaussParams()
	if paired {
		cs.cellFrac = cs.cellFracD
		cs.bafDev2 = cs.bafDev2D
		cs.setGaussParams()
	}

	tprob, err := cnvInitTprob(nstates, opts.xyProb(), opts.sameProb())
	if err != nil {
		// A misconfigured -x/-P leaves no valid transition matrix.
		// Upstream's init_tprob_matrix calls error() — a fatal exit;
		// fabricating plausible CN2 output for a user misconfiguration
		// would silently mask the mistake, so propagate the error.
		return nil, err
	}
	iprobs := cnvInitIProbs(nstates, opts.sameProb())
	h := newHMM(nstates, tprob, 10000)
	h.initStates(iprobs)

	// --optimize: iterate forward-backward to estimate the CN3
	// aberrant-cell fraction. Upstream caps the loop at 20 iterations.
	if opts.Optimize > 0 && opts.Optimize < 1 {
		niter := 0
		for {
			eprob := cnvEmissionProbs(qs, cs, nsites, nstates, paired, opts, c.nonrefAF)
			h.runFwdBwd(nsites, eprob, c.pos)
			done := cnvUpdateSample(qs, h.fwd, nstates, nsites, 0, paired, opts.Optimize)
			if paired {
				done2 := cnvUpdateSample(cs, h.fwd, nstates, nsites, 1, paired, opts.Optimize)
				done = done && done2
			}
			qs.setGaussParams()
			if paired {
				cs.setGaussParams()
			}
			niter++
			if done || niter >= 20 {
				if niter >= 20 {
					qs.cellFrac = qs.cellFracD
					qs.bafDev2 = qs.bafDev2D
					qs.setGaussParams()
					if paired {
						cs.cellFrac = cs.cellFracD
						cs.bafDev2 = cs.bafDev2D
						cs.setGaussParams()
					}
				}
				break
			}
		}
	}

	eprob := cnvEmissionProbs(qs, cs, nsites, nstates, paired, opts, c.nonrefAF)

	// --baum-welch: re-estimate the transition matrix until the mean
	// self-transition probability stabilises.
	if opts.BaumWelch != 0 {
		for iter := 0; iter < 100; iter++ {
			oriII := cnvAvgIIProb(nstates, h.tprobMatrix())
			bw := h.runBaumWelch(nsites, eprob, c.pos)
			newII := cnvAvgIIProb(nstates, bw)
			nt, err := cnvInitTprob(nstates, 1-newII, opts.sameProb())
			if err != nil {
				break
			}
			h.setTprobMatrix(nt, 10000)
			if math.Abs(newII-oriII) < opts.BaumWelch {
				break
			}
		}
	}

	h.runViterbi(nsites, eprob, c.pos)
	h.runFwdBwd(nsites, eprob, c.pos)

	return cnvSummarisePath(c, query, control, paired, nstates, h.vpath, h.fwd, qs, cs), nil
}

// cnvSummarisePath turns a decoded Viterbi path into region rows,
// mirroring the output loop of vcfcnv.c's cnv_flush_viterbi.
func cnvSummarisePath(c *cnvContig, query, control string, paired bool, nstates int, vpath []uint8, fwd []float64, qs, cs *cnvSample) []CNVRow {
	nsites := len(c.pos)
	var rows []CNVRow

	startCN := int(vpath[0])
	startPos := c.pos[0]
	istart := 0
	qual := 0.0
	smplNtot, smplNhet := 0, 0
	ctrlNtot, ctrlNhet := 0, 0

	emit := func(endPos uint32, endIdx int) {
		q := phredScore(1 - qual/float64(endIdx-istart))
		row := CNVRow{
			Sample:  query,
			Chrom:   c.chrom,
			Start:   int(startPos) + 1,
			End:     int(endPos),
			NSites:  smplNtot,
			NHets:   smplNhet,
			Quality: q,
			CNCall:  cnvCopyNumberState(paired, startCN, 0),
		}
		if paired {
			row.ControlCNCall = cnvCopyNumberState(paired, startCN, 1)
		}
		rows = append(rows, row)
	}

	var isite int
	for isite = 0; isite < nsites; isite++ {
		state := int(vpath[isite])
		row := fwd[isite*nstates : (isite+1)*nstates]
		qual += row[startCN]

		if qs.baf[isite] >= 0 {
			if bafLikelyHet(qs.baf[isite]) {
				smplNhet++
			}
			smplNtot++
		}
		if paired && cs.baf[isite] >= 0 {
			if bafLikelyHet(cs.baf[isite]) {
				ctrlNhet++
			}
			ctrlNtot++
		}

		if startCN != state {
			emit(c.pos[isite], isite)
			istart = isite
			startPos = c.pos[isite]
			startCN = state
			qual = 0
			smplNtot, smplNhet, ctrlNtot, ctrlNhet = 0, 0, 0, 0
		}
	}
	emit(c.pos[nsites-1]+1, isite)
	return rows
}
