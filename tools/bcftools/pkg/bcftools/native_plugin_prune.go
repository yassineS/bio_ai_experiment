// Native port of the upstream `prune` plugin (plugins/prune.c) and the slice of
// vcfbuf.c it drives. It implements every upstream mode:
//
//   - -n/--nsites-per-win with the "1st", "maxAF" (incl. the default maxAF
//     without --AF-tag, computed from INFO/AC+AN or the genotypes) and "rand"
//     selection modes,
//   - -m count=N cluster removal,
//   - -m R2=/LD=/RD= linkage-disequilibrium thresholding (hard drop or, with
//     -f LABEL, soft FILTER),
//   - -a count|r2|LD|RD annotation,
//   - --keep-sites (-k), -i/-e filtering, -r/-R/-t/-T region/target selection,
//   - --randomize-missing and --random-seed.
//
// The LD measures are computed by calcR2LD (native_plugin_prune_ld.go), a
// byte-exact port of _calc_r2_ld. The "rand" selection and randomize-missing
// reuse the deterministic drand48 generator (native_drand48.go) seeded exactly
// as hts_srand48(rseed); upstream seeds it once in init_data, so a single
// shared generator reproduces the draw order byte-for-byte.
package bcftools

import (
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("prune", func() NativePlugin { return &prunePlugin{} }) }

// prune nsites selection modes, matching PRUNE_MODE_* in vcfbuf.c.
const (
	pruneModeMaxAF = 1
	pruneMode1st   = 2
	pruneModeRand  = 3
)

// prunePlugin implements every pruning/annotation mode of `prune`.
type prunePlugin struct {
	hdr        *vcf.Header
	win        int // vcfbuf window: <0 for bp, >0 for number of sites
	nsites     int // -n: keep at most this many sites per window (0 = off)
	nsitesMode int
	afTag      string
	maxCluster int  // -m count=N (0 = off)
	clusterSet bool // -m count= was given
	clusterCnt bool // -a count (CLUSTER_MODE_SIZE annotation)

	// LD thresholding (-m R2=/LD=/RD=) and annotation (-a r2/LD/RD).
	ldMaxSet [ldN]bool
	ldMax    [ldN]float64
	ldAnnot  [ldN]bool // -a requested this measure
	ldMask   bool      // any -m or -a (LD_SET_MAX | LD_ANNOTATE)
	ldSetMax bool      // any -m R2/LD/RD threshold
	ldFilter string    // -f LABEL ("" if none, "." means annotate-only no filter)
	randMiss bool      // --randomize-missing
	rseed    int64     // --random-seed (default 0)
	rseedSet bool      // an explicit --random-seed was given
	rng      *drand48  // shared deterministic generator (rand selection / rand-missing)
	keepSite bool      // -k/--keep-sites

	filter      *Filter // compiled -i/-e expression (nil if none)
	filterLogic int     // pruneFilterInclude / pruneFilterExclude
	filterExpr  string  // raw expression text

	clusterHdrWin int // window value used in the CLUSTER_SIZE header (pre bp-conversion)
}

// prune filter-logic values, mirroring FLT_INCLUDE / FLT_EXCLUDE.
const (
	pruneFilterInclude = 1 // -i
	pruneFilterExclude = 2 // -e
)

// INFO tag names upstream uses for the -a annotation, per measure index.
var (
	pruneAnnotTag = [ldN]string{"R2", "LD", "RD"}
	pruneAnnotPos = [ldN]string{"POS_R2", "POS_LD", "POS_RD"}
)

// Name returns the plugin name.
func (p *prunePlugin) Name() string { return "prune" }

// RegionTargetCaps opts prune into the shared -r/-R/-t/-T region/target filter,
// applied to the records before the window pruning.
func (p *prunePlugin) RegionTargetCaps() regionTargetCaps { return allRegionTargetCaps }

// About returns the one-line description, matching prune.c about().
func (p *prunePlugin) About() string {
	return "Annotate sites with or prune sites by linkage disequilibrium or number of sites within a window"
}

// RunStyle reports that prune is a run()-style plugin: its options precede the
// input file with no `--` separator.
func (p *prunePlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of prune's own value-taking flags consumes
// the following CLI token, used by the host to split the input-file positional.
func (p *prunePlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-n", "--nsites-per-win", "-N", "--nsites-per-win-mode", "-w", "--window",
		"-m", "--max", "-a", "--annotate", "-f", "--set-filter", "--AF-tag",
		"--random-seed", "-i", "--include", "-e", "--exclude",
		"-r", "--regions", "-R", "--regions-file", "-t", "--targets", "-T", "--targets-file",
		"-o", "--output", "-O", "--output-type", "-v", "--verbosity":
		return true
	}
	return false
}

// Init parses the plugin options and prepares the header, supporting every
// upstream mode.
func (p *prunePlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.hdr = hdr
	p.win = -100000 // -100e3, the upstream default
	p.nsitesMode = pruneModeMaxAF

	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("prune: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-n", "--nsites-per-win":
			v, err := next()
			if err != nil {
				return nil, err
			}
			n, err := strconv.Atoi(v)
			if err != nil {
				return nil, fmt.Errorf("prune: could not parse: --nsites-per-win %s", v)
			}
			p.nsites = n
		case "-N", "--nsites-per-win-mode":
			v, err := next()
			if err != nil {
				return nil, err
			}
			switch {
			case eqFoldStr(v, "maxAF"):
				p.nsitesMode = pruneModeMaxAF
			case eqFoldStr(v, "1st"):
				p.nsitesMode = pruneMode1st
			case eqFoldStr(v, "rand"):
				p.nsitesMode = pruneModeRand
			default:
				return nil, fmt.Errorf("prune: the mode %q is not recognised", v)
			}
		case "--AF-tag":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.afTag = v
		case "-w", "--window":
			v, err := next()
			if err != nil {
				return nil, err
			}
			if err := p.parseWindow(v); err != nil {
				return nil, err
			}
		case "-m", "--max":
			v, err := next()
			if err != nil {
				return nil, err
			}
			if err := p.parseMax(v); err != nil {
				return nil, err
			}
			p.ldMask = true
		case "-a", "--annotate":
			v, err := next()
			if err != nil {
				return nil, err
			}
			if err := p.parseAnnotate(v); err != nil {
				return nil, err
			}
			p.ldMask = true
		case "-f", "--set-filter":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.ldFilter = v
		case "-i", "--include", "-e", "--exclude":
			v, err := next()
			if err != nil {
				return nil, err
			}
			if p.filterExpr != "" {
				return nil, fmt.Errorf("prune: only one of -i or -e can be given")
			}
			p.filterExpr = v
			if a == "-e" || a == "--exclude" {
				p.filterLogic = pruneFilterExclude
			} else {
				p.filterLogic = pruneFilterInclude
			}
		case "--random-seed":
			v, err := next()
			if err != nil {
				return nil, err
			}
			n, err := strconv.ParseInt(v, 10, 64)
			if err != nil {
				return nil, fmt.Errorf("prune: could not parse: --random-seed %s", v)
			}
			p.rseed = n
			p.rseedSet = true
		case "--randomize-missing":
			p.randMiss = true
		case "-k", "--keep-sites":
			p.keepSite = true
		case "--no-version":
			// no-op
		case "-o", "--output", "-O", "--output-type", "-v", "--verbosity":
			if _, err := next(); err != nil {
				return nil, err
			}
		default:
			return nil, fmt.Errorf("prune: unsupported option %q", a)
		}
	}

	if !p.ldMask && p.nsites == 0 {
		return nil, fmt.Errorf("prune: expected pruning (--max, --nsites-per-win) or annotation (--annotate) options")
	}
	if p.ldFilter != "" && p.ldFilter != "." && !p.ldSetMax {
		return nil, fmt.Errorf("prune: the --set-filter option requires --max")
	}
	if p.keepSite && p.nsites != 0 {
		return nil, fmt.Errorf("prune: the --keep-sites option cannot be combined with --nsites-per-win")
	}

	// The CLUSTER_SIZE header line (printed by init_data BEFORE the window
	// conversion below) uses the raw, pre-conversion window value.
	p.clusterHdrWin = p.win

	// Upstream warns and converts a positive window to bp when used with -m
	// count or -a count; mirror the conversion (the warning is stderr-only).
	if p.win > 0 && (p.maxCluster != 0 || p.clusterCnt) {
		p.win *= -1
		if p.maxCluster != 0 && -p.win <= p.maxCluster {
			return nil, fmt.Errorf("prune: -w must be bigger than -m")
		}
	}

	// Seed the shared RNG exactly as init_data does when rand-missing or the
	// rand selection mode is active. An unset --random-seed defaults to 0; the
	// emitted ##bcftools_plugin_prune_RandomSeed header line is provenance and
	// is stripped from parity comparisons.
	if p.randMiss || p.nsitesMode == pruneModeRand {
		p.rng = newDrand48(p.rseed)
	}

	if p.filterExpr != "" {
		f, err := CompileFilterWithHeader(p.filterExpr, p.hdr)
		if err != nil {
			return nil, fmt.Errorf("prune: %w", err)
		}
		p.filter = f
	}

	return p.buildHeader(), nil
}

// parseMax ports the -m/--max parsing: R2=/LD=/RD=/HD= set an LD threshold,
// count= sets the cluster size, and a bare number is an r2 threshold.
func (p *prunePlugin) parseMax(v string) error {
	parseFloat := func(s string) (float64, error) {
		f, err := strconv.ParseFloat(strings.TrimSpace(s), 64)
		if err != nil {
			return 0, fmt.Errorf("prune: could not parse: --max %s", v)
		}
		return f, nil
	}
	lc := strings.ToLower(v)
	switch {
	case strings.HasPrefix(lc, "r2="):
		f, err := parseFloat(v[3:])
		if err != nil {
			return err
		}
		p.ldMaxSet[ldIdxR2] = true
		p.ldMax[ldIdxR2] = f
		p.ldSetMax = true
	case strings.HasPrefix(lc, "ld="):
		f, err := parseFloat(v[3:])
		if err != nil {
			return err
		}
		p.ldMaxSet[ldIdxLD] = true
		p.ldMax[ldIdxLD] = f
		p.ldSetMax = true
	case strings.HasPrefix(lc, "rd=") || strings.HasPrefix(lc, "hd="):
		f, err := parseFloat(v[3:])
		if err != nil {
			return err
		}
		p.ldMaxSet[ldIdxRD] = true
		p.ldMax[ldIdxRD] = f
		p.ldSetMax = true
	case strings.HasPrefix(lc, "count="):
		n, err := strconv.Atoi(strings.TrimSpace(v[6:]))
		if err != nil {
			return fmt.Errorf("prune: could not parse: --max %s", v)
		}
		p.maxCluster = n
		p.clusterSet = true
	default:
		f, err := parseFloat(v)
		if err != nil {
			return err
		}
		p.ldMaxSet[ldIdxR2] = true
		p.ldMax[ldIdxR2] = f
		p.ldSetMax = true
	}
	return nil
}

// parseAnnotate ports the -a/--annotate parsing: a comma-separated list of
// count|r2|LD|RD (HD is an alias for RD), each enabling the matching annotation.
func (p *prunePlugin) parseAnnotate(v string) error {
	for _, tag := range strings.Split(v, ",") {
		switch {
		case eqFoldStr(tag, "r2"):
			p.ldAnnot[ldIdxR2] = true
		case eqFoldStr(tag, "ld"):
			p.ldAnnot[ldIdxLD] = true
		case eqFoldStr(tag, "rd") || eqFoldStr(tag, "hd"):
			p.ldAnnot[ldIdxRD] = true
		case eqFoldStr(tag, "count"):
			p.clusterCnt = true
		default:
			return fmt.Errorf("prune: the tag %q is not supported", tag)
		}
	}
	return nil
}

// buildHeader appends the ##FILTER / ##INFO lines exactly as init_data does
// (and in the same order), and returns the augmented header.
func (p *prunePlugin) buildHeader() *vcf.Header {
	out := &vcf.Header{Samples: p.hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, p.hdr.MetaInfo...)
	add := func(line string) { out.MetaInfo = appendInfoHeader(out.MetaInfo, line) }

	if p.ldFilter != "" && p.ldFilter != "." {
		var sb strings.Builder
		if p.ldMaxSet[ldIdxR2] {
			sb.WriteString("R2 bigger than ")
			sb.WriteString(formatVCFFloat(p.ldMax[ldIdxR2]))
		}
		if p.ldMaxSet[ldIdxLD] {
			if sb.Len() > 0 {
				sb.WriteString(" or ")
			}
			sb.WriteString("LD bigger than ")
			sb.WriteString(formatVCFFloat(p.ldMax[ldIdxLD]))
		}
		if p.ldMaxSet[ldIdxRD] {
			if sb.Len() > 0 {
				sb.WriteString(" or ")
			}
			sb.WriteString("RD bigger than ")
			sb.WriteString(formatVCFFloat(p.ldMax[ldIdxRD]))
		}
		winNum := p.win
		winUnit := " sites"
		if p.win < 0 {
			winNum = -p.win / 1000 // upstream: -ld_win/1000 (integer truncation)
			winUnit = "kb"
		}
		add(fmt.Sprintf(`##FILTER=<ID=%s,Description="An upstream site within %d%s with %s">`,
			p.ldFilter, winNum, winUnit, sb.String()))
	}
	// Annotation INFO lines, positions paired after their value line per
	// measure (matching init_data's order).
	if p.ldAnnot[ldIdxR2] {
		add(`##INFO=<ID=R2,Number=1,Type=Float,Description="Pairwise r2 with the POS_R2 site">`)
		add(`##INFO=<ID=POS_R2,Number=1,Type=Integer,Description="The position of the site for which R2 was calculated">`)
	}
	if p.ldAnnot[ldIdxLD] {
		add(`##INFO=<ID=LD,Number=1,Type=Float,Description="Pairwise Lewontin's D' (PMID:19433632) with the POS_LD site">`)
		add(`##INFO=<ID=POS_LD,Number=1,Type=Integer,Description="The position of the site for which LD was calculated">`)
	}
	if p.ldAnnot[ldIdxRD] {
		add(`##INFO=<ID=RD,Number=1,Type=Float,Description="Pairwise Ragsdale's \hat{D} (PMID:31697386) with the POS_RD site">`)
		add(`##INFO=<ID=POS_RD,Number=1,Type=Integer,Description="The position of the site for which RD was calculated">`)
	}
	if p.clusterCnt {
		add(fmt.Sprintf(`##INFO=<ID=CLUSTER_SIZE,Number=1,Type=Integer,Description="The number of variants within %d bp of the site">`, p.clusterHdrWin))
	}
	return out
}

// parseWindow ports the -w/--window parsing: a bare integer is a site count
// (positive), while a bp/kb/Mb suffix yields a negative base-pair window.
func (p *prunePlugin) parseWindow(s string) error {
	num := s
	suffix := ""
	for i := 0; i < len(s); i++ {
		if (s[i] < '0' || s[i] > '9') && s[i] != '-' && s[i] != '+' {
			num = s[:i]
			suffix = s[i:]
			break
		}
	}
	n, err := strconv.Atoi(num)
	if err != nil {
		return fmt.Errorf("prune: could not parse: --window %s", s)
	}
	switch {
	case suffix == "":
		p.win = n
	case eqFoldStr(suffix, "bp"):
		p.win = -n
	case eqFoldStr(suffix, "kb"):
		p.win = -n * 1000
	case eqFoldStr(suffix, "Mb"):
		p.win = -n * 1000000
	default:
		return fmt.Errorf("prune: could not parse: --window %s", s)
	}
	return nil
}

// ProcessAll runs the chosen pruning/annotation mode over the whole ordered
// stream, faithfully simulating the vcfbuf push + flush(0) per record then
// flush(1) cycle.
func (p *prunePlugin) ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error) {
	if p.ldMask && !p.clusterSet && !p.clusterCnt {
		// LD thresholding / annotation path (process() calls vcfbuf_ld before
		// pushing, then the pruning-window flush slides the buffer).
		return p.runLD(variants), nil
	}
	if p.clusterSet || p.clusterCnt {
		return p.runCluster(variants), nil
	}
	// -n nsites path: a non-matching record is simply dropped (keep-sites is
	// disallowed with -n).
	if p.filter != nil {
		kept := make([]*vcf.Variant, 0, len(variants))
		for _, v := range variants {
			if p.passFilter(v) {
				kept = append(kept, v)
			}
		}
		variants = kept
	}
	return p.runNsites(variants), nil
}

// passFilter returns whether the record survives the -i/-e expression (true
// when no filter is set). It does NOT apply keep-sites semantics.
func (p *prunePlugin) passFilter(v *vcf.Variant) bool {
	if p.filter == nil {
		return true
	}
	pass := p.filter.Eval(v)
	if p.filterLogic == pruneFilterInclude {
		return pass
	}
	return !pass
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *prunePlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *prunePlugin) Destroy() error { return nil }

// addPruneFilter ports bcf_add_filter: a record that is PASS/"." has its
// FILTER replaced by the label; otherwise the label is appended (deduplicated).
func addPruneFilter(v *vcf.Variant, label string) {
	if len(v.Filter) == 0 || (len(v.Filter) == 1 && (v.Filter[0] == "PASS" || v.Filter[0] == ".")) {
		v.Filter = []string{label}
		return
	}
	for _, f := range v.Filter {
		if f == label {
			return
		}
	}
	v.Filter = append(v.Filter, label)
}

// nsiteRec pairs a record with its (lazily computed) allele frequency for the
// maxAF selection sort.
type nsiteRec struct {
	src   int
	v     *vcf.Variant
	af    float64
	afSet bool
}

// runNsites simulates the vcfbuf window buffer for the -n path, porting
// vcfbuf_flush's `if (buf->win)` block and _prune_sites.
func (p *prunePlugin) runNsites(variants []*vcf.Variant) []*vcf.Variant {
	out := make([]*vcf.Variant, 0, len(variants))
	var buf []nsiteRec
	emitFront := func() {
		out = append(out, buf[0].v)
		buf = buf[1:]
	}
	flush := func(flushAll bool) {
		for len(buf) > 0 {
			if !p.windowCanFlush(buf, flushAll) {
				return
			}
			if p.nsites != 0 && p.nsites < len(buf) {
				p.pruneSites(&buf, flushAll)
			}
			if len(buf) == 0 {
				return
			}
			emitFront()
		}
	}
	for i, v := range variants {
		buf = append(buf, nsiteRec{src: i, v: v})
		flush(false)
	}
	flush(true)
	return out
}

// windowCanFlush ports the `if (buf->win)` can_flush test in vcfbuf_flush.
func (p *prunePlugin) windowCanFlush(buf []nsiteRec, flushAll bool) bool {
	if flushAll {
		return true
	}
	first := buf[0].v
	last := buf[len(buf)-1].v
	if first.Chrom != last.Chrom {
		return true
	}
	if p.win > 0 {
		return len(buf) > p.win
	}
	// win < 0: flushable once the span reaches the window size.
	return !(first.Pos-last.Pos > p.win)
}

// pruneSites ports _prune_sites for the 1st, rand and maxAF modes, removing
// nprune lowest-priority records from the buffer. nbuf excludes the just-added
// last record unless flushing all.
func (p *prunePlugin) pruneSites(bufp *[]nsiteRec, flushAll bool) {
	buf := *bufp
	nbuf := len(buf)
	if !flushAll {
		nbuf = len(buf) - 1
	}
	nprune := nbuf - p.nsites
	if nprune <= 0 {
		return
	}

	switch p.nsitesMode {
	case pruneMode1st:
		eoff := 2
		if flushAll {
			eoff = 1
		}
		for i := 0; i < nprune; i++ {
			k := len(buf) - eoff
			buf = append(buf[:k], buf[k+1:]...)
		}
		*bufp = buf
		return
	case pruneModeRand:
		eoff := 1
		if flushAll {
			eoff = 0
		}
		for i := 0; i < nprune; i++ {
			j := int(float64(len(buf)-eoff) * p.rng.float64())
			buf = append(buf[:j], buf[j+1:]...)
		}
		*bufp = buf
		return
	}

	// maxAF: compute AF for the first nbuf records, sort ascending (low AF
	// removed preferentially), then remove the nprune lowest by descending
	// rbuf index so earlier removals do not shift later ones. The sort is
	// stable (glibc qsort behaves stably for these tiny arrays), so the tie
	// order matches upstream — see UPSTREAM_BUGS.md.
	recs := make([]*nsiteRec, nbuf)
	for i := 0; i < nbuf; i++ {
		if !buf[i].afSet {
			buf[i].af = p.computeAF(buf[i].v)
			buf[i].afSet = true
		}
		recs[i] = &buf[i]
	}
	idxOf := make(map[*nsiteRec]int, nbuf)
	for i := 0; i < nbuf; i++ {
		idxOf[recs[i]] = i
	}
	sort.SliceStable(recs, func(a, b int) bool { return recs[a].af < recs[b].af })
	rm := make([]int, nprune)
	for i := 0; i < nprune; i++ {
		rm[i] = idxOf[recs[i]]
	}
	sort.Sort(sort.Reverse(sort.IntSlice(rm)))
	for _, k := range rm {
		buf = append(buf[:k], buf[k+1:]...)
	}
	*bufp = buf
}

// computeAF mirrors the AF used by _prune_sites: INFO/<af_tag>[0] when an AF tag
// is given, otherwise nalt/ntot from INFO/AC+AN or, failing that, the genotypes
// (bcf_calc_ac). The result is narrowed to float32 to match upstream exactly.
func (p *prunePlugin) computeAF(v *vcf.Variant) float64 {
	if p.afTag != "" {
		if s, ok := v.Info[p.afTag]; ok {
			if c := strings.IndexByte(s, ','); c >= 0 {
				s = s[:c]
			}
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				return float64(float32(f))
			}
		}
		return 0
	}
	// bcf_calc_ac fills ac[0] with the REFERENCE allele count and ac[1..] with
	// the ALT counts. _prune_sites then uses ntot=ac[0] (ref count) as the
	// denominator: af = nalt/nref, or 0 when nref==0. This is an upstream
	// quirk (af is alt/ref, not a true frequency) reproduced for parity — see
	// UPSTREAM_BUGS.md.
	nref, nalt, ok := calcAC(v)
	if !ok || nref == 0 {
		return 0
	}
	return float64(float32(float64(nalt) / float64(nref)))
}

// clusterBuf is a faithful re-implementation of the vcfbuf ring buffer driving
// the cluster paths. The parallel size ring accumulates the monotonic (max)
// cluster size as the window slides.
type clusterBuf struct {
	v      []*vcf.Variant
	size   []int
	filter []bool
	dirty  bool
}

// runCluster simulates the CLUSTER_MODE_PRUNE / CLUSTER_MODE_SIZE path: clusters
// of more than max_cluster sites within a -win window are removed (-m count=),
// or each site is annotated with its cluster size (-a count). It ports
// cluster_can_flush_ and the cluster branch of vcfbuf_flush, including the
// per-record `filter` flag set by an -i/-e expression.
func (p *prunePlugin) runCluster(variants []*vcf.Variant) []*vcf.Variant {
	out := make([]*vcf.Variant, 0, len(variants))
	win := -p.win // positive bp window
	b := &clusterBuf{}

	push := func(v *vcf.Variant) {
		// process(): with -i/-e, a non-matching record is dropped before push
		// unless --keep-sites, in which case it is pushed with filter=1.
		filtered := false
		if p.filter != nil && !p.passFilter(v) {
			if !p.keepSite {
				return
			}
			filtered = true
		}
		b.v = append(b.v, v)
		b.size = append(b.size, 0)
		b.filter = append(b.filter, filtered)
		b.dirty = true
		p.clusterFlush(b, win, false, &out)
	}
	for _, v := range variants {
		push(v)
	}
	p.clusterFlush(b, win, true, &out)
	return out
}

// clusterFlush drains all flushable records, porting the caller's
// `while ((rec=vcfbuf_flush(...)))` loop together with cluster_can_flush_ and
// the CLUSTER_MODE_PRUNE pruning / CLUSTER_MODE_SIZE annotation.
func (p *prunePlugin) clusterFlush(b *clusterBuf, win int, flushAll bool, out *[]*vcf.Variant) {
	for len(b.v) > 0 {
		size, ok := p.clusterCanFlush(b, win, flushAll)
		if !ok {
			return
		}
		front := b.v[0]
		if p.clusterCnt && size > 0 {
			setInfo(front, "CLUSTER_SIZE", strconv.Itoa(size))
		}
		*out = append(*out, front)
		b.v = b.v[1:]
		b.size = b.size[1:]
		b.filter = b.filter[1:]
	}
}

// clusterCanFlush ports cluster_can_flush_: it recomputes the monotonic cluster
// sizes when dirty, checks the window flush condition, and (in PRUNE mode)
// prunes leading oversized-cluster records before allowing the front to be
// emitted. It returns the front record's cluster size and whether it flushes.
func (p *prunePlugin) clusterCanFlush(b *clusterBuf, win int, flushAll bool) (int, bool) {
	if b.dirty {
		for ib := 0; ib < len(b.v); ib++ {
			ie := ib + 1
			for ie < len(b.v) {
				if b.v[ie].Chrom != b.v[ib].Chrom {
					break
				}
				if b.v[ie].Pos-b.v[ib].Pos+1 > win {
					break
				}
				ie++
			}
			// Count only the unfiltered sites that contribute to the cluster.
			n := 0
			for ix := ib; ix < ie; ix++ {
				if b.filter[ix] {
					continue
				}
				n++
			}
			for ix := ib; ix < ie; ix++ {
				if b.size[ix] < n {
					b.size[ix] = n
				}
			}
		}
		b.dirty = false
	}

	first := b.v[0]
	last := b.v[len(b.v)-1]
	canFlush := flushAll
	if first.Chrom != last.Chrom {
		canFlush = true
	}
	if last.Pos-first.Pos+1 > win {
		canFlush = true
	}
	if !canFlush {
		return 0, false
	}

	if p.clusterSet { // CLUSTER_MODE_PRUNE
		flush := false
		for len(b.v) > 0 {
			flush = false
			first = b.v[0]
			last = b.v[len(b.v)-1]
			if b.filter[0] {
				// not to be pruned, not counted as part of the cluster
				flush = true
				break
			}
			switch {
			case flushAll:
				flush = true
			case first.Chrom != last.Chrom:
				flush = true
			case last.Pos-first.Pos+1 > win:
				flush = true
			}
			if !flush {
				break
			}
			if b.size[0] <= p.maxCluster {
				break // front not to be pruned
			}
			b.v = b.v[1:]
			b.size = b.size[1:]
			b.filter = b.filter[1:]
		}
		if !flush {
			return 0, false
		}
	}

	if len(b.v) == 0 {
		return 0, false
	}
	size := b.size[0]
	if b.filter[0] {
		size = 0
	}
	return size, true
}

// ldRec is a buffered record for the LD path, carrying its -i/-e filter flag
// and a cached rand-missing allele frequency.
type ldRec struct {
	v        *vcf.Variant
	filter   bool
	randMiss float64
	afSet    bool
}

// runLD ports process()/vcfbuf_ld/vcfbuf_flush for the LD threshold/annotate
// path. For each record it computes the maximum LD against the buffered sites
// within the window, applies the hard/soft filter and annotations, then pushes
// the record and slides the window.
func (p *prunePlugin) runLD(variants []*vcf.Variant) []*vcf.Variant {
	out := make([]*vcf.Variant, 0, len(variants))
	var buf []ldRec

	flushWindow := func(flushAll bool) {
		for len(buf) > 0 {
			if !flushAll {
				first := buf[0].v
				last := buf[len(buf)-1].v
				canFlush := false
				if first.Chrom != last.Chrom {
					canFlush = true
				} else if p.win > 0 {
					if len(buf) > p.win {
						canFlush = true
					}
				} else if p.win < 0 {
					if !(first.Pos-last.Pos > p.win) {
						canFlush = true
					}
				}
				if !canFlush {
					return
				}
			}
			out = append(out, buf[0].v)
			buf = buf[1:]
		}
	}

	for _, rec := range variants {
		filtered := false
		if p.filter != nil && !p.passFilter(rec) {
			filtered = true
			if !p.keepSite {
				continue
			}
		}
		drop := false
		if !filtered {
			drop = p.applyLD(buf, rec)
		}
		if drop {
			// hard filter: the record is pruned and never pushed.
			continue
		}
		buf = append(buf, ldRec{v: rec, filter: filtered})
		flushWindow(false)
	}
	flushWindow(true)
	return out
}

// applyLD ports the vcfbuf_ld scan plus the hard/soft-filter and annotation
// logic in process() for the current record rec against the buffer. It returns
// true when the record is hard-filtered (pruned, not pushed).
func (p *prunePlugin) applyLD(buf []ldRec, rec *vcf.Variant) bool {
	if len(buf) == 0 {
		return false
	}
	// vcfbuf_ld bails immediately if the buffer's front is on another chrom
	// (the buffer is always flushed so all sites share rec's chromosome).
	if buf[0].v.Chrom != rec.Chrom {
		return false
	}

	var best ldResult
	for j := 0; j < ldN; j++ {
		best.val[j] = -1e308 // -HUGE_VAL
	}
	any := false

	baf := 0.0
	if p.randMiss {
		baf = pruneEstimateAF(rec)
	}

	n := len(buf)
	for i := 0; i < n; i++ {
		if buf[i].filter {
			continue
		}
		if !p.ldInsideWindow(buf, i, rec, n) {
			continue
		}
		aaf := 0.0
		if p.randMiss {
			if !buf[i].afSet {
				buf[i].randMiss = pruneEstimateAF(buf[i].v)
				buf[i].afSet = true
			}
			aaf = buf[i].randMiss
		}
		tmp, ok := calcR2LD(buf[i].v, rec, p.randMiss, aaf, baf, p.rng)
		if !ok {
			continue
		}
		done := false
		for j := 0; j < ldN; j++ {
			if best.val[j] < tmp.val[j] {
				best.val[j] = tmp.val[j]
				best.pos[j] = buf[i].v.Pos
				best.set[j] = true
			}
			if p.ldMaxSet[j] && p.ldMax[j] < tmp.val[j] {
				done = true
			}
			any = true
		}
		if done {
			break
		}
	}
	if !any {
		return false
	}

	// Hard/soft filter: a record is pruned (or soft-filtered) when any
	// thresholded measure exceeds its max.
	pass := true
	for j := 0; j < ldN; j++ {
		if !p.ldMaxSet[j] {
			continue
		}
		if best.val[j] > p.ldMax[j] {
			pass = false
			break
		}
	}
	dropHard := false
	if !pass {
		if p.ldFilter == "" {
			dropHard = true // hard filter: drop
		} else if p.ldFilter != "." {
			addPruneFilter(rec, p.ldFilter)
		}
	}

	// Annotation: positions first, then values (matching process()).
	for j := 0; j < ldN; j++ {
		if p.ldAnnot[j] && best.set[j] {
			setInfo(rec, pruneAnnotPos[j], strconv.Itoa(best.pos[j]))
		}
	}
	for j := 0; j < ldN; j++ {
		if p.ldAnnot[j] && best.set[j] {
			setInfo(rec, pruneAnnotTag[j], formatVCFFloat(best.val[j]))
		}
	}
	return dropHard
}

// ldInsideWindow ports the inside_win test in vcfbuf_ld for buffered record i
// relative to the current record rec (n = current buffer length).
func (p *prunePlugin) ldInsideWindow(buf []ldRec, i int, rec *vcf.Variant, n int) bool {
	if buf[i].v.Chrom != rec.Chrom {
		return false
	}
	if p.win > 0 {
		// rbuf_ridx2l: 0-based linear index from the front; n - linidx is the
		// distance (in sites) from the back. Outside if that exceeds win.
		if n-i > p.win {
			return false
		}
	} else if p.win < 0 {
		// inside_win=0 when !(pos_i - pos_rec > win); win is negative.
		if !(buf[i].v.Pos-rec.Pos > p.win) {
			return false
		}
	}
	return true
}
