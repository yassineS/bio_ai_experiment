// Native port of the upstream `prune` plugin (plugins/prune.c) for its window
// based pruning modes: -n/--nsites-per-win (keeping at most N sites per window,
// in "1st" and "maxAF" selection modes) and -m count=N (removing clusters of
// more than N sites within a window). These drive the vcfbuf window/cluster
// machinery, which is re-implemented here as a faithful streaming simulation.
//
// The linkage-disequilibrium modes (-a LD/RD/r2 annotation and -m R2=/LD=/RD=
// thresholding) require genotype-correlation calculations (calc_ld in
// vcfbuf.c) that the native pipeline does not provide, and the "rand" selection
// mode depends on htslib's hts_drand48 RNG; all of these are reported as
// unsupported rather than produced incorrectly.
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
	pruneModeMaxAF = iota
	pruneMode1st
)

// prunePlugin implements the supported window/count pruning modes of `prune`.
type prunePlugin struct {
	hdr        *vcf.Header
	win        int // vcfbuf window: <0 for bp, >0 for number of sites
	nsites     int // -n: keep at most this many sites per window (0 = off)
	nsitesMode int
	afTag      string
	maxCluster int // -m count=N (0 = off)
	clusterSet bool

	filter      *Filter // compiled -i/-e expression (nil if none)
	filterLogic int     // pruneFilterInclude / pruneFilterExclude
	filterExpr  string  // raw expression text
}

// prune filter-logic values, mirroring FLT_INCLUDE / FLT_EXCLUDE.
const (
	pruneFilterInclude = 1 // -i
	pruneFilterExclude = 2 // -e
)

// Name returns the plugin name.
func (p *prunePlugin) Name() string { return "prune" }

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

// Init parses the plugin options, supporting only the window/count pruning
// modes and reporting the LD/annotation/rand paths as unsupported.
func (p *prunePlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.hdr = hdr
	p.win = -100000 // -100e3, the upstream default
	p.nsitesMode = pruneModeMaxAF
	ldMask := false

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
				return nil, fmt.Errorf("prune: the \"rand\" selection mode is not supported by the native plugin (RNG cannot be matched byte-for-byte)")
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
			lc := strings.ToLower(v)
			if strings.HasPrefix(lc, "count=") {
				n, err := strconv.Atoi(v[6:])
				if err != nil {
					return nil, fmt.Errorf("prune: could not parse: --max %s", v)
				}
				p.maxCluster = n
				p.clusterSet = true
			} else {
				return nil, fmt.Errorf("prune: LD/r2/RD thresholding (--max %s) is not supported by the native plugin", v)
			}
			ldMask = true
		case "-a", "--annotate":
			return nil, fmt.Errorf("prune: LD/r2/RD/count annotation (-a) is not supported by the native plugin")
		case "-f", "--set-filter":
			return nil, fmt.Errorf("prune: soft-filter pruning (-f) is not supported by the native plugin")
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
		case "-r", "--regions", "-R", "--regions-file", "-t", "--targets", "-T", "--targets-file":
			return nil, fmt.Errorf("prune: index/stream region selection is not supported by the native plugin")
		case "--random-seed", "--randomize-missing":
			return nil, fmt.Errorf("prune: randomization options are not supported by the native plugin")
		case "-k", "--keep-sites":
			return nil, fmt.Errorf("prune: --keep-sites is not supported by the native plugin")
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

	if !ldMask && p.nsites == 0 {
		return nil, fmt.Errorf("prune: expected pruning options (--max count= or --nsites-per-win)")
	}

	// The default maxAF selection without an explicit --AF-tag relies on
	// htslib's bcf_calc_ac plus an unstable qsort whose tie ordering is
	// implementation-defined; it cannot be reproduced byte-for-byte. The
	// maxAF mode is therefore supported only when an --AF-tag supplies an
	// unambiguous, deterministic allele frequency.
	if p.nsites != 0 && p.nsitesMode == pruneModeMaxAF && p.afTag == "" {
		return nil, fmt.Errorf("prune: the default maxAF selection (without --AF-tag) is not supported by the native plugin; use -N 1st or supply --AF-tag for deterministic pruning")
	}

	// Upstream warns and converts a positive window to bp when used with -m
	// count; mirror the conversion (the warning goes to stderr and is not part
	// of stdout parity).
	if p.win > 0 && p.maxCluster != 0 {
		p.win *= -1
		if -p.win <= p.maxCluster {
			return nil, fmt.Errorf("prune: -w must be bigger than -m")
		}
	}

	if p.filterExpr != "" {
		f, err := CompileFilter(p.filterExpr)
		if err != nil {
			return nil, fmt.Errorf("prune: %w", err)
		}
		p.filter = f
	}
	return hdr, nil
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

// ProcessAll runs the chosen pruning mode over the whole ordered stream and
// returns the surviving records, faithfully simulating the vcfbuf window/cluster
// flush cycle (push + flush(0) per record, then flush(1)).
func (p *prunePlugin) ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error) {
	// A -i/-e expression drops non-matching records before they enter the
	// window/cluster buffer (prune.c process(): with !keep_sites a filtered
	// record returns early and is never pushed). --keep-sites is unsupported,
	// so the pre-filter is an unconditional drop.
	if p.filter != nil {
		kept := make([]*vcf.Variant, 0, len(variants))
		for _, v := range variants {
			pass := p.filter.Eval(v)
			if p.filterLogic == pruneFilterInclude {
				if !pass {
					continue
				}
			} else if pass {
				continue
			}
			kept = append(kept, v)
		}
		variants = kept
	}
	if p.clusterSet {
		return p.runCluster(variants), nil
	}
	return p.runNsites(variants), nil
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *prunePlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *prunePlugin) Destroy() error { return nil }

// nsiteRec pairs a record with its (lazily computed) allele frequency for the
// maxAF selection sort.
type nsiteRec struct {
	src   int
	v     *vcf.Variant
	af    float32
	afSet bool
}

// runNsites simulates the vcfbuf window buffer for the -n path. Records are
// pushed in order; whenever the buffer can flush (window exceeded, chromosome
// change, or end of stream), excess sites beyond -n are pruned per the selection
// mode, then the front record is emitted.
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
				p.pruneSites(buf, flushAll, &buf)
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
	// upstream: can_flush if !(first.pos - last.pos > win), i.e. span >= -win.
	return !(first.Pos-last.Pos > p.win)
}

// pruneSites ports _prune_sites for the 1st and maxAF modes, removing nprune
// lowest-priority records from the buffer. nbuf excludes the just-added last
// record unless flushing all.
func (p *prunePlugin) pruneSites(buf []nsiteRec, flushAll bool, out *[]nsiteRec) {
	nbuf := len(buf)
	if !flushAll {
		nbuf = len(buf) - 1
	}
	nprune := nbuf - p.nsites
	if nprune <= 0 {
		return
	}

	if p.nsitesMode == pruneMode1st {
		eoff := 2
		if flushAll {
			eoff = 1
		}
		for i := 0; i < nprune; i++ {
			k := len(buf) - eoff
			buf = append(buf[:k], buf[k+1:]...)
		}
		*out = buf
		return
	}

	// maxAF: compute AF for the first nbuf records, sort ascending (low AF
	// removed preferentially), then remove the nprune lowest by descending index
	// so earlier removals do not shift later ones.
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
	*out = buf
}

// computeAF mirrors the AF used by _prune_sites: INFO/<af_tag>[0] when an AF tag
// is given, otherwise nalt/ntot from INFO/AC+AN or, failing that, from the
// genotypes (bcf_calc_ac).
func (p *prunePlugin) computeAF(v *vcf.Variant) float32 {
	if p.afTag != "" {
		if s, ok := v.Info[p.afTag]; ok {
			if c := strings.IndexByte(s, ','); c >= 0 {
				s = s[:c]
			}
			if f, err := strconv.ParseFloat(s, 64); err == nil {
				return float32(f)
			}
		}
		return 0
	}
	ntot, nalt, ok := calcAC(v)
	if !ok || ntot == 0 {
		return 0
	}
	return float32(nalt) / float32(ntot)
}

// clusterBuf is a faithful re-implementation of the vcfbuf ring buffer driving
// the CLUSTER_MODE_PRUNE path. The main ring and the parallel cluster-size ring
// are simulated together so the monotonic (max) size accumulation across flush
// cycles matches cluster_can_flush_ exactly.
type clusterBuf struct {
	v     []*vcf.Variant
	size  []int
	dirty bool
}

// runCluster simulates the CLUSTER_MODE_PRUNE path: clusters of more than
// max_cluster sites within a -win window are removed. It ports
// cluster_can_flush_ and the CLUSTER_MODE_PRUNE branch of vcfbuf_flush,
// including the persistent (max) cluster sizes that accumulate as the window
// slides. Filters are not supported by the native plugin, so every buffered
// record contributes to its cluster.
func (p *prunePlugin) runCluster(variants []*vcf.Variant) []*vcf.Variant {
	out := make([]*vcf.Variant, 0, len(variants))
	win := -p.win // positive bp window
	b := &clusterBuf{}

	for _, v := range variants {
		b.v = append(b.v, v)
		b.size = append(b.size, 0)
		b.dirty = true
		p.clusterFlush(b, win, false, &out)
	}
	p.clusterFlush(b, win, true, &out)
	return out
}

// clusterFlush drains all flushable records for the current buffer state,
// porting the caller's `while ((rec=vcfbuf_flush(...)))` loop together with
// cluster_can_flush_ and the CLUSTER_MODE_PRUNE pruning.
func (p *prunePlugin) clusterFlush(b *clusterBuf, win int, flushAll bool, out *[]*vcf.Variant) {
	for len(b.v) > 0 {
		if !p.clusterCanFlush(b, win, flushAll) {
			return
		}
		// vcfbuf_flush "ret": shift the front record and emit it.
		*out = append(*out, b.v[0])
		b.v = b.v[1:]
		b.size = b.size[1:]
	}
}

// clusterCanFlush ports cluster_can_flush_: it (re)computes the monotonic
// cluster sizes when dirty, checks the window flush condition, and then prunes
// leading oversized-cluster records before allowing the front to be emitted.
func (p *prunePlugin) clusterCanFlush(b *clusterBuf, win int, flushAll bool) bool {
	if b.dirty {
		// For each window start ib, find ie (first index outside the window),
		// count the sites, and raise each member's size to that count (max).
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
			n := ie - ib
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
		return false
	}

	// CLUSTER_MODE_PRUNE: remove leading records belonging to an oversized
	// cluster, one at a time, until the front is keepable or no longer ready.
	flush := false
	for len(b.v) > 0 {
		flush = false
		first = b.v[0]
		last = b.v[len(b.v)-1]
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
	}
	if !flush {
		return false
	}
	return len(b.v) > 0
}

// calcAC computes (total alleles, alt alleles) for a record, preferring
// INFO/AN+AC and falling back to the genotypes, mirroring bcf_calc_ac for the
// INFO|FMT request used by _prune_sites.
func calcAC(v *vcf.Variant) (ntot, nalt int, ok bool) {
	an, anOK := v.Info["AN"]
	ac, acOK := v.Info["AC"]
	if anOK && acOK {
		if anv, err := strconv.Atoi(an); err == nil {
			sum := 0
			good := true
			for _, s := range strings.Split(ac, ",") {
				n, err := strconv.Atoi(s)
				if err != nil {
					good = false
					break
				}
				sum += n
			}
			if good {
				return anv, sum, true
			}
		}
	}
	// Fall back to genotypes.
	tot, alt := 0, 0
	any := false
	for i := range v.Samples {
		gt, gok := sampleGT(v, i)
		if !gok {
			continue
		}
		for _, a := range gt.alleles {
			if a == missingAllele {
				continue
			}
			any = true
			tot++
			if a > 0 {
				alt++
			}
		}
	}
	if !any {
		return 0, 0, false
	}
	return tot, alt, true
}
