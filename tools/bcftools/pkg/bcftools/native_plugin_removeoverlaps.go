// Native port of the upstream `remove-overlaps` plugin
// (plugins/remove-overlaps.c) together with the vcfbuf MARK machinery it drives
// (the mark_overlap/mark_dup/mark_expr logic in vcfbuf.c). It removes, lists, or
// marks overlapping or duplicate variants within a streaming window and prints a
// "Processed/Removed" summary to stderr. The plugin needs cross-record state, so
// it runs as a serial bufferedPlugin.
package bcftools

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("remove-overlaps", func() NativePlugin { return &removeOverlapsPlugin{} })
}

// remove-overlaps mark modes, matching MARK_OVERLAP/MARK_DUP/MARK_EXPR.
const (
	markOverlap = iota
	markDup
	markExpr
)

// removeOverlapsPlugin implements the `remove-overlaps` plugin.
type removeOverlapsPlugin struct {
	hdr      *vcf.Header
	markExpr string // "overlap" (default), "dup", or "min(QUAL)"
	markMode int
	markTag  string // -M: mark instead of remove, via this INFO flag
	reverse  bool   // --reverse: invert the keep/remove decision
	ntot     int
	nrm      int
	stderr   io.Writer
}

// Name returns the plugin name.
func (p *removeOverlapsPlugin) Name() string { return "remove-overlaps" }

// About returns the one-line description, matching remove-overlaps.c about().
func (p *removeOverlapsPlugin) About() string {
	return "Remove, list or mark overlapping variants"
}

// RunStyle reports that remove-overlaps is a run()-style plugin: upstream
// exports a `run` symbol, so its options precede the input file with no `--`
// separator.
func (p *removeOverlapsPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of remove-overlaps' own value-taking flags
// consumes the following CLI token, so the host can split the input file out of
// the run-style options.
func (p *removeOverlapsPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-m", "--mark", "-M", "--mark-tag", "--missing",
		"-i", "--include", "-e", "--exclude",
		"-r", "--regions", "-R", "--regions-file",
		"-t", "--targets", "-T", "--targets-file",
		"-o", "--output", "-O", "--output-type", "-v", "--verbosity":
		return true
	}
	return false
}

// SetStderr wires the host stderr writer the summary line is printed to.
func (p *removeOverlapsPlugin) SetStderr(w io.Writer) { p.stderr = w }

// Init parses the plugin options and rejects modes that depend on machinery the
// native pipeline does not provide (filter expressions, index regions, the
// min(QUAL) --missing DP heuristic, and text-list output). The supported set is
// the default overlap removal/marking, dup removal/marking, the min(QUAL)
// resolution, the -M mark tag, and --reverse.
func (p *removeOverlapsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.hdr = hdr
	p.markExpr = "overlap"
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("remove-overlaps: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-m", "--mark":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.markExpr = v
		case "-M", "--mark-tag":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.markTag = v
		case "--reverse":
			p.reverse = true
		case "--no-version":
			// no-op for the native path (provenance is handled elsewhere)
		case "-i", "--include", "-e", "--exclude":
			return nil, fmt.Errorf("remove-overlaps: filter expressions (-i/-e) are not supported by the native plugin")
		case "--missing":
			return nil, fmt.Errorf("remove-overlaps: the --missing DP heuristic is not supported by the native plugin")
		case "-r", "--regions", "-R", "--regions-file", "-t", "--targets", "-T", "--targets-file":
			return nil, fmt.Errorf("remove-overlaps: index/stream region selection (-r/-R/-t/-T) is not supported by the native plugin")
		case "-O", "--output-type":
			v, err := next()
			if err != nil {
				return nil, err
			}
			if len(v) > 0 && (v[0] == 't' || v[0] == 'T') {
				return nil, fmt.Errorf("remove-overlaps: text list output (-Ot) is not supported by the native plugin")
			}
			// other -O forms are handled by the host pipeline; ignore here
		case "-o", "--output", "-v", "--verbosity":
			if _, err := next(); err != nil {
				return nil, err
			}
		default:
			return nil, fmt.Errorf("remove-overlaps: unsupported option %q", a)
		}
	}

	switch {
	case eqFoldStr(p.markExpr, "overlap"):
		p.markMode = markOverlap
	case eqFoldStr(p.markExpr, "dup"):
		p.markMode = markDup
	case eqFoldStr(p.markExpr, "min(QUAL)"):
		p.markMode = markExpr
	default:
		return nil, fmt.Errorf("remove-overlaps: the --mark expression %q is not supported by the native plugin", p.markExpr)
	}

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	if p.markTag != "" {
		line := fmt.Sprintf(`##INFO=<ID=%s,Type=Flag,Number=0,Description="Marked by +remove-overlaps">`, p.markTag)
		out.MetaInfo = appendInfoHeader(out.MetaInfo, line)
	}
	return out, nil
}

// ProcessAll marks the whole ordered stream according to the chosen mode, then
// keeps or drops (or flags, with -M) each record exactly as the upstream flush()
// does. The "Processed/Removed" summary is recorded for Destroy.
func (p *removeOverlapsPlugin) ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error) {
	p.ntot = len(variants)

	var marked []bool
	switch p.markMode {
	case markOverlap:
		marked = markOverlapStream(variants)
	case markDup:
		marked = markDupStream(variants)
	case markExpr:
		marked = markMinQualStream(variants)
	}

	out := make([]*vcf.Variant, 0, len(variants))
	for i, v := range variants {
		keep := !marked[i]
		if p.reverse {
			keep = !keep
		}
		if !keep {
			p.nrm++
			if p.markTag == "" {
				continue // remove
			}
			setInfoFlag(v, p.markTag) // mark, keep
		}
		out = append(out, v)
	}
	return out, nil
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *removeOverlapsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy prints the "Processed/Removed" summary to stderr, matching upstream.
func (p *removeOverlapsPlugin) Destroy() error {
	if p.stderr != nil {
		fmt.Fprintf(p.stderr, "Processed/Removed\t%d\t%d\n", p.ntot, p.nrm)
	}
	return nil
}

// recordSpanBeg computes the left-aligned begin offset (imin) for a record,
// porting the per-allele common-prefix trimming in mark_overlap_helper_.
func recordSpanImin(v *vcf.Variant) int {
	rlen := refLen(v)
	imin := rlen
	for ai := 0; ai <= len(v.Alt); ai++ {
		var alt string
		if ai == 0 {
			alt = v.Ref
		} else {
			alt = v.Alt[ai-1]
		}
		if len(alt) > 0 && alt[0] == '<' {
			continue // symbolic allele
		}
		r, a := 0, 0
		for r < len(v.Ref) && a < len(alt) && upperByte(v.Ref[r]) == upperByte(alt[a]) {
			r++
			a++
		}
		if imin > r {
			imin = r
		}
	}
	return imin
}

// overlapBuf is a faithful re-implementation of the vcfbuf ring buffer driving
// the mark_overlap path. It simulates the upstream push/flush cycle so the
// marking (including the subtle buffer-size==1 handling) is byte-equivalent.
type overlapBuf struct {
	idx        []int  // main ring: indices into the source slice, in order
	mk         []bool // parallel mark ring
	dirty      bool
	overlapRid string // "" mirrors overlap_rid==-1 (no chromosome yet)
	overlapEnd int
}

// markOverlapStream marks overlapping records by simulating the streaming
// push/flush cycle of vcfbuf.c (process()/flush()), exactly as upstream does.
func markOverlapStream(variants []*vcf.Variant) []bool {
	marked := make([]bool, len(variants))
	b := &overlapBuf{}
	emit := func(i int, mk bool) { marked[i] = mk }
	for i := range variants {
		b.push(i)
		b.flush(variants, false, emit)
	}
	b.flush(variants, true, emit)
	return marked
}

// push appends a record index to the buffer and marks it dirty, mirroring
// vcfbuf_push.
func (b *overlapBuf) push(i int) {
	b.idx = append(b.idx, i)
	b.dirty = true
}

// flush drains all flushable records, invoking emit(srcIdx, marked) for each,
// mirroring the caller's `while ((rec=vcfbuf_flush(...)))` loop.
func (b *overlapBuf) flush(variants []*vcf.Variant, flushAll bool, emit func(int, bool)) {
	for len(b.idx) > 0 {
		ok, mk := b.canFlush(variants, flushAll)
		if !ok {
			return
		}
		src := b.idx[0]
		b.idx = b.idx[1:]
		emit(src, mk)
	}
}

// canFlush ports mark_overlap_can_flush_ + the shift in vcfbuf_flush, returning
// whether a record can be flushed and its mark.
func (b *overlapBuf) canFlush(variants []*vcf.Variant, flushAll bool) (bool, bool) {
	flush := flushAll
	if b.dirty {
		flush = b.helper(variants, flushAll)
	} else if len(b.idx) > 1 {
		flush = true
	}
	if !flush {
		return false, false
	}
	mk := b.mk[0]
	b.mk = b.mk[1:]
	return true, mk
}

// helper ports mark_overlap_helper_, updating the running overlap span and
// marking the last-added record (and its predecessor) on an overlap.
func (b *overlapBuf) helper(variants []*vcf.Variant, flushAll bool) bool {
	if !b.dirty {
		return flushAll
	}
	flush := flushAll
	b.dirty = false
	b.mk = append(b.mk, false)

	last := variants[b.idx[len(b.idx)-1]]
	if b.overlapRid != last.Chrom {
		b.overlapEnd = 0
	}
	begPos := last.Pos
	endPos := last.Pos + refLen(last) - 1
	imin := recordSpanImin(last)
	if begPos <= b.overlapEnd {
		begPos += imin
		if begPos > endPos {
			endPos = begPos
		}
	}
	if len(b.idx) == 1 {
		b.overlapRid = last.Chrom
		b.overlapEnd = endPos
		return flush
	}
	if begPos <= b.overlapEnd {
		if b.overlapEnd < endPos {
			b.overlapEnd = endPos
		}
		b.mk[len(b.mk)-1] = true
		b.mk[len(b.mk)-2] = true
	} else {
		if b.overlapEnd < endPos {
			b.overlapEnd = endPos
		}
		flush = true
	}
	return flush
}

// markDupStream marks records that share chromosome+position with the
// immediately preceding record, porting mark_dup_can_flush_.
func markDupStream(variants []*vcf.Variant) []bool {
	marked := make([]bool, len(variants))
	for i := 1; i < len(variants); i++ {
		if variants[i].Chrom == variants[i-1].Chrom && variants[i].Pos == variants[i-1].Pos {
			marked[i] = true
			marked[i-1] = true
		}
	}
	return marked
}

// markMinQualStream resolves overlaps by iteratively removing the lowest-QUAL
// site until no overlaps remain, porting mark_expr_can_flush_ for the
// "min(QUAL)" expression. Records are grouped into overlap-connected components
// over contiguous runs; within each run the overlap graph is built (records
// overlap when on the same chromosome and their reference spans touch), then the
// lowest-value record is repeatedly marked and its edges removed.
func markMinQualStream(variants []*vcf.Variant) []bool {
	marked := make([]bool, len(variants))
	n := len(variants)
	// Partition into maximal runs where each record overlaps the previous span,
	// matching the buffer that mark_overlap_helper_ accumulates before a flush.
	i := 0
	for i < n {
		j := i + 1
		runEnd := variants[i].Pos + refLen(variants[i]) - 1
		for j < n {
			if variants[j].Chrom != variants[i].Chrom {
				break
			}
			if variants[j].Pos > runEnd {
				break
			}
			e := variants[j].Pos + refLen(variants[j]) - 1
			if e > runEnd {
				runEnd = e
			}
			j++
		}
		if j-i > 1 {
			resolveMinQual(variants[i:j], marked[i:j])
		}
		i = j
	}
	return marked
}

// resolveMinQual marks the minimal set of records (lowest QUAL first) so that no
// two remaining records overlap, building the symmetric overlap graph as
// mark_expr_can_flush_ does and breaking ties by original order (stable sort).
func resolveMinQual(run []*vcf.Variant, marked []bool) {
	m := len(run)
	value := make([]float64, m)
	for i, v := range run {
		if v.Qual < 0 {
			value[i] = 0 // missing QUAL -> default missing value 0
		} else {
			value[i] = v.Qual
		}
	}
	adj := make([]map[int]bool, m)
	for i := range adj {
		adj[i] = map[int]bool{}
	}
	nolap := 0
	for i := 0; i < m; i++ {
		for j := i + 1; j < m; j++ {
			if recordsOverlap(run[i], run[j]) {
				adj[i][j] = true
				adj[j][i] = true
				nolap++
			}
		}
	}
	order := make([]int, m)
	for i := range order {
		order[i] = i
	}
	sort.SliceStable(order, func(a, b int) bool { return value[order[a]] < value[order[b]] })
	for _, oi := range order {
		if nolap == 0 {
			break
		}
		for oj := range adj[oi] {
			delete(adj[oi], oj)
			delete(adj[oj], oi)
			nolap--
		}
		marked[oi] = true
	}
}

// recordsOverlap reports whether two records share a chromosome and their
// reference spans touch, porting records_overlap() from vcfbuf.c.
func recordsOverlap(a, b *vcf.Variant) bool {
	if a.Chrom != b.Chrom {
		return false
	}
	if a.Pos+refLen(a)-1 < b.Pos {
		return false
	}
	return true
}

// setInfoFlag adds a value-less INFO flag (Number=0,Type=Flag) to v.
func setInfoFlag(v *vcf.Variant, key string) {
	if v.Info == nil {
		v.Info = make(map[string]string)
	}
	if _, ok := v.Info[key]; !ok {
		v.InfoOrder = append(v.InfoOrder, key)
	}
	v.Info[key] = ""
}

// eqFoldStr reports ASCII case-insensitive equality, used to match the --mark
// expression keywords the way strcasecmp does.
func eqFoldStr(a, b string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		if upperByte(a[i]) != upperByte(b[i]) {
			return false
		}
	}
	return true
}
