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
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
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

// remove-overlaps filter-logic values, mirroring FLT_INCLUDE / FLT_EXCLUDE.
const (
	roFilterInclude = 1 // -i
	roFilterExclude = 2 // -e
)

// remove-overlaps --missing modes, mirroring MARK_MISSING_SCALAR /
// MARK_MISSING_MAX_DP in vcfbuf.c.
const (
	markMissingScalar = iota // --missing 0 (the default): use a scalar value
	markMissingMaxDP         // --missing DP: scale max QUAL by INFO/DP
)

// removeOverlapsPlugin implements the `remove-overlaps` plugin.
type removeOverlapsPlugin struct {
	hdr          *vcf.Header
	markExpr     string // "overlap" (default), "dup", or "min(QUAL)"
	markMode     int
	markTag      string // -M: mark instead of remove, via this INFO flag
	reverse      bool   // --reverse: invert the keep/remove decision
	missingExpr  string // raw --missing value ("0" or "DP")
	missingMode  int    // markMissingScalar / markMissingMaxDP
	missingValue float64
	textOutput   bool           // -Ot / -Otz: emit a plain "chr\tpos" list instead of VCF
	textGz       bool           // -Otz: bgzip-compress the text list
	kept         []*vcf.Variant // records emitted in text mode (collected for Destroy)
	ntot         int
	nrm          int
	stderr       io.Writer
	stdout       io.Writer

	filter      *Filter // compiled -i/-e expression (nil if none)
	filterLogic int     // roFilterInclude / roFilterExclude
	filterExpr  string  // raw expression text
}

// Name returns the plugin name.
func (p *removeOverlapsPlugin) Name() string { return "remove-overlaps" }

// RegionTargetCaps opts remove-overlaps into the shared -r/-R/-t/-T
// region/target filter, applied to the records before the overlap/dup removal.
func (p *removeOverlapsPlugin) RegionTargetCaps() regionTargetCaps { return allRegionTargetCaps }

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

// stdout is the host stdout writer; only used for -Ot text-list output.
// (SetStdout / SuppressVCF are part of the outputSuppressor contract and only
// engage when textOutput is set; see SuppressVCF.)
func (p *removeOverlapsPlugin) SetStdout(w io.Writer) { p.stdout = w }

// SuppressVCF reports whether the framework should skip the VCF/BCF re-emit.
// remove-overlaps only suppresses it for -Ot text output, where the plugin
// writes its own "chr\tpos" list to stdout; in every other mode the VCF stream
// flows through the normal stage-3 writer.
func (p *removeOverlapsPlugin) SuppressVCF() bool { return p.textOutput }

// Init parses the plugin options. The supported set is the default overlap
// removal/marking, dup removal/marking, the min(QUAL) resolution (optionally
// with the --missing DP heuristic), the -M mark tag, --reverse, the -i/-e
// filter, region/target selection (handled by the host), and the -Ot/-Otz
// text-list output container.
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
			v, err := next()
			if err != nil {
				return nil, err
			}
			if p.filterExpr != "" {
				return nil, fmt.Errorf("remove-overlaps: only one of -i or -e can be given")
			}
			p.filterExpr = v
			if a == "-e" || a == "--exclude" {
				p.filterLogic = roFilterExclude
			} else {
				p.filterLogic = roFilterInclude
			}
		case "--missing":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.missingExpr = v
		case "-O", "--output-type":
			v, err := next()
			if err != nil {
				return nil, err
			}
			// -Ot (plain text list of chr,pos) and -Otz (bgzip-compressed list)
			// are handled by this plugin; the BCF/VCF -O forms (b/u/z/v) are
			// handled by the host pipeline and ignored here.
			if len(v) > 0 && v[0] == 't' {
				p.textOutput = true
				if len(v) > 1 && v[1] == 'z' {
					p.textGz = true
				}
			}
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
		// vcfbuf.c maps any other expression to MARK_EXPR, but
		// mark_expr_can_flush_ then errors "Todo; at this time only min(QUAL) is
		// supported". Reject early with that same restriction.
		return nil, fmt.Errorf("remove-overlaps: the --mark expression %q is not supported (only overlap, dup, or min(QUAL))", p.markExpr)
	}

	// --missing accepts only "0" (scalar 0, the default) or "DP" (scale the
	// window's maximum QUAL proportionally to INFO/DP), mirroring vcfbuf.c
	// MARK_MISSING_EXPR. DP requires -m 'min(QUAL)'.
	p.missingMode = markMissingScalar
	p.missingValue = 0
	if p.missingExpr != "" {
		switch {
		case p.missingExpr == "0":
			p.missingMode = markMissingScalar
			p.missingValue = 0
		case eqFoldStr(p.missingExpr, "DP"):
			if p.markMode != markExpr {
				return nil, fmt.Errorf("remove-overlaps: only the combination of --mark 'min(QUAL)' with --missing DP is currently supported")
			}
			p.missingMode = markMissingMaxDP
		default:
			return nil, fmt.Errorf("remove-overlaps: todo: MARK_MISSING_EXPR=%s", p.missingExpr)
		}
	}

	if p.filterExpr != "" {
		f, err := CompileFilterWithHeader(p.filterExpr, p.hdr)
		if err != nil {
			return nil, fmt.Errorf("remove-overlaps: %w", err)
		}
		p.filter = f
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

	// A -i/-e expression pre-filters the stream: upstream process() counts every
	// input record in ntot but only pushes records that pass into the overlap
	// buffer (remove-overlaps.c:199-211). Filtered-out records never reach the
	// window and are not emitted.
	if p.filter != nil {
		kept := make([]*vcf.Variant, 0, len(variants))
		for _, v := range variants {
			pass := p.filter.Eval(v)
			if p.filterLogic == roFilterInclude {
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

	var marked []bool
	switch p.markMode {
	case markOverlap:
		marked = markOverlapStream(variants)
	case markDup:
		marked = markDupStream(variants)
	case markExpr:
		marked = markMinQualStream(variants, p.missingMode, p.missingValue)
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
	if p.textOutput {
		// In text mode the framework suppresses the VCF re-emit (SuppressVCF);
		// the kept records are written as a "chr\tpos" list in Destroy, after
		// which the Processed/Removed summary goes to stderr — matching the order
		// upstream's flush() then run() produce.
		p.kept = out
		return nil, nil
	}
	return out, nil
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *removeOverlapsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy writes the -Ot text list (if any) to stdout, then prints the
// "Processed/Removed" summary to stderr, matching upstream's flush()-then-run()
// ordering.
func (p *removeOverlapsPlugin) Destroy() error {
	if p.textOutput && p.stdout != nil {
		// Each emitted record becomes one "chr\tpos\n" line, exactly as
		// remove-overlaps.c flush() ksprintf's it for FT_TAB_TEXT output. With
		// -Otz the list is bgzip-framed (upstream bgzf_open mode "wg").
		w := p.stdout
		var bw *bgzf.Writer
		if p.textGz {
			bw = bgzf.NewWriter(p.stdout)
			w = bw
		}
		for _, v := range p.kept {
			if v == nil {
				continue
			}
			fmt.Fprint(w, textListLine(v.Chrom, v.Pos))
		}
		if bw != nil {
			if err := bw.Close(); err != nil {
				return err
			}
		}
	}
	if p.stderr != nil {
		fmt.Fprintf(p.stderr, "Processed/Removed\t%d\t%d\n", p.ntot, p.nrm)
	}
	return nil
}

// textListLine formats one record as upstream's FT_TAB_TEXT line ("chr\tpos\n").
// It is exported as a pure helper for the binary-free unit test.
func textListLine(chrom string, pos int) string {
	return fmt.Sprintf("%s\t%d\n", chrom, pos)
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
func markMinQualStream(variants []*vcf.Variant, missingMode int, missingValue float64) []bool {
	marked := make([]bool, len(variants))
	b := &minQualBuf{
		variants:     variants,
		missingMode:  missingMode,
		missingValue: missingValue,
		emit:         func(src int, mk bool) { marked[src] = mk },
	}
	for i := range variants {
		b.push(i)
		b.flush(false)
	}
	b.flush(true)
	return marked
}

// minQualBuf is a faithful port of the vcfbuf MARK_EXPR state machine for the
// "min(QUAL)" expression (vcfbuf.c mark_expr_can_flush_ + mark_overlap_helper_).
// It reproduces upstream's subtle push/flush ordering: each pushed record sets
// the buffer dirty; a flush over a dirty buffer re-runs the full value-ordered
// overlap resolution over the WHOLE current buffer (the accumulated overlapping
// records plus the first non-overlapping one that triggered the flush), records
// the resulting marks, then drains records one at a time using those marks until
// a single record is left in the buffer (or, on flush_all, all of them). Because
// the marks are recomputed only on a dirty flush, this state machine — not a
// simple per-window resolution — is what matches upstream byte-for-byte.
type minQualBuf struct {
	variants     []*vcf.Variant
	missingMode  int
	missingValue float64
	emit         func(src int, mk bool)

	buf   []int  // FIFO of source indices currently buffered (buf->rbuf)
	mk    []bool // parallel mark FIFO (mark->mark / mark->rbuf)
	dirty bool

	overlapRid string
	overlapEnd int
}

// push appends a source index and marks the buffer dirty (vcfbuf_push).
func (b *minQualBuf) push(src int) {
	b.buf = append(b.buf, src)
	b.dirty = true
}

// flush drains every flushable record, calling emit(src, mark) for each,
// mirroring the caller's `while ((rec=vcfbuf_flush(buf,flush_all)))` loop.
func (b *minQualBuf) flush(flushAll bool) {
	for len(b.buf) > 0 {
		ok := b.canFlush(flushAll)
		if !ok {
			return
		}
		src := b.buf[0]
		mk := b.mk[0]
		b.buf = b.buf[1:]
		b.mk = b.mk[1:]
		b.emit(src, mk)
	}
}

// canFlush ports mark_expr_can_flush_: on a dirty buffer it runs the overlap
// helper to decide whether to flush and, if so, recomputes the marks over the
// whole buffer; on a clean buffer it flushes when more than one record remains
// (draining using the marks already computed).
func (b *minQualBuf) canFlush(flushAll bool) bool {
	flush := flushAll
	if b.dirty {
		flush = b.overlapHelper(flushAll)
		if !flush {
			return false
		}
		b.resolve()
	} else if len(b.buf) > 1 {
		flush = true
	}
	return flush
}

// overlapHelper ports mark_overlap_helper_: it updates the running overlap span
// for the just-added record and returns whether a flush is due (a new
// non-overlapping record arrived, or flush_all). It also clears the dirty flag
// and grows the mark FIFO for the new record.
func (b *minQualBuf) overlapHelper(flushAll bool) bool {
	if !b.dirty {
		return flushAll
	}
	flush := flushAll
	b.dirty = false
	b.mk = append(b.mk, false)

	last := b.variants[b.buf[len(b.buf)-1]]
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
	if len(b.buf) == 1 {
		b.overlapRid = last.Chrom
		b.overlapEnd = endPos
		return flush
	}
	if begPos <= b.overlapEnd {
		if b.overlapEnd < endPos {
			b.overlapEnd = endPos
		}
		// The just-added record overlaps its predecessor: mark BOTH, exactly as
		// mark_overlap_helper_ does (mark->mark[k1]=mark->mark[k2]=1). These
		// overlap marks persist unless a later flush re-runs resolve() (which
		// resets them); that persistence is what removes a trailing overlapping
		// pair at end-of-file, where no resolution is triggered.
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

// resolve recomputes mark->mark over the whole buffer, porting the overlap-graph
// construction and value-ordered marking in mark_expr_can_flush_.
func (b *minQualBuf) resolve() {
	m := len(b.buf)
	run := make([]*vcf.Variant, m)
	for i, src := range b.buf {
		run[i] = b.variants[src]
	}
	value := minQualValues(run, b.missingMode, b.missingValue)
	for i := range b.mk {
		b.mk[i] = false
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
	sort.SliceStable(order, func(a, c int) bool { return value[order[a]] < value[order[c]] })
	for _, oi := range order {
		if nolap == 0 {
			break
		}
		for oj := range adj[oi] {
			delete(adj[oi], oj)
			delete(adj[oj], oi)
			nolap--
		}
		b.mk[oi] = true
	}
}

// minQualValues computes the per-record overlap-resolution value for one window,
// porting the value assignment in mark_expr_can_flush_ together with the
// --missing DP heuristic (mark_expr_missing_*). Present-QUAL records use QUAL;
// missing-QUAL records use missingValue (scalar mode) or, in DP mode, the
// window's maximum QUAL scaled proportionally to INFO/DP:
// value = maxQual * DP / maxQualDP (0 when no usable DP). All arithmetic is done
// in float32 to match htslib's `float` value field byte-for-byte. It is exported
// as a pure helper for the binary-free unit test.
func minQualValues(run []*vcf.Variant, missingMode int, missingValue float64) []float64 {
	m := len(run)
	value := make([]float64, m)
	// First pass: present-QUAL records get QUAL; missing get the scalar default.
	// In DP mode, also track the highest present QUAL and its INFO/DP.
	var maxQual float32
	var maxQualDP int
	dp := make([]int, m)
	for i, v := range run {
		if missingMode == markMissingMaxDP {
			if d, ok := infoInt0(v, "DP"); ok {
				dp[i] = d
			}
		}
		if v.Qual < 0 { // missing QUAL
			value[i] = missingValue
			continue
		}
		value[i] = v.Qual
		if missingMode == markMissingMaxDP && maxQual < float32(v.Qual) {
			maxQual = float32(v.Qual)
			maxQualDP = dp[i]
		}
	}
	// Second pass (DP mode only): rescale missing-QUAL records by coverage.
	if missingMode == markMissingMaxDP && maxQualDP != 0 {
		for i, v := range run {
			if v.Qual >= 0 { // QUAL present, already valued
				continue
			}
			value[i] = float64(maxQual * float32(dp[i]) / float32(maxQualDP))
		}
	}
	return value
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

// infoInt0 returns the first integer value of an INFO tag, mirroring the
// bcf_get_info_int32(...,"DP")==1 use in mark_expr_missing_prep_ (a single
// scalar INFO/DP). It returns ok=false when the tag is absent, missing, or not
// a single integer.
func infoInt0(v *vcf.Variant, tag string) (int, bool) {
	s, ok := v.Info[tag]
	if !ok || s == "" || s == "." {
		return 0, false
	}
	if strings.IndexByte(s, ',') >= 0 {
		// More than one value: upstream's nval!=1 guard skips it.
		return 0, false
	}
	n, err := strconv.Atoi(s)
	if err != nil {
		return 0, false
	}
	return n, true
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
