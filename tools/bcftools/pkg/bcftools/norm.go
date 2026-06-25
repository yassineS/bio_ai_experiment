package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// CheckRefMode is a bitmask controlling how `bcftools norm` reacts when the
// REF field differs from the reference FASTA. It mirrors upstream
// vcfnorm.c's CHECK_REF_* flags, where the `-c`/`--check-ref` option accepts
// any combination of the letters w/x/s and the single letter e:
//
//	w  warn   (CheckRefWarn)  — print "REF_MISMATCH ..." to stderr, keep
//	x  skip   (CheckRefSkip)  — drop records whose REF cannot be realigned
//	s  set    (CheckRefFix)   — fix the REF (and swap/insert alleles) to match
//	e  error  (CheckRefExit)  — abort on any mismatch (the default; exclusive)
//
// As upstream does, `e` overrides every other bit; `w`, `x`, and `s` may be
// combined (e.g. `-c ws`).
type CheckRefMode int

const (
	// CheckRefExit aborts the run on the first REF/FASTA mismatch. This is
	// the upstream default. It is exclusive: combining it with any other
	// bit is not meaningful.
	CheckRefExit CheckRefMode = 0
	// CheckRefWarn prints one "REF_MISMATCH" line per mismatched record to
	// stderr and keeps the record unchanged.
	CheckRefWarn CheckRefMode = 1 << iota
	// CheckRefSkip drops records whose REF cannot be reconciled with the
	// FASTA (the upstream `x` letter).
	CheckRefSkip
	// CheckRefFix rewrites the REF (and, when possible, swaps/inserts
	// alleles) so it matches the FASTA (the upstream `s` letter).
	CheckRefFix
)

// CheckRefError is retained as a deprecated alias for CheckRefExit so that
// older callers keep compiling. New code should use CheckRefExit.
const CheckRefError = CheckRefExit

// ParseCheckRefMode turns the `-c`/`--check-ref` flag value into the bitmask.
// The recognised letters mirror upstream: any combination of w/x/s, or e
// (which overrides the others). The empty string is the default (exit).
func ParseCheckRefMode(s string) (CheckRefMode, error) {
	if s == "" {
		return CheckRefExit, nil
	}
	var mode CheckRefMode
	for _, c := range strings.ToLower(s) {
		switch c {
		case 'w':
			mode |= CheckRefWarn
		case 'x':
			mode |= CheckRefSkip
		case 's':
			mode |= CheckRefFix
		case 'e':
			return CheckRefExit, nil // overrides the above
		default:
			return 0, fmt.Errorf("bcftools norm: unknown --check-ref value %q (expect a combination of w/x/s, or e)", s)
		}
	}
	return mode, nil
}

// MultiallelicMode encodes the body of the `-m` flag (`-snps`, `+indels` etc.).
type MultiallelicMode struct {
	// Active is true when the user passed `-m`. When false the splitter
	// and joiner are both off and the `-m` switch was not supplied.
	Active bool
	// Split is true for `-` modes (split multiallelics into biallelics).
	// When false and Active is true the mode is `+` (join biallelics into
	// multiallelics).
	Split bool
	// Snps controls whether SNP records are affected.
	Snps bool
	// Indels controls whether indel records are affected.
	Indels bool
	// Any mirrors upstream COLLAPSE_ANY (the `any` type for -m+/-m-). When
	// joining, every biallelic record at the same position is merged into a
	// single multiallelic regardless of variant type; the default ("both")
	// instead buckets records by type category before merging.
	Any bool
}

// ParseMultiallelicMode turns the literal flag body (e.g. "-both") into the
// typed structure. An empty string yields an inactive value. A bare "-" or
// "+" with no suffix is treated as "any" (both SNPs and indels), matching
// upstream bcftools.
func ParseMultiallelicMode(s string) (MultiallelicMode, error) {
	if s == "" {
		return MultiallelicMode{}, nil
	}
	var m MultiallelicMode
	m.Active = true
	switch s[0] {
	case '-':
		m.Split = true
	case '+':
		m.Split = false
	default:
		return m, fmt.Errorf("bcftools norm: -m must start with + or -")
	}
	rest := strings.ToLower(s[1:])
	switch rest {
	case "", "both":
		m.Snps = true
		m.Indels = true
	case "any":
		m.Snps = true
		m.Indels = true
		m.Any = true
	case "snps":
		m.Snps = true
	case "indels":
		m.Indels = true
	default:
		return m, fmt.Errorf("bcftools norm: unknown -m type %q (expect snps|indels|both|any)", rest)
	}
	return m, nil
}

// RmDupMode covers the `-d` / `--rm-dup` flag.
type RmDupMode int

const (
	// RmDupNone (the default) keeps all records.
	RmDupNone RmDupMode = iota
	// RmDupSnps drops duplicate SNP records sharing chrom/pos.
	RmDupSnps
	// RmDupIndels drops duplicate indel records sharing chrom/pos.
	RmDupIndels
	// RmDupBoth drops duplicate SNP or indel records sharing chrom/pos.
	RmDupBoth
	// RmDupAll drops any duplicate sharing chrom/pos regardless of type.
	RmDupAll
	// RmDupExact drops records that are byte-for-byte identical at the
	// CHROM / POS / REF / ALT level.
	RmDupExact
)

// ParseRmDupMode parses the flag value.
func ParseRmDupMode(s string) (RmDupMode, error) {
	switch strings.ToLower(s) {
	case "", "none":
		return RmDupNone, nil
	case "snps":
		return RmDupSnps, nil
	case "indels":
		return RmDupIndels, nil
	case "both":
		return RmDupBoth, nil
	case "all":
		return RmDupAll, nil
	case "exact":
		return RmDupExact, nil
	}
	return 0, fmt.Errorf("bcftools norm: unknown --rm-dup value %q", s)
}

// NormOptions controls Norm / NormFile behaviour.
type NormOptions struct {
	// FastaRef is the path to the reference FASTA. Required for
	// left-alignment and for REF checking.
	FastaRef string
	// CheckRef controls the response to REF/FASTA disagreement.
	CheckRef CheckRefMode
	// Multiallelics enables splitting / joining.
	Multiallelics MultiallelicMode
	// RmDup enables duplicate-record removal.
	RmDup RmDupMode
	// Atomize decomposes complex variants into single-base atomic events.
	Atomize bool
	// DoNotNormalize skips left-alignment. Useful in `-m` only pipelines.
	DoNotNormalize bool
	// StrictFilter applies -f filters before splitting; when false they
	// run after splitting (matching upstream's default ordering).
	StrictFilter bool
	// Regions / Targets filter the input on the fly. Region / target
	// semantics mirror the `view` subcommand.
	Regions      []string
	RegionsFile  string
	Targets      []string
	TargetsFile  string
	ApplyFilters []string
	// OutputFormat / CompressLevel mirror view's writer wiring.
	OutputFormat  OutputFormat
	CompressLevel int
	// Threads is the -@/--threads value; >1 enables parallel BGZF compression
	// of -O z and -O b output via bgzf.MultiWriter (see ViewOptions.Threads).
	Threads int
}

// NormFile is the high-level entry point matching ViewFile's signature.
// It opens path through iohelper.OpenReader, dispatches on the magic bytes,
// and writes the normalized output to out.
func NormFile(path string, out io.Writer, opts NormOptions, stderr io.Writer) (int, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, err
	}
	defer in.Close()
	return Norm(in, out, opts, stderr)
}

// Norm runs the normalize pipeline on the supplied reader.
//
// The pipeline streams: rather than buffering the whole VCF/BCF into memory and
// globally sorting it (which made peak RSS O(file) — ~9 GiB on a 168 MB GIAB
// VCF), it mirrors upstream vcfnorm.c's bounded buffering. Records flow through
// the per-record transforms (region/strict-filter, atomize, split, left-align +
// REF-check) and into a small reorder window; only records whose position can
// no longer be overtaken by a left-shift are flushed onward to the stateful
// rmdup / multiallelic-join / lax-filter sink. Memory is O(window), not
// O(file). See normStreamRun.
func Norm(in io.Reader, out io.Writer, opts NormOptions, stderr io.Writer) (int, error) {
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return 0, err
	}
	src, hdr, err := openVariantSource(br, head)
	if err != nil {
		return 0, err
	}
	return normStreamRun(hdr, src, out, opts, stderr)
}

// variantSource yields the input variants one at a time. next returns
// (nil, io.EOF) when the stream is exhausted. It abstracts the VCF and BCF
// readers so the streaming driver is format-agnostic.
type variantSource interface {
	next() (*vcf.Variant, error)
}

// openVariantSource sniffs the magic bytes and returns the matching streaming
// source plus the parsed header.
func openVariantSource(in io.Reader, head []byte) (variantSource, *vcf.Header, error) {
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		br, err := bcf.NewReader(in)
		if err != nil {
			return nil, nil, err
		}
		return &bcfSource{r: br, hdr: br.Header()}, br.Header().VCF, nil
	}
	r := vcf.NewReader(in)
	hdr, err := r.ReadHeader()
	if err != nil {
		return nil, nil, err
	}
	return &vcfSource{r: r}, hdr, nil
}

// vcfSource streams *vcf.Variant from a text VCF reader.
type vcfSource struct{ r *vcf.Reader }

func (s *vcfSource) next() (*vcf.Variant, error) { return s.r.Read() }

// bcfSource streams *vcf.Variant from a BCF reader, converting each record.
type bcfSource struct {
	r   *bcf.Reader
	hdr *bcf.Header
}

func (s *bcfSource) next() (*vcf.Variant, error) {
	rec, err := s.r.Read()
	if err != nil {
		return nil, err
	}
	return rec.ToVariant(s.hdr), nil
}

// normBufWin is the reorder-window width in base pairs, mirroring upstream
// vcfnorm.c's buf_win default (the -w/--site-win option, default 1000). A
// record can only be flushed once every still-buffered record on the same
// contig sits at least this many bases ahead of it, because left-alignment can
// shift a variant's POS left by at most a bounded amount. Holding the window
// (rather than the whole file) is what keeps norm's memory O(window).
const normBufWin = 1000

// normStreamRun is the streaming entry point. It mirrors upstream
// vcfnorm.c::normalize_vcf: each record flows through the per-record transforms
// (region/strict-filter, atomize, split, left-align + REF-check) and into a
// bounded reorder window. When the window's span exceeds normBufWin (or a new
// contig appears) the settled prefix is flushed in (rid, pos) order to a
// stateful sink that applies duplicate removal, multiallelic join, and
// lax-mode filtering before emitting. Peak memory is O(window), not O(file).
//
// The transform order matches the former slice pipeline exactly:
//
//  1. region / target filter
//  2. strict-filter (if -s)
//  3. atomize
//  4. multiallelic split
//  5. left-align + REF check
//  6. duplicate removal
//  7. multiallelic join
//  8. non-strict filter
//  9. emit
func normStreamRun(hdr *vcf.Header, src variantSource, out io.Writer, opts NormOptions, stderr io.Writer) (int, error) {
	regions, err := parseRegions(opts.Regions)
	if err != nil {
		return 0, err
	}
	targets, err := parseRegions(opts.Targets)
	if err != nil {
		return 0, err
	}
	regionFilter := append([]region{}, regions...)
	regionFilter = append(regionFilter, targets...)

	// Open the reference once if we need it for normalize or REF check.
	var ref *fasta.RandomAccess
	if opts.FastaRef != "" {
		ref, err = fasta.OpenRandomAccess(opts.FastaRef)
		if err != nil {
			return 0, fmt.Errorf("bcftools norm: open reference: %w", err)
		}
		defer ref.Close()
	}

	// Per-field Number= lookup, used to re-index Number=A/R/G INFO and FORMAT
	// vectors when splitting or joining multiallelic sites.
	numbers := headerNumberMapsFrom(hdr.MetaInfo)

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}

	sink := newNormSink(w, opts, numbers)
	win := &normWindow{}

	// flushSettled drains every settled prefix of the reorder window, feeding
	// each to the stateful sink. It loops because flushing one contig's prefix
	// can expose a second contig that is now also fully settled (the window may
	// hold records from more than two contigs after a contig switch). Records
	// dropped earlier (region/filter/skip) never reach the window.
	flushSettled := func() error {
		for {
			settled := win.popSettled()
			if len(settled) == 0 {
				return nil
			}
			for _, v := range settled {
				if err := sink.push(v); err != nil {
					return err
				}
			}
		}
	}

	for {
		v, rerr := src.next()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return 0, rerr
		}

		// 1. region/target filtering.
		if len(regionFilter) > 0 && !overlapsAny(v, regionFilter) {
			continue
		}
		// 2. strict-filter mode runs the FILTER list before any splitting.
		if opts.StrictFilter && len(opts.ApplyFilters) > 0 && !filterMatches(v, opts.ApplyFilters) {
			continue
		}

		// 3. atomize, 4. split — each can fan one input record into several.
		expanded := []*vcf.Variant{v}
		if opts.Atomize {
			expanded = atomizeVariants(expanded)
		}
		if opts.Multiallelics.Active && opts.Multiallelics.Split {
			expanded = splitMultiallelics(expanded, opts.Multiallelics, numbers)
		}

		// 5. left-align + REF-check, then 6/7/8/9 happen at flush time.
		for _, e := range expanded {
			if ref != nil {
				kept, nerr := normalizeOne(e, ref, opts, stderr)
				if nerr != nil {
					return 0, nerr
				}
				if !kept {
					continue
				}
			}
			win.push(e)
		}
		if err := flushSettled(); err != nil {
			return 0, err
		}
	}

	// Drain the window, then the sink's pending join buffer.
	for _, v := range win.drain() {
		if err := sink.push(v); err != nil {
			return 0, err
		}
	}
	if err := sink.flush(); err != nil {
		return 0, err
	}
	if err := w.Flush(); err != nil {
		return 0, err
	}
	return sink.written, nil
}

// filterMatches reports whether v's FILTER column names any of the requested
// filters (the per-record predicate behind applyFilterList).
func filterMatches(v *vcf.Variant, filters []string) bool {
	for _, want := range filters {
		for _, f := range v.Filter {
			if f == want {
				return true
			}
		}
	}
	return false
}

// normalizeOne left-aligns one variant and applies the --check-ref policy. It
// returns ok=false when the record is dropped (the skip bit on a REF mismatch).
//
// The order mirrors upstream vcfnorm.c::normalize_line: when the `s` (fix) bit
// is set the REF is repaired first, then the realignment REF check runs; a
// residual mismatch is reported per the remaining bits (error / warn / skip).
// This is why `-c s` (or `-c ws`) emits no warning for a record it could fix:
// by the time the check runs, the REF already matches.
func normalizeOne(v *vcf.Variant, ref *fasta.RandomAccess, opts NormOptions, stderr io.Writer) (bool, error) {
	if opts.CheckRef&CheckRefFix != 0 {
		if err := fixRef(v, ref); err != nil {
			return false, err
		}
	}
	ok, err := checkRef(v, ref, opts.CheckRef, stderr)
	if err != nil {
		return false, err
	}
	if !ok {
		return false, nil
	}
	if !opts.DoNotNormalize && needsLeftAlign(v) {
		if err := leftAlignInPlace(v, ref); err != nil {
			return false, err
		}
	}
	return true, nil
}

// normWindow is the bounded reorder buffer (the port of upstream's rbuf). It
// holds variants in (rid, pos) order via stable insertion — a left-aligned
// indel that shifts upstream is moved back to its sorted slot among same-contig
// neighbours, while records on different contigs keep their arrival order
// (upstream streams across contigs rather than globally reordering).
type normWindow struct {
	buf []*vcf.Variant
}

// push inserts v into the window, restoring (rid, pos) order within v's contig.
// The insertion is stable: records sharing a contig and position keep arrival
// order, and records on a different contig are never reordered past it.
func (w *normWindow) push(v *vcf.Variant) {
	i := len(w.buf)
	w.buf = append(w.buf, v)
	for i > 0 {
		prev := w.buf[i-1]
		if prev.Chrom != v.Chrom {
			break // do not reorder across contigs
		}
		if prev.Pos <= v.Pos {
			break
		}
		w.buf[i-1], w.buf[i] = w.buf[i], w.buf[i-1]
		i--
	}
}

// popSettled removes and returns the prefix of the window that can no longer be
// overtaken by a later record. It mirrors upstream normalize_vcf's flush count:
// when two contigs are buffered, every record on the first contig is settled;
// otherwise a record is settled once the last buffered record sits at least
// normBufWin bases ahead of it on the same contig.
func (w *normWindow) popSettled() []*vcf.Variant {
	if len(w.buf) == 0 {
		return nil
	}
	last := w.buf[len(w.buf)-1]
	first := w.buf[0]
	n := 0
	if first.Chrom != last.Chrom {
		// Flush everything on the first contig.
		for n < len(w.buf) && w.buf[n].Chrom == first.Chrom {
			n++
		}
	} else {
		for n < len(w.buf) && last.Pos-w.buf[n].Pos >= normBufWin {
			n++
		}
	}
	if n == 0 {
		return nil
	}
	// Copy the settled prefix out before compacting: reslicing the shared
	// backing array and then shifting the tail over it would corrupt the
	// returned records.
	out := make([]*vcf.Variant, n)
	copy(out, w.buf[:n])
	w.buf = append(w.buf[:0], w.buf[n:]...)
	return out
}

// drain returns and clears every remaining buffered record (the final flush).
func (w *normWindow) drain() []*vcf.Variant {
	out := w.buf
	w.buf = nil
	return out
}

// normSink is the stateful output stage. It receives records in final
// (rid, pos) order from the reorder window and applies the position-local
// transforms — duplicate removal, multiallelic join, and lax-mode filtering —
// before writing them. Because rmdup and join are both confined to a single
// CHROM+POS group upstream (the rmdup dup test resets at each new position and
// the mrows join buffer flushes when the position changes), buffering exactly
// one position group at a time reproduces the slice pipeline byte-for-byte
// while keeping memory bounded.
type normSink struct {
	w       variantWriter
	opts    NormOptions
	numbers headerNumberMaps

	group   []*vcf.Variant // records buffered at the current CHROM+POS
	written int
}

// newNormSink builds a sink writing through w.
func newNormSink(w variantWriter, opts NormOptions, numbers headerNumberMaps) *normSink {
	return &normSink{w: w, opts: opts, numbers: numbers}
}

// push adds v to the sink. Records are accumulated into the current CHROM+POS
// group; when v opens a new group the previous one is processed and emitted.
func (s *normSink) push(v *vcf.Variant) error {
	if len(s.group) > 0 {
		g := s.group[0]
		if g.Chrom != v.Chrom || g.Pos != v.Pos {
			if err := s.flushGroup(); err != nil {
				return err
			}
		}
	}
	s.group = append(s.group, v)
	return nil
}

// flush processes any pending group (the final drain).
func (s *normSink) flush() error {
	return s.flushGroup()
}

// flushGroup applies duplicate removal, multiallelic join, and lax-mode
// filtering to the buffered position group, then writes the survivors.
func (s *normSink) flushGroup() error {
	if len(s.group) == 0 {
		return nil
	}
	group := s.group
	s.group = nil

	// 6. duplicate removal.
	if s.opts.RmDup != RmDupNone {
		group = removeDuplicates(group, s.opts.RmDup)
	}
	// 7. join biallelics back into a multiallelic.
	if s.opts.Multiallelics.Active && !s.opts.Multiallelics.Split {
		group = joinMultiallelics(group, s.opts.Multiallelics, s.numbers)
	}
	// 8. lax-mode filtering happens after splitting.
	if !s.opts.StrictFilter && len(s.opts.ApplyFilters) > 0 {
		group = applyFilterList(group, s.opts.ApplyFilters)
	}
	for _, v := range group {
		if err := s.w.Write(v); err != nil {
			return err
		}
		s.written++
	}
	return nil
}

// applyFilterList keeps variants whose FILTER column matches one of the
// requested names (mirrors view's behaviour).
func applyFilterList(variants []*vcf.Variant, filters []string) []*vcf.Variant {
	out := variants[:0:0]
	for _, v := range variants {
		for _, want := range filters {
			for _, f := range v.Filter {
				if f == want {
					out = append(out, v)
					goto next
				}
			}
		}
	next:
	}
	return out
}

// classifyAlt returns true if the (ref, alt) pair is a SNP (both length 1
// and different) and true if it is an indel (lengths differ).
func classifyAlt(ref, alt string) (isSNP, isIndel bool) {
	if len(ref) == 1 && len(alt) == 1 && ref != alt {
		return true, false
	}
	if len(ref) != len(alt) {
		return false, true
	}
	return false, false
}

// variantIsSnp returns true if every ALT against REF is a SNP.
func variantIsSnp(v *vcf.Variant) bool {
	if len(v.Alt) == 0 {
		return false
	}
	for _, a := range v.Alt {
		s, _ := classifyAlt(v.Ref, a)
		if !s {
			return false
		}
	}
	return true
}

// variantIsIndel returns true if any ALT against REF is an indel.
func variantIsIndel(v *vcf.Variant) bool {
	for _, a := range v.Alt {
		_, i := classifyAlt(v.Ref, a)
		if i {
			return true
		}
	}
	return false
}

// splitMultiallelics expands records carrying multiple ALTs into one
// record per ALT. INFO/AC, INFO/AF and FORMAT/GT are adjusted per allele:
//
//   - INFO/AC and INFO/AF are split into their k-th component (R-style).
//   - FORMAT/GT alleles are remapped: the chosen ALT (index k+1) becomes
//     "1" and every other ALT becomes "0" so the new record looks like
//     a clean biallelic.
//   - INFO/AN, FORMAT/DP, etc. carry over unchanged — they describe the
//     site rather than an individual allele.
func splitMultiallelics(variants []*vcf.Variant, m MultiallelicMode, numbers headerNumberMaps) []*vcf.Variant {
	if !m.Split {
		return variants
	}
	out := make([]*vcf.Variant, 0, len(variants))
	for _, v := range variants {
		if !shouldAffect(v, m) || len(v.Alt) <= 1 {
			out = append(out, v)
			continue
		}
		for i, alt := range v.Alt {
			child := cloneVariant(v)
			child.Alt = []string{alt}
			child.Info = perAlleleInfo(v.Info, i, len(v.Alt))
			child.Samples = perAlleleSamples(v.Samples, v.Format, i)
			reindexSplitFields(child, v, i, numbers)
			out = append(out, child)
		}
	}
	return out
}

// reindexSplitFields re-indexes every Number=A/R/G INFO and FORMAT vector in the
// freshly-split biallelic record child (the i-th ALT, 0-based, of the original
// multiallelic src). It mirrors upstream vcfnorm.c split_info_* /
// split_format_* which keep, for the chosen ALT index ialt=i:
//
//   - Number=A: the single value vals[i];
//   - Number=R: [vals[0], vals[i+1]] (REF + chosen ALT);
//   - Number=G: the diploid genotypes (0/0),(0/i+1),(i+1/i+1), i.e. the values
//     at gt-indices 0, gt(0,i+1) and gt(i+1,i+1) (haploid: [vals[0], vals[i+1]]).
//
// The work is expressed as "drop every ALT except i" and delegated to
// subsetNumberedList / subsetGenotypeList (the bcf_remove_allele_set port), which
// already implement the A/R/G layouts identically. INFO/AC and INFO/AF are not
// touched here — perAlleleInfo already narrows those — but re-indexing them again
// would be a no-op (a single remaining value), so the two paths stay consistent.
func reindexSplitFields(child, src *vcf.Variant, i int, numbers headerNumberMaps) {
	nROri := 1 + len(src.Alt) // REF + all original ALTs
	if nROri <= 2 {
		return // already biallelic, nothing to re-index
	}
	// rm[j] flags the j-th allele (0=REF) for removal: keep REF and ALT i only.
	rm := make([]bool, nROri)
	for j := 1; j < nROri; j++ {
		rm[j] = j != i+1
	}

	// INFO fields.
	for tag, val := range child.Info {
		if tag == "AC" || tag == "AF" {
			continue // already narrowed by perAlleleInfo
		}
		num := numbers.info[tag]
		if num != "A" && num != "R" && num != "G" {
			continue
		}
		if nv, changed := subsetNumberedList(val, num, rm, nil, nROri); changed {
			child.Info[tag] = nv
		}
	}

	// FORMAT fields, per sample.
	for _, tag := range child.Format {
		num := numbers.format[tag]
		if num != "A" && num != "R" && num != "G" {
			continue
		}
		for s := range child.Samples {
			val, ok := child.Samples[s].Data[tag]
			if !ok {
				continue
			}
			if nv, changed := subsetNumberedList(val, num, rm, nil, nROri); changed {
				child.Samples[s].Data[tag] = nv
			}
		}
	}
}

// shouldAffect returns true when the multiallelic switch applies to v.
func shouldAffect(v *vcf.Variant, m MultiallelicMode) bool {
	if !m.Active {
		return false
	}
	isSnp := variantIsSnp(v)
	isIndel := variantIsIndel(v)
	if isSnp && m.Snps {
		return true
	}
	if isIndel && m.Indels {
		return true
	}
	return false
}

// perAlleleInfo returns a copy of info with AC / AF narrowed to the i-th
// allele (a la R-format). Unknown / unrelated tags pass through unchanged.
func perAlleleInfo(info map[string]string, i, n int) map[string]string {
	out := make(map[string]string, len(info))
	for k, v := range info {
		switch k {
		case "AC", "AF":
			parts := strings.Split(v, ",")
			if len(parts) == n && i < len(parts) {
				out[k] = parts[i]
				continue
			}
			// If the field doesn't match the allele count we leave it
			// alone — better than dropping data on a malformed record.
			out[k] = v
		default:
			out[k] = v
		}
	}
	return out
}

// perAlleleSamples returns a copy of samples with each GT remapped so the
// (i+1)-th allele of the original ALT becomes "1" and every other ALT is
// re-coded as "0". Heterozygous calls that include the chosen ALT become
// "0/1"; calls that don't reference it become "0/0".
func perAlleleSamples(samples []vcf.Sample, format []string, i int) []vcf.Sample {
	if len(samples) == 0 {
		return nil
	}
	wantedAllele := strconv.Itoa(i + 1)
	out := make([]vcf.Sample, len(samples))
	for j, s := range samples {
		ns := vcf.Sample{
			Name: s.Name,
			Data: make(map[string]string, len(s.Data)),
		}
		for k, v := range s.Data {
			if k == "GT" {
				ns.Data[k] = remapGT(v, wantedAllele)
				continue
			}
			ns.Data[k] = v
		}
		out[j] = ns
	}
	_ = format
	return out
}

// remapGT rewrites a genotype string so the wantedAllele (e.g. "2") is
// reported as "1" and any other non-zero allele becomes "0".
func remapGT(gt, wanted string) string {
	if gt == "" || gt == "." {
		return gt
	}
	var b strings.Builder
	b.Grow(len(gt))
	i := 0
	for i < len(gt) {
		// Skip allele separators verbatim.
		c := gt[i]
		if c == '/' || c == '|' {
			b.WriteByte(c)
			i++
			continue
		}
		// Capture the next allele token.
		j := i
		for j < len(gt) && gt[j] != '/' && gt[j] != '|' {
			j++
		}
		token := gt[i:j]
		switch token {
		case ".":
			b.WriteString(".")
		case "0":
			b.WriteString("0")
		case wanted:
			b.WriteString("1")
		default:
			// Any other non-zero allele becomes 0 in the per-allele view.
			b.WriteString("0")
		}
		i = j
	}
	return b.String()
}

// cloneVariant returns a deep enough copy of v that downstream mutation of
// Alt / Info / Samples won't leak back into the source slice.
func cloneVariant(v *vcf.Variant) *vcf.Variant {
	c := *v
	c.Alt = append([]string{}, v.Alt...)
	c.Filter = append([]string{}, v.Filter...)
	c.Format = append([]string{}, v.Format...)
	c.Info = make(map[string]string, len(v.Info))
	for k, val := range v.Info {
		c.Info[k] = val
	}
	c.InfoOrder = append([]string{}, v.InfoOrder...)
	c.Samples = make([]vcf.Sample, len(v.Samples))
	for i, s := range v.Samples {
		ns := vcf.Sample{Name: s.Name, Data: make(map[string]string, len(s.Data))}
		for k, val := range s.Data {
			ns.Data[k] = val
		}
		c.Samples[i] = ns
	}
	return &c
}

// joinMultiallelics joins biallelic records at the same position into
// multiallelics, mirroring upstream vcfnorm.c's mrows buffer (mrows_push /
// mrows_can_flush / mrows_flush). All records sharing CHROM+POS are buffered;
// within that window they are bucketed by variant-type category (or merged
// wholesale for the "any" mode). Records of a category that occurs only once,
// and records the mode does not target, pass through unchanged in their
// original relative order. A bucket of two or more records is merged via
// mergeBiallelicsToMultiallelic, which computes a common REF, remaps allele
// indices, and combines FORMAT/GT exactly as merge_format_genotype does.
func joinMultiallelics(variants []*vcf.Variant, m MultiallelicMode, numbers headerNumberMaps) []*vcf.Variant {
	out := make([]*vcf.Variant, 0, len(variants))
	i := 0
	for i < len(variants) {
		// Collect the window of records at this CHROM+POS (the flush window).
		j := i + 1
		for j < len(variants) && variants[j].Chrom == variants[i].Chrom && variants[j].Pos == variants[i].Pos {
			j++
		}
		out = append(out, joinWindow(variants[i:j], m, numbers)...)
		i = j
	}
	return out
}

// joinTypeBit returns the upstream bcf_get_variant_types category bit for a
// record participating in the join (a biallelic record contributes its single
// ALT's type). Records that are not biallelic, or whose type the mode does not
// target, are reported with ok=false so they pass through unmerged.
func joinTypeBit(v *vcf.Variant, m MultiallelicMode) (bit int, ok bool) {
	if len(v.Alt) != 1 {
		return 0, false
	}
	t := variantTypeBit(v.Ref, v.Alt[0])
	switch {
	case m.Any:
		// Any non-reference biallelic joins regardless of type.
		return t, true
	case m.Snps && m.Indels:
		// Default "both": every category is eligible, but kept in separate
		// buckets (SNP, MNP, INDEL, OTHER are distinct).
		return t, true
	case m.Snps:
		if t == vtSNP {
			return t, true
		}
	case m.Indels:
		if t == vtINDEL {
			return t, true
		}
	}
	return 0, false
}

// joinWindow joins the records at a single CHROM+POS. For COLLAPSE_ANY all
// eligible records are merged into one multiallelic; otherwise records are
// grouped by type category and only same-category runs of two or more are
// merged, preserving the input order of everything else (matching upstream's
// mrows_flush, which emits a single-record bucket verbatim).
func joinWindow(window []*vcf.Variant, m MultiallelicMode, numbers headerNumberMaps) []*vcf.Variant {
	if !m.Active || m.Split || len(window) == 1 {
		return window
	}
	if m.Any {
		eligible := make([]*vcf.Variant, 0, len(window))
		var passthrough []*vcf.Variant
		var firstIdx = -1
		result := make([]*vcf.Variant, 0, len(window))
		for idx, v := range window {
			if _, ok := joinTypeBit(v, m); ok {
				if firstIdx < 0 {
					firstIdx = len(result)
					result = append(result, nil) // placeholder for merged record
				}
				eligible = append(eligible, v)
			} else {
				passthrough = append(passthrough, v)
				result = append(result, v)
			}
			_ = idx
		}
		if len(eligible) == 0 {
			return window
		}
		if len(eligible) == 1 {
			result[firstIdx] = eligible[0]
		} else {
			result[firstIdx] = mergeBiallelicsToMultiallelic(eligible, m, numbers)
		}
		_ = passthrough
		return result
	}
	// Bucket by type category, preserving the first-seen order of categories
	// and of pass-through records.
	type bucket struct {
		bit     int
		records []*vcf.Variant
		anchor  int // position of this bucket's output slot
	}
	var buckets []*bucket
	byBit := map[int]*bucket{}
	result := make([]*vcf.Variant, 0, len(window))
	for _, v := range window {
		bit, ok := joinTypeBit(v, m)
		if !ok {
			result = append(result, v)
			continue
		}
		b := byBit[bit]
		if b == nil {
			b = &bucket{bit: bit, anchor: len(result)}
			byBit[bit] = b
			buckets = append(buckets, b)
			result = append(result, nil) // placeholder, filled below
		}
		b.records = append(b.records, v)
	}
	for _, b := range buckets {
		if len(b.records) == 1 {
			result[b.anchor] = b.records[0]
		} else {
			result[b.anchor] = mergeBiallelicsToMultiallelic(b.records, m, numbers)
		}
	}
	// Upstream emits the records at a position ordered by variant-type bit
	// ascending (VCF_SNP=1 < VCF_MNP=2 < VCF_INDEL=4 < VCF_OTHER=8), not in
	// input order — e.g. a SNP and an indel sharing a POS come out SNP-first
	// regardless of which appeared first in the input (vcfnorm.c mrows_flush).
	// A stable sort by the (first-ALT) type bit reproduces that while keeping
	// the within-type order the bucketing already fixed.
	sort.SliceStable(result, func(i, j int) bool {
		return joinSortBit(result[i]) < joinSortBit(result[j])
	})
	return result
}

// joinSortBit returns the variant-type bit used to order a record within a
// CHROM+POS group on norm-join output, mirroring upstream's type ordering. For
// a merged multiallelic record the first ALT's type is representative (every
// ALT in a bucket shares the type category).
func joinSortBit(v *vcf.Variant) int {
	if len(v.Alt) == 0 {
		return vtSNP
	}
	return variantTypeBit(v.Ref, v.Alt[0])
}

// mergeBiallelicsToMultiallelic merges a run of biallelic records (already
// determined to belong to one type bucket) into a single multiallelic record,
// porting upstream vcfnorm.c's merge_biallelics_to_multiallelic. It computes a
// common REF via mergeAlleles, takes the maximum QUAL, accumulates FILTERs,
// joins per-allele INFO/AC and INFO/AF, and combines FORMAT/GT through
// mergeFormatGenotype.
func mergeBiallelicsToMultiallelic(records []*vcf.Variant, m MultiallelicMode, numbers headerNumberMaps) *vcf.Variant {
	if len(records) == 1 {
		return records[0]
	}
	dst := cloneVariant(records[0])

	// ID: upstream concatenates the records' IDs with ';', preserving first-seen
	// order and skipping missing ("." / empty) and duplicate IDs (vcfnorm.c
	// merge_lines). If every record is missing an ID the result stays ".".
	var ids []string
	seenID := map[string]bool{}
	for _, rec := range records {
		if rec.ID == "" || rec.ID == "." {
			continue
		}
		for _, part := range strings.Split(rec.ID, ";") {
			if part == "" || part == "." || seenID[part] {
				continue
			}
			seenID[part] = true
			ids = append(ids, part)
		}
	}
	if len(ids) > 0 {
		dst.ID = strings.Join(ids, ";")
	}

	// Build the merged allele list and per-record allele-index maps.
	als := []string{records[0].Ref}
	als = append(als, records[0].Alt...)
	maps := make([][]int, len(records))
	maps[0] = make([]int, len(als))
	for i := range maps[0] {
		maps[0][i] = i
	}
	for i := 1; i < len(records); i++ {
		rcAlleles := append([]string{records[i].Ref}, records[i].Alt...)
		merged, mp, ok := mergeAlleles(rcAlleles, als)
		if !ok {
			// Cannot merge (incompatible REF prefixes): fall back to keeping
			// the records unmerged is not possible here, so leave dst as-is.
			return dst
		}
		als = merged
		maps[i] = mp
		// merging may have right-padded earlier maps' alleles by lengthening
		// the common REF; the allele *indices* are unchanged, so earlier maps
		// remain valid.
	}
	dst.Ref = als[0]
	dst.Alt = als[1:]

	// FILTER: accumulate the filters across all joined records, mirroring
	// vcfnorm.c merge_lines + htslib bcf_add_filter. The merged site starts
	// from the first record's FILTER and then unions in every later record's
	// non-PASS filters (deduplicated, first-seen order). A "PASS" entry on a
	// later record is skipped in the default (non-strict) mode, and adding any
	// real filter drops a lone "PASS" — so a PASS+q10 join becomes q10, while
	// a PASS+PASS join stays PASS.
	dst.Filter = mergeJoinFilters(records)

	// QUAL: take the maximum across records (missing QUAL is -1 in our model).
	dst.Qual = maxQual(records)

	// INFO/AC and INFO/AF: join per-allele values in the merged allele order
	// when present in every record (Number=A semantics). Other INFO tags keep
	// the first record's value (already copied by cloneVariant).
	joinPerAlleleInfo(dst, records, maps, len(dst.Alt))

	// Every other Number=A/R/G INFO and FORMAT field is re-indexed into the
	// merged allele order (AC/AF already handled above).
	joinReindexFields(dst, records, maps, numbers)

	// FORMAT/GT: combine genotypes via the upstream algorithm.
	if len(dst.Samples) > 0 && hasGT(dst.Format) {
		mergeFormatGenotype(dst, records, maps)
	}
	return dst
}

// joinReindexFields re-indexes the Number=A/R/G INFO and FORMAT vectors of the
// joined records into the merged allele order of dst, porting the per-tag
// scatter in upstream vcfnorm.c merge_info_numeric / merge_format_numeric (and
// their string variants). For each record i with allele map maps[i] (source
// allele index -> merged allele index):
//
//   - Number=A: source value k (0-based over the record's ALTs) lands at merged
//     ALT slot maps[i][k+1]-1.
//   - Number=R: source value k (0-based over REF+ALTs) lands at merged allele
//     slot maps[i][k].
//   - Number=G: the diploid value for source genotype (ia,ib) lands at merged
//     genotype index gt(maps[i][ia], maps[i][ib]); a haploid (per-allele) source
//     scatters like Number=R.
//
// AC and AF are skipped (joinPerAlleleInfo owns them). Slots no record fills
// stay ".", matching upstream's missing-fill.
func joinReindexFields(dst *vcf.Variant, records []*vcf.Variant, maps [][]int, numbers headerNumberMaps) {
	nAlt := len(dst.Alt)
	nR := nAlt + 1

	// INFO fields.
	infoTags := collectTags(records, true)
	for _, tag := range infoTags {
		if tag == "AC" || tag == "AF" {
			continue
		}
		num := numbers.info[tag]
		if num != "A" && num != "R" && num != "G" {
			continue
		}
		// Re-index only when every record carries the tag (matches upstream,
		// which errors otherwise; we conservatively leave dst's first-record
		// value alone on a partial set).
		if !allHaveInfo(records, tag) {
			continue
		}
		vals := make([][]string, len(records))
		for i, r := range records {
			vals[i] = strings.Split(r.Info[tag], ",")
		}
		merged, ok := scatterNumbered(num, vals, maps, nR)
		if ok {
			setInfo(dst, tag, strings.Join(merged, ","))
		}
	}

	// FORMAT fields, per sample.
	formatTags := collectFormatTags(records)
	for _, tag := range formatTags {
		if tag == "GT" {
			continue
		}
		num := numbers.format[tag]
		if num != "A" && num != "R" && num != "G" {
			continue
		}
		for s := range dst.Samples {
			vals := make([][]string, len(records))
			haveAll := true
			for i, r := range records {
				if s >= len(r.Samples) {
					haveAll = false
					break
				}
				val, ok := r.Samples[s].Data[tag]
				if !ok || val == "" || val == "." {
					haveAll = false
					break
				}
				vals[i] = strings.Split(val, ",")
			}
			if !haveAll {
				continue
			}
			merged, ok := scatterNumbered(num, vals, maps, nR)
			if ok {
				dst.Samples[s].Data[tag] = strings.Join(merged, ",")
			}
		}
	}
}

// scatterNumbered scatters the per-record value lists vals (one slice per joined
// record, comma-split) into the merged allele order described by maps, returning
// the merged comma-fields. number is "A", "R", or "G"; nR is the merged allele
// count (REF + ALTs). It returns ok=false when any record's element count does
// not match its declared cardinality (mirroring upstream's "could not merge"
// guard — we leave the field unchanged rather than emit a malformed vector).
func scatterNumbered(number string, vals [][]string, maps [][]int, nR int) (merged []string, ok bool) {
	switch number {
	case "A":
		merged = fillMissing(nR - 1)
		for i, v := range vals {
			nAlleles := len(maps[i])
			if len(v) != nAlleles-1 {
				return nil, false
			}
			for k := 0; k < len(v); k++ {
				dstIdx := maps[i][k+1] - 1
				if dstIdx >= 0 && dstIdx < len(merged) {
					merged[dstIdx] = v[k]
				}
			}
		}
		return merged, true
	case "R":
		merged = fillMissing(nR)
		for i, v := range vals {
			nAlleles := len(maps[i])
			if len(v) != nAlleles {
				return nil, false
			}
			for k := 0; k < len(v); k++ {
				dstIdx := maps[i][k]
				if dstIdx >= 0 && dstIdx < len(merged) {
					merged[dstIdx] = v[k]
				}
			}
		}
		return merged, true
	default: // "G"
		nG := nR * (nR + 1) / 2
		merged = fillMissing(nG)
		for i, v := range vals {
			nAlleles := len(maps[i])
			nGsrc := nAlleles * (nAlleles + 1) / 2
			if len(v) == nGsrc {
				// Diploid: source genotype (ia,ib), ib<=ia.
				k := 0
				for ia := 0; ia < nAlleles; ia++ {
					for ib := 0; ib <= ia; ib++ {
						if k >= len(v) {
							break
						}
						l := gtIndex(maps[i][ia], maps[i][ib])
						if l >= 0 && l < len(merged) {
							merged[l] = v[k]
						}
						k++
					}
				}
			} else if len(v) == nAlleles {
				// Haploid: per-allele, scatter like R.
				for k := 0; k < len(v); k++ {
					dstIdx := maps[i][k]
					if dstIdx >= 0 && dstIdx < len(merged) {
						merged[dstIdx] = v[k]
					}
				}
			} else {
				return nil, false
			}
		}
		return merged, true
	}
}

// fillMissing returns a slice of n "." placeholders.
func fillMissing(n int) []string {
	if n < 0 {
		n = 0
	}
	out := make([]string, n)
	for i := range out {
		out[i] = "."
	}
	return out
}

// allHaveInfo reports whether every record carries the named INFO tag.
func allHaveInfo(records []*vcf.Variant, tag string) bool {
	for _, r := range records {
		if v, ok := r.Info[tag]; !ok || v == "" || v == "." {
			return false
		}
	}
	return true
}

// collectTags returns the INFO tags present on records, in first-seen order.
// The info flag is reserved for future use; it is true when collecting INFO.
func collectTags(records []*vcf.Variant, info bool) []string {
	_ = info
	var out []string
	seen := map[string]bool{}
	for _, r := range records {
		for _, k := range r.InfoOrder {
			if !seen[k] {
				seen[k] = true
				out = append(out, k)
			}
		}
		// Fall back to map keys for records without a recorded order.
		for k := range r.Info {
			if !seen[k] {
				seen[k] = true
				out = append(out, k)
			}
		}
	}
	return out
}

// collectFormatTags returns the FORMAT tags present on records, in first-seen
// order across the records' Format slices.
func collectFormatTags(records []*vcf.Variant) []string {
	var out []string
	seen := map[string]bool{}
	for _, r := range records {
		for _, k := range r.Format {
			if !seen[k] {
				seen[k] = true
				out = append(out, k)
			}
		}
	}
	return out
}

// mergeJoinFilters combines the FILTER columns of the records joined into one
// multiallelic site, reproducing vcfnorm.c merge_lines together with htslib's
// bcf_add_filter set semantics (default, non-strict -m+ behaviour):
//
//   - The merged set starts from the first record's FILTER.
//   - For every later record, each of its filters is added in turn: a "PASS"
//     entry is skipped, and any non-PASS filter is unioned in (deduplicated,
//     first-seen order). Adding a real filter to a set that is exactly {PASS}
//     replaces the PASS — PASS and real filters are mutually exclusive.
//   - The missing value "." behaves like a normal filter token (it is not
//     special-cased by bcf_add_filter), so it is unioned like any other.
//
// The net effect: the join is "PASS" only when every record is PASS; otherwise
// it is the ordered union of all non-PASS filters.
func mergeJoinFilters(records []*vcf.Variant) []string {
	out := append([]string{}, records[0].Filter...)
	seen := map[string]bool{}
	for _, f := range out {
		seen[f] = true
	}
	addFilter := func(f string) {
		if seen[f] {
			return
		}
		if f == "PASS" {
			// Setting PASS clears everything else.
			out = []string{"PASS"}
			seen = map[string]bool{"PASS": true}
			return
		}
		if len(out) == 1 && out[0] == "PASS" {
			// Replace a lone PASS with the first real filter.
			out[0] = f
			delete(seen, "PASS")
			seen[f] = true
			return
		}
		out = append(out, f)
		seen[f] = true
	}
	for i := 1; i < len(records); i++ {
		for _, f := range records[i].Filter {
			if f == "PASS" {
				// Non-strict: a later PASS is ignored.
				continue
			}
			addFilter(f)
		}
	}
	return out
}

// maxQual returns the maximum QUAL across records, treating -1 as missing.
// When all are missing the result is -1 (rendered as "." downstream), matching
// upstream which sets a missing float when no record has a QUAL.
func maxQual(records []*vcf.Variant) float64 {
	best := -1.0
	haveBest := false
	for _, r := range records {
		if r.Qual < 0 {
			continue
		}
		if !haveBest || r.Qual > best {
			best = r.Qual
			haveBest = true
		}
	}
	if !haveBest {
		return -1
	}
	return best
}

// joinPerAlleleInfo joins INFO/AC and INFO/AF (Number=A) across records using
// the merged allele maps, mirroring upstream merge_info_numeric for VL_A tags:
// a value at source ALT index k lands at the merged ALT index maps[i][k+1]-1.
// The join is applied only when the tag is present in every record so that a
// partially-annotated set is left to the first record's value.
func joinPerAlleleInfo(dst *vcf.Variant, records []*vcf.Variant, maps [][]int, nAlt int) {
	for _, tag := range []string{"AC", "AF"} {
		present := true
		for _, r := range records {
			if _, ok := r.Info[tag]; !ok {
				present = false
				break
			}
		}
		if !present {
			continue
		}
		merged := make([]string, nAlt)
		for i := range merged {
			merged[i] = "."
		}
		for i, r := range records {
			parts := strings.Split(r.Info[tag], ",")
			for k, p := range parts {
				dstIdx := maps[i][k+1] - 1
				if dstIdx >= 0 && dstIdx < nAlt {
					merged[dstIdx] = p
				}
			}
		}
		setInfo(dst, tag, strings.Join(merged, ","))
	}
}

// mergeAlleles merges the new record's allele list a into the accumulated
// allele list b, porting htslib vcfmerge.c::merge_alleles verbatim (a is the
// incoming line, b is the running accumulator). It returns the new accumulator,
// a map from a's allele indices to accumulator indices, and ok=false when the
// REF prefixes are incompatible. The reference allele is always index 0.
//
// When the new record's REF (a[0]) is longer than the accumulator's REF, every
// accumulator allele (including the REF at index 0) is right-padded with a[0]'s
// suffix so the merged REF becomes the longer one. Conversely, when the
// accumulator's REF is longer, the new record's ALT alleles are padded before
// being matched/added.
func mergeAlleles(a, b []string) (merged []string, mp []int, ok bool) {
	mp = make([]int, len(a))
	mp[0] = 0
	rla := len(a[0])
	rlb := len(b[0])

	// Most common case: identical single-base SNPs.
	if len(a) == 2 && len(b) == 2 && rla == 1 && rlb == 1 && a[1] == b[1] && len(a[1]) == 1 && len(b[1]) == 1 {
		mp[1] = 1
		return b, mp, true
	}

	// REF prefixes must agree (case-insensitively).
	minRL := rla
	if rlb < minRL {
		minRL = rlb
	}
	if !strings.EqualFold(a[0][:minRL], b[0][:minRL]) {
		return nil, nil, false
	}

	out := append([]string{}, b...)

	// b alleles need expanding (incl. REF) when a's REF is longer.
	if rla > rlb {
		suffix := a[0][rlb:]
		for i := range out {
			if len(out[i]) > 0 && (out[i][0] == '<' || out[i][0] == '*') {
				continue // symbolic / overlapping-deletion alleles unchanged
			}
			out[i] = out[i] + suffix
		}
	}

	// Add a's ALT alleles, expanding them when b's REF is longer.
	for i := 1; i < len(a); i++ {
		ai := a[i]
		if rlb > rla && len(a[i]) > 0 && a[i][0] != '<' && a[i][0] != '*' {
			ai = a[i] + b[0][rla:]
		}
		found := -1
		for j := 1; j < len(out); j++ {
			if strings.EqualFold(ai, out[j]) {
				found = j
				break
			}
		}
		if found >= 0 {
			mp[i] = found
			continue
		}
		mp[i] = len(out)
		out = append(out, ai)
	}
	return out, mp, true
}

// mergeFormatGenotype combines the FORMAT/GT of a run of biallelic records
// into dst, porting upstream vcfnorm.c::merge_format_genotype. For each sample
// it starts from the first record's GT and folds in each subsequent record's
// alleles using the per-record allele map: a non-missing, non-reference source
// allele is written into the matching destination strand, never overwriting a
// non-reference allele already there — a conflict instead falls back to the
// first free (reference/missing) strand, which is why "0/1" + "0/1" joins to
// "2/1" rather than "0/2".
func mergeFormatGenotype(dst *vcf.Variant, records []*vcf.Variant, maps [][]int) {
	for s := range dst.Samples {
		gt := splitStrands(records[0].Samples[s].Data["GT"])
		for i := 1; i < len(records); i++ {
			if s >= len(records[i].Samples) {
				continue
			}
			gt2 := splitStrands(records[i].Samples[s].Data["GT"])
			for k2 := range gt2 {
				src := gt2[k2].allele
				if src == "." || src == "" {
					continue // don't overwrite with missing
				}
				ial2, err := strconv.Atoi(src)
				if err != nil {
					continue
				}
				// Don't overwrite with ref unless the destination strand is
				// missing.
				if ial2 == 0 {
					if k2 < len(gt) && gt[k2].allele != "." && gt[k2].allele != "" {
						continue
					}
				}
				ial := ial2
				if ial2 < len(maps[i]) {
					ial = maps[i][ial2]
				}
				dstStr := strconv.Itoa(ial)
				if k2 < len(gt) {
					cur := gt[k2].allele
					if cur == "." || cur == "" || cur == "0" {
						// Preserve phasing if the strand was phased.
						gt[k2].allele = dstStr
						continue
					}
				}
				// Conflict: find the first free (ref/missing/end) strand.
				placed := false
				for k := range gt {
					if gt[k].allele == "." || gt[k].allele == "" || gt[k].allele == "0" {
						gt[k].allele = dstStr
						placed = true
						break
					}
				}
				if !placed && k2 >= len(gt) {
					// Ploidy grew: append a new strand.
					gt = append(gt, strand{sep: '/', allele: dstStr})
				}
			}
		}
		dst.Samples[s].Data["GT"] = joinStrands(gt)
	}
}

// strand pairs an allele with the separator that preceded it (so we can
// round-trip phased and unphased genotypes verbatim).
type strand struct {
	sep    byte
	allele string
}

func splitStrands(gt string) []strand {
	if gt == "" {
		return nil
	}
	var out []strand
	i := 0
	for i < len(gt) {
		var sep byte
		if i > 0 {
			sep = gt[i]
			i++
		}
		j := i
		for j < len(gt) && gt[j] != '/' && gt[j] != '|' {
			j++
		}
		out = append(out, strand{sep: sep, allele: gt[i:j]})
		i = j
	}
	return out
}

func joinStrands(s []strand) string {
	var b strings.Builder
	for i, st := range s {
		if i > 0 {
			b.WriteByte(st.sep)
		}
		b.WriteString(st.allele)
	}
	return b.String()
}

func hasGT(format []string) bool {
	for _, f := range format {
		if f == "GT" {
			return true
		}
	}
	return false
}

// atomizeVariants decomposes complex variants (REF and ALT both >1bp and
// equal length, e.g. "ACG" → "AGT") into a sequence of single-base
// substitutions. Length-changing complex variants are left intact because
// atomizing them safely requires a true alignment which is beyond the
// scope of this slice.
func atomizeVariants(variants []*vcf.Variant) []*vcf.Variant {
	out := make([]*vcf.Variant, 0, len(variants))
	for _, v := range variants {
		if len(v.Alt) != 1 || len(v.Ref) <= 1 || len(v.Ref) != len(v.Alt[0]) {
			out = append(out, v)
			continue
		}
		alt := v.Alt[0]
		for k := 0; k < len(v.Ref); k++ {
			if v.Ref[k] == alt[k] {
				continue
			}
			child := cloneVariant(v)
			child.Pos = v.Pos + k
			child.Ref = string(v.Ref[k])
			child.Alt = []string{string(alt[k])}
			out = append(out, child)
		}
	}
	return out
}

// fixRef rewrites v.Ref (and, for a simple allele swap, the ALT list and the
// per-sample GT) so the record agrees with the FASTA. It implements the
// common branches of upstream vcfnorm.c::fix_ref:
//
//   - REF already matches the reference: no change.
//   - REF is the missing marker ".": set it to the reference base.
//   - exactly one ALT equals the reference span: swap REF <-> that ALT and
//     re-index the genotypes (a REF/ALT swap).
//   - otherwise: replace REF with the reference span verbatim.
//
// Symbolic ALTs (<DEL> etc.) and the N-fill / IUPAC branches are out of
// scope here; such records are left untouched (the subsequent check-ref pass
// still reports them per the active policy).
func fixRef(v *vcf.Variant, ref *fasta.RandomAccess) error {
	if len(v.Ref) == 0 {
		return nil
	}
	for _, a := range v.Alt {
		if a != "" && a[0] == '<' {
			return nil // symbolic; leave for the check-ref pass
		}
	}
	refSeq, err := ref.Fetch(v.Chrom, int64(v.Pos-1), int64(v.Pos-1+len(v.Ref)))
	if err != nil {
		return nil // fetch issues are surfaced by checkRef
	}
	want := strings.ToUpper(string(refSeq))
	if want == strings.ToUpper(v.Ref) {
		return nil
	}
	// Missing REF marker: fill in the reference base.
	if v.Ref == "." {
		v.Ref = want
		return nil
	}
	// Simple swap: one ALT equals the reference span.
	for i, alt := range v.Alt {
		if strings.EqualFold(alt, want) {
			swapRefAlt(v, i)
			return nil
		}
	}
	// Otherwise set REF to the reference span.
	v.Ref = want
	return nil
}

// swapRefAlt swaps v.Ref with v.Alt[idx] and re-indexes every sample GT so
// allele 0 becomes (idx+1) and (idx+1) becomes 0, matching upstream's
// fix_ref genotype swap for a REF/ALT inversion.
func swapRefAlt(v *vcf.Variant, idx int) {
	oldRef := v.Ref
	v.Ref = v.Alt[idx]
	v.Alt[idx] = oldRef
	swapTo := idx + 1
	for si := range v.Samples {
		gt, ok := v.Samples[si].Data["GT"]
		if !ok {
			continue
		}
		v.Samples[si].Data["GT"] = swapGTAlleles(gt, swapTo)
	}
}

// swapGTAlleles rewrites a GT string, swapping allele index 0 with index n
// (and vice versa) while preserving phasing separators and missing alleles.
func swapGTAlleles(gt string, n int) string {
	var b strings.Builder
	cur := strings.Builder{}
	flush := func() {
		s := cur.String()
		cur.Reset()
		switch s {
		case "0":
			b.WriteString(strconv.Itoa(n))
		case strconv.Itoa(n):
			b.WriteString("0")
		default:
			b.WriteString(s)
		}
	}
	for i := 0; i < len(gt); i++ {
		c := gt[i]
		if c == '/' || c == '|' {
			flush()
			b.WriteByte(c)
			continue
		}
		cur.WriteByte(c)
	}
	flush()
	return b.String()
}

// needsLeftAlign returns true when the variant has at least one ALT whose
// length differs from REF (a candidate for shifting left).
func needsLeftAlign(v *vcf.Variant) bool {
	if len(v.Ref) > 1 {
		return true
	}
	for _, a := range v.Alt {
		if len(a) > 1 {
			return true
		}
	}
	return false
}

// checkRef compares v.Ref to the FASTA and applies the --check-ref policy
// bitmask. It returns (true, nil) when the record is kept (it matches, or
// only the warn bit is set) and (false, nil) when the skip bit drops it. An
// error is returned for the default exit policy on a mismatch.
//
// The stderr warning and the fatal error reproduce upstream's wording
// (vcfnorm.c): the warning is "REF_MISMATCH\t<chr>\t<pos>\t<vcf-ref>\t
// <ref-seq>" and the error is "Reference allele mismatch at <chr>:<pos> ..
// REF_SEQ:'<ref-seq>' vs VCF:'<vcf-ref>'".
func checkRef(v *vcf.Variant, ref *fasta.RandomAccess, mode CheckRefMode, stderr io.Writer) (bool, error) {
	refSeq, err := ref.Fetch(v.Chrom, int64(v.Pos-1), int64(v.Pos-1+len(v.Ref)))
	if err != nil {
		// Missing contig / out-of-range fetch: surface in exit mode,
		// drop in skip mode, warn otherwise.
		switch {
		case mode == CheckRefExit:
			return false, fmt.Errorf("bcftools norm: REF lookup %s:%d failed: %w", v.Chrom, v.Pos, err)
		case mode&CheckRefSkip != 0:
			return false, nil
		default:
			if mode&CheckRefWarn != 0 && stderr != nil {
				fmt.Fprintf(stderr, "bcftools norm: REF lookup %s:%d failed: %v (kept)\n", v.Chrom, v.Pos, err)
			}
			return true, nil
		}
	}
	if strings.EqualFold(string(refSeq), v.Ref) {
		return true, nil
	}
	// Mismatch. The default (exit) is exclusive and overrides the rest.
	if mode == CheckRefExit {
		return false, fmt.Errorf("Reference allele mismatch at %s:%d .. REF_SEQ:'%s' vs VCF:'%s'",
			v.Chrom, v.Pos, refSeq, v.Ref)
	}
	if mode&CheckRefWarn != 0 && stderr != nil {
		fmt.Fprintf(stderr, "REF_MISMATCH\t%s\t%d\t%s\t%s\n", v.Chrom, v.Pos, v.Ref, refSeq)
	}
	if mode&CheckRefSkip != 0 {
		return false, nil
	}
	return true, nil
}

// leftAlignInPlace shifts a variant left using the Tan-Abecasis-Durbin
// algorithm that bcftools / vt / GATK use:
//
//	repeat:
//	  if alleles all end in the same byte AND not all are length 1:
//	    trim the trailing byte from every allele
//	  if any allele is now empty:
//	    prepend the upstream reference base
//	  continue until no change is made.
//
// The result is guaranteed to be parsimonious and left-aligned: no shared
// suffix bases, no shared leading bases beyond the single anchor required
// for VCF representation, position as small as possible.
func leftAlignInPlace(v *vcf.Variant, ref *fasta.RandomAccess) error {
	alleles := make([]string, 0, 1+len(v.Alt))
	alleles = append(alleles, v.Ref)
	alleles = append(alleles, v.Alt...)
	pos := v.Pos
	for {
		changed := false
		// 1) Trim a shared trailing base when all alleles match in the
		// last position and they aren't all single-byte (we never trim
		// the only base out of an allele on this pass alone).
		if last, ok := commonLastByte(alleles); ok && !allLengthOne(alleles) {
			_ = last
			for i := range alleles {
				alleles[i] = alleles[i][:len(alleles[i])-1]
			}
			changed = true
		}
		// 2) If trimming produced an empty allele, prepend the
		// upstream reference base — this is the actual "shift left"
		// step.
		needPrepend := false
		for _, a := range alleles {
			if len(a) == 0 {
				needPrepend = true
				break
			}
		}
		if needPrepend {
			if pos <= 1 {
				return fmt.Errorf("bcftools norm: cannot left-align past chrom start at %s:%d", v.Chrom, v.Pos)
			}
			upstream, err := ref.Fetch(v.Chrom, int64(pos-2), int64(pos-1))
			if err != nil {
				return err
			}
			if len(upstream) != 1 {
				return fmt.Errorf("bcftools norm: bad upstream fetch on %s:%d", v.Chrom, pos)
			}
			for i := range alleles {
				alleles[i] = string(upstream) + alleles[i]
			}
			pos--
			changed = true
		}
		if !changed {
			break
		}
	}
	// 3) Trim a shared leading base when every allele has length > 1.
	// This collapses cases like (AACC, AAC) → (ACC, AC) into the
	// minimum representation. We always keep at least one base in every
	// allele so the VCF representation stays well-formed.
	for {
		if !allLongerThanOne(alleles) {
			break
		}
		first := alleles[0][0]
		same := true
		for _, a := range alleles[1:] {
			if !equalASCII(a[0], first) {
				same = false
				break
			}
		}
		if !same {
			break
		}
		for i := range alleles {
			alleles[i] = alleles[i][1:]
		}
		pos++
	}
	v.Pos = pos
	v.Ref = alleles[0]
	v.Alt = alleles[1:]
	return nil
}

// allLengthOne returns true when every allele is exactly one base. This is
// the terminating condition for the trim-trailing step.
func allLengthOne(alleles []string) bool {
	for _, a := range alleles {
		if len(a) != 1 {
			return false
		}
	}
	return true
}

// allLongerThanOne returns true when every allele has length >= 2.
func allLongerThanOne(alleles []string) bool {
	for _, a := range alleles {
		if len(a) < 2 {
			return false
		}
	}
	return true
}

// commonLastByte returns the last byte if every allele ends in the same
// letter, and false otherwise. Case-insensitive comparison keeps mixed-
// case FASTAs working.
func commonLastByte(alleles []string) (byte, bool) {
	last := alleles[0][len(alleles[0])-1]
	for _, a := range alleles[1:] {
		if !equalASCII(a[len(a)-1], last) {
			return 0, false
		}
	}
	return last, true
}

// equalASCII returns true when two ASCII bytes match ignoring case.
func equalASCII(a, b byte) bool {
	if a >= 'a' && a <= 'z' {
		a -= 'a' - 'A'
	}
	if b >= 'a' && b <= 'z' {
		b -= 'a' - 'A'
	}
	return a == b
}

// removeDuplicates filters variants per the --rm-dup policy. Apart from
// "exact", entries are matched on chrom+pos and (optionally) type.
func removeDuplicates(variants []*vcf.Variant, mode RmDupMode) []*vcf.Variant {
	if mode == RmDupNone {
		return variants
	}
	out := make([]*vcf.Variant, 0, len(variants))
	if mode == RmDupExact {
		seen := make(map[string]struct{})
		for _, v := range variants {
			key := exactKey(v)
			if _, dup := seen[key]; dup {
				continue
			}
			seen[key] = struct{}{}
			out = append(out, v)
		}
		return out
	}
	// For the non-exact modes, only the first record of a (chrom,pos)
	// pair (optionally narrowed by type) is retained.
	type seenKey struct {
		chrom string
		pos   int
		kind  string
	}
	seen := make(map[seenKey]struct{})
	for _, v := range variants {
		kind := "other"
		if variantIsSnp(v) {
			kind = "snp"
		}
		if variantIsIndel(v) {
			kind = "indel"
		}
		affect := false
		groupKind := kind
		switch mode {
		case RmDupSnps:
			affect = kind == "snp"
		case RmDupIndels:
			affect = kind == "indel"
		case RmDupBoth:
			affect = kind == "snp" || kind == "indel"
		case RmDupAll:
			affect = true
			groupKind = "all" // collapse every type into one bucket
		}
		if !affect {
			out = append(out, v)
			continue
		}
		key := seenKey{v.Chrom, v.Pos, groupKind}
		if _, dup := seen[key]; dup {
			continue
		}
		seen[key] = struct{}{}
		out = append(out, v)
	}
	return out
}

// exactKey returns a string used to detect byte-for-byte duplicate records
// in the "exact" mode.
func exactKey(v *vcf.Variant) string {
	return strings.Join([]string{
		v.Chrom,
		strconv.Itoa(v.Pos),
		v.Ref,
		strings.Join(v.Alt, ","),
	}, "\x00")
}
