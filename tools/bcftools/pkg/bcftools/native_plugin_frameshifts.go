// Native port of the upstream `frameshifts` plugin (plugins/frameshifts.c):
// annotate indels that overlap a list of exons with INFO/OOF.
//
// The plugin reads exons from -e/--exons (a BED or region-list file, optionally
// bgzipped/tabixed) into htslib's bcf_sr_regions_t cursor and, for every indel
// record that overlaps an exon, adds INFO/OOF with one Integer per ALT allele.
//
// IMPORTANT — upstream behaviour reproduced byte-for-byte (see
// docs/UPSTREAM_BUGS.md#bcftools-frameshifts-oof-dead-code): the per-allele
// in-frame/out-of-frame computation is DEAD CODE in the shipped binary. It
// guards on `rec->d.var[i].type != VCF_INDEL`, but htslib's bcf_set_variant_type
// sets the per-allele type to VCF_INDEL|VCF_INS or VCF_INDEL|VCF_DEL (never the
// bare VCF_INDEL bit that flag was when frameshifts was written in 2014). So the
// guard is ALWAYS true and every indel allele that reaches the loop is annotated
// OOF=-1 ("not applicable"). The default native port matches this exactly so the
// CLI-to-CLI oracle against the real 1.23.1 binary passes. The corrected
// in-frame/out-of-frame computation (the exon-trim + length-mod-3 logic the
// upstream source intended) is implemented and unit-tested as a pure helper and
// is available via the non-default --fix-oof flag.
//
// The cursor is a faithful port of bcf_sr_regions_overlap: a monotonic forward
// cursor over per-chromosome sorted+merged 0-based regions, re-seeking when the
// query position goes backwards or to a new chromosome. The plugin is run as a
// serial bufferedPlugin because the cursor state is shared across records (and a
// backwards record forces a re-seek exactly as upstream).
package bcftools

import (
	"bufio"
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("frameshifts", func() NativePlugin { return &frameshiftsPlugin{} })
}

// frameshiftsPlugin implements the frameshifts plugin.
type frameshiftsPlugin struct {
	exonsFile string
	exons     *exonCursor
	fixOOF    bool // --fix-oof: compute the real in-frame/out-of-frame value
	hdr       *vcf.Header
}

// Name returns the plugin name.
func (p *frameshiftsPlugin) Name() string { return "frameshifts" }

// About returns the one-line description, matching frameshifts.c about().
func (p *frameshiftsPlugin) About() string { return "Annotate frameshift indels." }

// Init parses the plugin options (the generic init/process form: options after
// the `--`), loads the exon cursor, and appends the OOF INFO definition.
func (p *frameshiftsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.hdr = hdr
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("frameshifts: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-e", "--exons":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.exonsFile = v
		case "--fix-oof":
			p.fixOOF = true
		case "--no-version":
			// provenance handled elsewhere in the native path
		default:
			return nil, fmt.Errorf("frameshifts: unsupported option %q", a)
		}
	}
	if p.exonsFile == "" {
		return nil, fmt.Errorf("frameshifts: missing the -e option")
	}

	cur, err := loadExonCursor(p.exonsFile)
	if err != nil {
		return nil, fmt.Errorf("frameshifts: error occurred while reading %q: %w", p.exonsFile, err)
	}
	p.exons = cur

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	out.MetaInfo = appendInfoHeader(out.MetaInfo,
		`##INFO=<ID=OOF,Number=A,Type=Integer,Description="Frameshift Indels: out-of-frame (1), in-frame (0), not-applicable (-1 or missing)">`)
	return out, nil
}

// ProcessAll annotates each indel record in order, porting process(). It is
// buffered (serial) because the exon cursor carries state between records.
func (p *frameshiftsPlugin) ProcessAll(variants []*vcf.Variant) ([]*vcf.Variant, error) {
	for _, rec := range variants {
		p.annotate(rec)
	}
	return variants, nil
}

// annotate applies the OOF computation to a single record, mutating it in place.
// It mirrors process(): non-variants and non-indels are passed through; an indel
// that does not overlap any exon is passed through; otherwise OOF (one value per
// ALT) is added.
func (p *frameshiftsPlugin) annotate(rec *vcf.Variant) {
	if len(rec.Alt) < 1 {
		return // not a variant
	}
	types := recordVariantTypes(rec)
	if types&vtINDEL == 0 {
		return // not an indel
	}

	// len = the most-negative per-allele n (the largest deletion), porting
	// `for (i...) if (len > rec->d.var[i].n) len = rec->d.var[i].n;`.
	minN := 0
	for i := range rec.Alt {
		n := alleleVariant(rec.Ref, rec.Alt[i]).n
		if minN > n {
			minN = n
		}
	}
	pos0 := rec.Pos - 1 // 0-based, matching rec->pos
	posTo := pos0
	if minN == 0 {
		posTo = pos0 - minN // == pos0 (minN==0)
	}

	if !p.exons.overlap(rec.Chrom, pos0, posTo) {
		return // no overlap
	}

	frm := make([]int, len(rec.Alt))
	for i := range rec.Alt {
		frm[i] = oofForAllele(rec.Ref, rec.Alt[i], pos0, p.exons.start, p.exons.end, p.fixOOF)
	}

	parts := make([]string, len(frm))
	for i, f := range frm {
		parts[i] = strconv.Itoa(f)
	}
	setInfo(rec, "OOF", strings.Join(parts, ","))
}

// Process is never called for a bufferedPlugin but satisfies NativePlugin.
func (p *frameshiftsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *frameshiftsPlugin) Destroy() error { return nil }

// --- pure helpers (also exercised by the binary-free unit tests) -----------

// oofForAllele computes the OOF value for one ALT allele, porting the per-allele
// body of process(). When fix is false it reproduces upstream's shipped (buggy)
// behaviour: the per-allele type guard `var[i].type != VCF_INDEL` is always true
// (htslib sets VCF_INDEL|VCF_INS / VCF_INDEL|VCF_DEL), so every allele yields -1.
// When fix is true it performs the intended exon-trim + length-mod-3 computation:
// the number of inserted/deleted bases that fall inside the exon is taken mod 3
// (out-of-frame => 1, in-frame => 0, none inside => -1).
func oofForAllele(ref, alt string, pos0, exStart, exEnd int, fix bool) int {
	av := alleleVariant(ref, alt)
	if !fix {
		// Shipped upstream behaviour: var[i].type is VCF_INDEL|VCF_INS or
		// VCF_INDEL|VCF_DEL, never the bare VCF_INDEL the source compares against,
		// so the guard `type != VCF_INDEL` is always true and frm[i-1] = -1.
		return -1
	}
	if !av.isIndel {
		return -1
	}
	length := av.n
	tlen := 0
	if length > 0 {
		// insertion
		if exStart <= pos0 && exEnd > pos0 {
			tlen = abs(length)
		}
	} else if exStart <= pos0+abs(length) {
		// deletion
		tlen = abs(length)
		if pos0 < exStart { // trim the beginning
			tlen -= exStart - pos0 + 1
		}
		if exEnd < pos0+abs(length) { // trim the end
			tlen -= pos0 + abs(length) - exEnd
		}
	}
	if tlen != 0 {
		if tlen%3 != 0 {
			return 1 // out-of-frame
		}
		return 0 // in-frame
	}
	return -1 // not applicable (outside)
}

// alleleInfo holds the htslib-equivalent (n, type-is-indel) classification of a
// single REF/ALT pair.
type alleleInfo struct {
	n       int  // signed indel length, matching bcf_variant_t.n
	isIndel bool // whether the per-allele type carries VCF_INDEL
}

// alleleVariant ports bcf_set_variant_type for a single REF/ALT pair, returning
// the signed length n and whether the allele is an indel. Only the fields
// frameshifts consumes (n and the VCF_INDEL bit) are computed; the case-folding
// and prefix/suffix-trimming match htslib byte-for-byte.
func alleleVariant(ref, alt string) alleleInfo {
	if alt == "*" {
		return alleleInfo{n: 0} // VCF_OVERLAP
	}
	// The most frequent case: single-base REF and ALT.
	if len(ref) == 1 && len(alt) == 1 {
		if alt == "." || ref[0] == alt[0] || alt == "X" {
			return alleleInfo{n: 0} // VCF_REF
		}
		return alleleInfo{n: 1} // VCF_SNP
	}
	if len(alt) > 0 && alt[0] == '<' {
		if alt == "<X>" || alt == "<*>" || alt == "<NON_REF>" {
			return alleleInfo{n: 0} // VCF_REF
		}
		return alleleInfo{} // VCF_OTHER
	}
	if len(alt) > 0 && (alt[0] == ']' || alt[0] == '[') {
		return alleleInfo{} // VCF_BND
	}

	// Iterate through alt characters that match the reference (case-insensitive).
	r, a := 0, 0
	for r < len(ref) && a < len(alt) && upperByte(ref[r]) == upperByte(alt[a]) {
		r++
		a++
	}
	if a < len(alt) && r == len(ref) {
		// consume to end of alt
		ae := len(alt)
		if alt[ae-1] == ']' || alt[ae-1] == '[' {
			return alleleInfo{} // VCF_BND
		}
		return alleleInfo{n: ae - a, isIndel: true} // VCF_INDEL|VCF_INS
	}
	if r < len(ref) && a == len(alt) {
		re := len(ref)
		return alleleInfo{n: (a - len(alt)) - (re - r), isIndel: true} // VCF_INDEL|VCF_DEL
	}
	if r == len(ref) && a == len(alt) {
		return alleleInfo{n: 0} // VCF_REF
	}

	// Trim matching suffix.
	re, ae := len(ref)-1, len(alt)-1
	if alt[ae] == ']' || alt[ae] == '[' {
		return alleleInfo{} // VCF_BND
	}
	for re > r && ae > a && upperByte(ref[re]) == upperByte(alt[ae]) {
		re--
		ae--
	}
	if ae == a {
		if re == r {
			return alleleInfo{n: 1} // VCF_SNP
		}
		n := -(re - r)
		if upperByte(ref[re]) == upperByte(alt[ae]) {
			return alleleInfo{n: n, isIndel: true} // VCF_INDEL|VCF_DEL
		}
		return alleleInfo{n: n} // VCF_OTHER
	}
	if re == r {
		n := ae - a
		if upperByte(ref[re]) == upperByte(alt[ae]) {
			return alleleInfo{n: n, isIndel: true} // VCF_INDEL|VCF_INS
		}
		return alleleInfo{n: n} // VCF_OTHER
	}

	// MNP / OTHER
	var n int
	if re-r > ae-a {
		n = -(re - r + 1)
	} else {
		n = ae - a + 1
	}
	return alleleInfo{n: n} // VCF_MNP or VCF_OTHER (no VCF_INDEL bit)
}

// recordVariantTypes returns the bitwise-OR of the per-allele types for a
// record, restricted to the VCF_INDEL bit that frameshifts tests against (it
// only ever asks `type & VCF_INDEL`). It ports bcf_get_variant_types' OR.
func recordVariantTypes(rec *vcf.Variant) int {
	mask := 0
	for i := range rec.Alt {
		if alleleVariant(rec.Ref, rec.Alt[i]).isIndel {
			mask |= vtINDEL
		}
	}
	return mask
}

// --- exon cursor (port of htslib bcf_sr_regions_t) --------------------------

// exonRegion is one 0-based, inclusive-end exon interval (start = from-1,
// end = to-1 in htslib terms).
type exonRegion struct {
	start int
	end   int
}

// exonCursor is a monotonic forward cursor over per-chromosome sorted+merged
// exon regions, porting the subset of bcf_sr_regions_t that frameshifts uses
// (bcf_sr_regions_overlap and its cursor state start/end/prev_seq/prev_start).
type exonCursor struct {
	byChrom map[string][]exonRegion
	order   []string // chromosome order as first seen (for the iseq monotonic test)
	iseqOf  map[string]int

	// cursor state, mirroring reg->start/end/iseq/creg/prev_seq/prev_start.
	start    int // current region start (0-based), -1 when before first
	end      int // current region end (0-based inclusive), -1 when before first
	iseq     int // current chromosome index, -1 when exhausted on the chromosome
	creg     int // index into byChrom[curSeq] of the current region
	curSeq   string
	prevSeq  int
	prevStrt int
}

// loadExonCursor reads the exon file into a cursor, porting the BED/region-list
// autodetection (.bed/.bed.gz => 0-based half-open) and the
// _regions_add/_regions_sort_and_merge pipeline (store 0-based, sort, merge
// overlapping/adjacent regions).
func loadExonCursor(path string) (*exonCursor, error) {
	isBED := strings.HasSuffix(path, ".bed") || strings.HasSuffix(path, ".bed.gz")
	f, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	cur := &exonCursor{
		byChrom: map[string][]exonRegion{},
		iseqOf:  map[string]int{},
	}
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	any := false
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("could not parse line %q, expected columns chr,beg[,end]", line)
		}
		chr := fields[0]
		from, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("bad start in %q", line)
		}
		to := from
		if len(fields) >= 3 {
			to, err = strconv.Atoi(fields[2])
			if err != nil {
				return nil, fmt.Errorf("bad end in %q", line)
			}
		}
		if isBED {
			from++ // _regions_next: `if (is_bed) from++`
		}
		// _regions_add stores 0-based: start = from-1, end = to-1.
		cur.add(chr, from-1, to-1)
		any = true
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if !any {
		return nil, fmt.Errorf("no regions in %q", path)
	}
	cur.sortAndMerge()
	cur.reset()
	return cur, nil
}

// add appends a 0-based region to a chromosome, recording the chromosome's
// first-seen order (the iseq the synced reader assigns).
func (c *exonCursor) add(chr string, start, end int) {
	if _, ok := c.iseqOf[chr]; !ok {
		c.iseqOf[chr] = len(c.order)
		c.order = append(c.order, chr)
	}
	c.byChrom[chr] = append(c.byChrom[chr], exonRegion{start: start, end: end})
}

// sortAndMerge ports _regions_sort_and_merge: per chromosome, sort by (start,end)
// then merge overlapping/adjacent (end >= next.start) regions.
func (c *exonCursor) sortAndMerge() {
	for chr, regs := range c.byChrom {
		sort.SliceStable(regs, func(i, j int) bool {
			if regs[i].start != regs[j].start {
				return regs[i].start < regs[j].start
			}
			return regs[i].end < regs[j].end
		})
		merged := regs[:0]
		for _, r := range regs {
			if len(merged) > 0 && merged[len(merged)-1].end >= r.start {
				if merged[len(merged)-1].end < r.end {
					merged[len(merged)-1].end = r.end
				}
				continue
			}
			merged = append(merged, r)
		}
		c.byChrom[chr] = merged
	}
}

// reset returns the cursor to the pre-first state, matching a freshly
// initialised bcf_sr_regions_t.
func (c *exonCursor) reset() {
	c.start, c.end = -1, -1
	c.iseq = -1
	c.creg = -1
	c.curSeq = ""
	c.prevSeq = -1
	c.prevStrt = -1
}

// seek positions the cursor at the start of chromosome seq, porting
// bcf_sr_regions_seek. It returns false if the chromosome has no regions.
func (c *exonCursor) seek(seq string) bool {
	c.iseq, c.start, c.end, c.creg = -1, -1, -1, -1
	idx, ok := c.iseqOf[seq]
	if !ok {
		return false
	}
	c.iseq = idx
	c.curSeq = seq
	return true
}

// next advances to the next region on the current chromosome, porting
// bcf_sr_regions_next. It returns false when no regions remain.
func (c *exonCursor) next() bool {
	c.start, c.end = -1, -1
	if c.iseq < 0 || c.curSeq == "" {
		return false
	}
	regs := c.byChrom[c.curSeq]
	if c.creg+1 >= len(regs) {
		c.iseq = -1 // exhausted this chromosome
		return false
	}
	c.creg++
	c.start = regs[c.creg].start
	c.end = regs[c.creg].end
	return true
}

// overlap ports bcf_sr_regions_overlap: it advances the cursor forward to the
// first region whose end >= start and reports whether that region's start <= end
// (i.e. the query [start,end] overlaps an exon). start/end are 0-based; after a
// true result reg.start/reg.end hold the overlapping exon, exactly as the
// frameshifts per-allele body reads exons->start/exons->end.
func (c *exonCursor) overlap(seq string, start, end int) bool {
	iseq, ok := c.iseqOf[seq]
	if !ok {
		return false // no such sequence
	}
	if c.prevSeq == -1 || iseq != c.prevSeq || c.prevStrt > start {
		// new chromosome or after a seek
		c.seek(seq)
		c.start, c.end = -1, -1
	}
	if c.prevSeq == iseq && c.iseq != iseq {
		return false // no more regions on this chromosome
	}
	c.prevSeq = c.iseq
	c.prevStrt = start

	for iseq == c.iseq && c.end < start {
		if !c.next() {
			return false // no more regions left
		}
		if c.iseq != iseq {
			return false // moved off the chromosome
		}
	}
	return c.start <= end // region overlap
}
