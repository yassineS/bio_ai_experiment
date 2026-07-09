// Package bcftools — `bcftools merge` subcommand.
//
// `bcftools merge` combines multiple per-sample VCF/BCF files into a single
// multi-sample VCF/BCF. The inputs MUST be sorted, and they MUST contain
// disjoint sample sets (the upstream binary errors out otherwise; we mirror
// that behaviour).
//
// Merge rules:
//
//   - Sites with identical (CHROM, POS, REF) are collected; their ALTs are
//     unioned into a single ALT list and the per-sample FORMAT fields are
//     fanned out across the union sample set.
//   - The `-m`/`--merge` flag controls which kinds of records collapse:
//     `none` (no merging beyond identical records), `snps`, `indels`,
//     `both`, `all`, or `id` (collapse on ID only). The default is
//     `both`, matching upstream.
//   - For each sample missing in a given input, the output emits a
//     "missing" sample column (FORMAT values rendered as `.`).
//   - INFO fields are unioned by tag; conflicting numeric Number=A/R lines
//     are remapped onto the merged ALT order. Tags whose Number is fixed
//     (1, 0, ...) take the value from the first input that defines them.
//
// This is a streaming sort-merge: the cursor advances on whichever input has
// the lowest (contig-order, POS) and the records at that key are merged
// before being emitted.
package bcftools

import (
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// MergeMode controls which kinds of records are eligible for collapsing
// when their (CHROM, POS) coincide.
type MergeMode int

const (
	// MergeBoth collapses both SNPs and indels (upstream's default).
	MergeBoth MergeMode = iota
	// MergeNone disables collapsing — only byte-identical records merge.
	MergeNone
	// MergeSNPs collapses SNPs only; indels are emitted as separate rows.
	MergeSNPs
	// MergeIndels collapses indels only.
	MergeIndels
	// MergeAll collapses every record at the same (CHROM, POS).
	MergeAll
	// MergeID collapses records sharing the same ID column.
	MergeID
)

// ParseMergeMode parses the -m/--merge flag. Empty string returns the
// default MergeBoth.
func ParseMergeMode(s string) (MergeMode, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "", "both":
		return MergeBoth, nil
	case "none":
		return MergeNone, nil
	case "snps":
		return MergeSNPs, nil
	case "indels":
		return MergeIndels, nil
	case "all":
		return MergeAll, nil
	case "id":
		return MergeID, nil
	}
	return MergeBoth, fmt.Errorf("bcftools merge: unknown -m mode %q (expect none|snps|indels|both|all|id)", s)
}

// MergeOptions controls the `bcftools merge` subcommand.
type MergeOptions struct {
	// FileList, when non-empty, is appended to the positional file list
	// (one path per line, comments with `#` allowed).
	FileList string
	// MergeMode controls which records collapse — see MergeMode.
	MergeMode MergeMode
	// Regions, when non-empty, restricts output to records overlapping at
	// least one window. v1 applies regions as a post-filter (no index seek).
	Regions []string
	// RegionsFile is the BED-like sidecar parsed via LoadRegionsFile.
	RegionsFile string
	// OutputFormat selects the output encoding. Defaults to OutputVCF.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
	// Threads is the -@/--threads value; >1 enables parallel BGZF compression
	// of -O z and -O b output via bgzf.MultiWriter (see ViewOptions.Threads).
	Threads int
	// ForceSamples (--force-samples) resolves duplicate sample names across
	// inputs by prefixing the clashing name from input i (0-based) with
	// "<i+1>:" instead of erroring out, matching upstream vcfmerge.c
	// merge_headers (e.g. A + A -> A, 2:A).
	ForceSamples bool
	// InfoRules is the -i/--info-rules spec controlling how INFO fields combine
	// across merged records ("TAG:method,..."; method one of sum/avg/min/max/
	// join). "-" disables all rules. Empty selects the upstream default
	// (DP:sum,DP4:sum, plus AN:sum,AC:sum when the output has no samples).
	InfoRules string
}

// MergeFiles is the file-aware entry point. It opens each path through
// iohelper, fully reads them into memory, and runs Merge.
func MergeFiles(paths []string, out io.Writer, opts MergeOptions) (int, error) {
	all := append([]string{}, paths...)
	if opts.FileList != "" {
		extra, err := ReadFileList(opts.FileList)
		if err != nil {
			return 0, fmt.Errorf("bcftools merge: %w", err)
		}
		all = append(all, extra...)
	}
	if len(all) < 2 {
		return 0, fmt.Errorf("bcftools merge: need at least two input files (got %d)", len(all))
	}

	headers := make([]*vcf.Header, len(all))
	groups := make([][]*vcf.Variant, len(all))
	for i, p := range all {
		in, err := iohelper.OpenReader(p)
		if err != nil {
			return 0, fmt.Errorf("bcftools merge: open %s: %w", p, err)
		}
		hdr, recs, err := readAllVariants(in)
		_ = in.Close()
		if err != nil {
			return 0, fmt.Errorf("bcftools merge: %s: %w", p, err)
		}
		headers[i] = hdr
		groups[i] = recs
	}
	return Merge(headers, groups, out, opts)
}

// Merge takes pre-loaded per-input headers and variant slices, merges them,
// and writes the result to out. Inputs MUST be in sorted order (CHROM by
// header contig order, then POS). Sample sets across inputs MUST be disjoint.
func Merge(headers []*vcf.Header, groups [][]*vcf.Variant, out io.Writer, opts MergeOptions) (int, error) {
	if len(headers) != len(groups) {
		return 0, fmt.Errorf("bcftools merge: header / variant count mismatch")
	}
	if len(headers) == 0 {
		return 0, fmt.Errorf("bcftools merge: no inputs")
	}

	mergedHdr, renames, err := mergeMergeHeaders(headers, opts.ForceSamples)
	if err != nil {
		return 0, err
	}
	// Apply --force-samples renames to each input's per-record sample names so
	// the bucket fan-out (which keys by sample name) lines up with the de-duped
	// merged sample list.
	for gi, rn := range renames {
		if len(rn) == 0 {
			continue
		}
		for _, v := range groups[gi] {
			for si := range v.Samples {
				if nn, ok := rn[v.Samples[si].Name]; ok {
					v.Samples[si].Name = nn
				}
			}
		}
	}
	order := contigOrder(mergedHdr)
	infoRules, err := resolveInfoRules(opts.InfoRules, mergedHdr)
	if err != nil {
		return 0, err
	}
	regions, err := parseRegions(append(opts.Regions, []string{}...))
	if err != nil {
		return 0, err
	}
	if opts.RegionsFile != "" {
		regs, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return 0, fmt.Errorf("bcftools merge: -R %s: %w", opts.RegionsFile, err)
		}
		more, err := parseRegions(regs)
		if err != nil {
			return 0, err
		}
		regions = append(regions, more...)
	}

	// Sort each input group: upstream requires sorted input but the in-memory
	// representation might pre-sort by raw text order rather than the
	// merged-header contig order. Stable-sort by (contig-order, POS, REF)
	// to be safe.
	for i := range groups {
		sortVariantsByMergeKey(groups[i], order)
	}

	cursors := make([]int, len(groups))
	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, mergedHdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}

	count := 0
	for {
		// Find the lowest (contig, POS) across all live cursors.
		bestKey := mergeKey{}
		first := true
		any := false
		for gi, g := range groups {
			if cursors[gi] >= len(g) {
				continue
			}
			k := keyFor(g[cursors[gi]], order)
			if first || k.less(bestKey) {
				bestKey = k
				first = false
			}
			any = true
		}
		if !any {
			break
		}
		// Collect every record at the lowest key.
		participants := make([]variantRef, 0, len(groups))
		for gi, g := range groups {
			for cursors[gi] < len(g) {
				k := keyFor(g[cursors[gi]], order)
				if !k.equal(bestKey) {
					break
				}
				participants = append(participants, variantRef{groupIdx: gi, variant: g[cursors[gi]]})
				cursors[gi]++
			}
		}
		// Walk the participants and form merge buckets according to opts.
		// MergeMode. Each bucket becomes one output record.
		buckets := bucketize(participants, opts.MergeMode)
		for _, bk := range buckets {
			rec := mergeBucket(bk, headers, mergedHdr, infoRules)
			if len(regions) > 0 && !overlapsAny(rec, regions) {
				continue
			}
			if err := w.Write(rec); err != nil {
				return count, err
			}
			count++
		}
	}
	return count, w.Flush()
}

// mergeKey is the (contig-index, POS) tuple used to advance cursors.
type mergeKey struct {
	contigIdx int
	pos       int
}

func (k mergeKey) less(o mergeKey) bool {
	if k.contigIdx != o.contigIdx {
		return k.contigIdx < o.contigIdx
	}
	return k.pos < o.pos
}

func (k mergeKey) equal(o mergeKey) bool { return k.contigIdx == o.contigIdx && k.pos == o.pos }

func keyFor(v *vcf.Variant, order map[string]int) mergeKey {
	idx, ok := order[v.Chrom]
	if !ok {
		idx = 1<<30 + sortFallback(v.Chrom)
	}
	return mergeKey{contigIdx: idx, pos: v.Pos}
}

// sortVariantsByMergeKey stable-sorts in-place by (contig-order, POS, REF).
func sortVariantsByMergeKey(recs []*vcf.Variant, order map[string]int) {
	sort.SliceStable(recs, func(i, j int) bool {
		ki := keyFor(recs[i], order)
		kj := keyFor(recs[j], order)
		if !ki.equal(kj) {
			return ki.less(kj)
		}
		return recs[i].Ref < recs[j].Ref
	})
}

// variantRef carries the input-group index alongside the variant so we can
// trace each merged sample back to its source.
type variantRef struct {
	groupIdx int
	variant  *vcf.Variant
}

// pairByOccurrence splits a set of records that are eligible to merge into
// buckets by per-input occurrence index: the k-th record from each input lands
// in bucket k. This mirrors bcftools' maux, where each output record draws at
// most one line from any single input — so a file that holds two compatible
// records at a position (an intra-position duplicate or a split multiallelic)
// contributes to two distinct output records rather than collapsing them.
func pairByOccurrence(group []variantRef) [][]variantRef {
	occ := map[int]int{}
	var buckets [][]variantRef
	for _, p := range group {
		k := occ[p.groupIdx]
		occ[p.groupIdx]++
		for len(buckets) <= k {
			buckets = append(buckets, nil)
		}
		buckets[k] = append(buckets[k], p)
	}
	return buckets
}

// bucketize groups co-located records according to MergeMode. The result
// is a slice of buckets (each bucket becomes one output record).
func bucketize(parts []variantRef, mode MergeMode) [][]variantRef {
	if len(parts) == 0 {
		return nil
	}
	switch mode {
	case MergeAll:
		return pairByOccurrence(parts)
	case MergeNone:
		// Each input record stands on its own — except byte-identical
		// (REF + ALT) records collapse trivially.
		groupedByREFALT := map[string][]variantRef{}
		var keys []string
		for _, p := range parts {
			k := p.variant.Ref + "|" + strings.Join(p.variant.Alt, ",")
			if _, ok := groupedByREFALT[k]; !ok {
				keys = append(keys, k)
			}
			groupedByREFALT[k] = append(groupedByREFALT[k], p)
		}
		out := make([][]variantRef, 0, len(keys))
		for _, k := range keys {
			out = append(out, groupedByREFALT[k])
		}
		return out
	case MergeID:
		groupedByID := map[string][]variantRef{}
		var keys []string
		for _, p := range parts {
			k := p.variant.ID
			if k == "" {
				k = "."
			}
			if _, ok := groupedByID[k]; !ok {
				keys = append(keys, k)
			}
			groupedByID[k] = append(groupedByID[k], p)
		}
		out := make([][]variantRef, 0, len(keys))
		for _, k := range keys {
			out = append(out, groupedByID[k])
		}
		return out
	case MergeSNPs:
		var snp, other []variantRef
		for _, p := range parts {
			if isSNPRecord(p.variant) {
				snp = append(snp, p)
			} else {
				other = append(other, p)
			}
		}
		out := [][]variantRef{}
		// SNPs merge, but pair by per-input occurrence so duplicates split.
		out = append(out, pairByOccurrence(snp)...)
		// Non-SNP records are emitted one bucket per record.
		for _, p := range other {
			out = append(out, []variantRef{p})
		}
		return out
	case MergeIndels:
		var indel, other []variantRef
		for _, p := range parts {
			if isIndelRecord(p.variant) {
				indel = append(indel, p)
			} else {
				other = append(other, p)
			}
		}
		out := [][]variantRef{}
		out = append(out, pairByOccurrence(indel)...)
		for _, p := range other {
			out = append(out, []variantRef{p})
		}
		return out
	case MergeBoth:
		fallthrough
	default:
		var snp, indel, other []variantRef
		for _, p := range parts {
			switch {
			case isSNPRecord(p.variant):
				snp = append(snp, p)
			case isIndelRecord(p.variant):
				indel = append(indel, p)
			default:
				other = append(other, p)
			}
		}
		out := [][]variantRef{}
		out = append(out, pairByOccurrence(snp)...)
		out = append(out, pairByOccurrence(indel)...)
		for _, p := range other {
			out = append(out, []variantRef{p})
		}
		return orderBucketsByInput(out, parts)
	}
}

// orderBucketsByInput stable-sorts buckets by the original input position of
// their first-appearing record. bcftools emits the merged records at a site in
// the order the records appear in the inputs (first input first), not grouped
// by variant type, so a SNP and an indel at one position keep their file order.
func orderBucketsByInput(buckets [][]variantRef, parts []variantRef) [][]variantRef {
	idx := make(map[*vcf.Variant]int, len(parts))
	for i, p := range parts {
		if _, ok := idx[p.variant]; !ok {
			idx[p.variant] = i
		}
	}
	minIdx := func(bk []variantRef) int {
		m := 1 << 30
		for _, p := range bk {
			if i, ok := idx[p.variant]; ok && i < m {
				m = i
			}
		}
		return m
	}
	sort.SliceStable(buckets, func(i, j int) bool {
		return minIdx(buckets[i]) < minIdx(buckets[j])
	})
	return buckets
}

// isSNPRecord returns true when REF is one base and every ALT is one base.
func isSNPRecord(v *vcf.Variant) bool {
	if len(v.Ref) != 1 {
		return false
	}
	for _, a := range v.Alt {
		if len(a) != 1 {
			return false
		}
	}
	return true
}

// isIndelRecord returns true when at least one allele has a different length
// from REF.
func isIndelRecord(v *vcf.Variant) bool {
	for _, a := range v.Alt {
		if len(a) != len(v.Ref) {
			return true
		}
	}
	return false
}

// mergeBucket collapses a bucket into a single output record.
// Allele numbering is unified across the bucket and each input sample's
// FORMAT values are remapped to the new ALT indices. Samples present in the
// merged header but not in this bucket are emitted with `.` placeholders.
func mergeBucket(bk []variantRef, srcHeaders []*vcf.Header, mergedHdr *vcf.Header, infoRules map[string]infoRule) *vcf.Variant {
	mergedSamples := mergedHdr.Samples
	first := bk[0].variant

	// Build the union ALT list, in first-seen order.
	altIndex := map[string]int{}    // ALT string -> output index (1-based)
	altList := make([]string, 0, 4) // output order
	perSrcMap := make([]map[int]int, len(bk))
	for i, ref := range bk {
		m := map[int]int{0: 0} // 0 (REF) maps to 0
		for j, a := range ref.variant.Alt {
			if _, ok := altIndex[a]; !ok {
				altIndex[a] = len(altList) + 1
				altList = append(altList, a)
			}
			m[j+1] = altIndex[a]
		}
		perSrcMap[i] = m
	}

	out := &vcf.Variant{
		Chrom:  first.Chrom,
		Pos:    first.Pos,
		Ref:    first.Ref,
		ID:     mergeIDs(bk),
		Alt:    altList,
		Qual:   mergeQuals(bk),
		Filter: mergeFilters(bk),
		Info:   map[string]string{},
	}

	// FORMAT order: union, with GT first if any input has GT.
	out.Format = mergeFormat(bk)

	// INFO union — first input wins on collisions for tags whose Number
	// isn't A or R; tags with A/R are remapped to the merged allele order.
	// INFO rules (e.g. DP:sum) then combine values across the bucket.
	out.Info, out.InfoOrder = mergeInfo(bk, perSrcMap, len(altList), infoRules, mergedHdr)

	// Per-sample fan-out.
	out.Samples = make([]vcf.Sample, len(mergedSamples))
	// Build a fast lookup: sample name -> (which bucket index, which sample
	// index within that variant).
	type srcLoc struct {
		bk      int
		sampleI int
	}
	loc := map[string]srcLoc{}
	for bi, ref := range bk {
		for si, s := range ref.variant.Samples {
			loc[s.Name] = srcLoc{bk: bi, sampleI: si}
		}
	}
	for i, name := range mergedSamples {
		s := vcf.Sample{Name: name, Data: map[string]string{}}
		l, ok := loc[name]
		if !ok {
			// Samples not present in this bucket get all-missing data.
			for _, f := range out.Format {
				if f == "GT" {
					s.Data[f] = "./."
				} else {
					s.Data[f] = "."
				}
			}
			out.Samples[i] = s
			continue
		}
		srcVar := bk[l.bk].variant
		srcSample := srcVar.Samples[l.sampleI]
		alleleMap := perSrcMap[l.bk]
		for _, f := range out.Format {
			val, ok := srcSample.Data[f]
			if !ok || val == "" {
				if f == "GT" {
					s.Data[f] = "./."
				} else {
					s.Data[f] = "."
				}
				continue
			}
			if f == "GT" {
				s.Data[f] = remapGTByMap(val, alleleMap)
			} else {
				s.Data[f] = val
			}
		}
		out.Samples[i] = s
	}
	return out
}

// mergeIDs picks the first non-`.` ID across the bucket, joining unique
// values with `;` (matching upstream behaviour).
func mergeIDs(bk []variantRef) string {
	seen := map[string]bool{}
	var ids []string
	for _, ref := range bk {
		if ref.variant.ID == "" || ref.variant.ID == "." {
			continue
		}
		for _, p := range strings.Split(ref.variant.ID, ";") {
			if seen[p] {
				continue
			}
			seen[p] = true
			ids = append(ids, p)
		}
	}
	if len(ids) == 0 {
		return "."
	}
	return strings.Join(ids, ";")
}

// mergeQuals takes the maximum QUAL across the bucket. Records with `.`
// (missing) QUAL contribute nothing.
func mergeQuals(bk []variantRef) float64 {
	best := -1.0
	for _, ref := range bk {
		if ref.variant.Qual >= 0 && ref.variant.Qual > best {
			best = ref.variant.Qual
		}
	}
	return best
}

// mergeFilters returns the union of FILTER values across the bucket. PASS is
// only kept if every record was PASS — otherwise upstream drops it.
func mergeFilters(bk []variantRef) []string {
	allPass := true
	seen := map[string]bool{}
	var out []string
	for _, ref := range bk {
		if len(ref.variant.Filter) == 0 {
			continue
		}
		isPass := false
		for _, f := range ref.variant.Filter {
			if f == "PASS" {
				isPass = true
			}
		}
		if !isPass {
			allPass = false
		}
		for _, f := range ref.variant.Filter {
			if f == "PASS" {
				continue
			}
			if seen[f] {
				continue
			}
			seen[f] = true
			out = append(out, f)
		}
	}
	if len(out) == 0 {
		if allPass {
			return []string{"PASS"}
		}
		return []string{"."}
	}
	if allPass && len(out) == 0 {
		return []string{"PASS"}
	}
	return out
}

// mergeFormat returns the union of FORMAT keys, with GT first when present
// in any input.
func mergeFormat(bk []variantRef) []string {
	seen := map[string]bool{}
	var out []string
	for _, ref := range bk {
		for _, f := range ref.variant.Format {
			if seen[f] {
				continue
			}
			seen[f] = true
			out = append(out, f)
		}
	}
	// Move GT to the front if present.
	for i, f := range out {
		if f == "GT" && i > 0 {
			out = append([]string{"GT"}, append(out[:i], out[i+1:]...)...)
			break
		}
	}
	return out
}

// mergeInfo returns the union of INFO tags together with the output key order.
//
// Output order mirrors upstream vcfmerge.c merge_info(), which emits INFO tags
// in three classes:
//
//  1. plain scalar / fixed-Number tags with no rule — in first-seen order;
//  2. rule-combined tags (e.g. DP:sum) — in alphabetical tag order (upstream
//     applies the rules via bcf_update_info after sorting them by tag);
//  3. Number=A/R/G tags with no rule (e.g. AF) — in first-seen order, last.
//
// A tag that carries a rule is always treated as class (2), even when its
// Number is A/R/G (upstream's info_rules_add_values path pre-empts the AGR
// path). Class (3) AGR-no-rule tags combine per-allele with LAST-non-missing
// wins across every record in the bucket, remapping each source record's
// per-ALT values through its allele map into the output allele slots.
func mergeInfo(bk []variantRef, perSrc []map[int]int, nAlt int, rules map[string]infoRule, hdr *vcf.Header) (map[string]string, []string) {
	out := map[string]string{}

	// isAGR reports whether the tag's declared Number is A, R or G.
	isAGR := func(tag string) bool {
		switch infoNumber(hdr, tag) {
		case "A", "R", "G":
			return true
		}
		return false
	}

	var scalarOrder []string // class 1: first-seen
	var agrOrder []string    // class 3: first-seen
	var ruleOrder []string   // class 2: gathered, then sorted alphabetically
	scalarSeen := map[string]bool{}
	agrSeen := map[string]bool{}
	ruleSeen := map[string]bool{}

	// agrVals accumulates the per-allele output slots for each class-3 tag,
	// applying last-non-missing-wins as records are visited in input order.
	agrVals := map[string][]string{}

	for bi, ref := range bk {
		for _, k := range ref.variant.InfoOrder {
			v := ref.variant.Info[k]
			if _, ok := rules[k]; ok {
				// class 2 — rule tag. Combined below via applyInfoRule.
				if !ruleSeen[k] {
					ruleSeen[k] = true
					ruleOrder = append(ruleOrder, k)
				}
				continue
			}
			if isAGR(k) {
				// class 3 — Number=A/R/G, no rule: per-allele last-wins.
				if !agrSeen[k] {
					agrSeen[k] = true
					agrOrder = append(agrOrder, k)
					slots := make([]string, agrSlots(hdr, k, nAlt))
					for i := range slots {
						slots[i] = "."
					}
					agrVals[k] = slots
				}
				remapAGRInto(agrVals[k], v, infoNumber(hdr, k), infoType(hdr, k), perSrc[bi])
				continue
			}
			// class 1 — plain scalar / fixed-Number, first-wins.
			if !scalarSeen[k] {
				scalarSeen[k] = true
				scalarOrder = append(scalarOrder, k)
				out[k] = v
			}
		}
	}

	// class 2 — apply rules, emit alphabetically.
	sort.Strings(ruleOrder)
	for _, k := range ruleOrder {
		out[k] = applyInfoRule(bk, k, rules[k], nAlt, infoNumber(hdr, k), infoType(hdr, k), perSrc)
	}

	// class 3 — materialise the accumulated per-allele slots.
	for _, k := range agrOrder {
		out[k] = strings.Join(agrVals[k], ",")
	}

	order := make([]string, 0, len(scalarOrder)+len(ruleOrder)+len(agrOrder))
	order = append(order, scalarOrder...)
	order = append(order, ruleOrder...)
	order = append(order, agrOrder...)
	return out, order
}

// infoNumber returns the declared Number= field ("A", "R", "G", "0", "1",
// ".", a positive integer, ...) for INFO tag from the header, or "" when the
// tag is not declared. It reads the authoritative ##INFO header line rather
// than guessing from the value's arity.
func infoNumber(hdr *vcf.Header, tag string) string {
	if hdr == nil {
		return ""
	}
	for _, ln := range hdr.MetaInfo {
		if !strings.HasPrefix(ln, "##INFO=<") {
			continue
		}
		if headerLineID(ln) != tag {
			continue
		}
		return headerLineField(ln, "Number")
	}
	return ""
}

// headerLineField extracts a key=value field (e.g. Number=A) from a structured
// header line, returning "" when the key is absent.
func headerLineField(ln, key string) string {
	i := strings.Index(ln, key+"=")
	if i < 0 {
		return ""
	}
	s := ln[i+len(key)+1:]
	for j := 0; j < len(s); j++ {
		if s[j] == ',' || s[j] == '>' {
			return s[:j]
		}
	}
	return s
}

// agrSlots returns the number of output value slots a Number=A/R/G tag occupies
// given nAlt output ALT alleles. G is the number of diploid genotypes over
// nAlt+1 alleles; R is nAlt+1; A is nAlt.
func agrSlots(hdr *vcf.Header, tag string, nAlt int) int {
	switch infoNumber(hdr, tag) {
	case "R":
		return nAlt + 1
	case "G":
		n := nAlt + 1 // total alleles incl. REF
		return n * (n + 1) / 2
	default: // "A"
		return nAlt
	}
}

// remapAGRInto writes a source record's comma-separated Number=A/R/G value into
// the output slots, remapping each source allele index through alleleMap and
// overwriting with the last non-missing value seen per slot. Float values are
// re-normalised (e.g. 0.50 -> 0.5) to match upstream, which stores AGR values
// numerically and re-renders them. G-typed tags are copied positionally
// (genotype remapping across differing ALT sets is not attempted; this matches
// the common single-ALT case byte-exactly).
func remapAGRInto(slots []string, val, number, typ string, alleleMap map[int]int) {
	parts := strings.Split(val, ",")
	norm := func(p string) string {
		if typ != "Float" {
			return p
		}
		if f, err := strconv.ParseFloat(p, 64); err == nil {
			return formatInfoFloat(f)
		}
		return p
	}
	switch number {
	case "A":
		// parts[j] is for source ALT j+1 -> output allele alleleMap[j+1].
		for j, p := range parts {
			if p == "." || p == "" {
				continue
			}
			dst, ok := alleleMap[j+1]
			if !ok {
				continue
			}
			if idx := dst - 1; idx >= 0 && idx < len(slots) {
				slots[idx] = norm(p)
			}
		}
	case "R":
		// parts[j] is for source allele j (0=REF) -> output allele alleleMap[j].
		for j, p := range parts {
			if p == "." || p == "" {
				continue
			}
			dst, ok := alleleMap[j]
			if !ok {
				continue
			}
			if dst >= 0 && dst < len(slots) {
				slots[dst] = norm(p)
			}
		}
	default: // "G" — positional copy.
		for j, p := range parts {
			if p == "." || p == "" || j >= len(slots) {
				continue
			}
			slots[j] = norm(p)
		}
	}
}

// infoRule is one INFO-combine method for the -i/--info-rules machinery.
type infoRule int

const (
	infoRuleSum  infoRule = iota // element-wise sum
	infoRuleAvg                  // element-wise mean
	infoRuleMin                  // element-wise minimum
	infoRuleMax                  // element-wise maximum
	infoRuleJoin                 // comma-join every value across records
)

// resolveInfoRules turns the -i/--info-rules spec into a tag->rule map. "-"
// disables all rules; "" selects upstream's default (DP:sum, DP4:sum, and
// AN:sum,AC:sum when the merged header carries no samples), restricted to tags
// that are actually declared in the merged header.
func resolveInfoRules(spec string, hdr *vcf.Header) (map[string]infoRule, error) {
	if spec == "-" {
		return nil, nil
	}
	if spec == "" {
		declared := infoTagSet(hdr)
		var parts []string
		if declared["DP"] {
			parts = append(parts, "DP:sum")
		}
		if declared["DP4"] {
			parts = append(parts, "DP4:sum")
		}
		if len(hdr.Samples) == 0 {
			if declared["AN"] {
				parts = append(parts, "AN:sum")
			}
			if declared["AC"] {
				parts = append(parts, "AC:sum")
			}
		}
		spec = strings.Join(parts, ",")
	}
	return parseInfoRules(spec, hdr)
}

// parseInfoRules parses a "TAG:method,TAG:method" spec against the merged
// header. It mirrors upstream's info_rules_init validation:
//
//   - the tag must be declared in the merged header (else an error);
//   - a numeric method (sum/avg/min/max) on a String/Type-flag tag is rejected;
//   - join on an ALREADY variable-length A/R/G tag rewrites that tag's declared
//     Number to '.' in the merged header (upstream only rewrites when the tag's
//     Number was variable, i.e. 0xfffff; a fixed Number=1 is left untouched).
func parseInfoRules(spec string, hdr *vcf.Header) (map[string]infoRule, error) {
	if spec == "" {
		return nil, nil
	}
	out := map[string]infoRule{}
	for _, item := range strings.Split(spec, ",") {
		kv := strings.SplitN(item, ":", 2)
		if len(kv) != 2 || kv[0] == "" {
			return nil, fmt.Errorf("bcftools merge: could not parse INFO rule %q", item)
		}
		tag := kv[0]
		var r infoRule
		switch kv[1] {
		case "sum":
			r = infoRuleSum
		case "avg":
			r = infoRuleAvg
		case "min":
			r = infoRuleMin
		case "max":
			r = infoRuleMax
		case "join":
			r = infoRuleJoin
		default:
			return nil, fmt.Errorf("bcftools merge: unknown INFO rule method %q", kv[1])
		}
		if hdr != nil {
			if !infoTagSet(hdr)[tag] {
				return nil, fmt.Errorf("The INFO tag is not defined in the header: \"%s\"", tag)
			}
			typ := infoType(hdr, tag)
			if r != infoRuleJoin && (typ == "String" || typ == "Flag") {
				return nil, fmt.Errorf("Numeric operation \"%s\" requested on non-numeric field: %s", kv[1], tag)
			}
			if r == infoRuleJoin {
				// Upstream only rewrites Number when it was ALREADY variable
				// (A/R/G) and not the plain "." variable form; a fixed Number
				// (1, 4, ...) is left untouched.
				switch infoNumber(hdr, tag) {
				case "A", "R", "G":
					setInfoNumber(hdr, tag, ".")
				}
			}
		}
		out[tag] = r
	}
	return out, nil
}

// infoType returns the declared Type= field ("Integer", "Float", "String",
// "Flag", ...) for INFO tag from the header, or "" when undeclared.
func infoType(hdr *vcf.Header, tag string) string {
	if hdr == nil {
		return ""
	}
	for _, ln := range hdr.MetaInfo {
		if !strings.HasPrefix(ln, "##INFO=<") || headerLineID(ln) != tag {
			continue
		}
		return headerLineField(ln, "Type")
	}
	return ""
}

// setInfoNumber rewrites the Number= field of INFO tag's header line to n,
// mutating the merged header in place (used by join on an A/R/G tag). It
// mirrors upstream's bcf_hdr_remove + bcf_hdr_add_hrec: the rewritten line is
// removed from its original position and re-appended to the very end of the
// meta-information block (after any ##FORMAT lines), so it moves to the end.
func setInfoNumber(hdr *vcf.Header, tag, n string) {
	if hdr == nil {
		return
	}
	idx := -1
	for i, ln := range hdr.MetaInfo {
		if strings.HasPrefix(ln, "##INFO=<") && headerLineID(ln) == tag {
			idx = i
			break
		}
	}
	if idx < 0 {
		return
	}
	old := headerLineField(hdr.MetaInfo[idx], "Number")
	rewritten := strings.Replace(hdr.MetaInfo[idx], "Number="+old, "Number="+n, 1)
	hdr.MetaInfo = append(hdr.MetaInfo[:idx], hdr.MetaInfo[idx+1:]...)
	hdr.MetaInfo = append(hdr.MetaInfo, rewritten)
}

// infoTagSet returns the set of INFO tag IDs declared in the header.
func infoTagSet(hdr *vcf.Header) map[string]bool {
	out := map[string]bool{}
	for _, ln := range hdr.MetaInfo {
		if !strings.HasPrefix(ln, "##INFO=<") {
			continue
		}
		if id := headerLineID(ln); id != "" {
			out[id] = true
		}
	}
	return out
}

// headerLineID extracts the ID=<value> field from a structured header line.
func headerLineID(ln string) string {
	i := strings.Index(ln, "ID=")
	if i < 0 {
		return ""
	}
	s := ln[i+3:]
	for j := 0; j < len(s); j++ {
		if s[j] == ',' || s[j] == '>' {
			return s[:j]
		}
	}
	return s
}

// applyInfoRule combines a tag's value across every record in the bucket per
// the rule, mirroring upstream's info_rules_* mergers:
//
//   - join   — concatenates every record's raw value with commas;
//   - sum/avg — element-wise; a per-record missing element counts as 0. avg
//     divides each slot by the number of records that CARRIED the tag (not by
//     the number of missing-aware contributions);
//   - min/max — element-wise, skipping missing elements; a slot missing in
//     every record renders as '.'.
//
// When number is A/R/G the source per-allele values are remapped through each
// record's allele map (perSrc) into the output allele slots before combining,
// so inputs with differing ALT sets combine correctly.
func applyInfoRule(bk []variantRef, tag string, r infoRule, nAlt int, number, typ string, perSrc []map[int]int) string {
	if r == infoRuleJoin {
		var vals []string
		for _, ref := range bk {
			v, ok := ref.variant.Info[tag]
			if !ok {
				continue
			}
			// Numeric joins re-render each element (upstream stores the joined
			// vector numerically); Float normalises 0.50 -> 0.5.
			if typ == "Float" {
				elems := strings.Split(v, ",")
				for i, e := range elems {
					if e == "." || e == "" {
						continue
					}
					if f, err := strconv.ParseFloat(e, 64); err == nil {
						elems[i] = formatInfoFloat(f)
					}
				}
				v = strings.Join(elems, ",")
			}
			vals = append(vals, v)
		}
		return strings.Join(vals, ",")
	}

	isAGR := number == "A" || number == "R" || number == "G"

	var acc []float64 // running accumulator per output slot
	var present []int // count of records contributing a non-missing value per slot
	count := 0        // records that carried the tag (avg divisor)
	anyFloat := false // render as float when any input looked float-like
	blockSize := -1   // fixed-Number vector length (non-AGR)
	for bi, ref := range bk {
		v, ok := ref.variant.Info[tag]
		if !ok {
			continue
		}
		count++
		parts := strings.Split(v, ",")

		// Build this record's per-slot value list. For AGR tags remap through
		// the allele map; for fixed-Number tags use the values positionally.
		var slots []string
		if isAGR {
			slots = make([]string, agrSlotsN(number, nAlt))
			for i := range slots {
				slots[i] = "."
			}
			remapAGRInto(slots, v, number, "", perSrc[bi])
		} else {
			if blockSize < 0 {
				blockSize = len(parts)
			} else if len(parts) != blockSize {
				// Vector-length mismatch: fall back to the raw value.
				return v
			}
			slots = parts
		}

		if acc == nil {
			acc = make([]float64, len(slots))
			present = make([]int, len(slots))
		}
		if len(slots) != len(acc) {
			return v
		}
		for i, p := range slots {
			missing := p == "." || p == ""
			var f float64
			if !missing {
				var err error
				f, err = strconv.ParseFloat(p, 64)
				if err != nil {
					return v
				}
				if strings.ContainsAny(p, ".eE") {
					anyFloat = true
				}
			}
			switch r {
			case infoRuleSum, infoRuleAvg:
				// Missing counts as 0 (upstream sets missing->0 pre-sum).
				acc[i] += f
				present[i]++
			case infoRuleMin:
				if missing {
					continue
				}
				if present[i] == 0 || f < acc[i] {
					acc[i] = f
				}
				present[i]++
			case infoRuleMax:
				if missing {
					continue
				}
				if present[i] == 0 || f > acc[i] {
					acc[i] = f
				}
				present[i]++
			}
		}
	}
	if count == 0 {
		return ""
	}
	if r == infoRuleAvg {
		for i := range acc {
			acc[i] /= float64(count)
		}
		anyFloat = true
	}
	strs := make([]string, len(acc))
	for i, a := range acc {
		if (r == infoRuleMin || r == infoRuleMax) && present[i] == 0 {
			strs[i] = "."
			continue
		}
		if anyFloat {
			strs[i] = formatInfoFloat(a)
		} else {
			strs[i] = strconv.FormatInt(int64(a), 10)
		}
	}
	return strings.Join(strs, ",")
}

// agrSlotsN returns the output slot count for a Number=A/R/G tag given nAlt
// output ALT alleles (the header-free counterpart of agrSlots).
func agrSlotsN(number string, nAlt int) int {
	switch number {
	case "R":
		return nAlt + 1
	case "G":
		n := nAlt + 1
		return n * (n + 1) / 2
	default: // "A"
		return nAlt
	}
}

// formatInfoFloat renders a combined float the way bcftools prints INFO floats
// (shortest round-trippable form).
func formatInfoFloat(f float64) string {
	return strconv.FormatFloat(f, 'g', -1, 64)
}

// remapGTByMap rewrites a GT string ("0/1", "1|0", "./.", ...) using a
// per-input allele-index map so the result references the merged ALT
// positions.
func remapGTByMap(gt string, alleleMap map[int]int) string {
	if gt == "" || gt == "." {
		return gt
	}
	// Preserve the separator (/, |) used in the input. We split, remap,
	// then re-join. Mixed separators are uncommon and we treat the whole
	// string as if it had a single separator — same as upstream.
	sep := "/"
	if strings.IndexByte(gt, '|') >= 0 {
		sep = "|"
	}
	parts := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
	if len(parts) == 0 {
		return gt
	}
	out := make([]string, len(parts))
	for i, p := range parts {
		if p == "." {
			out[i] = "."
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			out[i] = p
			continue
		}
		mapped, ok := alleleMap[n]
		if !ok {
			out[i] = p
		} else {
			out[i] = strconv.Itoa(mapped)
		}
	}
	return strings.Join(out, sep)
}

// mergeMergeHeaders returns a *vcf.Header that's the union across the inputs,
// with the sample list being the *concatenation* of each input's samples in
// input order. Conflicting structured definitions are an error.
func mergeMergeHeaders(heads []*vcf.Header, forceSamples bool) (*vcf.Header, []map[string]string, error) {
	if len(heads) == 0 {
		return nil, nil, fmt.Errorf("bcftools merge: no headers")
	}
	out := &vcf.Header{}
	seenLine := map[string]bool{}
	definitions := map[string]string{}

	// fileformat: take the first one.
	for _, h := range heads {
		for _, m := range h.MetaInfo {
			if strings.HasPrefix(m, "##fileformat=") {
				if !seenLine[m] {
					out.MetaInfo = append(out.MetaInfo, m)
					seenLine[m] = true
				}
				break
			}
		}
	}

	for _, h := range heads {
		for _, m := range h.MetaInfo {
			if strings.HasPrefix(m, "##fileformat=") {
				continue
			}
			key, structID := structuredID(m)
			if key != "" {
				dkey := key + "/" + structID
				if prev, ok := definitions[dkey]; ok {
					if prev != m {
						return nil, nil, fmt.Errorf("bcftools merge: conflicting %s ID %q definitions:\n  %s\n  %s", key, structID, prev, m)
					}
					continue
				}
				definitions[dkey] = m
				out.MetaInfo = append(out.MetaInfo, m)
				continue
			}
			if seenLine[m] {
				continue
			}
			seenLine[m] = true
			out.MetaInfo = append(out.MetaInfo, m)
		}
	}

	// Concatenate sample lists. Samples must normally be disjoint across
	// inputs; with --force-samples a clashing name from input i (0-based) is
	// prefixed with "<i+1>:" (repeatedly, if the prefixed form also clashes),
	// mirroring upstream vcfmerge.c merge_headers. renames[i] maps an original
	// sample name in input i to its de-duped output name (only entries that
	// actually changed are recorded).
	seenSample := map[string]bool{}
	renames := make([]map[string]string, len(heads))
	for i, h := range heads {
		renames[i] = map[string]string{}
		prefix := strconv.Itoa(i + 1)
		for _, s := range h.Samples {
			name := s
			if seenSample[name] {
				if !forceSamples {
					return nil, nil, fmt.Errorf("Duplicate sample names (%s), use --force-samples to proceed anyway", s)
				}
				for seenSample[name] {
					name = prefix + ":" + name
				}
				renames[i][s] = name
			}
			seenSample[name] = true
			out.Samples = append(out.Samples, name)
		}
	}
	return out, renames, nil
}
