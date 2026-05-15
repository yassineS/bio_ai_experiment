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

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
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

	mergedHdr, err := mergeMergeHeaders(headers)
	if err != nil {
		return 0, err
	}
	order := contigOrder(mergedHdr)
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
			rec := mergeBucket(bk, headers, mergedHdr.Samples, order)
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

// bucketize groups co-located records according to MergeMode. The result
// is a slice of buckets (each bucket becomes one output record).
func bucketize(parts []variantRef, mode MergeMode) [][]variantRef {
	if len(parts) == 0 {
		return nil
	}
	switch mode {
	case MergeAll:
		return [][]variantRef{parts}
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
		if len(snp) > 0 {
			out = append(out, snp)
		}
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
		if len(indel) > 0 {
			out = append(out, indel)
		}
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
		if len(snp) > 0 {
			out = append(out, snp)
		}
		if len(indel) > 0 {
			out = append(out, indel)
		}
		for _, p := range other {
			out = append(out, []variantRef{p})
		}
		return out
	}
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
func mergeBucket(bk []variantRef, srcHeaders []*vcf.Header, mergedSamples []string, _ map[string]int) *vcf.Variant {
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
	out.Info, out.InfoOrder = mergeInfo(bk, perSrcMap, len(altList))

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

// mergeInfo returns the union of INFO tags. Tags with comma-separated
// per-allele values (one per ALT) are remapped to the new allele order. For
// fixed-Number tags the first input that defines the tag wins.
func mergeInfo(bk []variantRef, perSrc []map[int]int, nAlt int) (map[string]string, []string) {
	out := map[string]string{}
	var order []string
	tagSeen := map[string]bool{}
	for bi, ref := range bk {
		for _, k := range ref.variant.InfoOrder {
			v := ref.variant.Info[k]
			if !tagSeen[k] {
				tagSeen[k] = true
				order = append(order, k)
			}
			if _, exists := out[k]; exists {
				continue
			}
			// Heuristic: a comma-separated value with len(parts) == #ALT
			// in this input is a per-allele tag.
			parts := strings.Split(v, ",")
			if len(parts) == len(ref.variant.Alt) && nAlt > 0 {
				remapped := make([]string, nAlt)
				for i := range remapped {
					remapped[i] = "."
				}
				for j, p := range parts {
					mappedIdx, ok := perSrc[bi][j+1]
					if !ok {
						continue
					}
					if mappedIdx-1 >= 0 && mappedIdx-1 < len(remapped) {
						remapped[mappedIdx-1] = p
					}
				}
				out[k] = strings.Join(remapped, ",")
			} else {
				out[k] = v
			}
		}
	}
	return out, order
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
func mergeMergeHeaders(heads []*vcf.Header) (*vcf.Header, error) {
	if len(heads) == 0 {
		return nil, fmt.Errorf("bcftools merge: no headers")
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
						return nil, fmt.Errorf("bcftools merge: conflicting %s ID %q definitions:\n  %s\n  %s", key, structID, prev, m)
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

	// Concatenate sample lists; samples must be disjoint across inputs.
	seenSample := map[string]bool{}
	for _, h := range heads {
		for _, s := range h.Samples {
			if seenSample[s] {
				return nil, fmt.Errorf("bcftools merge: sample %q appears in more than one input", s)
			}
			seenSample[s] = true
			out.Samples = append(out.Samples, s)
		}
	}
	return out, nil
}
