package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// ConcatOptions controls the behaviour of Concat / ConcatFiles.
type ConcatOptions struct {
	// OutputFormat selects the output format. The default OutputVCF emits
	// plain VCF text.
	OutputFormat OutputFormat
	// AllowOverlaps requests a sort-merge concatenation: records are
	// returned in (contig-order, POS) order across all inputs. When false
	// the inputs are concatenated in the order given.
	AllowOverlaps bool
	// RemoveDuplicates is the `-D/--remove-duplicates` flag: it is an alias
	// for `-d exact`. Upstream only honours duplicate removal together with
	// AllowOverlaps (`-a`); standalone use is rejected (see ConcatFiles).
	RemoveDuplicates bool
	// RmDupMode is the `-d/--rm-dup` collapse mode. The empty string means
	// no de-duplication. Recognised modes mirror upstream vcfconcat.c:
	// "exact"/"none" (drop only records with identical REF+ALT seen earlier
	// across files), "snps", "indels", "both", "any"/"all". As with `-D`,
	// dedup requires AllowOverlaps.
	RmDupMode string
	// FileList, when non-empty, is read line-by-line and prepended to the
	// list of inputs supplied to ConcatFiles.
	FileList string
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
	// Threads is the -@/--threads value; >1 enables parallel BGZF compression
	// of -O z and -O b output via bgzf.MultiWriter (see ViewOptions.Threads).
	Threads int
	// MinPQ is the `-q/--min-PQ` value: when a sample's phase quality at a
	// chunk boundary is below MinPQ a new phase set starts for that sample.
	// Upstream's default is 30.
	MinPQ int
	// Ligate enables phased concatenation (`-l/--ligate`): overlapping phased
	// chunks are ligated, the overlap is emitted once, phase is reconciled
	// across chunks, and FORMAT/PS and FORMAT/PQ are added.
	Ligate bool
	// LigateForce (`--ligate-force`) downgrades the "chunks do not overlap"
	// error to a clean join, keeping all sites.
	LigateForce bool
	// LigateWarn (`--ligate-warn`) drops sites in imperfect overlaps instead
	// of erroring.
	LigateWarn bool
}

// Concat reads VCF input from readers in order, merges their headers, and
// emits a single VCF/BCF stream to out. The slice of (path, reader) pairs is
// passed so error messages can reference filenames; ConcatFiles is the more
// common entry point.
func Concat(inputs []NamedReader, out io.Writer, opts ConcatOptions) (int, error) {
	if len(inputs) == 0 {
		return 0, fmt.Errorf("bcftools concat: no input files")
	}
	// Upstream vcfconcat.c rejects -D/-d unless -a is in effect:
	//   "The -D option is supported only with -a".
	if (opts.RemoveDuplicates || opts.RmDupMode != "") && !opts.AllowOverlaps {
		return 0, fmt.Errorf("The -D option is supported only with -a")
	}
	// Read each input fully into memory: bcftools concat already requires
	// inputs to be sorted, and we operate on the in-memory variants to
	// implement sort-merge and de-duplication cleanly. For v1 this is the
	// right trade-off; a streaming sort-merge can land later.
	heads := make([]*vcf.Header, 0, len(inputs))
	groups := make([][]*vcf.Variant, 0, len(inputs))
	for _, in := range inputs {
		h, vs, err := readAllVariants(in.Reader)
		if err != nil {
			return 0, fmt.Errorf("bcftools concat: %s: %w", in.Name, err)
		}
		heads = append(heads, h)
		groups = append(groups, vs)
	}

	merged, err := MergeHeaders(heads)
	if err != nil {
		return 0, err
	}

	// Build the merged record stream.
	var records []*vcf.Variant
	switch {
	case opts.Ligate:
		// Phased concat: ligate overlapping phased chunks and inject PS/PQ.
		ensureLigateHeaders(merged)
		records, err = ligateConcat(merged, groups, opts)
		if err != nil {
			return 0, fmt.Errorf("bcftools concat: %w", err)
		}
	case opts.AllowOverlaps:
		// Upstream's `-a` traverses the genome through a synced reader
		// whose contig order is the first-seen-in-data order across the
		// inputs (the index/header seqnames of each file, appended in
		// file order), NOT the merged-header ##contig declaration order.
		// See vcfconcat.c (init_data builds out_hdr via bcf_hdr_merge but
		// the synced reader region list is populated from each reader's
		// tbx_seqnames in synced_bcf_reader.c add_reader).
		order := dataContigOrder(groups)
		mode := opts.RmDupMode
		if opts.RemoveDuplicates && mode == "" {
			mode = "exact"
		}
		records = mergeSortedDedup(groups, order, mode)
	default:
		for _, g := range groups {
			records = append(records, g...)
		}
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
	}, merged)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}
	for _, v := range records {
		if err := w.Write(v); err != nil {
			return 0, err
		}
	}
	return len(records), w.Flush()
}

// NamedReader pairs an io.Reader with a human-friendly name used in error
// messages.
type NamedReader struct {
	Name   string
	Reader io.Reader
}

// ConcatFiles is the file-aware entry point used by the CLI. It opens each
// path through iohelper, threads them into Concat, and adds support for the
// `-f file-list` option.
func ConcatFiles(paths []string, out io.Writer, opts ConcatOptions) (int, error) {
	all := paths
	if opts.FileList != "" {
		extra, err := ReadFileList(opts.FileList)
		if err != nil {
			return 0, fmt.Errorf("bcftools concat: %w", err)
		}
		all = append(extra, paths...)
	}
	if len(all) == 0 {
		return 0, fmt.Errorf("bcftools concat: no input files")
	}
	readers := make([]NamedReader, 0, len(all))
	closes := make([]io.Closer, 0, len(all))
	defer func() {
		for _, c := range closes {
			_ = c.Close()
		}
	}()
	for _, p := range all {
		f, err := iohelper.OpenReader(p)
		if err != nil {
			return 0, fmt.Errorf("bcftools concat: open %s: %w", p, err)
		}
		closes = append(closes, f)
		readers = append(readers, NamedReader{Name: p, Reader: f})
	}
	return Concat(readers, out, opts)
}

// ReadFileList reads the inputs listed in a `-f file-list.txt` file.
// Blank lines and lines beginning with '#' are skipped.
func ReadFileList(path string) ([]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	var out []string
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		out = append(out, line)
	}
	return out, sc.Err()
}

// readAllVariants returns the header and complete record list from a VCF or
// BCF reader. BCF inputs are converted to vcf.Variant via the BCF decoder.
func readAllVariants(in io.Reader) (*vcf.Header, []*vcf.Variant, error) {
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return nil, nil, err
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		return readAllBCF(br)
	}
	return readAllVCF(br)
}

// readAllVCF reads every record from a VCF stream.
func readAllVCF(in io.Reader) (*vcf.Header, []*vcf.Variant, error) {
	r := vcf.NewReader(in)
	hdr, err := r.ReadHeader()
	if err != nil {
		return nil, nil, err
	}
	var out []*vcf.Variant
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return hdr, out, err
		}
		out = append(out, v)
	}
	return hdr, out, nil
}

// readAllBCF reads every record from a BCF stream and returns the embedded
// vcf.Header plus the variants.
func readAllBCF(in io.Reader) (*vcf.Header, []*vcf.Variant, error) {
	br, err := bcf.NewReader(in)
	if err != nil {
		return nil, nil, err
	}
	hdr := br.Header()
	var out []*vcf.Variant
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return hdr.VCF, out, err
		}
		out = append(out, rec.ToVariant(hdr))
	}
	return hdr.VCF, out, nil
}

// MergeHeaders returns a single vcf.Header that is the union of the inputs.
// The rules implemented:
//
//   - The first header's fileformat / verbatim ordering is preserved when
//     possible.
//   - ##contig lines are union-merged in first-seen order.
//   - ##INFO / ##FORMAT / ##FILTER lines are union-merged by ID; if the same
//     ID appears with a different definition the function returns an error.
//   - Other meta lines are de-duplicated by exact-string equality.
//   - The sample list is the intersection-preserving union of all inputs;
//     headers with differing sample sets yield an error.
func MergeHeaders(heads []*vcf.Header) (*vcf.Header, error) {
	if len(heads) == 0 {
		return nil, fmt.Errorf("bcftools concat: no headers to merge")
	}
	// Track sample sets across all headers; first non-empty list wins as
	// the canonical order.
	var canonicalSamples []string
	for _, h := range heads {
		if len(h.Samples) == 0 {
			continue
		}
		if canonicalSamples == nil {
			canonicalSamples = append([]string{}, h.Samples...)
			continue
		}
		if !sameStringSlice(canonicalSamples, h.Samples) {
			return nil, fmt.Errorf("bcftools concat: input sample sets differ: %v vs %v", canonicalSamples, h.Samples)
		}
	}

	out := &vcf.Header{Samples: canonicalSamples}
	seenLine := make(map[string]bool)
	definitions := make(map[string]string) // "INFO/ID" -> raw line
	// fileformat must be the first line; track and emit it once.
	var fileformat string
	for _, h := range heads {
		for _, m := range h.MetaInfo {
			if strings.HasPrefix(m, "##fileformat=") && fileformat == "" {
				fileformat = m
				continue
			}
		}
	}
	if fileformat != "" {
		out.MetaInfo = append(out.MetaInfo, fileformat)
		seenLine[fileformat] = true
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
						return nil, fmt.Errorf("bcftools concat: conflicting %s ID %q definitions:\n  %s\n  %s", key, structID, prev, m)
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
	return out, nil
}

// structuredID returns the category ("INFO", "FORMAT", "FILTER", "contig",
// "ALT", ...) and the structured-line ID for meta lines of the form
// `##KEY=<ID=name,...>`. For non-structured lines it returns "".
func structuredID(line string) (kind, id string) {
	if !strings.HasPrefix(line, "##") {
		return "", ""
	}
	eq := strings.IndexByte(line, '=')
	if eq < 0 {
		return "", ""
	}
	key := line[2:eq]
	rest := line[eq+1:]
	if len(rest) < 2 || rest[0] != '<' {
		return "", ""
	}
	body := strings.TrimSuffix(rest[1:], ">")
	parts := splitTopLevel(body)
	for _, p := range parts {
		eq := strings.IndexByte(p, '=')
		if eq < 0 {
			continue
		}
		if p[:eq] == "ID" {
			return key, strings.TrimSpace(p[eq+1:])
		}
	}
	return "", ""
}

// splitTopLevel splits a structured meta body on commas while honouring
// double-quoted strings so descriptions can contain commas.
func splitTopLevel(s string) []string {
	var out []string
	var cur strings.Builder
	inQ := false
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c == '"' {
			inQ = !inQ
			cur.WriteByte(c)
			continue
		}
		if c == ',' && !inQ {
			out = append(out, cur.String())
			cur.Reset()
			continue
		}
		cur.WriteByte(c)
	}
	if cur.Len() > 0 {
		out = append(out, cur.String())
	}
	return out
}

// sameStringSlice reports whether a and b contain the same strings in the
// same order.
func sameStringSlice(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// contigOrder returns a map of contig name -> index from the ##contig lines
// in hdr. Contigs not declared in the header sort after declared ones, in
// lexical order (matching upstream behaviour for missing entries).
func contigOrder(hdr *vcf.Header) map[string]int {
	out := make(map[string]int)
	idx := 0
	for _, m := range hdr.MetaInfo {
		if !strings.HasPrefix(m, "##contig=") {
			continue
		}
		_, id := structuredID(m)
		if id == "" {
			continue
		}
		if _, ok := out[id]; ok {
			continue
		}
		out[id] = idx
		idx++
	}
	return out
}

// dataContigOrder returns a contig name -> rank map reproducing the
// first-seen-in-data ordering that upstream's `-a`/`--allow-overlaps`
// synced reader uses. Contigs are ranked in the order they first appear in
// the record stream, scanning the input groups in file order and, within a
// file, in record order. This matches synced_bcf_reader.c, which seeds its
// region list from each reader's index (or header) seqnames appended in
// file order, NOT from the merged-header ##contig declaration order.
func dataContigOrder(groups [][]*vcf.Variant) map[string]int {
	order := make(map[string]int)
	idx := 0
	for _, g := range groups {
		for _, v := range g {
			if _, seen := order[v.Chrom]; !seen {
				order[v.Chrom] = idx
				idx++
			}
		}
	}
	return order
}

// mergeSortedDedup performs the stable n-way merge that `concat -a` emits.
//
// Records are ordered by (data-contig rank, POS). At an equal (contig, POS)
// the file index is the tiebreaker: every record of the lower-numbered file
// is emitted before any record of a higher-numbered file, and within a file
// the original record order is preserved. This reproduces the synced
// reader's reader-by-reader emission at a shared position.
//
// When mode is non-empty (`-D`/`-d`) duplicate removal is applied across
// files only: at a given (contig, POS) a record from file i>0 is dropped
// when an already-emitted record from an earlier file at the same position
// collapses with it under the mode's rule. Records from the same file are
// never de-duplicated against each other, matching upstream (the synced
// reader's per-position `break` skips only the lines of later readers).
func mergeSortedDedup(groups [][]*vcf.Variant, order map[string]int, mode string) []*vcf.Variant {
	// rank ranks a record by (contig rank, POS) for the merge.
	rank := func(v *vcf.Variant) (int, int) {
		ci, ok := order[v.Chrom]
		if !ok {
			ci = 1<<30 + sortFallback(v.Chrom)
		}
		return ci, v.Pos
	}

	// Flatten all inputs, tagging each record with its source file index.
	// Upstream's synced reader visits each file through its index in the
	// region-traversal order, so a file whose records are not physically in
	// coordinate order is still consumed in coordinate order. Flattening and
	// stable-sorting reproduces that without an index.
	total := 0
	for _, g := range groups {
		total += len(g)
	}
	all := make([]taggedRec, 0, total)
	for fi, g := range groups {
		for _, v := range g {
			ci, pos := rank(v)
			all = append(all, taggedRec{v: v, file: fi, ci: ci, pos: pos})
		}
	}
	// Stable sort by (contig rank, POS). Stability preserves, at a shared
	// position, the lower-numbered file's records first and each file's
	// original record order — matching the synced reader's reader-by-reader
	// emission.
	sort.SliceStable(all, func(i, j int) bool {
		if all[i].ci != all[j].ci {
			return all[i].ci < all[j].ci
		}
		if all[i].pos != all[j].pos {
			return all[i].pos < all[j].pos
		}
		return all[i].file < all[j].file
	})

	out := make([]*vcf.Variant, 0, total)
	for i := 0; i < len(all); {
		// Walk one (contig, POS) group.
		j := i
		for j < len(all) && all[j].ci == all[i].ci && all[j].pos == all[i].pos {
			j++
		}
		// Within the group, keep all records of the first file, then for
		// each later file keep only the records that do not collapse with a
		// record already kept by a strictly-earlier file. This mirrors the
		// synced reader, which at a shared position writes only the lowest
		// reader and skips matching lines of later readers. Records from the
		// same file are never de-duplicated against each other.
		var earlier []*vcf.Variant  // kept records from files < curFile
		var sameFile []*vcf.Variant // kept records from curFile
		var kept []*vcf.Variant     // kept records in file order
		curFile := -1
		for k := i; k < j; k++ {
			rec := all[k]
			if rec.file != curFile {
				earlier = append(earlier, sameFile...)
				sameFile = sameFile[:0]
				curFile = rec.file
			}
			if mode != "" && len(earlier) > 0 && collapsesWithAny(rec.v, earlier, mode) {
				continue
			}
			kept = append(kept, rec.v)
			sameFile = append(sameFile, rec.v)
		}
		// Reorder the kept records into the synced-reader emission order: the
		// records are grouped by their unique REF>ALT variant string and the
		// groups are emitted by descending pre-dedup record count, ties broken
		// by first-appearance order. See orderPositionGroup.
		out = append(out, orderPositionGroup(all[i:j], kept)...)
		i = j
	}
	return out
}

// orderPositionGroup reorders the kept records at a single (contig, POS) into
// the order htslib's synced reader (bcf_sr_sort.c) emits them under
// concat -a. The reader buckets the records of every active file by their
// variant string ("REF>ALT,REF>ALT,…", with per-file duplicate disambiguation)
// and emits the buckets by descending count of records that share the variant
// (var->nvcf, computed across all files BEFORE de-duplication), breaking ties
// by the order in which the variant was first encountered (files in order,
// records within a file in order). Within a bucket the kept records keep their
// file/record order. The `all` slice is the full (pre-dedup) set of tagged
// records at this position, in (file, record) order; `kept` is the subset that
// survived de-duplication, also in (file, record) order.
func orderPositionGroup(all []taggedRec, kept []*vcf.Variant) []*vcf.Variant {
	if len(kept) <= 1 {
		return kept
	}
	// Build variant groups over ALL records (pre-dedup) to get first-appearance
	// order and the count used for ordering. Duplicate variant strings within
	// the same file are disambiguated with a trailing counter, matching
	// upstream's var_str2int handling.
	type vgroup struct {
		key   string
		order int // first-appearance index
		count int // pre-dedup record count (var->nvcf)
	}
	groupByKey := map[string]*vgroup{}
	var order []*vgroup
	perFileSeen := map[string]int{} // (file:key) -> times seen, for disambiguation
	keyOf := func(rec taggedRec) string {
		base := concatVariantKey(rec.v)
		fk := strconv.Itoa(rec.file) + "\x00" + base
		n := perFileSeen[fk]
		perFileSeen[fk] = n + 1
		if n > 0 {
			return base + "\x01" + strconv.Itoa(n)
		}
		return base
	}
	keyByRecPtr := map[*vcf.Variant]string{}
	for _, rec := range all {
		k := keyOf(rec)
		keyByRecPtr[rec.v] = k
		g := groupByKey[k]
		if g == nil {
			g = &vgroup{key: k, order: len(order)}
			groupByKey[k] = g
			order = append(order, g)
		}
		g.count++
	}
	// Bucket the kept records by their variant key, preserving file/record
	// order within a bucket.
	bucket := map[string][]*vcf.Variant{}
	for _, v := range kept {
		k := keyByRecPtr[v]
		bucket[k] = append(bucket[k], v)
	}
	// Emit buckets by descending count, ties broken by first-appearance order
	// — a stable sort over the first-appearance-ordered group list.
	sort.SliceStable(order, func(i, j int) bool {
		return order[i].count > order[j].count
	})
	out := make([]*vcf.Variant, 0, len(kept))
	for _, g := range order {
		out = append(out, bucket[g.key]...)
	}
	return out
}

// taggedRec is one record tagged with its source file index, used by
// mergeSortedDedup and orderPositionGroup.
type taggedRec struct {
	v    *vcf.Variant
	file int
	ci   int
	pos  int
}

// concatVariantKey builds the "REF>ALT,REF>ALT,…" string htslib's synced reader uses
// to group records at a shared position (bcf_sr_sort.c). A no-ALT record is
// keyed as "REF>.".
func concatVariantKey(v *vcf.Variant) string {
	if len(v.Alt) == 0 {
		return v.Ref + ">."
	}
	var b strings.Builder
	for i, a := range v.Alt {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(v.Ref)
		b.WriteByte('>')
		b.WriteString(a)
	}
	return b.String()
}

// collapsesWithAny reports whether v collapses with any record in prev
// under the given `-d`/`-D` mode.
func collapsesWithAny(v *vcf.Variant, prev []*vcf.Variant, mode string) bool {
	for _, p := range prev {
		if collapses(v, p, mode) {
			return true
		}
	}
	return false
}

// collapses implements the per-pair duplicate test for the `-d`/`-D`
// collapse modes recognised by upstream vcfconcat.c, which maps them to the
// synced reader's BCF_SR_PAIR_* logic (bcf_sr_sort.c::pairing_score). Both
// records are assumed to share (CHROM, POS); only the allele content is
// compared.
//
// In every non-exact mode an identical REF+ALT pair still collapses
// (pairing_score returns the best score for an exact match regardless of the
// type flags). The type-specific modes add further pairings on top:
//
//   - "exact"/"none": only identical REF and ALT (the `-D` default).
//   - "snps":  identical, or both records are SNVs/MNVs.
//   - "indels": identical, or both records are indels.
//   - "both":  identical, or both SNV/MNV, or both indels.
//   - "any"/"all": always collapse (any record at the position).
func collapses(a, b *vcf.Variant, mode string) bool {
	switch mode {
	case "any", "all":
		return true
	case "snps":
		return sameVariant(a, b) || (isSNPRecord(a) && isSNPRecord(b))
	case "indels":
		return sameVariant(a, b) || (isIndelRecord(a) && isIndelRecord(b))
	case "both":
		return sameVariant(a, b) ||
			(isSNPRecord(a) && isSNPRecord(b)) ||
			(isIndelRecord(a) && isIndelRecord(b))
	default: // "exact", "none"
		return sameVariant(a, b)
	}
}

// lessVariant orders two variants by (contig-index, POS, REF, ALT).
func lessVariant(a, b *vcf.Variant, order map[string]int) bool {
	ai, aok := order[a.Chrom]
	bi, bok := order[b.Chrom]
	if !aok {
		ai = 1<<30 + sortFallback(a.Chrom)
	}
	if !bok {
		bi = 1<<30 + sortFallback(b.Chrom)
	}
	if ai != bi {
		return ai < bi
	}
	if a.Pos != b.Pos {
		return a.Pos < b.Pos
	}
	if a.Ref != b.Ref {
		return a.Ref < b.Ref
	}
	aAlt := strings.Join(a.Alt, ",")
	bAlt := strings.Join(b.Alt, ",")
	return aAlt < bAlt
}

// sortFallback gives a deterministic but pretty arbitrary ordering for
// contigs that are not declared in the merged header. We use a hash of the
// name so the order is stable across runs.
func sortFallback(name string) int {
	// FNV-1a hash compressed to a small positive int. Good enough for
	// deterministic ordering of unknown contigs.
	const (
		offset = 2166136261
		prime  = 16777619
	)
	h := uint32(offset)
	for i := 0; i < len(name); i++ {
		h ^= uint32(name[i])
		h *= prime
	}
	return int(h & 0x0fffffff)
}

// sameVariant returns true if a and b share (CHROM, POS, REF, ALT).
func sameVariant(a, b *vcf.Variant) bool {
	if a.Chrom != b.Chrom || a.Pos != b.Pos || a.Ref != b.Ref {
		return false
	}
	if len(a.Alt) != len(b.Alt) {
		return false
	}
	// Compare ALT sets order-insensitively by sorting copies first; concat
	// can see records emitted by tools that re-order alleles.
	aa := append([]string{}, a.Alt...)
	bb := append([]string{}, b.Alt...)
	sort.Strings(aa)
	sort.Strings(bb)
	for i := range aa {
		if aa[i] != bb[i] {
			return false
		}
	}
	return true
}
