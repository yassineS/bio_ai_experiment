package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
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
	// RemoveDuplicates collapses adjacent records sharing the same
	// (CHROM, POS, REF, ALT) tuple. Only kept records contribute to the
	// returned count.
	RemoveDuplicates bool
	// FileList, when non-empty, is read line-by-line and prepended to the
	// list of inputs supplied to ConcatFiles.
	FileList string
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
	// MinPQ is accepted but ignored in v1 (matches upstream's `-q`).
	MinPQ int
	// Ligate is accepted but ignored in v1 (matches upstream's `-l`).
	Ligate bool
}

// Concat reads VCF input from readers in order, merges their headers, and
// emits a single VCF/BCF stream to out. The slice of (path, reader) pairs is
// passed so error messages can reference filenames; ConcatFiles is the more
// common entry point.
func Concat(inputs []NamedReader, out io.Writer, opts ConcatOptions) (int, error) {
	if len(inputs) == 0 {
		return 0, fmt.Errorf("bcftools concat: no input files")
	}
	// Match upstream vcfconcat.c:1120 — `-D` is meaningful only when the
	// inputs are first sort-merged with `-a`, so reject `-D` standalone
	// rather than silently doing a stream-level adjacency dedupe (which
	// upstream simply never offers).
	if opts.RemoveDuplicates && !opts.AllowOverlaps {
		return 0, fmt.Errorf("bcftools concat: The -D option is supported only with -a")
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
	if opts.AllowOverlaps {
		// Match upstream's synced-reader behaviour: the contig sort key
		// is the order in which contigs first appear in any input's
		// *data* (mirroring how htslib's tabix/BCF index seq-names are
		// added to the synced reader's regions table), not the order
		// the contigs are declared in the merged header. See
		// reference_code/htslib/synced_bcf_reader.c:392-414 and
		// _regions_add (1028+): the first reader's index seq-names are
		// added in order, then later readers add only the previously
		// unseen ones at the end. We approximate that here without an
		// index by walking each input's record stream in order and
		// recording first-seen chromosomes.
		records = mergeSorted(groups, firstSeenContigOrder(groups))
	} else {
		for _, g := range groups {
			records = append(records, g...)
		}
	}
	if opts.RemoveDuplicates {
		records = dedupAdjacent(records)
	}

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
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

// firstSeenContigOrder returns a contig-name -> rank map matching
// upstream's synced-reader contig ordering. The synced reader walks the
// inputs in command-line order, asks each one's tabix/BCF index for its
// seq-names list (which, for a sorted-then-indexed file, reflects the
// order chromosomes first appear in the data), and adds them to its
// regions table; previously-seen names keep their original rank. We
// approximate the same ordering from the in-memory record stream: walk
// every group in order, and record the rank at which each new chromosome
// first appears across all of them. This matches upstream's observable
// concat -a contig order for the common cases the synced reader was
// designed for and, crucially, does NOT depend on the merged header's
// ##contig declaration order.
func firstSeenContigOrder(groups [][]*vcf.Variant) map[string]int {
	out := make(map[string]int)
	idx := 0
	for _, g := range groups {
		for _, v := range g {
			if _, ok := out[v.Chrom]; ok {
				continue
			}
			out[v.Chrom] = idx
			idx++
		}
	}
	return out
}

// mergeSorted mirrors upstream's synced-reader concat behaviour: emit one
// contig at a time in `order` rank order, and within each contig do a POS
// n-way merge across whichever input groups contain records for that
// contig.
//
// Unlike a single global n-way merge, this matches upstream even when an
// input file's records are not globally sorted by (contig, POS) — e.g.
// the upstream `concat.2.a.vcf` fixture lists chr2 records first, then
// chr1 records. The synced reader handles that by seeking each input to
// the current contig (via the tabix index); without an index we simulate
// the same effect by partitioning records by contig first.
func mergeSorted(groups [][]*vcf.Variant, order map[string]int) []*vcf.Variant {
	total := 0
	for _, g := range groups {
		total += len(g)
	}
	out := make([]*vcf.Variant, 0, total)

	// Partition each group by chromosome, preserving the original
	// per-chromosome record order (upstream relies on each input being
	// sorted within a chromosome, which is the contract of `concat -a`).
	type bucket struct {
		records []*vcf.Variant
	}
	perGroup := make([]map[string]*bucket, len(groups))
	for gi, g := range groups {
		perGroup[gi] = make(map[string]*bucket)
		for _, v := range g {
			b, ok := perGroup[gi][v.Chrom]
			if !ok {
				b = &bucket{}
				perGroup[gi][v.Chrom] = b
			}
			b.records = append(b.records, v)
		}
	}

	// Build contig emission order: ranked names from `order`, followed
	// by any contigs we saw in the data but aren't in `order` (rare:
	// happens when order was built from a different source than the
	// records; the firstSeenContigOrder caller covers every chromosome
	// it walked, so this is mostly defensive).
	allChroms := make(map[string]struct{})
	for gi := range groups {
		for c := range perGroup[gi] {
			allChroms[c] = struct{}{}
		}
	}
	ranked := make([]string, 0, len(allChroms))
	for c := range allChroms {
		ranked = append(ranked, c)
	}
	sort.Slice(ranked, func(i, j int) bool {
		ri, ok := order[ranked[i]]
		if !ok {
			ri = 1<<30 + sortFallback(ranked[i])
		}
		rj, ok := order[ranked[j]]
		if !ok {
			rj = 1<<30 + sortFallback(ranked[j])
		}
		if ri != rj {
			return ri < rj
		}
		return ranked[i] < ranked[j]
	})

	// Emit one contig at a time. Within a contig advance position-min
	// across the per-group buckets; at each tied position, emit records
	// grouped by their (REF,ALT) signature, with signatures ordered by
	// the order in which they first appear across the inputs. This
	// mirrors upstream's bcf_sr_sort behaviour (see
	// reference_code/htslib/bcf_sr_sort.c:bcf_sr_sort_next).
	cursors := make([]int, len(groups))
	for _, chrom := range ranked {
		slices := make([][]*vcf.Variant, len(groups))
		for gi := range groups {
			cursors[gi] = 0
			if b, ok := perGroup[gi][chrom]; ok {
				slices[gi] = b.records
			}
		}
		for {
			// Find the minimum position across all groups at their
			// current cursor.
			var minPos int = -1
			haveAny := false
			for gi, s := range slices {
				if cursors[gi] >= len(s) {
					continue
				}
				p := s[cursors[gi]].Pos
				if !haveAny || p < minPos {
					minPos = p
					haveAny = true
				}
			}
			if !haveAny {
				break
			}
			// Collect all records at minPos from each group, ordered
			// by signature-first-seen.
			type entry struct {
				gi  int
				rec *vcf.Variant
			}
			var entries []entry
			for gi, s := range slices {
				for cursors[gi] < len(s) && s[cursors[gi]].Pos == minPos {
					entries = append(entries, entry{gi: gi, rec: s[cursors[gi]]})
					cursors[gi]++
				}
			}
			// Build signature order: first appearance in the
			// `entries` slice (which is itself ordered by group
			// then by in-file order — matching upstream's
			// command-line-order traversal).
			sigOrder := make(map[string]int)
			for _, e := range entries {
				sig := e.rec.Ref + "\x00" + strings.Join(e.rec.Alt, ",")
				if _, ok := sigOrder[sig]; !ok {
					sigOrder[sig] = len(sigOrder)
				}
			}
			// Emit: for each signature in first-seen order, emit
			// all entries with that signature in group order
			// (entries is already group-ordered).
			emitted := make([]bool, len(entries))
			for sig := 0; sig < len(sigOrder); sig++ {
				for i, e := range entries {
					if emitted[i] {
						continue
					}
					thisSig := e.rec.Ref + "\x00" + strings.Join(e.rec.Alt, ",")
					if sigOrder[thisSig] != sig {
						continue
					}
					out = append(out, e.rec)
					emitted[i] = true
				}
			}
		}
	}
	return out
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

// dedupAdjacent removes adjacent variants whose (CHROM, POS, REF, ALT) tuple
// matches the previous record. Inputs must already be ordered for this to
// behave like `--remove-duplicates`.
func dedupAdjacent(in []*vcf.Variant) []*vcf.Variant {
	if len(in) <= 1 {
		return in
	}
	out := make([]*vcf.Variant, 0, len(in))
	out = append(out, in[0])
	for i := 1; i < len(in); i++ {
		if sameVariant(in[i], out[len(out)-1]) {
			continue
		}
		out = append(out, in[i])
	}
	return out
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
