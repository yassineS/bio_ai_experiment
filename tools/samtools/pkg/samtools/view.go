package samtools

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// ViewOptions configures the behaviour of View. Zero values disable the
// corresponding filter.
type ViewOptions struct {
	// OutputBAM forces BAM output. When false, output is text SAM unless the
	// caller pre-wraps the output writer in a BAM writer.
	OutputBAM bool
	// OutputCRAM forces CRAM output (samtools view -C / --output-fmt cram).
	// When set it takes precedence over OutputBAM. The CRAM the writer
	// emits is reference-free, so it decodes without an external FASTA.
	OutputCRAM bool
	// WithHeader emits the header alongside records.
	WithHeader bool
	// HeaderOnly emits the header and skips all records.
	HeaderOnly bool
	// Count returns just the record count (no record output).
	Count bool
	// IncludeFlags keeps only records whose flag has ALL bits in IncludeFlags set.
	IncludeFlags uint16
	// ExcludeFlags drops records whose flag has ANY bits in ExcludeFlags set.
	ExcludeFlags uint16
	// ExcludeFlagsAll drops records only if ALL bits in ExcludeFlagsAll are set.
	ExcludeFlagsAll uint16
	UseExcludeAll   bool
	// MinMAPQ keeps records with MAPQ ≥ MinMAPQ. 0 disables.
	MinMAPQ uint8
	// ReadGroup, if non-empty, restricts output to records whose RG:Z tag
	// matches this value.
	ReadGroup string
	// ReadGroupSet, if non-nil, restricts output to records whose RG:Z tag is
	// in the set. Combined with ReadGroup via union.
	ReadGroupSet map[string]struct{}
	// Subsample is the keep-fraction in (0, 1]. 0 or 1 means no subsampling.
	Subsample float64
	// SubsampleSeed seeds the RNG for reproducible subsampling.
	SubsampleSeed int64
	// Regions are CLI-style "chr:start-end" specifiers (kept for diagnostics —
	// see RegionsEnabled).
	Regions []string
	// RegionsEnabled signals whether region-query filtering is requested.
	// In this first slice it is always rejected, with a clear error message.
	RegionsEnabled bool
	// BedPath, when non-empty, is the path of a BED file whose intervals
	// restrict the emitted records. A record is kept only when its
	// [Pos, Pos+refLen) range intersects at least one BED interval on the
	// record's RName. Combines (logical AND) with Regions: when both are
	// set, records must satisfy BOTH predicates. Matches upstream samtools
	// view's `-L/--regions-file`.
	BedPath string
	// MultiRegion is accept-and-ignore for upstream's
	// `-M/--use-multi-region-iterator`. The flag controls whether upstream
	// uses its multi-region iterator (an indexed-scan optimisation); the
	// resulting filtered record set is identical, so we always perform the
	// full intersection regardless.
	MultiRegion bool
	// NoPG suppresses the @PG line emission (placeholder; the view pipeline
	// does not currently inject an @PG line, so this is a no-op kept for
	// flag compatibility).
	NoPG bool
	// TagFilters is the conjunction of aux-tag predicates derived from
	// `-d TAG[:VAL]` and `-D TAG:FILE`. A record is kept only when every
	// filter matches. Empty means no tag filtering.
	TagFilters []TagFilter
	// QNameSet, if non-nil, restricts output to records whose QNAME is in
	// the set (or, when QNameInvert is true, NOT in the set). Populated
	// from `-N FILE` or `-N ^FILE`.
	QNameSet map[string]struct{}
	// QNameInvert flips the QName filter to exclude-mode, matching
	// upstream samtools' `-N ^FILE` syntax (sam_view.c:352
	// rnhash_discard=1).
	QNameInvert bool
	// Reference names a FASTA file (with a sibling .fai) used as the decode
	// reference for reference-backed CRAM input — upstream samtools view's
	// `-T/--reference`. It is ignored for SAM and BAM input. When empty a
	// reference-backed CRAM still decodes, with reference-derived bases
	// filled with 'N'.
	Reference string
}

// TagFilter is a single aux-tag predicate as derived from samtools view's
// `-d`/`-D` flags. ExistsOnly true means "record must carry tag Tag";
// otherwise the record must carry Tag AND its stringified value must
// appear in Values.
//
// The string comparison mirrors upstream samtools (sam_view.c
// process_aln): integer aux values are formatted as base-10 decimals,
// 'A'-type values use the single character, and Z/H types use the raw
// string. Float / array tags are not currently exposed for matching
// because upstream does not stringify them in process_aln either.
type TagFilter struct {
	Tag        string
	ExistsOnly bool
	Values     map[string]struct{}
}

// ErrRegionsUnsupported is retained for backward-compatibility with callers
// that match against it. Region filtering via linear scan is now supported
// directly by View; indexed seek is available via ViewFile when a .bai
// sibling file exists.
//
// Deprecated: regions are now handled.
var ErrRegionsUnsupported = errors.New("samtools view: region-query support requires .bai indexing — not yet implemented")

// View streams alignment records from r, applying the configured filters,
// and emits them to out as either SAM text or BAM. Returns the number of
// records that passed all filters.
//
// When opts.Count is true, the matching count is the only thing emitted to
// out. When opts.Regions is non-empty View does a linear scan and filters
// records to those overlapping any region — for indexed seek use ViewFile.
func View(in io.Reader, out io.Writer, opts ViewOptions) (int, error) {
	r, err := alnio.NewReaderWithReference(in, opts.Reference)
	if err != nil {
		return 0, err
	}
	hdr := r.Header()

	// Resolve any region specifiers up front so unknown chroms don't abort
	// the linear scan — they just yield zero matches.
	resolved, _, err := region.ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
	if err != nil {
		return 0, err
	}
	regionFilter := buildRegionFilter(resolved, hdr)
	if regionFilter == nil && len(opts.Regions) > 0 {
		// User asked for regions but none resolved — surface a predicate
		// that rejects every record so the result is correctly empty
		// instead of accidentally returning the full stream.
		regionFilter = func(*sam.Record) bool { return false }
	}

	bedFilter, err := loadBedFilter(opts.BedPath)
	if err != nil {
		return 0, err
	}

	w, err := openViewWriter(out, hdr, opts)
	if err != nil {
		return 0, err
	}
	if opts.HeaderOnly {
		return 0, closeViewWriter(w)
	}

	rng := newSubsampleRNG(opts)
	matched := 0
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return matched, err
		}
		if !keepRecord(rec, &opts, rng) {
			continue
		}
		if regionFilter != nil && !regionFilter(rec) {
			continue
		}
		if bedFilter != nil && !bedFilter(rec) {
			continue
		}
		matched++
		if opts.Count {
			continue
		}
		if err := w.Write(rec); err != nil {
			return matched, err
		}
	}
	if opts.Count {
		fmt.Fprintln(out, matched)
		return matched, nil
	}
	if err := closeViewWriter(w); err != nil {
		return matched, err
	}
	return matched, nil
}

// ViewFile is the indexed entry point for samtools view: it opens the BAM
// at inPath, looks for a sibling <inPath>.bai, and (when found and the
// caller has supplied regions) uses the BAI chunks to seek to relevant
// portions of the file. When no .bai is found ViewFile falls back to the
// streaming View() linear-scan path; a warning is written to warnW (which
// may be nil to silence it).
//
// If inPath is empty or "-" ViewFile delegates to View on os.Stdin.
func ViewFile(inPath string, out io.Writer, opts ViewOptions, warnW io.Writer) (int, error) {
	if inPath == "" || inPath == "-" {
		return View(os.Stdin, out, opts)
	}
	// No CLI regions requested: take the streaming path. -L/--regions-file
	// (opts.BedPath) is handled inside View() as a post-filter and doesn't
	// require an index.
	if len(opts.Regions) == 0 {
		f, err := os.Open(inPath)
		if err != nil {
			return 0, err
		}
		defer f.Close()
		return View(f, out, opts)
	}
	baiPath := inPath + ".bai"
	baiBytes, baiErr := os.ReadFile(baiPath)
	if baiErr != nil {
		if warnW != nil {
			fmt.Fprintf(warnW, "samtools view: no index at %s, falling back to linear scan\n", baiPath)
		}
		f, err := os.Open(inPath)
		if err != nil {
			return 0, err
		}
		defer f.Close()
		return View(f, out, opts)
	}
	idx, ierr := bam.ReadBAI(strings.NewReader(string(baiBytes)))
	if ierr != nil {
		return 0, fmt.Errorf("samtools view: read %s: %w", baiPath, ierr)
	}

	f, err := os.Open(inPath)
	if err != nil {
		return 0, err
	}
	defer f.Close()
	return viewIndexed(f, idx, out, opts)
}

// viewIndexed performs an indexed region scan: it parses regions, computes
// chunk unions, seeks into each chunk's compressed offset, decodes records
// until each chunk's end virtual offset, and emits records overlapping any
// requested region.
//
// This path is BAM-only: it does BGZF virtual-offset seeks against a .bai
// index. CRAM uses a .crai index and a different seek model; indexed CRAM
// region query is a separate roadmap item, so a CRAM file reaches the
// streaming path above, not here.
func viewIndexed(f *os.File, idx *bam.BAIIndex, out io.Writer, opts ViewOptions) (int, error) {
	// Need the header — open a BAM reader first.
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		return 0, err
	}
	hdrReader, err := sam.NewBAMReader(f)
	if err != nil {
		return 0, err
	}
	hdr := hdrReader.Header()

	resolved, _, perr := region.ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
	if perr != nil {
		return 0, perr
	}
	chunks := bam.UnionChunks(idx, resolved)
	regionFilter := buildRegionFilter(resolved, hdr)
	if regionFilter == nil && len(opts.Regions) > 0 {
		regionFilter = func(*sam.Record) bool { return false }
	}
	bedFilter, berr := loadBedFilter(opts.BedPath)
	if berr != nil {
		return 0, berr
	}

	w, werr := openViewWriter(out, hdr, opts)
	if werr != nil {
		return 0, werr
	}
	if opts.HeaderOnly {
		return 0, closeViewWriter(w)
	}

	rng := newSubsampleRNG(opts)
	matched := 0
	for _, c := range chunks {
		if c.Beg >= c.End {
			continue
		}
		startBlock := int64(c.Beg >> 16)
		if _, err := f.Seek(startBlock, io.SeekStart); err != nil {
			return matched, err
		}
		bgz, err := bgzip.NewReader(f)
		if err != nil {
			return matched, err
		}
		// Skip in-block bytes — they are *before* the first record we want.
		uoff := int(c.Beg & 0xFFFF)
		if uoff > 0 {
			if _, err := io.CopyN(io.Discard, bgz, int64(uoff)); err != nil {
				return matched, err
			}
		}
		// Build a chunk-bounded reader so the BAM parser stops when the
		// chunk ends; we use a custom wrapper that watches the underlying
		// virtual offset.
		boundedSrc := &chunkBoundedReader{r: bgz, end: uint64(c.End)}
		br := sam.NewBAMBodyReader(boundedSrc, hdr)
		for {
			rec, err := br.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return matched, err
			}
			if !keepRecord(rec, &opts, rng) {
				continue
			}
			if regionFilter != nil && !regionFilter(rec) {
				continue
			}
			if bedFilter != nil && !bedFilter(rec) {
				continue
			}
			matched++
			if opts.Count {
				continue
			}
			if err := w.Write(rec); err != nil {
				return matched, err
			}
		}
		_ = bgz.Close()
	}
	if opts.Count {
		fmt.Fprintln(out, matched)
		return matched, nil
	}
	if err := closeViewWriter(w); err != nil {
		return matched, err
	}
	return matched, nil
}

// chunkBoundedReader stops returning bytes once its underlying BGZF reader
// has advanced past the chunk-end virtual offset. We watch the virtual
// offset of the bgzip.Reader (cheap — a pair of fields) after every read
// and report io.EOF as soon as we cross the end boundary.
type chunkBoundedReader struct {
	r   *bgzip.Reader
	end uint64
}

func (c *chunkBoundedReader) Read(p []byte) (int, error) {
	if c.r.VirtualOffset() >= c.end {
		return 0, io.EOF
	}
	return c.r.Read(p)
}

// buildRegionFilter returns nil when no regions are configured; otherwise
// returns a predicate that keeps records overlapping any region's range
// on the matching reference.
func buildRegionFilter(regions []region.ResolvedRegion, hdr *sam.Header) func(*sam.Record) bool {
	if len(regions) == 0 {
		return nil
	}
	return func(rec *sam.Record) bool {
		if rec.RName == "" || rec.RName == "*" {
			return false
		}
		rid := hdr.RefIndex(rec.RName)
		if rid < 0 {
			return false
		}
		pos0 := int(rec.Pos) - 1
		if pos0 < 0 {
			pos0 = 0
		}
		refLen := rec.Cigar.ReferenceLength()
		if refLen <= 0 {
			refLen = 1
		}
		for _, r := range regions {
			if r.RefID != rid {
				continue
			}
			recEnd := pos0 + refLen
			if pos0 < r.End0 && recEnd > r.Beg0 {
				return true
			}
		}
		return false
	}
}

// openViewWriter sets up the output writer and emits the header when the
// configured options call for it. The returned writer is nil only when
// opts.Count is true (which suppresses all record-emitting writes).
func openViewWriter(out io.Writer, hdr *sam.Header, opts ViewOptions) (sam.Writer, error) {
	if opts.Count {
		return nil, nil
	}
	var w sam.Writer
	switch {
	case opts.OutputCRAM:
		w = alnio.NewCRAMWriter(out)
	case opts.OutputBAM:
		w = sam.NewBAMWriter(out)
	default:
		w = sam.NewSAMWriter(out)
	}
	// BAM and CRAM are binary container formats whose header is structural
	// — it is always emitted regardless of -h/-H. CRAM in particular
	// cannot be written without a header.
	emitHeader := opts.HeaderOnly || opts.WithHeader || opts.OutputBAM || opts.OutputCRAM
	if emitHeader {
		if err := w.WriteHeader(hdr); err != nil {
			return nil, err
		}
	} else {
		if err := w.WriteHeader(nil); err != nil {
			return nil, err
		}
	}
	return w, nil
}

// closeViewWriter flushes the writer if it is non-nil.
func closeViewWriter(w sam.Writer) error {
	if w == nil {
		return nil
	}
	return w.Close()
}

// keepRecord applies the per-record filters and returns true if the record
// passes them all.
func keepRecord(rec *sam.Record, opts *ViewOptions, rng *rand.Rand) bool {
	if opts.IncludeFlags != 0 && rec.Flag&opts.IncludeFlags != opts.IncludeFlags {
		return false
	}
	if opts.ExcludeFlags != 0 && rec.Flag&opts.ExcludeFlags != 0 {
		return false
	}
	if opts.UseExcludeAll && rec.Flag&opts.ExcludeFlagsAll == opts.ExcludeFlagsAll {
		return false
	}
	if opts.MinMAPQ > 0 && rec.MapQ < opts.MinMAPQ {
		return false
	}
	if opts.ReadGroup != "" || len(opts.ReadGroupSet) > 0 {
		rg, ok := rec.GetAux("RG")
		if !ok {
			return false
		}
		s, _ := rg.String()
		match := false
		if opts.ReadGroup != "" && s == opts.ReadGroup {
			match = true
		}
		if _, ok := opts.ReadGroupSet[s]; ok {
			match = true
		}
		if !match {
			return false
		}
	}
	if len(opts.QNameSet) > 0 {
		_, present := opts.QNameSet[rec.QName]
		// Upstream sam_view.c:221: `(kh_get(...) == kh_end(h)) != rnhash_discard`
		// — keep when (NOT present) XOR invert flips.
		if present == opts.QNameInvert {
			return false
		}
	}
	for i := range opts.TagFilters {
		if !matchTagFilter(rec, &opts.TagFilters[i]) {
			return false
		}
	}
	if rng != nil && opts.Subsample > 0 && opts.Subsample < 1 {
		if rng.Float64() >= opts.Subsample {
			return false
		}
	}
	return true
}

// matchTagFilter returns true when rec satisfies the supplied TagFilter
// (used by samtools view -d / -D). It mirrors upstream samtools'
// process_aln logic: missing tag → drop; existence-only → keep; value
// list → keep when the record's stringified aux value is in the list.
func matchTagFilter(rec *sam.Record, f *TagFilter) bool {
	a, ok := rec.GetAux(f.Tag)
	if !ok {
		return false
	}
	if f.ExistsOnly {
		return true
	}
	val, ok := auxValueAsString(a)
	if !ok {
		return false
	}
	_, hit := f.Values[val]
	return hit
}

// auxValueAsString returns the string form of an aux value as used for
// `-d TAG:VAL` comparison: integers go to base-10 decimal, 'A' uses the
// single character, Z/H use the raw string. Other types (f, B, ...) are
// not stringified — upstream samtools does not match against them in
// process_aln, so we report `false` and let the caller drop the record.
func auxValueAsString(a sam.Aux) (string, bool) {
	switch a.Type {
	case 'c', 'C', 's', 'S', 'i', 'I':
		if v, ok := a.Int(); ok {
			return strconv.FormatInt(v, 10), true
		}
	case 'A':
		if v, ok := a.Value.(string); ok && len(v) > 0 {
			return v[:1], true
		}
	case 'Z', 'H':
		if v, ok := a.String(); ok {
			return v, true
		}
	}
	return "", false
}

// LoadLinesFile reads a UTF-8 file of whitespace-separated tokens and
// returns them as a set. Mirrors upstream samtools' `fscanf("%1023s")`
// reader in `populate_lookup_from_file` (sam_view.c:294) which has no
// comment-line handling — tokens beginning with `#` are kept verbatim.
// Used by samtools view's `-D TAG:FILE` and `-N FILE` flags.
func LoadLinesFile(path string) (map[string]struct{}, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	out := map[string]struct{}{}
	scn := bufio.NewScanner(f)
	// Allow long lines (BAM qnames can run to a few hundred chars; default
	// 64 KiB is plenty but we bump anyway in case of pathological files).
	scn.Buffer(make([]byte, 0, 64*1024), 1024*1024)
	for scn.Scan() {
		for _, tok := range strings.Fields(scn.Text()) {
			out[tok] = struct{}{}
		}
	}
	return out, scn.Err()
}

// ParseTagFilterSpec parses a samtools view `-d TAG[:VAL]` argument and
// returns the corresponding TagFilter. The tag must be exactly two
// characters; an optional ":VAL" suffix becomes a single-element value
// set. Returns an error when the spec is malformed (too short, missing
// ":" after a 2-char tag prefix, etc.).
func ParseTagFilterSpec(spec string) (TagFilter, error) {
	switch {
	case len(spec) == 2:
		return TagFilter{Tag: spec, ExistsOnly: true}, nil
	case len(spec) < 4 || spec[2] != ':':
		return TagFilter{}, fmt.Errorf("samtools view: invalid -d tag:value spec %q", spec)
	}
	tag := spec[:2]
	val := spec[3:]
	return TagFilter{
		Tag:    tag,
		Values: map[string]struct{}{val: {}},
	}, nil
}

// ParseTagFileSpec parses a samtools view `-D TAG:FILE` argument and
// returns a TagFilter populated from the file (one value per line).
// Both `:` and `;` are accepted as separators, matching upstream's
// MinGW path-translation workaround in sam_view.c.
func ParseTagFileSpec(spec string) (TagFilter, error) {
	if len(spec) < 4 || (spec[2] != ':' && spec[2] != ';') {
		return TagFilter{}, fmt.Errorf("samtools view: invalid -D tag:file spec %q", spec)
	}
	tag := spec[:2]
	path := spec[3:]
	vals, err := LoadLinesFile(path)
	if err != nil {
		return TagFilter{}, err
	}
	return TagFilter{Tag: tag, Values: vals}, nil
}

// MergeTagFilter folds a new TagFilter into an existing slice. Upstream
// samtools rejects mixing different tags in the same command line; we do
// the same. When the new filter shares its tag with an existing entry
// the value sets are unioned. Once any value has been supplied (by any
// `-d TAG:VAL` or `-D TAG:FILE`) the merged filter becomes value-bound
// — a bare `-d TAG` does not relax that constraint, mirroring
// process_aln in sam_view.c which only ever consults tvhash once it is
// non-nil.
func MergeTagFilter(dst []TagFilter, add TagFilter) ([]TagFilter, error) {
	for i := range dst {
		if dst[i].Tag != add.Tag {
			return nil, fmt.Errorf("samtools view: different tag %q specified before %q", dst[i].Tag, add.Tag)
		}
		// Same tag — merge value sets.
		if add.ExistsOnly && dst[i].ExistsOnly {
			return dst, nil
		}
		if dst[i].Values == nil {
			dst[i].Values = map[string]struct{}{}
		}
		for v := range add.Values {
			dst[i].Values[v] = struct{}{}
		}
		// As soon as any side has constrained values, the filter is no
		// longer a pure existence check.
		dst[i].ExistsOnly = false
		return dst, nil
	}
	return append(dst, add), nil
}

// newSubsampleRNG seeds a *rand.Rand for the subsample filter. Returns nil
// when no subsampling is requested.
func newSubsampleRNG(opts ViewOptions) *rand.Rand {
	if opts.Subsample <= 0 || opts.Subsample >= 1 {
		return nil
	}
	seed := opts.SubsampleSeed
	if seed == 0 {
		seed = time.Now().UnixNano()
	}
	// #nosec G404 — non-cryptographic subsample by design.
	return rand.New(rand.NewSource(seed))
}

// LoadReadGroupsFile reads a file of read group IDs (one per line) and
// returns them as a set. Lines starting with '#' and blank lines are skipped.
func LoadReadGroupsFile(path string) (map[string]struct{}, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	out := map[string]struct{}{}
	scn := bufio.NewScanner(f)
	for scn.Scan() {
		line := strings.TrimSpace(scn.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		out[line] = struct{}{}
	}
	return out, scn.Err()
}

// loadBedFilter reads a BED file at path and returns a predicate that keeps
// records whose `[Pos, Pos+refLen)` half-open range intersects any BED
// interval on the record's RName. Returns (nil, nil) when path is empty.
//
// The implementation builds one bed.IntervalTree per chromosome so per-record
// queries cost O(log n + k). Records whose RName has no BED entries are
// dropped, matching upstream's `bed_overlap` behaviour.
func loadBedFilter(path string) (func(*sam.Record) bool, error) {
	if path == "" {
		return nil, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("samtools view: open BED %q: %w", path, err)
	}
	defer f.Close()
	byChrom := map[string][]*bed.Record{}
	rd := bed.NewReader(f)
	for {
		rec, rerr := rd.Read()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return nil, fmt.Errorf("samtools view: read BED %q: %w", path, rerr)
		}
		byChrom[rec.Chrom] = append(byChrom[rec.Chrom], rec)
	}
	if len(byChrom) == 0 {
		// Empty BED means "no records pass" — match upstream's behaviour
		// where bed_overlap returns 0 for every query.
		return func(*sam.Record) bool { return false }, nil
	}
	// Sort intervals by ChromStart so the tree is balanced (NewIntervalTree
	// requires its input to be sorted to produce a balanced tree).
	trees := make(map[string]*bed.IntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		sort.Slice(recs, func(i, j int) bool {
			return recs[i].ChromStart < recs[j].ChromStart
		})
		trees[chrom] = bed.NewIntervalTree(recs)
	}
	return func(rec *sam.Record) bool {
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" {
			return false
		}
		t, ok := trees[rec.RName]
		if !ok {
			return false
		}
		pos0 := int(rec.Pos) - 1
		if pos0 < 0 {
			pos0 = 0
		}
		refLen := rec.Cigar.ReferenceLength()
		if refLen <= 0 {
			// Zero-length footprint (CIGAR `*` or all-clip) cannot
			// overlap any BED interval per upstream bed_overlap's strict
			// half-open semantics. Drop the record.
			return false
		}
		q := &bed.Record{Chrom: rec.RName, ChromStart: pos0, ChromEnd: pos0 + refLen}
		return len(t.Query(q)) > 0
	}, nil
}

// ParseSubsample parses upstream samtools' "<seed>.<fraction>" composite form
// used by `-s`. A plain float (e.g. "0.5") is also accepted and returns a
// zero seed (time-based).
func ParseSubsample(s string) (frac float64, seed int64, err error) {
	if s == "" {
		return 0, 0, fmt.Errorf("samtools view: empty subsample value")
	}
	// First try a plain float.
	if f, ferr := strconv.ParseFloat(s, 64); ferr == nil && (f >= 0 && f <= 1) {
		return f, 0, nil
	}
	// Then try "<seed>.<fraction>".
	dot := strings.IndexByte(s, '.')
	if dot < 0 {
		return 0, 0, fmt.Errorf("samtools view: bad subsample %q", s)
	}
	seedPart := s[:dot]
	fracPart := "0." + s[dot+1:]
	if seedPart != "" {
		sv, perr := strconv.ParseInt(seedPart, 10, 64)
		if perr != nil {
			return 0, 0, fmt.Errorf("samtools view: bad subsample seed in %q: %w", s, perr)
		}
		seed = sv
	}
	f, perr := strconv.ParseFloat(fracPart, 64)
	if perr != nil {
		return 0, 0, fmt.Errorf("samtools view: bad subsample fraction in %q: %w", s, perr)
	}
	if f < 0 || f > 1 {
		return 0, 0, fmt.Errorf("samtools view: subsample fraction out of range: %f", f)
	}
	return f, seed, nil
}
