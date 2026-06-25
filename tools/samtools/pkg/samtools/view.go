package samtools

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// openSeekable opens path for indexed random access. A remote URL (http(s)://,
// s3://, gs://) is opened through hfile with read-ahead buffering so the many
// small sequential reads the BGZF/BAM decoder performs coalesce into a few
// large ranged GETs; a local path is opened from disk. The caller must Close
// the result.
func openSeekable(path string) (hfile.SeekHandle, error) {
	return hfile.OpenSeekable(path)
}

// ViewOptions configures the behaviour of View. Zero values disable the
// corresponding filter.
type ViewOptions struct {
	// OutputBAM forces BAM output. When false, output is text SAM unless the
	// caller pre-wraps the output writer in a BAM writer.
	OutputBAM bool
	// Uncompressed selects uncompressed (level-0) BAM output — samtools view's
	// `-u` flag. It implies BAM output (upstream's `-u` sets the output format
	// to BAM with compression level 0). The BAM is still BGZF-framed, just with
	// stored DEFLATE blocks. It has no effect on SAM or CRAM output.
	Uncompressed bool
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
	// MultiRegion selects upstream's `-M/--use-multi-region-iterator`. It
	// changes the *output*, not just performance: by default `samtools view
	// reg1 reg2` walks each region in command-line order and emits its
	// overlapping records, so a record overlapping two regions is emitted once
	// per region (and the regions need not be coordinate-ordered). With -M the
	// regions are walked as one deduplicated, coordinate-ordered set — each
	// record at most once. See buildRegionScanPasses.
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
	// CRAMQualityBinning selects a lossy quality-score binning scheme for
	// CRAM output, set from `--output-fmt-option qbin=...`. The empty
	// string (and "none") disables binning. It is a no-op for SAM and BAM
	// output, where qualities are always stored verbatim. Recognised
	// values: "8" (Illumina 8-level), "4" (4-level), "2" (NovaSeq-style
	// 2-level). See alnio.ParseQualityBinning for the full grammar.
	CRAMQualityBinning string
	// IndexPath is an explicit index-file path supplied via samtools view's
	// `-X/--customized-index` flag. When non-empty it overrides the
	// conventional sibling `<input>.csi`/`<input>.bai` lookup; the index
	// kind (CSI or BAI) is auto-detected from the file's magic bytes.
	IndexPath string
	// Threads is the BGZF compression worker count from `-@/--threads`. A
	// value above 1 enables the parallel BGZF back end for BAM output; the
	// decoded records are byte-identical to single-threaded output. It has no
	// effect on SAM output (uncompressed text) and is a no-op for CRAM, whose
	// container does not yet parallelise. 0 or 1 means single-threaded.
	Threads int
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
	r, err := alnio.NewReaderThreaded(in, opts.Reference, opts.Threads)
	if err != nil {
		return 0, err
	}
	if rc, ok := r.(io.Closer); ok {
		defer rc.Close()
	}
	hdr := r.Header()

	// Resolve any region specifiers up front so unknown chroms don't abort
	// the linear scan — they just yield zero matches.
	resolved, _, err := region.ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
	if err != nil {
		return 0, err
	}
	// Fast path: a BGZF-wrapped BAM stream serialised to plain SAM with only
	// fixed-prefix-decidable filters (flag/MAPQ/region/BED) is written straight
	// from the raw BAM bytes, skipping the per-record Record build. This is the
	// same direct serialiser the indexed path uses; it applies here to whole-
	// file and stdin BAM->SAM. Other formats (text SAM, CRAM) and filters that
	// need the read name or an aux tag fall through to the decode loop below.
	if br, ok := r.(*sam.BAMReader); ok && samFastPathEligible(&opts) {
		return viewStreamFast(br, out, &opts, hdr, resolved)
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

	sub := newSubsampler(opts)
	matched := 0
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return matched, err
		}
		if !keepRecord(rec, &opts, sub) {
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
// at inPath, looks for a sibling index, and (when found and the caller has
// supplied regions) uses the index chunks to seek to relevant portions of
// the file. A coordinate-sorted index (<inPath>.csi) is preferred when
// present — it addresses references beyond the BAI 2^29 bp ceiling — and a
// <inPath>.bai is used otherwise. When no index is found ViewFile falls
// back to the streaming View() linear-scan path; a warning is written to
// warnW (which may be nil to silence it).
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
		f, err := openSeekable(inPath)
		if err != nil {
			return 0, err
		}
		defer f.Close()
		return View(f, out, opts)
	}
	// CRAM input uses a .crai index and a container-seek model rather than
	// the BGZF virtual-offset model the .bai/.csi paths below assume. Detect
	// a CRAM file from its 4-byte magic (sniffed from the seekable handle, so
	// a remote object is read with a tiny ranged GET, not downloaded whole)
	// and route it to the dedicated .crai-indexed query path.
	isCRAM, sniffErr := pathIsCRAM(inPath)
	if sniffErr != nil {
		return 0, sniffErr
	}
	if isCRAM {
		return viewCRAMIndexed(inPath, out, opts, warnW)
	}
	// An explicit -X/--customized-index path overrides the sibling-file
	// lookup. The index kind is auto-detected from its 4-byte magic, so a
	// caller may point at either a .csi or a .bai regardless of extension.
	if opts.IndexPath != "" {
		idxBytes, ierr := hfile.ReadFile(opts.IndexPath)
		if ierr != nil {
			return 0, fmt.Errorf("samtools view: read index %s: %w", opts.IndexPath, ierr)
		}
		f, err := openSeekable(inPath)
		if err != nil {
			return 0, err
		}
		defer f.Close()
		// A .bai is an uncompressed file beginning with the magic "BAI\1";
		// a .csi is BGZF-compressed (its "CSI\1" magic is inside the
		// stream), so its first two bytes are the gzip magic 0x1f 0x8b.
		if len(idxBytes) >= 2 && idxBytes[0] == 0x1f && idxBytes[1] == 0x8b {
			idx, cerr := bam.ReadCSI(bytes.NewReader(idxBytes))
			if cerr != nil {
				return 0, fmt.Errorf("samtools view: read %s: %w", opts.IndexPath, cerr)
			}
			return viewIndexedCSI(f, idx, out, opts)
		}
		idx, berr := bam.ReadBAI(bytes.NewReader(idxBytes))
		if berr != nil {
			return 0, fmt.Errorf("samtools view: read %s: %w", opts.IndexPath, berr)
		}
		return viewIndexed(f, idx, out, opts)
	}
	// Prefer a coordinate-sorted index (.csi) over .bai — it covers the
	// larger coordinate range CSI supports.
	csiPath := inPath + ".csi"
	if csiBytes, csiErr := hfile.ReadFile(csiPath); csiErr == nil {
		idx, ierr := bam.ReadCSI(bytes.NewReader(csiBytes))
		if ierr != nil {
			return 0, fmt.Errorf("samtools view: read %s: %w", csiPath, ierr)
		}
		f, err := openSeekable(inPath)
		if err != nil {
			return 0, err
		}
		defer f.Close()
		return viewIndexedCSI(f, idx, out, opts)
	}
	baiPath := inPath + ".bai"
	baiBytes, baiErr := hfile.ReadFile(baiPath)
	if baiErr != nil {
		if warnW != nil {
			fmt.Fprintf(warnW, "samtools view: no index at %s, falling back to linear scan\n", baiPath)
		}
		f, err := openSeekable(inPath)
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

	f, err := openSeekable(inPath)
	if err != nil {
		return 0, err
	}
	defer f.Close()
	return viewIndexed(f, idx, out, opts)
}

// pathIsCRAM reports whether the file at path begins with the four-byte
// CRAM magic. It reads only the first four bytes through openSeekable, so a
// remote object is probed with a tiny ranged GET rather than downloaded
// whole. A path that cannot be opened (e.g. a missing file) is not a CRAM —
// the error is left for the downstream open to report with full context.
func pathIsCRAM(path string) (bool, error) {
	f, err := openSeekable(path)
	if err != nil {
		// Defer the error: the BAM path below re-opens path and will report
		// the failure with its own context. Treat an unopenable path as
		// not-CRAM so the existing diagnostics are preserved.
		return false, nil
	}
	defer f.Close()
	var magic [4]byte
	if _, err := io.ReadFull(f, magic[:]); err != nil {
		return false, nil
	}
	return string(magic[:]) == "CRAM", nil
}

// viewCRAMIndexed answers a region query against a CRAM file using its
// sibling .crai index (or the explicit -X/--customized-index path), seeking
// to the containers the regions overlap instead of streaming the whole file.
// When no .crai is found it falls back to the streaming View path, writing a
// warning to warnW (which may be nil).
//
// It honours opts.Reference and the REF_CACHE environment variable exactly
// as the streaming path does, so reference-backed CRAM reconstructs its
// bases. Emission reuses the shared keepRecord / subsample / region / BED
// pipeline so a CRAM region query filters identically to BAM.
func viewCRAMIndexed(inPath string, out io.Writer, opts ViewOptions, warnW io.Writer) (int, error) {
	craiPath := inPath + ".crai"
	if opts.IndexPath != "" {
		craiPath = opts.IndexPath
	}
	craiBytes, craiErr := hfile.ReadFile(craiPath)
	if craiErr != nil {
		if warnW != nil {
			fmt.Fprintf(warnW, "samtools view: no CRAM index at %s, falling back to linear scan\n", craiPath)
		}
		f, err := openSeekable(inPath)
		if err != nil {
			return 0, err
		}
		defer f.Close()
		return View(f, out, opts)
	}
	idx, err := cram.ReadCRAI(bytes.NewReader(craiBytes))
	if err != nil {
		return 0, fmt.Errorf("samtools view: read %s: %w", craiPath, err)
	}

	f, err := openSeekable(inPath)
	if err != nil {
		return 0, err
	}
	defer f.Close()

	rr, err := cram.NewRegionReader(f, idx)
	if err != nil {
		return 0, err
	}
	defer rr.Close()
	if opts.Reference != "" {
		if err := rr.SetReferenceFASTA(opts.Reference); err != nil {
			return 0, err
		}
	}
	rr.UseRefCacheFromEnv()
	rr.UseRefPathFromEnv()

	hdr := rr.Header()
	resolved, _, perr := region.ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
	if perr != nil {
		return 0, perr
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

	sub := newSubsampler(opts)
	matched := 0
	// QueryRegion already restricts each query's records to its reference and
	// coordinate range; iterate the resolved regions in command-line order and
	// emit each region's overlapping records, applying the same per-record
	// filters the BAM indexed path uses. Records are de-duplicated within a
	// region by the RegionReader's container-overlap test; across regions a
	// record appears once per overlapping region, matching upstream's default
	// multi-region behaviour. (NB: -M's cross-region de-duplication is honoured
	// for BAM/CSI but not yet for CRAM — a rare combination; see the parity
	// roadmap.)
	for _, reg := range resolved {
		recs, qerr := rr.Query(reg)
		if qerr != nil {
			return matched, qerr
		}
		for _, rec := range recs {
			if !keepRecord(rec, &opts, sub) {
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

// viewIndexed performs an indexed region scan against a .bai index: it
// parses regions, computes chunk unions, seeks into each chunk's
// compressed offset, decodes records until each chunk's end virtual
// offset, and emits records overlapping any requested region.
//
// This path is BAM-only: it does BGZF virtual-offset seeks against a .bai
// index. CRAM uses a .crai index and a different seek model, handled by
// viewCRAMIndexed; ViewFile sniffs the input's magic and routes a CRAM file
// there before reaching this function.
func viewIndexed(f io.ReadSeeker, idx *bam.BAIIndex, out io.Writer, opts ViewOptions) (int, error) {
	return viewIndexedChunks(f, out, opts, func(hdr *sam.Header, resolved []region.ResolvedRegion) []bam.BAIChunk {
		return bam.UnionChunks(idx, resolved)
	})
}

// viewIndexedCSI performs an indexed region scan against a BAM .csi index.
// CSI shares the BAI chunk model, so the seek-and-scan loop is identical;
// only the chunk-lookup index differs. CSI is preferred when present
// because it addresses references beyond the BAI 2^29 bp ceiling.
func viewIndexedCSI(f io.ReadSeeker, idx *bam.CSIIndex, out io.Writer, opts ViewOptions) (int, error) {
	return viewIndexedChunks(f, out, opts, func(hdr *sam.Header, resolved []region.ResolvedRegion) []bam.BAIChunk {
		return bam.UnionChunksCSI(idx, resolved)
	})
}

// viewIndexedChunks is the index-kind-agnostic seek-and-scan core shared by
// the .bai and .csi indexed-query paths. unionFn resolves the requested
// regions to a sorted, merged BAIChunk slice for the relevant index.
func viewIndexedChunks(f io.ReadSeeker, out io.Writer, opts ViewOptions, unionFn func(*sam.Header, []region.ResolvedRegion) []bam.BAIChunk) (int, error) {
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
	// Build the indexed-scan passes. Default `samtools view reg1 reg2 ...`
	// processes each region in command-line order and emits its overlapping
	// records, so a record overlapping two regions is emitted once per region;
	// -M / --use-multi-region-iterator instead walks all regions as one
	// deduplicated, coordinate-ordered set. A single region is one pass either
	// way, so the common case is unchanged.
	passes := buildRegionScanPasses(hdr, resolved, opts.MultiRegion, unionFn)

	// Fast path: plain-SAM output with only fixed-prefix-decidable filters
	// (flag/MAPQ/region/BED) serialises each survivor straight from the raw BAM
	// bytes, skipping the per-record Record build and its string allocations.
	// This is the dominant samtools-view workload (a region query to SAM) and
	// the biggest win; output is byte-identical to the decode path below.
	if samFastPathEligible(&opts) {
		return viewIndexedChunksFast(f, out, &opts, hdr, passes)
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

	sub := newSubsampler(opts)
	matched := 0
	for i := range passes {
		regionFilter := buildRegionFilter(passes[i].regions, hdr)
		if regionFilter == nil && len(opts.Regions) > 0 {
			regionFilter = func(*sam.Record) bool { return false }
		}
		n, serr := scanChunksDecode(f, passes[i].chunks, hdr, &opts, regionFilter, bedFilter, sub, w)
		matched += n
		if serr != nil {
			return matched, serr
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

// regionScanPass is one indexed-scan pass: the BAI/CSI chunks to walk and the
// regions a record must overlap to be emitted. See buildRegionScanPasses.
type regionScanPass struct {
	regions []region.ResolvedRegion
	chunks  []bam.BAIChunk
}

// buildRegionScanPasses splits the resolved regions into indexed-scan passes.
// In the default (non -M) mode each region is its own pass, walked in
// command-line order, reproducing upstream `samtools view reg1 reg2`'s
// once-per-region emission (a record overlapping two regions is emitted twice).
// With -M (or a single region) all regions collapse into one deduplicated,
// coordinate-ordered pass over the union of their chunks.
func buildRegionScanPasses(hdr *sam.Header, resolved []region.ResolvedRegion, multiRegion bool, unionFn func(*sam.Header, []region.ResolvedRegion) []bam.BAIChunk) []regionScanPass {
	if multiRegion || len(resolved) <= 1 {
		return []regionScanPass{{regions: resolved, chunks: unionFn(hdr, resolved)}}
	}
	passes := make([]regionScanPass, 0, len(resolved))
	for i := range resolved {
		one := resolved[i : i+1]
		passes = append(passes, regionScanPass{regions: one, chunks: unionFn(hdr, one)})
	}
	return passes
}

// scanChunksDecode walks chunks, decoding each record and emitting those that
// pass the flag/MAPQ filters, regionFilter, and bedFilter into w (or only
// counting when opts.Count). It is the per-pass body of the decode-path indexed
// scan: the caller opens and closes w (and the header) so several passes — one
// per region in default multi-region mode — share a single output stream.
func scanChunksDecode(f io.ReadSeeker, chunks []bam.BAIChunk, hdr *sam.Header, opts *ViewOptions, regionFilter, bedFilter func(*sam.Record) bool, sub *subsampler, w sam.Writer) (int, error) {
	matched := 0
	for _, c := range chunks {
		if c.Beg >= c.End {
			continue
		}
		startBlock := int64(c.Beg >> 16)
		if _, err := f.Seek(startBlock, io.SeekStart); err != nil {
			return matched, err
		}
		// NewReaderAt(startBlock) so the reader's VirtualOffset is absolute and
		// the chunk-bounded wrapper stops exactly at c.End (see bgzf.NewReaderAt).
		bgz, err := bgzip.NewReaderAt(f, startBlock)
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
		// Chunk-bounded reader so the BAM parser stops when the chunk ends.
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
			if !keepRecord(rec, opts, sub) {
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
	return matched, nil
}

// viewIndexedChunksFast is the BAM->SAM text fast path for the indexed
// seek-and-scan loop. It mirrors viewIndexedChunks' chunk iteration but, for
// each record that survives the flag/MAPQ/region/BED filters, serialises the
// SAM line straight from the raw BAM bytes via sam.BAMReader.WriteSAMBody —
// avoiding the intermediate Record and its per-field/per-aux string
// allocations. Callers must have verified samFastPathEligible(opts).
func viewIndexedChunksFast(f io.ReadSeeker, out io.Writer, opts *ViewOptions, hdr *sam.Header, passes []regionScanPass) (int, error) {
	bedPred, berr := fastBedFilter(opts.BedPath)
	if berr != nil {
		return 0, berr
	}

	// When counting (-c) we still iterate and filter via the cheap fixed-prefix
	// decode, but emit no header and write no record bytes; bw stays nil so
	// fastSAMScan only counts. Otherwise serialise to a buffered writer.
	var bw *bufio.Writer
	if !opts.Count {
		bw = bufio.NewWriter(out)
		if opts.HeaderOnly || opts.WithHeader {
			if _, err := hdr.WriteTo(bw); err != nil {
				return 0, err
			}
		}
		if opts.HeaderOnly {
			return 0, bw.Flush()
		}
	}

	matched := 0
	for pi := range passes {
		regionPred := fastRegionPredicate(passes[pi].regions)
		if regionPred == nil && len(opts.Regions) > 0 {
			regionPred = func(*sam.FastFields) bool { return false }
		}
		n, serr := scanChunksFast(f, passes[pi].chunks, hdr, opts, bw, regionPred, bedPred)
		matched += n
		if serr != nil {
			return matched, serr
		}
	}
	if opts.Count {
		fmt.Fprintln(out, matched)
		return matched, nil
	}
	if err := bw.Flush(); err != nil {
		return matched, err
	}
	return matched, nil
}

// scanChunksFast is the per-pass body of the fast indexed scan: it walks chunks,
// serialising each survivor straight from raw BAM bytes into bw (a nil bw counts
// only). The caller writes the header and flushes bw once across all passes, so
// default multi-region mode can run one pass per region into a single stream.
func scanChunksFast(f io.ReadSeeker, chunks []bam.BAIChunk, hdr *sam.Header, opts *ViewOptions, bw *bufio.Writer, regionPred, bedPred func(*sam.FastFields) bool) (int, error) {
	matched := 0
	for _, c := range chunks {
		if c.Beg >= c.End {
			continue
		}
		startBlock := int64(c.Beg >> 16)
		if _, err := f.Seek(startBlock, io.SeekStart); err != nil {
			return matched, err
		}
		// NewReaderAt(startBlock) so VirtualOffset is absolute and the
		// chunk-bounded wrapper stops exactly at c.End (see bgzf.NewReaderAt).
		bgz, err := bgzip.NewReaderAt(f, startBlock)
		if err != nil {
			return matched, err
		}
		uoff := int(c.Beg & 0xFFFF)
		if uoff > 0 {
			if _, err := io.CopyN(io.Discard, bgz, int64(uoff)); err != nil {
				return matched, err
			}
		}
		boundedSrc := &chunkBoundedReader{r: bgz, end: uint64(c.End)}
		br := sam.NewBAMBodyReader(boundedSrc, hdr)
		n, scanErr := fastSAMScan(br, bw, opts, regionPred, bedPred)
		matched += n
		_ = bgz.Close()
		if scanErr != nil {
			return matched, scanErr
		}
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
		// CRAMQualityBinning selects optional lossy quality-score binning;
		// the empty string maps to BinningNone, so the default CRAM output
		// stays losslessly exact.
		binning, err := alnio.ParseQualityBinning(opts.CRAMQualityBinning)
		if err != nil {
			return nil, err
		}
		// With -T/--reference, load the reference so the writer can encode
		// reads reference-based (only mismatches stored) — far smaller and
		// faster, matching upstream CRAM. Without it the writer falls back to
		// the self-contained reference-free encoding.
		var ref map[string][]byte
		if opts.Reference != "" {
			ref, err = loadReferenceMap(opts.Reference, hdr)
			if err != nil {
				return nil, err
			}
		}
		// opts.Reference is the -T path: pass it as the UR: tag source and the
		// loaded bases (ref) as the M5 / reference-based-encoding source, so the
		// CRAM @SQ lines carry M5+UR exactly as upstream samtools writes them.
		w = alnio.NewCRAMWriterOpts(out, alnio.CRAMWriteOptions{
			QualityBinning: binning,
			Reference:      ref,
			ReferencePath:  opts.Reference,
		})
	case opts.OutputBAM:
		bw, err := sam.NewBAMWriterOptions(out, sam.BAMWriterOptions{
			Uncompressed: opts.Uncompressed,
			Threads:      opts.Threads,
		})
		if err != nil {
			return nil, err
		}
		w = bw
	default:
		w = sam.NewSAMWriter(out)
	}
	// BAM and CRAM are binary container formats whose header is structural
	// — it is always emitted regardless of -h/-H. CRAM in particular
	// cannot be written without a header.
	emitHeader := opts.HeaderOnly || opts.WithHeader || opts.OutputBAM || opts.OutputCRAM
	hdrArg := hdr
	if !emitHeader {
		hdrArg = nil
	}
	if err := w.WriteHeader(hdrArg); err != nil {
		// Close the writer so a parallel BGZF back end's worker goroutines
		// drain rather than blocking forever on the unconsumed job channel.
		_ = w.Close()
		return nil, err
	}
	return w, nil
}

// loadReferenceMap loads the bases of every header @SQ contig present in the
// FASTA at refPath into a name->bases map for reference-based CRAM encoding.
// Contigs absent from the FASTA are skipped (reads on them fall back to the
// reference-free per-read encoding). The whole referenced sequence is held in
// memory, matching how a CRAM encode needs random access to it; a lazy
// faidx-backed provider would lower peak memory for very large genomes.
func loadReferenceMap(refPath string, hdr *sam.Header) (map[string][]byte, error) {
	ra, err := fasta.OpenRandomAccess(refPath)
	if err != nil {
		return nil, fmt.Errorf("samtools view: open reference %s: %w", refPath, err)
	}
	defer ra.Close()
	ref := make(map[string][]byte, len(hdr.Refs))
	for _, sq := range hdr.Refs {
		n := ra.Length(sq.Name)
		if n <= 0 {
			continue // contig not in the FASTA index
		}
		seq, ferr := ra.Fetch(sq.Name, 0, n)
		if ferr != nil {
			return nil, fmt.Errorf("samtools view: fetch reference %s: %w", sq.Name, ferr)
		}
		ref[sq.Name] = seq
	}
	return ref, nil
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
func keepRecord(rec *sam.Record, opts *ViewOptions, sub *subsampler) bool {
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
	if sub != nil {
		if !sub.keep(rec.QName) {
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

// subsampler reproduces upstream samtools view's deterministic per-read-name
// subsample decision (sam_view.c process_aln): a read is kept when
//
//	k = Wang(X31(qname) ^ seed); (k & 0xffffff) / 0x1000000 < frac
//
// where Wang is khash.h's __ac_Wang_hash, X31 is __ac_X31_hash_string, and
// seed is the user seed run through glibc's srand()/rand() transform when
// non-zero (sam_view.c:1307). Because the decision is a pure function of the
// read name, every alignment of a template/read-pair makes the same keep/drop
// choice — so mates stay together — and the kept set is independent of input
// order, exactly matching upstream (and replacing the old, non-portable
// *rand.Rand approach that drew a fresh random number per record).
type subsampler struct {
	frac float64
	seed uint32
}

// newSubsampler builds the subsample filter for the requested fraction/seed,
// or returns nil when no subsampling is requested (frac <= 0 or >= 1). The
// seed is the SEED component of `-s SEED.FRAC`; a non-zero seed is passed
// through glibcSrandRand to match upstream's entropy-spreading step.
func newSubsampler(opts ViewOptions) *subsampler {
	if opts.Subsample <= 0 || opts.Subsample >= 1 {
		return nil
	}
	seed := uint32(opts.SubsampleSeed)
	if seed != 0 {
		seed = glibcSrandRand(seed)
	}
	return &subsampler{frac: opts.Subsample, seed: seed}
}

// keep reports whether the named read survives subsampling.
func (s *subsampler) keep(qname string) bool {
	k := acWangHash(acX31HashString(qname) ^ s.seed)
	return float64(k&0xffffff)/0x1000000 < s.frac
}

// acX31HashString is a byte-for-byte port of htslib khash.h's
// __ac_X31_hash_string (h = h*31 + c, seeded with the first byte).
func acX31HashString(s string) uint32 {
	if len(s) == 0 {
		return 0
	}
	h := uint32(s[0])
	for i := 1; i < len(s); i++ {
		h = (h << 5) - h + uint32(s[i])
	}
	return h
}

// acWangHash is a byte-for-byte port of htslib khash.h's __ac_Wang_hash
// integer mixing function.
func acWangHash(key uint32) uint32 {
	key += ^(key << 15)
	key ^= key >> 10
	key += key << 3
	key ^= key >> 6
	key += ^(key << 11)
	key ^= key >> 16
	return key
}

// glibcSrandRand reproduces the result of `srand(seed); return rand();` under
// glibc's default TYPE_3 additive-feedback generator. samtools view runs the
// user seed through this transform (sam_view.c:1307-1311) to spread the
// entropy of small integer seeds before XOR-ing it into the name hash, so we
// must reproduce it exactly for byte-identical subsampling.
func glibcSrandRand(seed uint32) uint32 {
	if seed == 0 {
		seed = 1
	}
	var r [344]int32
	r[0] = int32(seed)
	for i := 1; i < 31; i++ {
		// r[i] = (16807 * r[i-1]) % 2147483647, via Schrage's method to
		// avoid 32-bit signed overflow (matches glibc's int arithmetic).
		hi := r[i-1] / 127773
		lo := r[i-1] % 127773
		w := 16807*lo - 2836*hi
		if w < 0 {
			w += 2147483647
		}
		r[i] = w
	}
	for i := 31; i < 34; i++ {
		r[i] = r[i-31]
	}
	for i := 34; i < 344; i++ {
		r[i] = r[i-31] + r[i-3]
	}
	val := r[344-31] + r[344-3]
	return (uint32(val) >> 1) & 0x7fffffff
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
