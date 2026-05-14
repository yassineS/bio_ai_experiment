package samtools

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
)

// ViewOptions configures the behaviour of View. Zero values disable the
// corresponding filter.
type ViewOptions struct {
	// OutputBAM forces BAM output. When false, output is text SAM unless the
	// caller pre-wraps the output writer in a BAM writer.
	OutputBAM bool
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
	// NoPG suppresses the @PG line emission (placeholder; the view pipeline
	// does not currently inject an @PG line, so this is a no-op kept for
	// flag compatibility).
	NoPG bool
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
	r, err := sam.NewReader(in)
	if err != nil {
		return 0, err
	}
	hdr := r.Header()

	// Resolve any region specifiers up front so unknown chroms don't abort
	// the linear scan — they just yield zero matches.
	resolved, _, err := ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
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
	// No regions requested: take the streaming path.
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
	idx, ierr := ReadBAI(strings.NewReader(string(baiBytes)))
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
func viewIndexed(f *os.File, idx *BAIIndex, out io.Writer, opts ViewOptions) (int, error) {
	// Need the header — open a BAM reader first.
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		return 0, err
	}
	hdrReader, err := sam.NewBAMReader(f)
	if err != nil {
		return 0, err
	}
	hdr := hdrReader.Header()

	resolved, _, perr := ResolveRegions(opts.Regions, func(name string) int { return hdr.RefIndex(name) })
	if perr != nil {
		return 0, perr
	}
	chunks := UnionChunks(idx, resolved)
	regionFilter := buildRegionFilter(resolved, hdr)
	if regionFilter == nil && len(opts.Regions) > 0 {
		regionFilter = func(*sam.Record) bool { return false }
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
func buildRegionFilter(regions []ResolvedRegion, hdr *sam.Header) func(*sam.Record) bool {
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
	if opts.OutputBAM {
		w = sam.NewBAMWriter(out)
	} else {
		w = sam.NewSAMWriter(out)
	}
	emitHeader := opts.HeaderOnly || opts.WithHeader || opts.OutputBAM
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
	if rng != nil && opts.Subsample > 0 && opts.Subsample < 1 {
		if rng.Float64() >= opts.Subsample {
			return false
		}
	}
	return true
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
