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

// ErrRegionsUnsupported is returned when a caller requests region-query
// filtering. Region queries depend on a BAI index that is not yet implemented.
var ErrRegionsUnsupported = errors.New("samtools view: region-query support requires .bai indexing — not yet implemented")

// View streams alignment records from r, applying the configured filters,
// and emits them to out as either SAM text or BAM. Returns the number of
// records that passed all filters.
//
// When opts.Count is true, the matching count is the only thing emitted to out.
func View(in io.Reader, out io.Writer, opts ViewOptions) (int, error) {
	if opts.RegionsEnabled {
		return 0, ErrRegionsUnsupported
	}

	r, err := sam.NewReader(in)
	if err != nil {
		return 0, err
	}
	hdr := r.Header()

	var w sam.Writer
	if !opts.Count {
		if opts.OutputBAM {
			bw := sam.NewBAMWriter(out)
			w = bw
		} else {
			w = sam.NewSAMWriter(out)
		}
		emitHeader := opts.HeaderOnly || opts.WithHeader || opts.OutputBAM
		if emitHeader {
			if err := w.WriteHeader(hdr); err != nil {
				return 0, err
			}
		} else {
			// Still call WriteHeader on the underlying writer with nil so the
			// internal state is correctly initialised, but emit nothing.
			if err := w.WriteHeader(nil); err != nil {
				return 0, err
			}
		}
		// HeaderOnly stops here.
		if opts.HeaderOnly {
			return 0, w.Close()
		}
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
	if w != nil {
		if err := w.Close(); err != nil {
			return matched, err
		}
	}
	return matched, nil
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
