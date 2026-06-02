package bcftools

// `bcftools filter --mask [^]REGION` and `-M/--mask-file [^]FILE` mark
// records inside (or, with the `^` prefix, outside) a region set with
// the configured `--soft-filter` ID. Behaviour mirrors upstream
// `vcffilter.c` lines 92-110 (region-set construction) and 700-709
// (per-record overlap test) — translated to in-memory regions because
// our streaming filter is index-free.

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// filterMask is the resolved mask: a per-contig sorted list of regions
// plus the negate / overlap-mode flags.
type filterMask struct {
	regions map[string][]maskRegion // CHROM -> sorted []maskRegion
	negate  bool
	overlap int
}

// maskRegion is a single 1-based inclusive interval.
type maskRegion struct {
	beg int
	end int
}

// buildFilterMask compiles MaskRegion / MaskFile / MaskOverlap into a
// filterMask. Returns (nil, nil) when no mask is configured. Returns a
// non-nil error when the user asked for a mask but did not also pass a
// `-s/--soft-filter` name (mirrors upstream's exact text from
// vcffilter.c:656).
func buildFilterMask(opts VCFFilterOptions) (*filterMask, error) {
	if opts.MaskRegion == "" && opts.MaskFile == "" {
		return nil, nil
	}
	if opts.SoftFilter == "" {
		return nil, fmt.Errorf("The option --soft-filter is required with --mask and --mask-file options")
	}
	m := &filterMask{
		regions: map[string][]maskRegion{},
		overlap: opts.MaskOverlap,
	}
	if m.overlap < 0 || m.overlap > 2 {
		m.overlap = 1
	}
	switch {
	case opts.MaskFile != "":
		spec := opts.MaskFile
		if strings.HasPrefix(spec, "^") {
			m.negate = true
			spec = spec[1:]
		}
		if err := loadMaskBED(spec, m); err != nil {
			return nil, err
		}
	case opts.MaskRegion != "":
		spec := opts.MaskRegion
		if strings.HasPrefix(spec, "^") {
			m.negate = true
			spec = spec[1:]
		}
		// Upstream accepts comma-separated regions in --mask too
		// (regidx_init_string splits on commas). Match that.
		for _, part := range strings.Split(spec, ",") {
			part = strings.TrimSpace(part)
			if part == "" {
				continue
			}
			chrom, beg, end, err := parseMaskRegionSpec(part)
			if err != nil {
				return nil, err
			}
			m.regions[chrom] = append(m.regions[chrom], maskRegion{beg: beg, end: end})
		}
	}
	for chrom := range m.regions {
		sort.Slice(m.regions[chrom], func(i, j int) bool {
			return m.regions[chrom][i].beg < m.regions[chrom][j].beg
		})
	}
	return m, nil
}

// parseMaskRegionSpec parses "chr", "chr:pos", or "chr:beg-end" into a
// 1-based inclusive triple.
func parseMaskRegionSpec(s string) (string, int, int, error) {
	colon := strings.IndexByte(s, ':')
	if colon < 0 {
		// Whole-contig mask.
		return s, 1, 1 << 30, nil
	}
	chrom := s[:colon]
	rest := s[colon+1:]
	dash := strings.IndexByte(rest, '-')
	if dash < 0 {
		pos, err := strconv.Atoi(rest)
		if err != nil {
			return "", 0, 0, fmt.Errorf("bad mask region %q", s)
		}
		return chrom, pos, pos, nil
	}
	beg, err := strconv.Atoi(rest[:dash])
	if err != nil {
		return "", 0, 0, fmt.Errorf("bad mask region %q", s)
	}
	endStr := rest[dash+1:]
	end := 1 << 30
	if endStr != "" {
		end, err = strconv.Atoi(endStr)
		if err != nil {
			return "", 0, 0, fmt.Errorf("bad mask region %q", s)
		}
	}
	return chrom, beg, end, nil
}

// loadMaskBED reads a BED-style mask file. BED is 0-based half-open, so
// `chrom\tstart\tend` becomes the 1-based inclusive range
// `[start+1, end]` (matching htslib's regidx_parse_bed). Blank lines
// and `#`-prefixed comments are skipped.
func loadMaskBED(path string, m *filterMask) error {
	f, err := os.Open(path)
	if err != nil {
		return fmt.Errorf("mask-file %s: %w", path, err)
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 3 {
			return fmt.Errorf("mask-file %s: bad line %q (need CHROM\\tBEG\\tEND)", path, line)
		}
		beg, err := strconv.Atoi(fields[1])
		if err != nil {
			return fmt.Errorf("mask-file %s: bad start %q", path, fields[1])
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return fmt.Errorf("mask-file %s: bad end %q", path, fields[2])
		}
		m.regions[fields[0]] = append(m.regions[fields[0]], maskRegion{beg: beg + 1, end: end})
	}
	return sc.Err()
}

// maskMatches reports whether v should be soft-filtered by the mask.
// The overlap mode determines the record span; "negate" inverts the
// final answer (so `^REGION` keeps only records outside the region).
func maskMatches(m *filterMask, v *vcf.Variant) bool {
	beg, end := variantSpan(v, m.overlap)
	regs := m.regions[v.Chrom]
	hit := false
	// Linear scan is fine for v1; the typical mask is dozens of
	// regions, not millions. If real workloads change we can swap in
	// a sorted-binary-search hop using the sort already in place.
	for _, r := range regs {
		if r.beg <= end && r.end >= beg {
			hit = true
			break
		}
	}
	if m.negate {
		return !hit
	}
	return hit
}

// variantSpan returns the 1-based inclusive interval to test against
// the mask for a given overlap mode.
//
//	0 - POS only (single point).
//	1 - the record's REF span: [POS, POS+len(REF)-1] (upstream default).
//	2 - the *variant* span: spans every ALT (longest of REF + each ALT
//	    boundary, to cover deletions / insertions / MNPs alike).
func variantSpan(v *vcf.Variant, overlap int) (int, int) {
	pos := v.Pos
	switch overlap {
	case 0:
		return pos, pos
	case 2:
		// Variant span: union of every allele's affected window.
		// Each ALT may extend the right edge to POS+max(len(REF),
		// len(ALT))-1; we take the max length.
		maxLen := len(v.Ref)
		for _, a := range v.Alt {
			if l := len(a); l > maxLen {
				maxLen = l
			}
		}
		if maxLen <= 0 {
			maxLen = 1
		}
		return pos, pos + maxLen - 1
	}
	// Default: REF span.
	rlen := len(v.Ref)
	if rlen <= 0 {
		rlen = 1
	}
	return pos, pos + rlen - 1
}

// (The header description for the mask soft-filter is set by
// ensureSoftFilterHeaderForMode in vcffilter.go — that path also owns
// the precedence rule "expression description wins over mask
// description" from upstream vcffilter.c:140-148.)
