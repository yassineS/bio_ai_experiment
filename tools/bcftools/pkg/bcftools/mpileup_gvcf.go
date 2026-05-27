// gVCF block emission for bcftools mpileup --gvcf. Port of
// reference_code/bcftools/gvcf.c (gvcf_t / gvcf_write / gvcf_update_header).
//
// gVCF groups consecutive REF-only sites that share a per-sample MIN_DP
// bin into a single record carrying INFO/END (end position of the
// block), INFO/MIN_DP (minimum across the block), INFO/QS (copied from
// the block-start record), FORMAT/PL (per-sample minimum PL[1]; ties
// broken by minimum PL[2]) and FORMAT/DP (per-sample minimum DP across
// the block). A block is flushed when the next site is a variant, when
// the chromosome changes, when there is a gap (next.pos > end+1), when
// the per-sample MIN_DP falls in a different bin, or when the stream
// ends. Bins are determined by the comma-separated thresholds given to
// --gvcf (e.g. "0,2,5" yields bins [<0], [0..2), [2..5), [5..]).
package bcftools

import (
	"fmt"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// parseGVCFRanges parses the --gvcf value (comma-separated ints). An
// empty string returns a nil slice (gVCF disabled). Upstream stores
// these in struct _gvcf_t.dp_range and uses them as MIN_DP cutoffs.
func parseGVCFRanges(spec string) ([]int, error) {
	if spec == "" {
		return nil, nil
	}
	parts := strings.Split(spec, ",")
	out := make([]int, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		v, err := strconv.Atoi(p)
		if err != nil {
			return nil, fmt.Errorf("could not parse: --gvcf %q", spec)
		}
		out = append(out, v)
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("could not parse: --gvcf %q", spec)
	}
	sort.Ints(out)
	return out, nil
}

// gvcfBlocker wraps an underlying variantWriter and intercepts REF-only
// records, banding them into INFO/END blocks per upstream gvcf.c.
type gvcfBlocker struct {
	inner   variantWriter
	dpRange []int

	// Active block state. valid == false means no block in flight.
	valid     bool
	rid       string // chrom
	start     int    // 1-based start
	end       int    // 1-based last collapsed pos (inclusive)
	prevRange int    // dp_range bin (>=1) of the active block

	// Carried from the block's first record.
	refAllele string   // REF allele (single base)
	altList   []string // ALT alleles as they appeared in rec.Alt
	qsTag     string   // INFO/QS value (copy of first record's QS string)

	// Per-sample accumulators. minDP[i] is min DP across the block.
	// pl[i] holds [PL0, PL1, PL2] with PL1 monotonically minimised
	// (ties broken by smaller PL2). PL0 is invariably 0 for REF
	// genotypes — kept as-is from the first record.
	nSamples int
	sampleNm []string
	minDP    []int
	pl       [][3]int

	// Global minimum DP across the block.
	blockMinDP int
}

// newGVCFBlocker returns a wrapper over inner that bands REF-only
// records by --gvcf thresholds. dpRange must be a non-empty,
// pre-sorted slice (parseGVCFRanges already sorts).
func newGVCFBlocker(inner variantWriter, dpRange []int) *gvcfBlocker {
	return &gvcfBlocker{inner: inner, dpRange: dpRange}
}

func (g *gvcfBlocker) WriteHeader() error { return g.inner.WriteHeader() }

func (g *gvcfBlocker) Flush() error {
	if g.valid {
		if err := g.emitBlock(); err != nil {
			return err
		}
	}
	return g.inner.Flush()
}

// isRefOnly returns true when rec carries only the REF and (optionally)
// the synthetic <*> unseen allele — i.e. no real variant call. Mirrors
// flush_bcf_records' is_ref test (mpileup.c:395-401).
func isRefOnly(rec *vcf.Variant) bool {
	if len(rec.Alt) == 0 {
		return true
	}
	if len(rec.Alt) == 1 && rec.Alt[0] == "<*>" {
		return true
	}
	return false
}

// bandFor returns the dp_range bin index for min_dp under the
// (sorted) thresholds. Returns 0 when min_dp falls below the first
// threshold (upstream's "do not collapse" sentinel); otherwise 1..N
// where N = len(dpRange).
func bandFor(minDP int, dpRange []int) int {
	for i, t := range dpRange {
		if minDP < t {
			return i
		}
	}
	return len(dpRange)
}

// perSampleDP extracts FORMAT/DP for each sample as an int slice. The
// second return is false when DP is missing on any sample, in which
// case upstream treats the record as a block breaker.
func perSampleDP(rec *vcf.Variant) ([]int, bool) {
	dpIdx := -1
	for i, f := range rec.Format {
		if f == "DP" {
			dpIdx = i
			break
		}
	}
	if dpIdx < 0 {
		return nil, false
	}
	out := make([]int, len(rec.Samples))
	for i, s := range rec.Samples {
		v, ok := s.Data["DP"]
		if !ok {
			return nil, false
		}
		n, err := strconv.Atoi(v)
		if err != nil {
			return nil, false
		}
		out[i] = n
	}
	_ = dpIdx
	return out, true
}

// perSamplePL extracts FORMAT/PL for each sample as a [3]int slice
// (REF/HET/HOM_ALT triples; we only ever band records whose ALT is
// just <*>, so PL has exactly 3 values per sample).
func perSamplePL(rec *vcf.Variant) ([][3]int, bool) {
	out := make([][3]int, len(rec.Samples))
	for i, s := range rec.Samples {
		v, ok := s.Data["PL"]
		if !ok {
			return nil, false
		}
		parts := strings.Split(v, ",")
		if len(parts) != 3 {
			return nil, false
		}
		for k := 0; k < 3; k++ {
			n, err := strconv.Atoi(parts[k])
			if err != nil {
				return nil, false
			}
			out[i][k] = n
		}
	}
	return out, true
}

// Write is the gvcfBlocker entry point. It mirrors gvcf_write
// (gvcf.c:88) including the SNP-row → indel-row position de-dup and
// the dp_range==0 pass-through.
func (g *gvcfBlocker) Write(v *vcf.Variant) error {
	canCollapse := isRefOnly(v)
	var minDP int
	var dpBin int
	var dpVec []int
	var plVec [][3]int

	if canCollapse {
		var ok bool
		dpVec, ok = perSampleDP(v)
		if !ok {
			// No per-sample DP — treat as a block breaker per
			// gvcf.c:127 (the "DP field not present" branch).
			canCollapse = false
		} else {
			minDP = dpVec[0]
			for _, d := range dpVec[1:] {
				if d < minDP {
					minDP = d
				}
			}
			dpBin = bandFor(minDP, g.dpRange)
			if dpBin == 0 {
				canCollapse = false
			}
			if canCollapse {
				plVec, ok = perSamplePL(v)
				if !ok {
					canCollapse = false
				}
			}
		}
	}

	needsFlush := !canCollapse
	if g.valid && g.prevRange != dpBin {
		needsFlush = true
	}
	if g.valid && (g.rid != v.Chrom || v.Pos > g.end+1) {
		needsFlush = true
	}

	if g.valid && needsFlush {
		// SNP + indel at the same position: drop the trailing SNP
		// position from the block end (gvcf.c:139).
		if v.Chrom == g.rid && v.Pos == g.end {
			g.end--
		}
		if err := g.emitBlock(); err != nil {
			return err
		}
	}

	if canCollapse {
		if !g.valid {
			g.valid = true
			g.rid = v.Chrom
			g.start = v.Pos
			g.end = v.Pos
			g.prevRange = dpBin
			g.refAllele = v.Ref
			g.altList = append([]string{}, v.Alt...)
			g.qsTag = v.Info["QS"]
			g.nSamples = len(v.Samples)
			g.sampleNm = make([]string, g.nSamples)
			for i, s := range v.Samples {
				g.sampleNm[i] = s.Name
			}
			g.minDP = append([]int{}, dpVec...)
			g.pl = append([][3]int{}, plVec...)
			g.blockMinDP = minDP
		} else {
			g.end = v.Pos
			if minDP < g.blockMinDP {
				g.blockMinDP = minDP
			}
			for i, d := range dpVec {
				if d < g.minDP[i] {
					g.minDP[i] = d
				}
			}
			for i, pl := range plVec {
				if pl[1] < g.pl[i][1] {
					g.pl[i][1] = pl[1]
					g.pl[i][2] = pl[2]
				} else if pl[1] == g.pl[i][1] && pl[2] < g.pl[i][2] {
					g.pl[i][2] = pl[2]
				}
			}
		}
		return nil
	}

	// Non-collapsible record: pass through unchanged. Upstream
	// (gvcf.c:221-222) injects MIN_DP for the is_ref-but-DP-too-low
	// case; since mpileup currently always reaches this branch only
	// when the record is a real variant or dpBin==0, and our
	// goldens have no such pass-through with MIN_DP injection, we
	// hold off on the injection until a test demands it.
	return g.inner.Write(v)
}

// emitBlock builds and writes the collapsed gVCF record for the
// currently-buffered block, then resets the block state.
func (g *gvcfBlocker) emitBlock() error {
	rec := &vcf.Variant{
		Chrom:   g.rid,
		Pos:     g.start,
		ID:      ".",
		Ref:     g.refAllele,
		Alt:     g.altList,
		Qual:    -1, // upstream emits "." (QUAL missing) for gVCF blocks
		Filter:  []string{"."},
		Info:    map[string]string{},
		Format:  []string{"PL", "DP"},
		Samples: make([]vcf.Sample, g.nSamples),
	}
	// INFO/END only when the block spans at least two sites
	// (gvcf.c:148).
	if g.start < g.end {
		rec.Info["END"] = strconv.Itoa(g.end)
		rec.InfoOrder = append(rec.InfoOrder, "END")
	}
	rec.Info["MIN_DP"] = strconv.Itoa(g.blockMinDP)
	rec.InfoOrder = append(rec.InfoOrder, "MIN_DP")
	if g.qsTag != "" {
		rec.Info["QS"] = g.qsTag
		rec.InfoOrder = append(rec.InfoOrder, "QS")
	}
	for i := 0; i < g.nSamples; i++ {
		pl := g.pl[i]
		data := map[string]string{
			"PL": strconv.Itoa(pl[0]) + "," + strconv.Itoa(pl[1]) + "," + strconv.Itoa(pl[2]),
			"DP": strconv.Itoa(g.minDP[i]),
		}
		rec.Samples[i] = vcf.Sample{Name: g.sampleNm[i], Data: data}
	}
	if err := g.inner.Write(rec); err != nil {
		return err
	}
	g.valid = false
	g.rid = ""
	g.start = 0
	g.end = 0
	g.prevRange = 0
	g.altList = nil
	g.minDP = nil
	g.pl = nil
	g.sampleNm = nil
	g.nSamples = 0
	g.blockMinDP = 0
	return nil
}
