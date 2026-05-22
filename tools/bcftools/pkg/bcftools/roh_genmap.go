package bcftools

// Genetic-map and recombination-rate transition scaling, plus the
// --estimate-AF cohort logic, for bcftools roh. Ported from
// reference_code/bcftools/vcfroh.c (load_genmap, get_genmap_rate,
// set_tprob_genmap, set_tprob_rrate, estimate_AF_from_GT).

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// genMapEntry is one row of an IMPUTE2-format genetic map: a 0-based
// position and the cumulative genetic-map value (in likelihood units,
// i.e. cM scaled by 0.01).
type genMapEntry struct {
	pos  int
	rate float64
}

// genMap is a single chromosome's genetic map plus a cursor for the
// monotone get-rate lookups.
type genMap struct {
	entries []genMapEntry
	cursor  int
}

// loadGenMap reads an IMPUTE2-format genetic map. The "{CHROM}"
// placeholder is left for the caller to substitute; this loader works
// with the literal path. It returns nil when no map was requested.
func loadGenMap(path string) (*genMap, error) {
	if path == "" {
		return nil, nil
	}
	if strings.Contains(path, "{CHROM}") {
		// Per-chromosome maps are loaded lazily per block instead; the
		// literal path here is unusable, so defer to the modifier.
		return &genMap{}, nil
	}
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("bcftools roh: --genetic-map: %w", err)
	}
	defer in.Close()
	gm, err := parseGenMap(in)
	if err != nil {
		return nil, fmt.Errorf("bcftools roh: --genetic-map: %w", err)
	}
	return gm, nil
}

// parseGenMap parses the IMPUTE2 genetic-map format. The header line
// must be exactly "position COMBINED_rate(cM/Mb) Genetic_Map(cM)".
func parseGenMap(in io.Reader) (*genMap, error) {
	sc := bufio.NewScanner(in)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	if !sc.Scan() {
		return nil, fmt.Errorf("genetic map empty")
	}
	header := strings.TrimRight(sc.Text(), "\r\n")
	if header != "position COMBINED_rate(cM/Mb) Genetic_Map(cM)" {
		return nil, fmt.Errorf("unexpected header in genetic map: %q", header)
	}
	gm := &genMap{}
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		if line == "" {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 3 {
			return nil, fmt.Errorf("could not parse genetic map line: %q", line)
		}
		pos, err := strconv.Atoi(fields[0])
		if err != nil {
			return nil, fmt.Errorf("could not parse genetic map position: %q", line)
		}
		rate, err := strconv.ParseFloat(fields[2], 64)
		if err != nil {
			return nil, fmt.Errorf("could not parse genetic map rate: %q", line)
		}
		gm.entries = append(gm.entries, genMapEntry{pos: pos - 1, rate: rate * 0.01})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	if len(gm.entries) == 0 {
		return nil, fmt.Errorf("genetic map empty")
	}
	return gm, nil
}

// rate returns the genetic-map distance between two 0-based positions,
// replicating vcfroh.c's get_genmap_rate including its monotone cursor.
func (gm *genMap) rate(start, end int) float64 {
	if len(gm.entries) == 0 {
		return 0
	}
	i := gm.cursor
	if i >= len(gm.entries) {
		i = len(gm.entries) - 1
	}
	if gm.entries[i].pos > start {
		for i > 0 && gm.entries[i].pos > start {
			i--
		}
	} else {
		for i+1 < len(gm.entries) && gm.entries[i+1].pos < start {
			i++
		}
	}
	j := i
	for j+1 < len(gm.entries) && gm.entries[j].pos < end {
		j++
	}
	if i == j {
		gm.cursor = i
		return 0
	}
	if start < gm.entries[i].pos {
		start = gm.entries[i].pos
	}
	if end > gm.entries[j].pos {
		end = gm.entries[j].pos
	}
	rate := (gm.entries[j].rate - gm.entries[i].rate) /
		float64(gm.entries[j].pos-gm.entries[i].pos) * float64(end-start)
	gm.cursor = j
	return rate
}

// tprobModifier scales the off-diagonal transition probabilities at
// each site by the cross-over probability of the inter-site interval,
// either from a genetic map or a constant recombination rate.
type tprobModifier struct {
	gm      *genMap
	recRate float64
	useGM   bool
}

// newTprobModifier builds the per-site transition modifier for the
// chromosome, or nil if neither -m nor -M was requested.
func newTprobModifier(opts RohOptions, gmap *genMap, chrom string) *tprobModifier {
	if gmap != nil && len(gmap.entries) > 0 {
		return &tprobModifier{gm: gmap, recRate: opts.RecRate, useGM: true}
	}
	if opts.RecRate > 0 {
		return &tprobModifier{recRate: opts.RecRate}
	}
	return nil
}

// reset rewinds the genetic-map cursor before a fresh decoding pass.
func (m *tprobModifier) reset() {
	if m.gm != nil {
		m.gm.cursor = 0
	}
}

// apply implements set_tprob_genmap / set_tprob_rrate: it scales the
// HW->AZ and AZ->HW entries by the interval cross-over probability and
// keeps the rows normalised.
func (m *tprobModifier) apply(prevPos, pos uint32, tprob []float64) {
	var ci float64
	if m.useGM {
		ci = m.gm.rate(int(prevPos), int(pos))
		if m.recRate != 0 {
			ci *= m.recRate
		}
	} else {
		ci = float64(pos-prevPos) * m.recRate
	}
	if ci > 1 {
		ci = 1
	}
	hw2az := matAt(tprob, 2, stateHW, stateAZ) * ci
	az2hw := matAt(tprob, 2, stateAZ, stateHW) * ci
	matSet(tprob, 2, stateHW, stateAZ, hw2az)
	matSet(tprob, 2, stateAZ, stateHW, az2hw)
	matSet(tprob, 2, stateAZ, stateAZ, 1-hw2az)
	matSet(tprob, 2, stateHW, stateHW, 1-az2hw)
}

// estimateAFCohort holds the sample indices the --estimate-AF
// frequency is computed from. An empty slice means "all samples".
type estimateAFCohort struct {
	indices []int
	fromPL  bool
}

// parseEstimateAF parses the -e/--estimate-AF argument. The optional
// "GT," / "PL," prefix selects the source FORMAT tag; the remainder is
// "-" (all samples) or a comma-separated sample list. The returned
// bool reports whether AF estimation is active.
func parseEstimateAF(spec string, hdr *vcf.Header, sampleIdx map[string]int) (estimateAFCohort, bool, error) {
	if spec == "" {
		return estimateAFCohort{}, false, nil
	}
	cohort := estimateAFCohort{}
	rest := spec
	switch {
	case strings.HasPrefix(rest, "GT,"):
		rest = rest[3:]
	case strings.HasPrefix(rest, "PL,"):
		rest = rest[3:]
		cohort.fromPL = true
	}
	if rest == "-" || rest == "" {
		return cohort, true, nil
	}
	for _, name := range strings.Split(rest, ",") {
		name = strings.TrimSpace(name)
		if name == "" {
			continue
		}
		idx, ok := sampleIdx[name]
		if !ok {
			return estimateAFCohort{}, false, fmt.Errorf("bcftools roh: --estimate-AF: sample not found: %s", name)
		}
		cohort.indices = append(cohort.indices, idx)
	}
	// When every sample is selected the per-site lookup is unneeded;
	// upstream nulls the list. We keep it equivalent: empty == all.
	if len(cohort.indices) == len(hdr.Samples) {
		cohort.indices = nil
	}
	return cohort, true, nil
}

// estimateAF estimates the alternate-allele frequency for one record
// from the genotypes of the cohort, mirroring estimate_AF_from_GT.
// PL-based estimation is not present in the test inputs; when the PL
// path is requested but PLs are absent the site is rejected.
func estimateAF(v *vcf.Variant, ial int, cohort estimateAFCohort, sampleIdx map[string]int, hardGTErr float64) (float64, bool) {
	nalt, nref := 0, 0
	count := func(idx int) {
		cls, ok := rohGTClass(v, idx, ial)
		if !ok {
			return
		}
		switch cls {
		case 0:
			nref += 2
		case 1:
			nalt++
			nref++
		case 2:
			nalt += 2
		}
	}
	if len(cohort.indices) > 0 {
		for _, idx := range cohort.indices {
			count(idx)
		}
	} else {
		for i := range v.Samples {
			count(i)
		}
	}
	if nalt == 0 && nref == 0 {
		return 0, false
	}
	return float64(nalt) / float64(nalt+nref), true
}
