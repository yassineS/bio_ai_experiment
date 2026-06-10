package bcftools

// --distinctive-sites support for bcftools gtcheck, mirroring the
// diff_sites_* / report_distinctive_sites() machinery in vcfgtcheck.c.
//
// The goal is to find a minimal set of sites that, taken together,
// distinguish at least NUM sample pairs. Upstream:
//
//  1. For every site, records the bitset of pairs that are discordant
//     there (diff_sites_push), tagging the site with its discordant-pair
//     count "ndiff" and a random tie-break key.
//  2. Sorts all recorded sites by ndiff descending (ties broken by the
//     random key) via an external sort.
//  3. Walks the sorted sites, greedily adding a site if it distinguishes
//     at least one not-yet-distinguished pair, emitting a "DS" line with
//     the running cumulative count and a block id. A block closes once
//     the cumulative count reaches the requested NUM, after which the
//     "already distinguished" set resets.
//
// Note on determinism: upstream calls hts_srand48(0) in init_data() and
// breaks ndiff ties with hts_lrand48(). On Linux htslib delegates to
// glibc's srand48/lrand48, so the generator is a 48-bit LCG seeded to
// X = (0<<16)|0x330E. We reproduce that exact generator (lrand48Gen) and
// assign a random key per pushed site in read order, matching upstream
// byte-for-byte including tie resolution.

import (
	"bufio"
	"io"
	"sort"
	"strconv"
)

// diffSite records, for one site, which pairs are discordant there.
type diffSite struct {
	chrom string
	pos   int
	ndiff int
	mask  []uint64 // bitset over pair indices
	rand  uint32   // hts_lrand48() tie-break key, in read order
	order int      // read order, secondary tie-break for determinism
}

// lrand48Gen reproduces glibc's srand48(0)+lrand48(): a 48-bit LCG
// X' = (a*X + c) mod 2^48 returning X' >> 17, seeded to (0<<16)|0x330E.
type lrand48Gen struct {
	state uint64
}

const (
	lrand48A    = 0x5DEECE66D
	lrand48C    = 0xB
	lrand48M    = (uint64(1) << 48) - 1
	lrand48Seed = 0x330E // srand48(0): (seed<<16) | 0x330E
)

// newLrand48 returns a generator matching glibc's srand48(0) state.
func newLrand48() lrand48Gen {
	return lrand48Gen{state: lrand48Seed}
}

// next advances the generator and returns the next 31-bit value, exactly
// as glibc's lrand48() would.
func (g *lrand48Gen) next() uint32 {
	g.state = (lrand48A*g.state + lrand48C) & lrand48M
	return uint32(g.state >> 17)
}

// distinctiveCollector accumulates per-site discordance bitsets when
// --distinctive-sites is active. It is a no-op when disabled.
type distinctiveCollector struct {
	enabled  bool
	npairs   int
	words    int
	nsites   int // resolved NUM threshold (absolute pair count)
	cur      []uint64
	curNdiff int
	sites    []diffSite
	readSeq  int
	rng      lrand48Gen
}

// resolveDistinctiveNsites resolves the --distinctive-sites NUM field to
// an absolute pair-count threshold, mirroring diff_sites_init: a value <=
// 1 is a fraction of npairs, a value > 1 is an absolute count.
func resolveDistinctiveNsites(num float64, npairs int) int {
	if num <= 1 {
		return int(float64(npairs) * num)
	}
	return int(num)
}

// newDistinctiveCollector builds a collector for npairs sample pairs. If
// --distinctive-sites is not set, the returned collector is inert.
func newDistinctiveCollector(opts GtcheckOptions, npairs int) *distinctiveCollector {
	dc := &distinctiveCollector{enabled: opts.HasDistinctiveSites, npairs: npairs}
	if !dc.enabled {
		return dc
	}
	dc.words = (npairs + 63) / 64
	nsites := resolveDistinctiveNsites(opts.DistinctiveSites, npairs)
	if nsites > npairs {
		nsites = npairs
	}
	dc.nsites = nsites
	dc.cur = make([]uint64, dc.words)
	dc.rng = newLrand48()
	return dc
}

// resetSite clears the per-site discordance bitset (diff_sites_reset).
func (dc *distinctiveCollector) resetSite() {
	if !dc.enabled {
		return
	}
	for i := range dc.cur {
		dc.cur[i] = 0
	}
	dc.curNdiff = 0
}

// markDiff records that pair index pi is discordant at the current site.
func (dc *distinctiveCollector) markDiff(pi int) {
	if !dc.enabled {
		return
	}
	dc.cur[pi/64] |= 1 << uint(pi%64)
	dc.curNdiff++
}

// pushSite stores the current site's bitset if it has any discordances
// (diff_sites_push), then resets the read-order counter for tie-breaks.
func (dc *distinctiveCollector) pushSite(chrom string, pos int) {
	if !dc.enabled || dc.curNdiff == 0 {
		return
	}
	mask := make([]uint64, dc.words)
	copy(mask, dc.cur)
	dc.sites = append(dc.sites, diffSite{
		chrom: chrom,
		pos:   pos,
		ndiff: dc.curNdiff,
		mask:  mask,
		rand:  dc.rng.next(),
		order: dc.readSeq,
	})
	dc.readSeq++
}

// report writes the "#DS" block after the discordance table, mirroring
// report_distinctive_sites().
func (dc *distinctiveCollector) report(out io.Writer) error {
	if !dc.enabled {
		return nil
	}
	w := bufio.NewWriter(out)

	// NOTE: upstream report_distinctive_sites() builds a "# DS" comment
	// header into its kstring but never writes it to the output stream
	// (only the DS data rows are emitted via bgzf_write). We reproduce
	// that exactly and emit no comment block — see dsHeader for the text
	// upstream constructs but drops.

	// Sort by ndiff descending, ties broken by the hts_lrand48() key
	// ascending (matching upstream diff_sites_cmp), with the read order
	// as a final deterministic tiebreak.
	sorted := make([]diffSite, len(dc.sites))
	copy(sorted, dc.sites)
	sort.SliceStable(sorted, func(i, j int) bool {
		if sorted[i].ndiff != sorted[j].ndiff {
			return sorted[i].ndiff > sorted[j].ndiff
		}
		if sorted[i].rand != sorted[j].rand {
			return sorted[i].rand < sorted[j].rand
		}
		return sorted[i].order < sorted[j].order
	})

	ndiffMin := dc.nsites
	if ndiffMin > dc.npairs {
		ndiffMin = dc.npairs
	}

	seen := make([]uint64, dc.words)
	ndiffTot := 0
	iblock := 0
	for _, s := range sorted {
		ndiffNew := 0
		for word := 0; word < dc.words; word++ {
			bits := s.mask[word] &^ seen[word]
			if bits == 0 {
				continue
			}
			ndiffNew += popcount(bits)
			seen[word] |= bits
		}
		if ndiffNew == 0 {
			continue
		}
		ndiffTot += ndiffNew
		w.WriteString("DS\t")
		w.WriteString(s.chrom)
		w.WriteByte('\t')
		w.WriteString(strconv.Itoa(s.pos))
		w.WriteByte('\t')
		w.WriteString(strconv.Itoa(ndiffTot))
		w.WriteByte('\t')
		w.WriteString(strconv.Itoa(iblock))
		w.WriteByte('\n')
		if ndiffTot < ndiffMin {
			continue
		}
		iblock++
		ndiffTot = 0
		for i := range seen {
			seen[i] = 0
		}
	}
	return w.Flush()
}

// dsHeader is the literal "#DS" comment block upstream writes.
const dsHeader = `# DS, distinctive sites:
#     - chromosome
#     - position
#     - cumulative number of pairs distinguished by this block
#     - block id
#DS	[2]Chromosome	[3]Position	[4]Cumulative number of distinct pairs	[5]Block id
`

// popcount returns the number of set bits in x.
func popcount(x uint64) int {
	n := 0
	for x != 0 {
		x &= x - 1
		n++
	}
	return n
}
