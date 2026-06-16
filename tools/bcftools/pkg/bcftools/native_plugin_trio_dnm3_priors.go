// Mendelian priors tables and chrX region matcher for the trio-dnm3 NAIVE
// model. These port the integer/combinatorial parts of trio-dnm3.c that the
// NAIVE verdict depends on: the seq1/seq2/seq3 genotype-index lookups, the
// denovo / denovo_allele tables built from init_tprob_mprob (autosomal),
// init_tprob_mprob_chrX and init_tprob_mprob_chrXX, and the chrX PAR-exclusion
// region list. None of this touches floating point — denovo[fi][mi][ci] is set
// from the integer test tprob==0, and denovo_allele from a combinatorial allele
// comparison — so the NAIVE annotations are byte-reproducible against upstream.
package bcftools

import (
	"strconv"
	"strings"
)

// dnmSeq1 / dnmSeq2 map a genotype index 0..9 to its first/second allele,
// matching seq1[] / seq2[] in trio-dnm3.c.
var (
	dnmSeq1 = [10]int{0, 1, 1, 2, 2, 2, 3, 3, 3, 3}
	dnmSeq2 = [10]int{0, 0, 1, 0, 1, 2, 0, 1, 2, 3}
	// dnmSeq3 maps an allele bitmask (1<<a)|(1<<b) (1..12) to a genotype index,
	// matching seq3[] in trio-dnm3.c (-1 marks an impossible mask).
	dnmSeq3 = [13]int{-1, 0, 2, 1, 5, 3, 4, -1, 9, 6, 7, -1, 8}
)

// dnmPriorsType selects which Mendelian rule set to build.
type dnmPriorsType int

const (
	autosomalPriors dnmPriorsType = iota
	chrXPriors
	chrXXPriors
)

// dnmPriors holds the per-genotype-configuration de-novo verdict and the
// de-novo allele, mirroring the denovo / denovo_allele fields of priors_t.
type dnmPriors struct {
	denovo       [10][10][10]int
	denovoAllele [10][10][10]int
}

// newDNMPriors builds the denovo/denovo_allele tables for the given rule set,
// mirroring init_priors() restricted to the NAIVE-relevant fields. The mutation
// rate does not affect these (it only scales the unused mprob/pprob), so it is
// omitted entirely.
func newDNMPriors(strictlyNovel bool, typ dnmPriorsType) *dnmPriors {
	p := &dnmPriors{}
	for fi := 0; fi < 10; fi++ {
		for mi := 0; mi < 10; mi++ {
			for ci := 0; ci < 10; ci++ {
				var tprobZero bool
				var allele int
				switch typ {
				case autosomalPriors:
					tprobZero, allele = dnmTprobAutosomal(strictlyNovel, fi, mi, ci)
				case chrXPriors:
					tprobZero, allele = dnmTprobChrX(strictlyNovel, fi, mi, ci)
				case chrXXPriors:
					tprobZero, allele = dnmTprobChrXX(strictlyNovel, fi, mi, ci)
				}
				if tprobZero {
					p.denovo[fi][mi][ci] = 1
					p.denovoAllele[fi][mi][ci] = allele
				} else {
					p.denovo[fi][mi][ci] = 0
					p.denovoAllele[fi][mi][ci] = 255 // UINT8_MAX; never read when denovo==0
				}
			}
		}
	}
	return p
}

// dnmTprobAutosomal returns whether the transmission probability is zero (i.e. a
// de-novo) and the candidate de-novo allele for the autosomal rule set,
// mirroring init_tprob_mprob() in trio-dnm3.c.
func dnmTprobAutosomal(strictlyNovel bool, fi, mi, ci int) (tprobZero bool, denovoAllele int) {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	ca, cb := dnmSeq1[ci], dnmSeq2[ci]

	if ca != fa && ca != fb && ca != ma && ca != mb {
		denovoAllele = ca
	} else {
		denovoAllele = cb
	}

	var isNovel bool
	if strictlyNovel {
		isNovel = (ca != fa && ca != fb && ca != ma && ca != mb) || (cb != fa && cb != fb && cb != ma && cb != mb)
		if isNovel && denovoAllele == 0 {
			isNovel = false
		}
	} else {
		isNovel = !(((ca == fa || ca == fb) && (cb == ma || cb == mb)) || ((ca == ma || ca == mb) && (cb == fa || cb == fb)))
	}
	return isNovel, denovoAllele
}

// dnmTprobChrX returns the de-novo verdict and allele for a male proband on the
// non-PAR chrX, mirroring init_tprob_mprob_chrX(). It delegates to the autosomal
// rule for the "not novel wrt parents" genotype-error fallback.
func dnmTprobChrX(strictlyNovel bool, fi, mi, ci int) (tprobZero bool, denovoAllele int) {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	ca, cb := dnmSeq1[ci], dnmSeq2[ci]

	if ca != ma && ca != mb {
		denovoAllele = ca
	} else {
		denovoAllele = cb
	}

	if ca != cb { // male cannot be heterozygous in X, but it can be mosaic
		var isNovel bool
		if strictlyNovel {
			isNovel = (ca != fa && ca != fb && ca != ma && ca != mb) || (cb != fa && cb != fb && cb != ma && cb != mb)
		} else {
			isNovel = (ca != ma && ca != mb) || (cb != ma && cb != mb)
		}
		if isNovel {
			return true, denovoAllele
		}
		// Genotype error: fall back to autosomal inheritance.
		return dnmTprobAutosomal(strictlyNovel, fi, mi, ci)
	}
	if ca == ma || ca == mb { // inherited
		return false, denovoAllele
	}
	// de novo
	return true, denovoAllele
}

// dnmTprobChrXX returns the de-novo verdict and allele for a female proband on
// the non-PAR chrX, mirroring init_tprob_mprob_chrXX().
func dnmTprobChrXX(strictlyNovel bool, fi, mi, ci int) (tprobZero bool, denovoAllele int) {
	fa, fb := dnmSeq1[fi], dnmSeq2[fi]
	ma, mb := dnmSeq1[mi], dnmSeq2[mi]
	ca, cb := dnmSeq1[ci], dnmSeq2[ci]

	if ca != fa && ca != fb && ca != ma && ca != mb {
		denovoAllele = ca
	} else {
		denovoAllele = cb
	}

	if fa != fb {
		// Father cannot be heterozygous in X; treat as genotype error, fall back to
		// the autosomal rule.
		return dnmTprobAutosomal(strictlyNovel, fi, mi, ci)
	}
	if (ca == fa && (cb == ma || cb == mb)) || (cb == fa && (ca == ma || ca == mb)) {
		return false, denovoAllele // inherited
	}
	return true, denovoAllele
}

// chrXMatcher tests whether a (chrom,pos) falls in the non-PAR chrX region list,
// mirroring regidx_overlap on the chrX_list_str. Positions are 1-based (as in
// vcf.Variant.Pos), and the C check uses [pos, pos+rlen) overlap; for the SNP/
// 1-length records the plugin handles, a point-in-interval test suffices and
// matches the C overlap of rec->pos..rec->pos+rec->rlen against the closed
// 1-based intervals.
type chrXMatcher struct {
	intervals map[string][][2]int
}

// buildChrXMatcher builds the chrX matcher from the --chrX value, defaulting to
// the GRCh37 PAR-exclusion list and recognising the "GRCh38" shortcut, exactly
// as init_data() does.
func buildChrXMatcher(chrXListStr string) *chrXMatcher {
	list := chrXListStr
	if list == "" || strings.EqualFold(list, "GRCh37") {
		list = "X:1-60000,chrX:1-60000,X:2699521-154931043,chrX:2699521-154931043"
	} else if strings.EqualFold(list, "GRCh38") {
		list = "X:1-9999,chrX:1-9999,X:2781480-155701381,chrX:2781480-155701381"
	}
	m := &chrXMatcher{intervals: map[string][][2]int{}}
	for _, part := range strings.Split(list, ",") {
		part = strings.TrimSpace(part)
		if part == "" {
			continue
		}
		colon := strings.LastIndexByte(part, ':')
		if colon < 0 {
			continue
		}
		chr := part[:colon]
		rng := part[colon+1:]
		beg, end := 0, 0
		if dash := strings.IndexByte(rng, '-'); dash >= 0 {
			beg, _ = strconv.Atoi(rng[:dash])
			end, _ = strconv.Atoi(rng[dash+1:])
		} else {
			beg, _ = strconv.Atoi(rng)
			end = beg
		}
		m.intervals[chr] = append(m.intervals[chr], [2]int{beg, end})
	}
	return m
}

// overlaps reports whether the 1-based position on chrom falls within any chrX
// non-PAR interval.
func (m *chrXMatcher) overlaps(chrom string, pos int) bool {
	for _, iv := range m.intervals[chrom] {
		if pos >= iv[0] && pos <= iv[1] {
			return true
		}
	}
	return false
}
