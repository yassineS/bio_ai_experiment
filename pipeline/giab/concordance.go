package giab

import (
	"fmt"
	"sort"
	"strings"
)

// DefaultQualULP is the absolute QUAL difference, in Phred units, at or below
// which a QUAL discrepancy is classified as a last-place (ULP) wobble rather
// than a substantive change. bcftools emits QUAL to two decimals, and the libm
// floor that motivates this experiment perturbs the last printed digit at most.
// A value of 0.5 comfortably brackets a one-printed-unit difference while
// staying well below any genotype-changing QUAL swing.
const DefaultQualULP = 0.5

// SiteDiff describes how a single locus differs between the two call sets.
type SiteDiff struct {
	Key        string `json:"key"`
	Chrom      string `json:"chrom"`
	Pos        int    `json:"pos"`
	OursQual   string `json:"ours_qual"`
	UpQual     string `json:"up_qual"`
	OursPL     string `json:"ours_pl,omitempty"`
	UpPL       string `json:"up_pl,omitempty"`
	OursGT     string `json:"ours_gt"`
	UpGT       string `json:"up_gt"`
	OursPass   bool   `json:"ours_pass"`
	UpPass     bool   `json:"up_pass"`
	QualULP    bool   `json:"qual_ulp_only"` // QUAL/PL differ only at the ULP/Phred floor
	GTFlip     bool   `json:"gt_flip"`       // the genotype changed
	FilterFlip bool   `json:"filter_flip"`   // PASS/FAIL changed
	Note       string `json:"note,omitempty"`
}

// Flips reports whether this difference flips a genotype or a PASS/FAIL verdict.
func (d SiteDiff) Flips() bool { return d.GTFlip || d.FilterFlip }

// Concordance is the result of comparing OUR call set against UPSTREAM's,
// restricted to a region set (the high-confidence BED).
type Concordance struct {
	// Sites considered (intersection of both call sets within the region).
	Common int `json:"common"`
	// Records present in ours but not upstream, and vice versa.
	OnlyOurs int `json:"only_ours"`
	OnlyUp   int `json:"only_upstream"`
	// Identical records (REF/ALT/QUAL/FILTER/GT/PL all equal).
	Identical int `json:"identical"`
	// Records that differ in any compared field.
	Differ int `json:"differ"`
	// Of the differing records, how many differ ONLY in QUAL (and/or PL) at the
	// ULP/Phred floor with identical GT and FILTER.
	QualULPOnly int `json:"qual_ulp_only"`
	// Of the differing records, how many flip a genotype or PASS/FAIL verdict.
	GenotypeOrFilterFlips int `json:"genotype_or_filter_flips"`
	// The differing sites (capped for report sanity by the caller if needed).
	Diffs []SiteDiff `json:"diffs,omitempty"`
}

// Headline renders the marquee sentence: "N sites differ in QUAL, 0 of which
// flip a genotype or PASS/FAIL".
func (c Concordance) Headline() string {
	return fmt.Sprintf("%d/%d records identical; %d differ (%d only at the QUAL/PL ULP floor), %d of which flip a genotype or PASS/FAIL",
		c.Identical, c.Common, c.Differ, c.QualULPOnly, c.GenotypeOrFilterFlips)
}

// CompareCallSets compares OUR records against UPSTREAM's, considering only
// records whose locus falls inside the supplied region set. Pass a nil or empty
// RegionSet to consider all records (the comparator still works; it just is not
// restricted to the high-confidence BED).
//
// Two records are paired by Key (CHROM:POS:REF:ALT). A record present in only
// one set is counted in OnlyOurs/OnlyUp. Paired records are compared on QUAL,
// PL, GT and FILTER; any difference makes them "differ". A differing pair is
// classified as QUAL/PL-ULP-only when GT and FILTER are unchanged and every
// numeric QUAL/PL discrepancy is within qualULP Phred units (use
// DefaultQualULP); otherwise the genotype/filter flip is recorded.
func CompareCallSets(ours, up []VCFRecord, region *RegionSet, qualULP float64) Concordance {
	inRegion := func(r VCFRecord) bool {
		if region == nil || region.Empty() {
			return true
		}
		return region.Contains(r.Chrom, r.Pos)
	}

	oursByKey := map[string]VCFRecord{}
	for _, r := range ours {
		if inRegion(r) {
			oursByKey[r.Key()] = r
		}
	}
	upByKey := map[string]VCFRecord{}
	for _, r := range up {
		if inRegion(r) {
			upByKey[r.Key()] = r
		}
	}

	var c Concordance
	for k, o := range oursByKey {
		u, ok := upByKey[k]
		if !ok {
			c.OnlyOurs++
			continue
		}
		c.Common++
		diff, identical := classify(o, u, qualULP)
		if identical {
			c.Identical++
			continue
		}
		c.Differ++
		if diff.QualULP && !diff.Flips() {
			c.QualULPOnly++
		}
		if diff.Flips() {
			c.GenotypeOrFilterFlips++
		}
		c.Diffs = append(c.Diffs, diff)
	}
	for k := range upByKey {
		if _, ok := oursByKey[k]; !ok {
			c.OnlyUp++
		}
	}

	sort.Slice(c.Diffs, func(i, j int) bool {
		if c.Diffs[i].Chrom != c.Diffs[j].Chrom {
			return c.Diffs[i].Chrom < c.Diffs[j].Chrom
		}
		return c.Diffs[i].Pos < c.Diffs[j].Pos
	})
	return c
}

// classify compares two paired records and returns the SiteDiff plus whether
// they are fully identical (in which case the SiteDiff is the zero value and
// should be ignored).
func classify(o, u VCFRecord, qualULP float64) (SiteDiff, bool) {
	oGT, uGT := o.GT(), u.GT()
	oPass, uPass := o.PassFail(), u.PassFail()
	oQual, uQual := o.Qual, u.Qual
	oPL, uPL := o.PL(), u.PL()

	gtSame := oGT == uGT
	passSame := oPass == uPass
	qualSame := strings.TrimSpace(oQual) == strings.TrimSpace(uQual)
	plSame := oPL == uPL

	if gtSame && passSame && qualSame && plSame {
		return SiteDiff{}, true
	}

	d := SiteDiff{
		Key:        o.Key(),
		Chrom:      o.Chrom,
		Pos:        o.Pos,
		OursQual:   oQual,
		UpQual:     uQual,
		OursPL:     oPL,
		UpPL:       uPL,
		OursGT:     oGT,
		UpGT:       uGT,
		OursPass:   oPass,
		UpPass:     uPass,
		GTFlip:     !gtSame,
		FilterFlip: !passSame,
	}

	// A "QUAL/PL ULP-only" difference is one where the genotype and the
	// pass/fail verdict are unchanged AND every numeric QUAL/PL discrepancy is
	// within the ULP tolerance. If GT or FILTER changed, it is a real flip and
	// not a ULP wobble regardless of the numbers.
	qualULPOnly := qualWithinULP(oQual, uQual, qualULP) && plWithinULP(oPL, uPL, qualULP)
	d.QualULP = gtSame && passSame && qualULPOnly && (!qualSame || !plSame)

	switch {
	case d.GTFlip && d.FilterFlip:
		d.Note = "genotype and FILTER both changed"
	case d.GTFlip:
		d.Note = "genotype changed"
	case d.FilterFlip:
		d.Note = "FILTER (PASS/FAIL) changed"
	case d.QualULP:
		d.Note = "QUAL/PL differ only at the ULP/Phred floor; GT and FILTER unchanged"
	default:
		d.Note = "QUAL/PL differ beyond the ULP floor; GT and FILTER unchanged"
	}
	return d, false
}

// qualWithinULP reports whether two QUAL tokens are equal or differ only within
// the ULP tolerance. Missing ("." or "") on both sides is within; missing on
// one side only is NOT within (a real change).
func qualWithinULP(a, b string, tol float64) bool {
	if strings.TrimSpace(a) == strings.TrimSpace(b) {
		return true
	}
	av, aok := qualToFloat(a)
	bv, bok := qualToFloat(b)
	if aok != bok {
		return false
	}
	if !aok && !bok {
		return true
	}
	return floatsClose(av, bv, tol)
}

// plWithinULP reports whether two PL strings (comma-separated per-genotype
// phred lists, e.g. "0,30,255") are equal or differ only within the ULP
// tolerance element-by-element. Mismatched element counts are a real change.
//
// PL is an integer Phred field, so its ULP/last-place unit is 1 — a difference
// of a single Phred point is the integer analogue of a QUAL last-decimal
// wobble. The effective per-element tolerance is therefore max(tol, 1.0), so
// the libm-floor case (PL entries differing by one) classifies as ULP-only.
func plWithinULP(a, b string, tol float64) bool {
	if a == b {
		return true
	}
	if a == "" || b == "" {
		return false
	}
	if tol < 1.0 {
		tol = 1.0
	}
	as := strings.Split(a, ",")
	bs := strings.Split(b, ",")
	if len(as) != len(bs) {
		return false
	}
	for i := range as {
		av, aok := qualToFloat(as[i])
		bv, bok := qualToFloat(bs[i])
		if aok != bok {
			return false
		}
		if !aok {
			if strings.TrimSpace(as[i]) != strings.TrimSpace(bs[i]) {
				return false
			}
			continue
		}
		if !floatsClose(av, bv, tol) {
			return false
		}
	}
	return true
}
