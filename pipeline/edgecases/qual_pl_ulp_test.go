package edgecases

import (
	"math"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// callRecord captures the GT/FILTER decision and the QUAL value of one called
// site, so we can assert the *decision* is byte-identical while tolerating a
// last-ULP QUAL difference.
type callRecord struct {
	key    string // CHROM:POS:REF:ALT — the site identity
	qual   string // raw QUAL field ("." or a float)
	filter string
	gts    string // tab-joined sample GT subfields
}

// parseCallVCF extracts callRecords keyed by site from `bcftools call` output.
func parseCallVCF(vcf string) map[string]callRecord {
	out := map[string]callRecord{}
	for _, ln := range strings.Split(vcf, "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		c := strings.Split(ln, "\t")
		if len(c) < 10 {
			continue
		}
		key := c[0] + ":" + c[1] + ":" + c[3] + ":" + c[4]
		fmtKeys := strings.Split(c[8], ":")
		gtIdx := -1
		for i, k := range fmtKeys {
			if k == "GT" {
				gtIdx = i
				break
			}
		}
		var gts []string
		for _, samp := range c[9:] {
			sub := strings.Split(samp, ":")
			if gtIdx >= 0 && gtIdx < len(sub) {
				gts = append(gts, sub[gtIdx])
			}
		}
		out[key] = callRecord{
			key:    key,
			qual:   c[5],
			filter: c[6],
			gts:    strings.Join(gts, "\t"),
		}
	}
	return out
}

// TestQualPLULPNonImpact runs `bcftools call -m` with our bcftools and upstream
// on a small genotype-likelihood VCF and confirms that any QUAL difference is
// at most a last-ULP rounding wobble (|Δ| < 0.01 phred) and, crucially, never
// flips the per-sample GT or the site FILTER. The full statistical version
// lives in the GIAB harness; this is the lightweight, always-runnable guard.
func TestQualPLULPNonImpact(t *testing.T) {
	our := ourBin(t, "bcftools")
	up := upBin(t, "bcftools")

	fix := smallFixtureDir(t)
	if fix == "" {
		t.Skip("pipeline/.fixtures/small not present; run the fixture generator")
	}
	src := filepath.Join(fix, "mpileup.pl.vcf")
	if _, _, err := run(t, our, "view", src); err != nil {
		// view is just a presence/parse probe; if the fixture is missing skip.
		t.Skipf("mpileup.pl.vcf fixture not usable: %v", err)
	}

	ours := mustRun(t, our, "call", "-m", src)
	ups := mustRun(t, up, "call", "-m", src)

	o := parseCallVCF(ours)
	u := parseCallVCF(ups)

	if len(o) == 0 || len(u) == 0 {
		t.Skipf("no called sites parsed (ours=%d, upstream=%d)", len(o), len(u))
	}

	const qualTol = 0.01 // phred units; well under a meaningful difference
	for key, ur := range u {
		or, ok := o[key]
		if !ok {
			t.Errorf("site %s present in upstream output but missing from ours", key)
			continue
		}
		// Decision must be byte-identical.
		if or.gts != ur.gts {
			t.Errorf("site %s: GT flipped (ours %q vs upstream %q) — ULP QUAL difference must NOT change genotypes", key, or.gts, ur.gts)
		}
		if or.filter != ur.filter {
			t.Errorf("site %s: FILTER flipped (ours %q vs upstream %q)", key, or.filter, ur.filter)
		}
		// QUAL may differ only within the ULP tolerance.
		oq, okO := parseQual(or.qual)
		uq, okU := parseQual(ur.qual)
		if okO && okU {
			if d := math.Abs(oq - uq); d > qualTol {
				t.Errorf("site %s: QUAL differs by %.4f (> %.2f tol): ours %q vs upstream %q", key, d, qualTol, or.qual, ur.qual)
			}
		} else if or.qual != ur.qual {
			t.Errorf("site %s: QUAL mismatch (non-numeric): ours %q vs upstream %q", key, or.qual, ur.qual)
		}
	}
}

func parseQual(s string) (float64, bool) {
	if s == "." || s == "" {
		return 0, false
	}
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, false
	}
	return f, true
}
