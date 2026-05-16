package bcftools

import (
	"bytes"
	"strings"
	"testing"
)

// fixtureGtcheck is a 3-sample, 4-site biallelic VCF. Site #2 has a
// missing GT for S2 (skip, not discordance).
const fixtureGtcheck = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	T	.	.	.	GT	0/0	0/0	1/1
chr1	200	.	C	G	.	.	.	GT	0/1	./.	1/1
chr1	300	.	G	A	.	.	.	GT	1/1	1/1	0/0
chr1	400	.	T	C	.	.	.	GT	0/0	0/1	1/1
`

// TestGtcheck_HeaderIsDCv2 pins the first output line to upstream's
// literal `#DCv2 + 6 column descriptors` shape, NOT the home-grown
// `#DC` header from the closed PR #106.
func TestGtcheck_HeaderIsDCv2(t *testing.T) {
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	first := strings.SplitN(out.String(), "\n", 2)[0]
	want := "#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[7]Number of matching genotypes"
	if first != want {
		t.Fatalf("first line:\n got %q\nwant %q", first, want)
	}
	// There must be NO trailing "# totals..." line.
	if strings.Contains(out.String(), "# totals") {
		t.Errorf("output contains forbidden '# totals' trailer:\n%s", out.String())
	}
}

// TestGtcheck_HammingScores: pin the score arithmetic. For the
// fixture above:
//
//	S1 dosages: 0, 1, 2, 0
//	S2 dosages: 0, MISS, 2, 1
//	S3 dosages: 2, 2, 0, 2
//
// Pair (S1,S2): sites where both non-missing = {100,300,400}; diffs = |0-0|+|2-2|+|0-1| = 1.
// Pair (S1,S3): all 4 sites; diffs = |0-2|+|1-2|+|2-0|+|0-2| = 7.
// Pair (S2,S3): sites where both non-missing = {100,300,400}; diffs = |0-2|+|2-0|+|1-2| = 5.
func TestGtcheck_HammingScores(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	want := map[[2]string]struct {
		disc, sites, match int
	}{
		{"S1", "S2"}: {1, 3, 2},
		{"S1", "S3"}: {7, 4, 0},
		{"S2", "S3"}: {5, 3, 0},
	}
	if len(r.Pairs) != len(want) {
		t.Fatalf("npairs got %d, want %d", len(r.Pairs), len(want))
	}
	for _, p := range r.Pairs {
		w, ok := want[[2]string{p.QuerySample, p.GenotypedSample}]
		if !ok {
			t.Errorf("unexpected pair (%s,%s)", p.QuerySample, p.GenotypedSample)
			continue
		}
		if p.Discordance != w.disc {
			t.Errorf("(%s,%s) discordance: got %d, want %d", p.QuerySample, p.GenotypedSample, p.Discordance, w.disc)
		}
		if p.NumSites != w.sites {
			t.Errorf("(%s,%s) sites: got %d, want %d", p.QuerySample, p.GenotypedSample, p.NumSites, w.sites)
		}
		if p.NumMatching != w.match {
			t.Errorf("(%s,%s) match: got %d, want %d", p.QuerySample, p.GenotypedSample, p.NumMatching, w.match)
		}
	}
}

// TestGtcheck_MissingGTIsSkipNotDiscordance directly tests the
// reviewer requirement: a `./.` GT must NOT count as a discordance,
// it must count as a skip.
func TestGtcheck_MissingGTIsSkipNotDiscordance(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	100	.	A	T	.	.	.	GT	0/0	./.
chr1	200	.	C	G	.	.	.	GT	0/0	0/0
`
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(src), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	if len(r.Pairs) != 1 {
		t.Fatalf("npairs=%d", len(r.Pairs))
	}
	p := r.Pairs[0]
	if p.Discordance != 0 {
		t.Errorf("discordance: got %d, want 0 (missing must skip)", p.Discordance)
	}
	if p.NumSites != 1 {
		t.Errorf("sites: got %d, want 1 (only 1 site has both GTs)", p.NumSites)
	}
}

// TestGtcheck_RejectsMultiAllelic mirrors upstream's "run `bcftools
// norm -m -` first" diagnostic.
func TestGtcheck_RejectsMultiAllelic(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	100	.	A	T,C	.	.	.	GT	0/1	1/2
`
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(src), &out, GtcheckOptions{})
	if err == nil {
		t.Fatalf("expected error for multi-allelic input, got nil")
	}
	if !strings.Contains(err.Error(), "bcftools norm -m -") {
		t.Errorf("error %q does not mention `bcftools norm -m -`", err)
	}
}

// TestGtcheck_PLDeferred ensures `-u PL` is rejected with the roadmap
// pointer rather than silently accepted.
func TestGtcheck_PLDeferred(t *testing.T) {
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{UseTag: "PL"})
	if err == nil {
		t.Fatalf("expected error for -u PL, got nil")
	}
	if !strings.Contains(err.Error(), "docs/PARITY_ROADMAP.md") {
		t.Errorf("PL error must mention PARITY_ROADMAP: got %q", err)
	}
}

// TestGtcheck_Symmetry verifies that swapping (q,g) yields the same
// discordance under cross-check mode; we only emit the i<j half by
// construction, but the score is symmetric in the dosage Hamming
// metric, so an independent swapped run must yield the same number.
func TestGtcheck_Symmetry(t *testing.T) {
	var out1, out2 bytes.Buffer
	r1, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out1, GtcheckOptions{})
	if err != nil {
		t.Fatal(err)
	}
	r2, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out2, GtcheckOptions{PairsSpec: "S3,S1,S2,S1,S3,S2"})
	if err != nil {
		t.Fatal(err)
	}
	get := func(rr GtcheckResult, q, g string) (int, int) {
		for _, p := range rr.Pairs {
			if p.QuerySample == q && p.GenotypedSample == g {
				return p.Discordance, p.NumSites
			}
		}
		return -1, -1
	}
	cases := []struct{ a, b string }{
		{"S1", "S3"}, {"S1", "S2"}, {"S2", "S3"},
	}
	for _, c := range cases {
		d1, s1 := get(r1, c.a, c.b)
		d2, s2 := get(r2, c.b, c.a)
		if d1 != d2 {
			t.Errorf("symmetry: (%s,%s)=%d vs (%s,%s)=%d", c.a, c.b, d1, c.b, c.a, d2)
		}
		if s1 != s2 {
			t.Errorf("symmetry sites: (%s,%s)=%d vs (%s,%s)=%d", c.a, c.b, s1, c.b, c.a, s2)
		}
	}
}

// TestGtcheck_PairsSpec validates the -p qry,gt[,...] parser.
func TestGtcheck_PairsSpec(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{PairsSpec: "S1,S2"})
	if err != nil {
		t.Fatal(err)
	}
	if len(r.Pairs) != 1 || r.Pairs[0].QuerySample != "S1" || r.Pairs[0].GenotypedSample != "S2" {
		t.Errorf("pairs: %+v", r.Pairs)
	}
}

// TestGtcheck_HomsOnly drops sites where the panel GT is heterozygous.
func TestGtcheck_HomsOnly(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{HomsOnly: true, PairsSpec: "S1,S2"})
	if err != nil {
		t.Fatal(err)
	}
	// In fixtureGtcheck pair (S1,S2) we had sites where both non-missing:
	// 100 (S2 GT=0/0 hom), 300 (S2 GT=1/1 hom), 400 (S2 GT=0/1 HET → drop).
	// So homs-only sites=2, discordance=0.
	if len(r.Pairs) != 1 {
		t.Fatalf("npairs=%d", len(r.Pairs))
	}
	p := r.Pairs[0]
	if p.NumSites != 2 || p.Discordance != 0 {
		t.Errorf("homs-only: got sites=%d disc=%d, want sites=2 disc=0", p.NumSites, p.Discordance)
	}
}

// TestGtcheck_NoHWEProb sets the HWE column to 0.0 across the board.
func TestGtcheck_NoHWEProb(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{NoHWEProb: true})
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range r.Pairs {
		if p.AvgLogPHWE != 0.0 {
			t.Errorf("--no-HWE-prob: pair (%s,%s) has nonzero HWE=%v", p.QuerySample, p.GenotypedSample, p.AvgLogPHWE)
		}
	}
}
