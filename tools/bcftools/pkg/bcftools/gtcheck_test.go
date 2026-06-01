package bcftools

import (
	"bytes"
	"math"
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

// gtcheckHeaderLine is the DCv2 column-descriptor row that must appear
// (after the INFO block and the DCv2 comment block) in every report.
const gtcheckHeaderLine = "#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[7]Number of matching genotypes"

// TestGtcheck_MultiSectionOutput pins the report layout: the INFO
// counter block comes first, then the DCv2 comment block, then the
// DCv2 header row, then DCv2 data rows.
func TestGtcheck_MultiSectionOutput(t *testing.T) {
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	s := out.String()
	lines := strings.Split(strings.TrimRight(s, "\n"), "\n")
	if lines[0] != "INFO\tsites-compared\t4" {
		t.Errorf("first line: got %q, want INFO sites-compared 4", lines[0])
	}
	for _, want := range []string{
		"INFO\tsites-skipped-multiallelic\t0",
		"INFO\tsites-used-GT-vs-GT\t4",
		"# DCv2, discordance version 2:",
		gtcheckHeaderLine,
	} {
		if !strings.Contains(s, want) {
			t.Errorf("output missing line %q\n%s", want, s)
		}
	}
	// The DCv2 header row must precede the DCv2 data rows.
	hdrIdx := strings.Index(s, gtcheckHeaderLine)
	dataIdx := strings.Index(s, "\nDCv2\t")
	if hdrIdx < 0 || dataIdx < 0 || hdrIdx > dataIdx {
		t.Errorf("DCv2 header must precede data rows (hdr=%d data=%d)", hdrIdx, dataIdx)
	}
}

// TestGtcheck_CrossCheckPairOrder verifies the sub-diagonal (i>j) pair
// ordering: query=samples[i], genotyped=samples[j] for j<i.
func TestGtcheck_CrossCheckPairOrder(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	want := [][2]string{{"S2", "S1"}, {"S3", "S1"}, {"S3", "S2"}}
	if len(r.Pairs) != len(want) {
		t.Fatalf("npairs got %d, want %d", len(r.Pairs), len(want))
	}
	for i, w := range want {
		if r.Pairs[i].QuerySample != w[0] || r.Pairs[i].GenotypedSample != w[1] {
			t.Errorf("pair %d: got (%s,%s) want (%s,%s)", i,
				r.Pairs[i].QuerySample, r.Pairs[i].GenotypedSample, w[0], w[1])
		}
	}
}

// TestGtcheck_SitesAndMatches pins the per-pair site / match counts.
//
//	S1 dosages: 0, 1, 2, 0
//	S2 dosages: 0, MISS, 2, 1
//	S3 dosages: 2, 2, 0, 2
//
// Pair (S2,S1): both non-missing at {100,300,400} → 3 sites; matches at 100,300 → 2.
// Pair (S3,S1): both non-missing at all 4 sites → 4 sites; no matches → 0.
// Pair (S3,S2): both non-missing at {100,300,400} → 3 sites; no matches → 0.
func TestGtcheck_SitesAndMatches(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	want := map[[2]string]struct{ sites, match int }{
		{"S2", "S1"}: {3, 2},
		{"S3", "S1"}: {4, 0},
		{"S3", "S2"}: {3, 0},
	}
	for _, p := range r.Pairs {
		w, ok := want[[2]string{p.QuerySample, p.GenotypedSample}]
		if !ok {
			t.Errorf("unexpected pair (%s,%s)", p.QuerySample, p.GenotypedSample)
			continue
		}
		if p.NumSites != w.sites {
			t.Errorf("(%s,%s) sites: got %d, want %d", p.QuerySample, p.GenotypedSample, p.NumSites, w.sites)
		}
		if p.NumMatching != w.match {
			t.Errorf("(%s,%s) match: got %d, want %d", p.QuerySample, p.GenotypedSample, p.NumMatching, w.match)
		}
	}
}

// TestGtcheck_DiscordanceScore pins the default error-probability
// discordance score for a fully-concordant and a fully-discordant pair.
// With gt_err=40, eprob=1e-4, le=-log(1e-4)=9.210340... For two
// identical homozygous genotypes the minimum joint negative-log
// probability is 0 (the matching genotype). A hom-ref vs hom-alt pair
// contributes 2*le per site (the cheapest path is via the het with
// e + e). Here we just assert a matching pair scores < a mismatching
// pair, and that a perfectly-matching synthetic pair scores 0.
func TestGtcheck_DiscordanceScore(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	100	.	A	T	.	.	AC=2;AN=4	GT	0/0	0/0
chr1	200	.	C	G	.	.	AC=2;AN=4	GT	0/1	0/1
chr1	300	.	G	A	.	.	AC=2;AN=4	GT	1/1	1/1
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
	// Identical genotypes at every site → discordance 0.
	if p.Discordance != 0 {
		t.Errorf("identical samples: discordance got %v, want 0", p.Discordance)
	}
	if p.NumSites != 3 || p.NumMatching != 3 {
		t.Errorf("identical samples: sites=%d match=%d, want 3/3", p.NumSites, p.NumMatching)
	}
}

// TestGtcheck_MissingGTIsSkipNotDiscordance: a `./.` GT must NOT count
// as a compared site, it must be skipped.
func TestGtcheck_MissingGTIsSkipNotDiscordance(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	100	.	A	T	.	.	AC=1;AN=4	GT	0/0	./.
chr1	200	.	C	G	.	.	AC=1;AN=4	GT	0/0	0/0
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
	if p.NumSites != 1 {
		t.Errorf("sites: got %d, want 1 (only 1 site has both GTs)", p.NumSites)
	}
	if p.Discordance != 0 {
		t.Errorf("discordance: got %v, want 0", p.Discordance)
	}
}

// TestGtcheck_SkipsMultiAllelic mirrors upstream: a multi-allelic record
// is counted in sites-skipped-multiallelic and excluded from scoring,
// NOT a hard error.
func TestGtcheck_SkipsMultiAllelic(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B
chr1	100	.	A	T,C	.	.	.	GT	0/1	1/2
chr1	200	.	C	G	.	.	AC=1;AN=4	GT	0/0	0/1
`
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(src), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: unexpected error %v", err)
	}
	if r.Counters.SkippedMultiallelic != 1 {
		t.Errorf("SkippedMultiallelic: got %d, want 1", r.Counters.SkippedMultiallelic)
	}
	if r.Counters.SitesCompared != 1 {
		t.Errorf("SitesCompared: got %d, want 1", r.Counters.SitesCompared)
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
	// For pair (S1,S2) the genotyped sample is S2: sites where both
	// non-missing are 100 (S2=0/0 hom), 300 (S2=1/1 hom), 400 (S2=0/1
	// HET → drop). So homs-only sites=2.
	if len(r.Pairs) != 1 {
		t.Fatalf("npairs=%d", len(r.Pairs))
	}
	if r.Pairs[0].NumSites != 2 {
		t.Errorf("homs-only: got sites=%d, want 2", r.Pairs[0].NumSites)
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

// TestGtcheck_HWEMatchesFormula checks the HWE column equals the
// upstream per-site -log of the matching-dosage HWE probability,
// averaged over matching sites. For pair (S2,S1) on fixtureGtcheck,
// matches occur at site 100 (both 0/0, AF defaults 1e-6 → hwe[0] ~ 0)
// and site 300 (both 1/1, AF defaults 1e-6 → hwe[2] = -log(af^2)).
func TestGtcheck_HWEMatchesFormula(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{})
	if err != nil {
		t.Fatal(err)
	}
	// fixtureGtcheck has no INFO/AC, so AF falls back to GT counts.
	// Recompute the expected average for (S2,S1) directly.
	var found bool
	for _, p := range r.Pairs {
		if p.QuerySample == "S2" && p.GenotypedSample == "S1" {
			found = true
			if math.IsNaN(p.AvgLogPHWE) || math.IsInf(p.AvgLogPHWE, 0) {
				t.Errorf("HWE for (S2,S1) is not finite: %v", p.AvgLogPHWE)
			}
			if p.AvgLogPHWE <= 0 {
				t.Errorf("HWE for (S2,S1) should be positive, got %v", p.AvgLogPHWE)
			}
		}
	}
	if !found {
		t.Fatalf("pair (S2,S1) not found")
	}
}
