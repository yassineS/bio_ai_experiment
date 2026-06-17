package bcftools

import (
	"bytes"
	"compress/gzip"
	"io"
	"math"
	"strings"
	"testing"
)

// fixtureGtcheck is a 3-sample, 4-site biallelic VCF. Site #2 (chr1:200)
// has a missing GT for S2 (skip, not a discordance). PL is attached so
// the -u PL path can be exercised on the same data.
const fixtureGtcheck = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	T	.	.	.	GT:PL	0/0:0,30,255	0/0:0,30,255	1/1:255,30,0
chr1	200	.	C	G	.	.	.	GT:PL	0/1:30,0,30	./.:0,0,0	1/1:255,30,0
chr1	300	.	G	A	.	.	.	GT:PL	1/1:255,30,0	1/1:255,30,0	0/0:0,30,255
chr1	400	.	T	C	.	.	.	GT:PL	0/0:0,30,255	0/1:30,0,30	1/1:255,30,0
`

// dcRow holds the parsed fields of one DCv2 output row.
type dcRow struct {
	qry, gt string
	disc    string // raw discordance field (int or %e)
	hwe     string
	nsites  string
	nmatch  string
}

// parseDCv2 extracts the DCv2 data rows from the gtcheck output.
func parseDCv2(out string) []dcRow {
	var rows []dcRow
	for _, line := range strings.Split(out, "\n") {
		if !strings.HasPrefix(line, "DCv2\t") {
			continue
		}
		f := strings.Split(line, "\t")
		if len(f) != 7 {
			continue
		}
		rows = append(rows, dcRow{f[1], f[2], f[3], f[4], f[5], f[6]})
	}
	return rows
}

// TestGtcheck_HeaderIsDCv2 pins the column header line to upstream's
// literal "#DCv2 + 6 column descriptors" shape.
func TestGtcheck_HeaderIsDCv2(t *testing.T) {
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	want := "#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[7]Number of matching genotypes"
	if !strings.Contains(out.String(), want) {
		t.Fatalf("output missing DCv2 header line:\n%s", out.String())
	}
	// The INFO stats block must precede the table.
	if !strings.Contains(out.String(), "INFO\tsites-compared\t4") {
		t.Errorf("output missing 'INFO\\tsites-compared\\t4':\n%s", out.String())
	}
}

// TestUnitGtcheckOutputTypeZ checks the -O z container selection at the
// library level with no upstream binary: OutputType "z" must produce a
// gzip/BGZF stream whose decompressed bytes are identical to the OutputType
// "t" text output, and the gzip magic must be present. This is the
// binary-free counterpart to the live BGZF parity test.
func TestUnitGtcheckOutputTypeZ(t *testing.T) {
	var textBuf bytes.Buffer
	if _, err := Gtcheck(strings.NewReader(fixtureGtcheck), &textBuf, GtcheckOptions{OutputType: "t"}); err != nil {
		t.Fatalf("Gtcheck -O t: %v", err)
	}
	for _, level := range []int{-1, 0, 6, 9} {
		var zBuf bytes.Buffer
		if _, err := Gtcheck(strings.NewReader(fixtureGtcheck), &zBuf, GtcheckOptions{OutputType: "z", CompressLevel: level}); err != nil {
			t.Fatalf("Gtcheck -O z (level %d): %v", level, err)
		}
		raw := zBuf.Bytes()
		if len(raw) < 2 || raw[0] != 0x1f || raw[1] != 0x8b {
			t.Fatalf("level %d: output is not gzip-framed: % x", level, raw[:min(2, len(raw))])
		}
		zr, err := gzip.NewReader(bytes.NewReader(raw))
		if err != nil {
			t.Fatalf("level %d: gzip.NewReader: %v", level, err)
		}
		zr.Multistream(true)
		got, err := io.ReadAll(zr)
		if err != nil {
			t.Fatalf("level %d: gzip read: %v", level, err)
		}
		if string(got) != textBuf.String() {
			t.Fatalf("level %d: -O z decompressed != -O t\n--- z ---\n%s\n--- t ---\n%s", level, got, textBuf.String())
		}
	}
}

// TestGtcheck_CrossCheckOrdering pins the lower-triangle pair order
// (query = larger-indexed sample): rows must be (S2,S1),(S3,S1),(S3,S2).
func TestGtcheck_CrossCheckOrdering(t *testing.T) {
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	rows := parseDCv2(out.String())
	want := [][2]string{{"S2", "S1"}, {"S3", "S1"}, {"S3", "S2"}}
	if len(rows) != 3 {
		t.Fatalf("got %d rows, want 3:\n%s", len(rows), out.String())
	}
	for i, w := range want {
		if rows[i].qry != w[0] || rows[i].gt != w[1] {
			t.Errorf("row %d: got (%s,%s), want (%s,%s)", i, rows[i].qry, rows[i].gt, w[0], w[1])
		}
	}
}

// TestGtcheck_IntegerDiscordance pins the -E 0 integer-mismatch path:
// discordance is the COUNT of discordant sites (a site is discordant
// when the dosage bitmasks share no bit), NOT the Hamming distance.
//
// Reference (upstream `bcftools gtcheck -u GT -E 0`):
//
//	DCv2 S2 S1 1 ... 3 2
//	DCv2 S3 S1 4 ... 4 0
//	DCv2 S3 S2 3 ... 3 0
func TestGtcheck_IntegerDiscordance(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{UseTag: "GT", ErrorProbabilityZero: true})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	want := map[[2]string]struct{ disc, sites, match int }{
		{"S2", "S1"}: {1, 3, 2},
		{"S3", "S1"}: {4, 4, 0},
		{"S3", "S2"}: {3, 3, 0},
	}
	for _, p := range r.Pairs {
		w, ok := want[[2]string{p.QuerySample, p.GenotypedSample}]
		if !ok {
			t.Fatalf("unexpected pair (%s,%s)", p.QuerySample, p.GenotypedSample)
		}
		if !p.IsInteger {
			t.Errorf("(%s,%s): expected integer scoring", p.QuerySample, p.GenotypedSample)
		}
		if p.DiscCount != w.disc || p.NumSites != w.sites || p.NumMatching != w.match {
			t.Errorf("(%s,%s): got disc=%d sites=%d match=%d, want disc=%d sites=%d match=%d",
				p.QuerySample, p.GenotypedSample, p.DiscCount, p.NumSites, p.NumMatching, w.disc, w.sites, w.match)
		}
	}
}

// TestGtcheck_ProbabilityDiscordance pins the default (-E 40) probability
// path. For pair (S2,S1) the only discordant site is chr1:400 (S2=0/1,
// S1=0/0), contributing -log(eprob) = 9.2103... ; matching sites add 0.
// The HWE average over the two matched sites is 0.8109302.
func TestGtcheck_ProbabilityDiscordance(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{UseTag: "GT"})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	var s2s1 *GtcheckPair
	for i := range r.Pairs {
		if r.Pairs[i].QuerySample == "S2" && r.Pairs[i].GenotypedSample == "S1" {
			s2s1 = &r.Pairs[i]
		}
	}
	if s2s1 == nil {
		t.Fatalf("pair (S2,S1) missing")
	}
	if s2s1.IsInteger {
		t.Fatalf("(S2,S1): expected probability scoring")
	}
	wantDisc := -math.Log(math.Pow(10, -0.1*40)) // 9.2103...
	if math.Abs(s2s1.DiscScore-wantDisc) > 1e-9 {
		t.Errorf("(S2,S1) disc score: got %v, want %v", s2s1.DiscScore, wantDisc)
	}
	if math.Abs(s2s1.AvgLogPHWE-0.8109302) > 1e-6 {
		t.Errorf("(S2,S1) HWE: got %v, want ~0.8109302", s2s1.AvgLogPHWE)
	}
	if s2s1.NumSites != 3 || s2s1.NumMatching != 2 {
		t.Errorf("(S2,S1): got sites=%d match=%d, want sites=3 match=2", s2s1.NumSites, s2s1.NumMatching)
	}
}

// TestGtcheck_PLMode exercises the -u PL path. With the fixture's PLs the
// minimum PL genotype is unambiguous everywhere except the missing S2
// site (PL=0,0,0 → ambiguous, but that site has GT ./. and PL is still
// read). The pair (S2,S1) compares all four sites under PL.
func TestGtcheck_PLMode(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{UseTag: "PL"})
	if err != nil {
		t.Fatalf("Gtcheck(-u PL): %v", err)
	}
	// All four sites contribute under PL (the S2 PL=0,0,0 at chr1:200 is
	// ambiguous-but-present), so (S2,S1) compares 4 sites, matching 3.
	var s2s1 *GtcheckPair
	for i := range r.Pairs {
		if r.Pairs[i].QuerySample == "S2" && r.Pairs[i].GenotypedSample == "S1" {
			s2s1 = &r.Pairs[i]
		}
	}
	if s2s1 == nil {
		t.Fatalf("pair (S2,S1) missing")
	}
	if s2s1.NumSites != 4 {
		t.Errorf("(S2,S1) PL sites: got %d, want 4 (GT mode would skip the ./. site)", s2s1.NumSites)
	}
	if s2s1.NumMatching != 3 {
		t.Errorf("(S2,S1) PL matches: got %d, want 3", s2s1.NumMatching)
	}
}

// TestGtcheck_NoHWEProb zeroes the HWE column.
func TestGtcheck_NoHWEProb(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{NoHWEProb: true})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	for _, p := range r.Pairs {
		if p.AvgLogPHWE != 0.0 {
			t.Errorf("--no-HWE-prob: pair (%s,%s) HWE=%v, want 0", p.QuerySample, p.GenotypedSample, p.AvgLogPHWE)
		}
		if p.NumMatching != 0 {
			t.Errorf("--no-HWE-prob: pair (%s,%s) nmatch=%d, want 0 (upstream prints 0)", p.QuerySample, p.GenotypedSample, p.NumMatching)
		}
	}
}

// TestGtcheck_HomsOnlyRequiresPanel mirrors upstream's rejection of
// --homs-only without -g.
func TestGtcheck_HomsOnlyRequiresPanel(t *testing.T) {
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{HomsOnly: true})
	if err == nil {
		t.Fatalf("expected error for --homs-only without -g")
	}
}

// TestGtcheck_Paired pins the -g panel mode: every query sample vs every
// panel sample, in query-major order. Here query==panel cohort so a 2x2
// fixture yields S1xS1, S1xS2, S2xS1, S2xS2.
func TestGtcheck_Paired(t *testing.T) {
	const panel = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	T	.	.	.	GT	0/0	1/1
chr1	300	.	G	A	.	.	.	GT	1/1	0/0
`
	const query = panel
	var out bytes.Buffer
	r, err := GtcheckPaired(strings.NewReader(query), strings.NewReader(panel), &out, GtcheckOptions{ErrorProbabilityZero: true})
	if err != nil {
		t.Fatalf("GtcheckPaired: %v", err)
	}
	if len(r.Pairs) != 4 {
		t.Fatalf("got %d pairs, want 4:\n%s", len(r.Pairs), out.String())
	}
	// Self pairs are concordant (disc 0); cross pairs differ at both sites.
	for _, p := range r.Pairs {
		if p.QuerySample == p.GenotypedSample {
			if p.DiscCount != 0 {
				t.Errorf("self pair (%s,%s) disc=%d, want 0", p.QuerySample, p.GenotypedSample, p.DiscCount)
			}
		} else {
			if p.DiscCount != 2 {
				t.Errorf("cross pair (%s,%s) disc=%d, want 2", p.QuerySample, p.GenotypedSample, p.DiscCount)
			}
		}
	}
}

// TestGtcheck_ExplicitPairs pins -p ordering: rows are emitted in the
// order given (qry,gt), not the cross-check triangle.
func TestGtcheck_ExplicitPairs(t *testing.T) {
	var out bytes.Buffer
	r, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{
		PairsSpec: "S1,S2,S1,S3",
	})
	if err != nil {
		t.Fatalf("Gtcheck(-p): %v", err)
	}
	if len(r.Pairs) != 2 {
		t.Fatalf("got %d pairs, want 2", len(r.Pairs))
	}
	if r.Pairs[0].QuerySample != "S1" || r.Pairs[0].GenotypedSample != "S2" {
		t.Errorf("pair 0: got (%s,%s), want (S1,S2)", r.Pairs[0].QuerySample, r.Pairs[0].GenotypedSample)
	}
	if r.Pairs[1].QuerySample != "S1" || r.Pairs[1].GenotypedSample != "S3" {
		t.Errorf("pair 1: got (%s,%s), want (S1,S3)", r.Pairs[1].QuerySample, r.Pairs[1].GenotypedSample)
	}
}

// TestGtcheck_DistinctiveSites pins the --distinctive-sites block.
func TestGtcheck_DistinctiveSites(t *testing.T) {
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{
		PairsSpec:           "S1,S2,S1,S3,S2,S3",
		DistinctiveSites:    1,
		HasDistinctiveSites: true,
	})
	if err != nil {
		t.Fatalf("Gtcheck(--distinctive-sites): %v", err)
	}
	// Upstream emits DS data rows but (a latent upstream quirk) NOT the
	// "#DS" comment header; we reproduce that.
	if strings.Contains(out.String(), "#DS\t[2]Chromosome") {
		t.Errorf("unexpected #DS comment header (upstream never writes it):\n%s", out.String())
	}
	// Site chr1:400 distinguishes the most pairs and must appear first.
	dsLines := []string{}
	for _, l := range strings.Split(out.String(), "\n") {
		if strings.HasPrefix(l, "DS\t") {
			dsLines = append(dsLines, l)
		}
	}
	if len(dsLines) == 0 {
		t.Fatalf("no DS data lines:\n%s", out.String())
	}
	if !strings.HasPrefix(dsLines[0], "DS\tchr1\t400\t") {
		t.Errorf("first DS line: got %q, want chr1:400 first", dsLines[0])
	}
}

// TestGtcheck_DistinctiveRequiresPairs mirrors upstream's requirement.
func TestGtcheck_DistinctiveRequiresPairs(t *testing.T) {
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{
		DistinctiveSites:    1,
		HasDistinctiveSites: true,
	})
	if err == nil {
		t.Fatalf("expected error: --distinctive-sites requires -p/-P")
	}
}

// TestGtcheck_RegionFilterIsolated locks in that a lone Regions filter
// gates the scored sites.
func TestGtcheck_RegionFilterIsolated(t *testing.T) {
	var outAll, outR bytes.Buffer
	rAll, err := Gtcheck(strings.NewReader(fixtureGtcheck), &outAll, GtcheckOptions{})
	if err != nil {
		t.Fatalf("Gtcheck(no regions): %v", err)
	}
	rR, err := Gtcheck(strings.NewReader(fixtureGtcheck), &outR, GtcheckOptions{
		Regions: []string{"chr1:150-250"},
	})
	if err != nil {
		t.Fatalf("Gtcheck(regions only): %v", err)
	}
	for _, p := range rAll.Pairs {
		if p.NumSites < 3 {
			t.Errorf("baseline: pair (%s,%s) scored %d sites, fixture has 3-4", p.QuerySample, p.GenotypedSample, p.NumSites)
		}
	}
	for _, p := range rR.Pairs {
		if p.NumSites > 1 {
			t.Errorf("regions filter: pair (%s,%s) scored %d sites, want <=1", p.QuerySample, p.GenotypedSample, p.NumSites)
		}
	}
}

// TestGtcheck_MultiAllelicRejected mirrors upstream's biallelic input
// requirement.
func TestGtcheck_MultiAllelicRejected(t *testing.T) {
	const multi = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	T,C	.	.	.	GT	0/0	1/2
`
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(multi), &out, GtcheckOptions{})
	if err == nil || !strings.Contains(err.Error(), "norm -m -") {
		t.Fatalf("expected multi-allelic rejection, got %v", err)
	}
}

// TestValidateUseTag pins the -u parser.
func TestValidateUseTag(t *testing.T) {
	for _, ok := range []string{"", "GT", "gt", "PL", "pl", "GT,PL", "PL,GT"} {
		if err := validateUseTag(ok); err != nil {
			t.Errorf("validateUseTag(%q) unexpected error: %v", ok, err)
		}
	}
	for _, bad := range []string{"XX", "GT,PL,GT", "foo"} {
		if err := validateUseTag(bad); err == nil {
			t.Errorf("validateUseTag(%q): expected error", bad)
		}
	}
}

// TestGtcheckCluster groups cross-check samples by pairwise discordance. The
// fixture has two replicate pairs (A1≈A2, B1≈B2) that differ between groups,
// so -c with a small max should yield exactly two clusters {A1,A2} and {B1,B2}.
func TestGtcheckCluster(t *testing.T) {
	const fixture = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A1	A2	B1	B2
chr1	10	.	A	T	.	.	.	GT	0/0	0/0	1/1	1/1
chr1	20	.	C	G	.	.	.	GT	0/1	0/1	0/0	0/0
chr1	30	.	G	A	.	.	.	GT	1/1	1/1	0/1	0/1
chr1	40	.	T	C	.	.	.	GT	0/0	0/0	1/1	1/1
chr1	50	.	A	G	.	.	.	GT	0/1	0/1	1/1	1/1
`
	var out bytes.Buffer
	_, err := Gtcheck(strings.NewReader(fixture), &out, GtcheckOptions{
		UseTag:               "GT",
		ErrorProbabilityZero: true, // integer mismatch scoring
		Cluster:              true,
		ClusterMin:           0.20,
		ClusterMax:           0.05,
	})
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	var clusters [][]string
	for _, line := range strings.Split(out.String(), "\n") {
		if !strings.HasPrefix(line, "CLUSTER\t") {
			continue
		}
		f := strings.Split(line, "\t")
		clusters = append(clusters, strings.Split(f[len(f)-1], ","))
	}
	if len(clusters) != 2 {
		t.Fatalf("got %d clusters, want 2:\n%s", len(clusters), out.String())
	}
	// Each cluster must be a single replicate pair (A* together, B* together).
	for _, c := range clusters {
		if len(c) != 2 {
			t.Fatalf("cluster %v has %d members, want 2", c, len(c))
		}
		if c[0][0] != c[1][0] { // same leading letter (A or B)
			t.Errorf("cluster %v mixes individuals", c)
		}
	}
}
