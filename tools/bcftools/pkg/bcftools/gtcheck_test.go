package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// gtcheckQueryFixture: 3 samples (Q1, Q2, Q3), 4 autosomal sites.
// Genotypes are chosen so that Q1 ≡ P1 exactly, Q2 differs at one
// site from P2, and Q3 is uncalled at one site.
func gtcheckQueryFixture() string {
	return `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000>
##INFO=<ID=AF,Number=1,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Q1	Q2	Q3
chr1	10	.	A	T	.	PASS	AF=0.1	GT	0/0	0/1	0/0
chr1	20	.	G	C	.	PASS	AF=0.2	GT	0/1	1/1	0/1
chr1	30	.	C	A	.	PASS	AF=0.3	GT	1/1	0/0	./.
chr1	40	.	T	G	.	PASS	AF=0.4	GT	0/1	0/1	0/1
`
}

func gtcheckPanelFixture() string {
	return `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000>
##INFO=<ID=AF,Number=1,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	P1	P2	P3
chr1	10	.	A	T	.	PASS	AF=0.1	GT	0/0	0/1	0/0
chr1	20	.	G	C	.	PASS	AF=0.2	GT	0/1	0/1	0/1
chr1	30	.	C	A	.	PASS	AF=0.3	GT	1/1	0/0	1/1
chr1	40	.	T	G	.	PASS	AF=0.4	GT	0/1	0/1	0/1
`
}

func TestParseGtcheckUseMode(t *testing.T) {
	cases := []struct {
		in   string
		want GtcheckUseMode
		err  bool
	}{
		{"", GtcheckUseGT, false},
		{"GT", GtcheckUseGT, false},
		{"gt", GtcheckUseGT, false},
		{"PL", GtcheckUsePL, false},
		{"GL", GtcheckUseGL, false},
		{"bogus", 0, true},
	}
	for _, c := range cases {
		got, err := ParseGtcheckUseMode(c.in)
		if c.err {
			if err == nil {
				t.Errorf("ParseGtcheckUseMode(%q) expected error", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseGtcheckUseMode(%q): %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseGtcheckUseMode(%q) = %v want %v", c.in, got, c.want)
		}
	}
}

func TestParseGtcheckPairsSpec(t *testing.T) {
	// Colon-separated form.
	got, err := ParseGtcheckPairsSpec("Q1:P1,Q2:P2")
	if err != nil {
		t.Fatalf("ParseGtcheckPairsSpec colon: %v", err)
	}
	want := []GtcheckPair{{"Q1", "P1"}, {"Q2", "P2"}}
	if len(got) != len(want) {
		t.Fatalf("got %d pairs want %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("pair %d: got %v want %v", i, got[i], want[i])
		}
	}

	// Flat alternating form.
	got, err = ParseGtcheckPairsSpec("Q1,P1,Q2,P2")
	if err != nil {
		t.Fatalf("ParseGtcheckPairsSpec flat: %v", err)
	}
	if len(got) != 2 || got[0] != (GtcheckPair{"Q1", "P1"}) || got[1] != (GtcheckPair{"Q2", "P2"}) {
		t.Errorf("flat parse: %+v", got)
	}

	// Odd-length flat form should error.
	if _, err := ParseGtcheckPairsSpec("Q1,P1,Q2"); err == nil {
		t.Error("expected error on odd-length flat spec")
	}

	// Empty colon side errors.
	if _, err := ParseGtcheckPairsSpec("Q1:"); err == nil {
		t.Error("expected error on empty panel side")
	}

	// Empty spec returns nil, no error.
	if pairs, err := ParseGtcheckPairsSpec(""); err != nil || pairs != nil {
		t.Errorf("empty spec: %+v %v", pairs, err)
	}
}

func TestGtcheck_HammingScoringExactValues(t *testing.T) {
	out := &bytes.Buffer{}
	res, err := Gtcheck(
		strings.NewReader(gtcheckQueryFixture()),
		strings.NewReader(gtcheckPanelFixture()),
		out,
		GtcheckOptions{
			Pairs: []GtcheckPair{
				{Query: "Q1", Panel: "P1"},
				{Query: "Q2", Panel: "P2"},
				{Query: "Q3", Panel: "P3"},
			},
		},
	)
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	// Expected per-pair behaviour:
	//   Q1 vs P1: identical at all 4 sites → score 0 / 4.
	//   Q2 vs P2: differ at chr1:20 (Q2=1/1, P2=0/1) → score 1 / 4.
	//   Q3 vs P3: differ at chr1:20 (Q3=0/1, P3=0/1 same!), wait recompute.
	//     site 10: Q3=0/0 P3=0/0 same.
	//     site 20: Q3=0/1 P3=0/1 same.
	//     site 30: Q3=./. (missing) → counted in NMissing.
	//     site 40: Q3=0/1 P3=0/1 same.
	//   So Q3 vs P3: score 0, NSites 3, NMissing 1.
	if len(res.Pairs) != 3 {
		t.Fatalf("Pairs len = %d want 3", len(res.Pairs))
	}
	if res.Pairs[0].Score != 0 || res.Pairs[0].NSites != 4 {
		t.Errorf("Q1/P1 = %+v want score=0 nsites=4", res.Pairs[0])
	}
	if res.Pairs[1].Score != 1 || res.Pairs[1].NSites != 4 {
		t.Errorf("Q2/P2 = %+v want score=1 nsites=4", res.Pairs[1])
	}
	if res.Pairs[2].Score != 0 || res.Pairs[2].NSites != 3 || res.Pairs[2].NMissing != 1 {
		t.Errorf("Q3/P3 = %+v want score=0 nsites=3 nmissing=1", res.Pairs[2])
	}

	// And the TSV should contain a "DC" row for each pair.
	body := out.String()
	for _, want := range []string{
		"# DC\tquery\tpanel\tscore\tn_sites\tn_missing\tdiscordance",
		"DC\tQ1\tP1\t0\t4\t0\t0",
		"DC\tQ2\tP2\t1\t4\t0\t",
		"DC\tQ3\tP3\t0\t3\t1\t0",
		"# totals\tn_sites_compared=4",
	} {
		if !strings.Contains(body, want) {
			t.Errorf("output missing %q\nfull:\n%s", want, body)
		}
	}
}

func TestGtcheck_CrossJoinWhenNoPairs(t *testing.T) {
	out := &bytes.Buffer{}
	res, err := Gtcheck(
		strings.NewReader(gtcheckQueryFixture()),
		strings.NewReader(gtcheckPanelFixture()),
		out,
		GtcheckOptions{}, // no Pairs → cross-join
	)
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	if len(res.Pairs) != 3*3 {
		t.Errorf("cross-join expected 9 pairs, got %d", len(res.Pairs))
	}
}

func TestGtcheck_SwapIsSymmetric(t *testing.T) {
	// Invariant: the score for (Q,P) computed against the panel must
	// equal the score for (P,Q) computed against the *swapped* files
	// (i.e. panel-as-query, query-as-panel). The metric is the
	// unordered-allele Hamming distance, so it's order-independent.
	pair := []GtcheckPair{{Query: "Q1", Panel: "P1"}}
	resQ, err := Gtcheck(
		strings.NewReader(gtcheckQueryFixture()),
		strings.NewReader(gtcheckPanelFixture()),
		&bytes.Buffer{},
		GtcheckOptions{Pairs: pair},
	)
	if err != nil {
		t.Fatalf("Gtcheck forward: %v", err)
	}
	pairSwapped := []GtcheckPair{{Query: "P1", Panel: "Q1"}}
	resP, err := Gtcheck(
		strings.NewReader(gtcheckPanelFixture()),
		strings.NewReader(gtcheckQueryFixture()),
		&bytes.Buffer{},
		GtcheckOptions{Pairs: pairSwapped},
	)
	if err != nil {
		t.Fatalf("Gtcheck swapped: %v", err)
	}
	if resQ.Pairs[0].Score != resP.Pairs[0].Score ||
		resQ.Pairs[0].NSites != resP.Pairs[0].NSites ||
		resQ.Pairs[0].NMissing != resP.Pairs[0].NMissing {
		t.Errorf("symmetric invariant broken: forward=%+v swapped=%+v",
			resQ.Pairs[0], resP.Pairs[0])
	}
}

func TestGtcheck_UnknownSampleError(t *testing.T) {
	_, err := Gtcheck(
		strings.NewReader(gtcheckQueryFixture()),
		strings.NewReader(gtcheckPanelFixture()),
		&bytes.Buffer{},
		GtcheckOptions{
			Pairs: []GtcheckPair{{Query: "NOPE", Panel: "P1"}},
		},
	)
	if err == nil {
		t.Fatal("expected error on missing query sample")
	}
}

func TestGtcheck_PLRejectedWithRoadmap(t *testing.T) {
	tmp := t.TempDir()
	queryPath := filepath.Join(tmp, "q.vcf")
	if err := os.WriteFile(queryPath, []byte(gtcheckQueryFixture()), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	panelPath := filepath.Join(tmp, "p.vcf")
	if err := os.WriteFile(panelPath, []byte(gtcheckPanelFixture()), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	_, err := GtcheckFile(queryPath, &bytes.Buffer{}, GtcheckOptions{
		PanelPath: panelPath,
		Use:       GtcheckUsePL,
	})
	if err == nil {
		t.Fatal("expected error on -u PL")
	}
	if !strings.Contains(err.Error(), "PARITY_ROADMAP") {
		t.Errorf("error should point at roadmap: %v", err)
	}
}

func TestGtcheck_PairsFileLoad(t *testing.T) {
	tmp := t.TempDir()
	queryPath := filepath.Join(tmp, "q.vcf")
	panelPath := filepath.Join(tmp, "p.vcf")
	pairsPath := filepath.Join(tmp, "pairs.txt")
	for _, w := range []struct{ p, c string }{
		{queryPath, gtcheckQueryFixture()},
		{panelPath, gtcheckPanelFixture()},
		{pairsPath, "# header\nQ1\tP1\nQ2,P2\n"},
	} {
		if err := os.WriteFile(w.p, []byte(w.c), 0o644); err != nil {
			t.Fatalf("write %s: %v", w.p, err)
		}
	}
	res, err := GtcheckFile(queryPath, &bytes.Buffer{}, GtcheckOptions{
		PanelPath: panelPath,
		PairsFile: pairsPath,
	})
	if err != nil {
		t.Fatalf("GtcheckFile: %v", err)
	}
	if len(res.Pairs) != 2 {
		t.Errorf("expected 2 pairs, got %d", len(res.Pairs))
	}
	// Q1 / P1 should be a perfect match (0 score, 4 sites).
	if res.Pairs[0].Score != 0 || res.Pairs[0].NSites != 4 {
		t.Errorf("Q1/P1 = %+v", res.Pairs[0])
	}
}

func TestGtcheck_RegionFilter(t *testing.T) {
	// chr1:5-15 keeps only the first record.
	out := &bytes.Buffer{}
	res, err := Gtcheck(
		strings.NewReader(gtcheckQueryFixture()),
		strings.NewReader(gtcheckPanelFixture()),
		out,
		GtcheckOptions{
			Pairs:   []GtcheckPair{{Query: "Q2", Panel: "P2"}},
			Regions: []string{"chr1:5-15"},
		},
	)
	if err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	if res.NSitesCompared != 1 {
		t.Errorf("expected 1 site compared, got %d", res.NSitesCompared)
	}
	// In this region Q2 = 0/1, P2 = 0/1 → 0 mismatches.
	if res.Pairs[0].Score != 0 || res.Pairs[0].NSites != 1 {
		t.Errorf("expected 0/1 score, got %+v", res.Pairs[0])
	}
}

func TestParseHardGT(t *testing.T) {
	cases := []struct {
		in     string
		want   [2]int
		wantOK bool
	}{
		{"0/0", [2]int{0, 0}, true},
		{"0|1", [2]int{0, 1}, true},
		{"1/0", [2]int{0, 1}, true}, // canonicalised
		{"1/1", [2]int{1, 1}, true},
		{"./.", [2]int{}, false},
		{".", [2]int{}, false},
		{"", [2]int{}, false},
		{"0", [2]int{}, false},
		{"0/1/2", [2]int{}, false},
	}
	for _, c := range cases {
		got, ok := parseHardGT(c.in)
		if ok != c.wantOK {
			t.Errorf("parseHardGT(%q) ok = %v want %v", c.in, ok, c.wantOK)
			continue
		}
		if ok && got != c.want {
			t.Errorf("parseHardGT(%q) = %v want %v", c.in, got, c.want)
		}
	}
}

func TestSortedGtcheckPairs(t *testing.T) {
	r := GtcheckResult{
		Pairs: []GtcheckPairResult{
			{Query: "B", Panel: "X", Discordance: 0.5},
			{Query: "A", Panel: "X", Discordance: 0.1},
			{Query: "C", Panel: "X", Discordance: 0.1},
		},
	}
	sorted := SortedGtcheckPairs(r)
	if sorted[0].Query != "A" || sorted[1].Query != "C" || sorted[2].Query != "B" {
		t.Errorf("sort order: %+v", sorted)
	}
}
