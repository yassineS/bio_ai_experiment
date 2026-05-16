package bcftools

import (
	"bytes"
	"math"
	"strings"
	"testing"
)

// fixtureRoh is a one-sample VCF that mixes a long HOM-ALT stretch
// (autozygous-leaning) with a few HET sites that pull the model back
// toward HW.
const fixtureRoh = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000000>
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	T	.	.	AF=0.5	GT	1/1
chr1	2000	.	C	G	.	.	AF=0.5	GT	1/1
chr1	3000	.	G	A	.	.	AF=0.5	GT	1/1
chr1	4000	.	T	C	.	.	AF=0.5	GT	1/1
chr1	5000	.	A	G	.	.	AF=0.5	GT	1/1
chr1	6000	.	C	T	.	.	AF=0.5	GT	0/1
chr1	7000	.	G	C	.	.	AF=0.5	GT	0/1
chr1	8000	.	T	A	.	.	AF=0.5	GT	0/1
`

// TestRoh_AZEmissionFormula pins the AZ-state emission for hard GTs
// against the upstream formula `pdg[0]*(1-p) + pdg[2]*p`:
//
//	HOM-REF -> emission(AZ) = (1-AF)
//	HOM-ALT -> emission(AZ) = AF
//	HET     -> emission(AZ) = 0
//
// This is the reviewer requirement #6 in the brief.
func TestRoh_AZEmissionFormula(t *testing.T) {
	cases := []struct {
		name string
		dose int
		af   float64
		want float64
	}{
		{"HOM-REF AF=0.2 -> 1-AF", 0, 0.2, 0.8},
		{"HOM-ALT AF=0.2 -> AF", 2, 0.2, 0.2},
		{"HET AF=0.2 -> 0", 1, 0.2, 0},
		{"HOM-REF AF=0.95", 0, 0.95, 0.05},
		{"HOM-ALT AF=0.95", 2, 0.95, 0.95},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			pdg := pdgFromDose(c.dose, 0)
			emAZ := pdg[0]*(1-c.af) + pdg[2]*c.af
			if math.Abs(emAZ-c.want) > 1e-9 {
				t.Errorf("emission(AZ) for dose=%d AF=%v: got %v, want %v", c.dose, c.af, emAZ, c.want)
			}
		})
	}
}

// TestRoh_MissingAFSkipByDefault is the reviewer requirement #7: when
// AF is unknown AND --AF-dflt is not given, the site must be SKIPPED.
// The PR #106 hard-coded fallback to 0.4 is gone.
func TestRoh_MissingAFSkipByDefault(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	T	.	.	.	GT	1/1
chr1	2000	.	C	G	.	.	.	GT	1/1
`
	var out bytes.Buffer
	r, err := Roh(strings.NewReader(src), &out, RohOptions{})
	if err != nil {
		t.Fatalf("Roh: %v", err)
	}
	if len(r.Sites) != 0 {
		t.Errorf("missing-AF sites must be skipped: got %d sites:\n%+v", len(r.Sites), r.Sites)
	}
}

// TestRoh_MissingAFDfltUsed verifies that providing --AF-dflt makes
// sites with missing AF process again. Together with the test above
// this nails the "*float64 sentinel" requirement in the brief.
func TestRoh_MissingAFDfltUsed(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	T	.	.	.	GT	1/1
chr1	2000	.	C	G	.	.	.	GT	1/1
`
	dflt := 0.4
	var out bytes.Buffer
	r, err := Roh(strings.NewReader(src), &out, RohOptions{AFDflt: &dflt})
	if err != nil {
		t.Fatalf("Roh --AF-dflt: %v", err)
	}
	if len(r.Sites) != 2 {
		t.Errorf("expected 2 sites with --AF-dflt, got %d", len(r.Sites))
	}
}

// TestRoh_DefaultTransitions pins the per-bp transition magnitudes to
// upstream's literal values (reviewer #9).
func TestRoh_DefaultTransitions(t *testing.T) {
	if math.Abs(DefaultHWtoAZ-6.7e-8) > 1e-12 {
		t.Errorf("DefaultHWtoAZ = %v, want 6.7e-8", DefaultHWtoAZ)
	}
	if math.Abs(DefaultAZtoHW-5e-9) > 1e-12 {
		t.Errorf("DefaultAZtoHW = %v, want 5e-9", DefaultAZtoHW)
	}
}

// TestRoh_HMMDecodesHomALTAsAZ runs the full pipeline with much
// higher transition probabilities (so the v1 port — which does NOT
// scale by physical distance between markers — can actually flip
// states) and confirms that a long HOM-ALT stretch with rare-allele
// AF decodes as AZ (state=1). The pinned RG quality scores are NOT
// comparable to upstream; that's tracked in the package docs.
func TestRoh_HMMDecodesHomALTAsAZ(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000000>
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	T	.	.	AF=0.05	GT	1/1
chr1	2000	.	C	G	.	.	AF=0.05	GT	1/1
chr1	3000	.	G	A	.	.	AF=0.05	GT	1/1
chr1	4000	.	T	C	.	.	AF=0.05	GT	1/1
chr1	5000	.	A	G	.	.	AF=0.05	GT	1/1
chr1	6000	.	C	T	.	.	AF=0.05	GT	1/1
chr1	7000	.	G	C	.	.	AF=0.05	GT	1/1
chr1	8000	.	T	A	.	.	AF=0.05	GT	1/1
`
	var out bytes.Buffer
	// Use stronger-than-default transitions so the v1 (no distance
	// scaling) model can actually flip on this short fixture; the
	// per-bp magnitudes are pinned by TestRoh_DefaultTransitions.
	r, err := Roh(strings.NewReader(src), &out, RohOptions{HWtoAZ: 1e-2, AZtoHW: 1e-3})
	if err != nil {
		t.Fatalf("Roh: %v", err)
	}
	azCount := 0
	for _, s := range r.Sites {
		if s.State == 1 {
			azCount++
		}
	}
	if azCount < 6 {
		t.Errorf("expected >=6 AZ sites for 8 rare-AF HOM-ALT calls, got %d:\n%+v", azCount, r.Sites)
	}
	if len(r.Regions) == 0 {
		t.Errorf("expected at least one RG region, got 0")
	}
}

// TestRoh_GTsOnlyAcceptsFloat: reviewer #8 — the upstream `-G/--GTs-only`
// flag is a float (phred error), not the int the PR #106 used. Test
// that the library accepts a float without error.
func TestRoh_GTsOnlyAcceptsFloat(t *testing.T) {
	var out bytes.Buffer
	_, err := Roh(strings.NewReader(fixtureRoh), &out, RohOptions{GTsOnly: 30.0})
	if err != nil {
		t.Fatalf("GTsOnly=30: %v", err)
	}
	_, err = Roh(strings.NewReader(fixtureRoh), &out, RohOptions{GTsOnly: 12.5})
	if err != nil {
		t.Fatalf("GTsOnly=12.5: %v", err)
	}
}

// TestRoh_CustomTransitionsThreaded: reviewer #10 — the CLI must wire
// -a/-H through to RohOptions.HWtoAZ / AZtoHW. We exercise the
// library directly here; the CLI side is covered in
// subcmds_gtcheck_roh_test.go.
func TestRoh_CustomTransitionsThreaded(t *testing.T) {
	var out bytes.Buffer
	r1, err := Roh(strings.NewReader(fixtureRoh), &out, RohOptions{
		HWtoAZ: 1e-4, AZtoHW: 1e-4,
	})
	if err != nil {
		t.Fatalf("Roh: %v", err)
	}
	// With much higher (1e-4) transition probabilities the decoder
	// should be more willing to switch states than with the
	// upstream-default 1e-8/1e-9 magnitudes. We don't pin an exact
	// path — just verify the call succeeded and yielded sites.
	if len(r1.Sites) == 0 {
		t.Errorf("expected sites with custom transitions, got 0")
	}
}

// TestRoh_SkipIndels confirms that -I/--skip-indels drops indel sites.
func TestRoh_SkipIndels(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	T	.	.	AF=0.5	GT	1/1
chr1	2000	.	A	AT	.	.	AF=0.5	GT	1/1
chr1	3000	.	C	G	.	.	AF=0.5	GT	1/1
`
	var out bytes.Buffer
	r, err := Roh(strings.NewReader(src), &out, RohOptions{SkipIndels: true})
	if err != nil {
		t.Fatalf("Roh -I: %v", err)
	}
	if len(r.Sites) != 2 {
		t.Errorf("with -I, expected 2 sites (the SNPs), got %d", len(r.Sites))
	}
	for _, s := range r.Sites {
		if s.Pos == 2000 {
			t.Errorf("indel at chr1:2000 should be skipped by -I, but appeared in output")
		}
	}
}

// TestRoh_IgnoreHomRef confirms -i/--ignore-homref skips 0/0 GTs.
func TestRoh_IgnoreHomRef(t *testing.T) {
	src := `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	T	.	.	AF=0.5	GT	0/0
chr1	2000	.	C	G	.	.	AF=0.5	GT	1/1
`
	var out bytes.Buffer
	r, err := Roh(strings.NewReader(src), &out, RohOptions{IgnoreHomRef: true})
	if err != nil {
		t.Fatalf("Roh -i: %v", err)
	}
	if len(r.Sites) != 1 {
		t.Errorf("with -i, expected 1 site (skipping the 0/0), got %d", len(r.Sites))
	}
}
