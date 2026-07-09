package bcftools

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// TestGtcheckProjectionSurvivesBufferReuse pins the projection RSS fix: reading
// via vcf.ReadInto reuses a single scratch record whose string fields alias the
// reader's line buffer, so projectVariant must copy CHROM/REF/ALT. This test
// reads the fixture through readGtSites and asserts every projected site's
// CHROM/REF/ALT matches an independent, fully-owned parse — a regression guard
// for the aliasing bug where a later read clobbered an earlier site's ALT and
// mis-flagged it multiallelic.
func TestGtcheckProjectionSurvivesBufferReuse(t *testing.T) {
	_, sites, dropped, err := readGtSites(strings.NewReader(fixtureGtcheck), nil, nil)
	if err != nil {
		t.Fatalf("readGtSites: %v", err)
	}
	if len(dropped) != 0 {
		t.Fatalf("no filter given, expected 0 dropped, got %d", len(dropped))
	}

	// Independent, owned parse of the same fixture for the reference key columns.
	r := vcf.NewReader(strings.NewReader(fixtureGtcheck))
	if _, herr := r.ReadHeader(); herr != nil {
		t.Fatalf("ReadHeader: %v", herr)
	}
	var want []*vcf.Variant
	for {
		v, rerr := r.Read()
		if rerr != nil {
			break
		}
		want = append(want, v)
	}
	if len(sites) != len(want) {
		t.Fatalf("site count: projection=%d owned=%d", len(sites), len(want))
	}
	for i := range want {
		if sites[i].chrom != want[i].Chrom {
			t.Errorf("site %d chrom: %q vs %q", i, sites[i].chrom, want[i].Chrom)
		}
		if sites[i].pos != want[i].Pos {
			t.Errorf("site %d pos: %d vs %d", i, sites[i].pos, want[i].Pos)
		}
		if sites[i].ref != want[i].Ref {
			t.Errorf("site %d ref: %q vs %q", i, sites[i].ref, want[i].Ref)
		}
		if strings.Join(sites[i].alt, ",") != strings.Join(want[i].Alt, ",") {
			t.Errorf("site %d alt: %q vs %q", i, sites[i].alt, want[i].Alt)
		}
		if isMultiAllelic(&sites[i]) {
			t.Errorf("site %d unexpectedly multiallelic (alt=%q)", i, sites[i].alt)
		}
	}
}

// TestGtcheckProjectionTagPresence pins the record-level tag-presence semantics
// the projection must preserve: a present-but-unparseable PL still counts as
// "PL tag present" (so recordTag does not drop the record to no-data), while a
// GT-only record has no PL key. This mirrors the upstream recordHasTag behaviour
// that the projection replaced with the anyGT / anyPLkey flags.
func TestGtcheckProjectionTagPresence(t *testing.T) {
	const fixture = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	T	.	.	.	GT:PL	0/0:0,30,255	0/1:.
chr1	200	.	C	G	.	.	.	GT	0/1	1/1
`
	_, sites, _, err := readGtSites(strings.NewReader(fixture), nil, nil)
	if err != nil {
		t.Fatalf("readGtSites: %v", err)
	}
	if len(sites) != 2 {
		t.Fatalf("expected 2 sites, got %d", len(sites))
	}
	// Site 0 carries both GT and a PL key (the S2 PL "." is present but not a
	// parseable diploid triple).
	if !sites[0].anyGT || !sites[0].anyPLkey {
		t.Errorf("site 0: anyGT=%v anyPLkey=%v, want both true", sites[0].anyGT, sites[0].anyPLkey)
	}
	if !sites[0].samples[1].hasPLkey || sites[0].samples[1].hasPL {
		t.Errorf("site 0 S2: hasPLkey=%v hasPL=%v, want key present but unparsed",
			sites[0].samples[1].hasPLkey, sites[0].samples[1].hasPL)
	}
	// Site 1 is GT-only: no PL key anywhere.
	if !sites[1].anyGT || sites[1].anyPLkey {
		t.Errorf("site 1: anyGT=%v anyPLkey=%v, want GT-only", sites[1].anyGT, sites[1].anyPLkey)
	}
}

// TestGtcheckProjectionAFFromInfo verifies the projected INFO/AC,AN reproduce
// siteAF's arithmetic (AF = AC/AN) exactly, so the HWE column is unchanged by
// the projection.
func TestGtcheckProjectionAFFromInfo(t *testing.T) {
	const fixture = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">
##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	T	.	.	AC=3;AN=4	GT	0/1	1/1
`
	_, sites, _, err := readGtSites(strings.NewReader(fixture), nil, nil)
	if err != nil {
		t.Fatalf("readGtSites: %v", err)
	}
	if len(sites) != 1 {
		t.Fatalf("expected 1 site, got %d", len(sites))
	}
	got := siteAF(&sites[0])
	want := 3.0 / 4.0
	if got != want {
		t.Errorf("siteAF from INFO AC=3,AN=4: got %v, want %v", got, want)
	}
}

// TestGtcheckProjectionOutputUnchanged confirms the whole projected scoring path
// still emits the same #DCv2 table it did before (a self-contained regression
// pin, independent of any upstream binary).
func TestGtcheckProjectionOutputUnchanged(t *testing.T) {
	var out bytes.Buffer
	if _, err := Gtcheck(strings.NewReader(fixtureGtcheck), &out, GtcheckOptions{}); err != nil {
		t.Fatalf("Gtcheck: %v", err)
	}
	rows := parseDCv2(out.String())
	// 3 samples cross-checked → lower triangle of 3 pairs.
	if len(rows) != 3 {
		t.Fatalf("expected 3 DCv2 rows, got %d:\n%s", len(rows), out.String())
	}
	if !strings.Contains(out.String(), "INFO\tsites-compared\t4") {
		t.Errorf("projected output missing sites-compared:\n%s", out.String())
	}
}
