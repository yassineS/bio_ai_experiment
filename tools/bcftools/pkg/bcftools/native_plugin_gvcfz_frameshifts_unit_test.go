package bcftools

import (
	"reflect"
	"strconv"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Binary-free unit tests for the pure helpers of the gvcfz and frameshifts
// native ports. These pass with the upstream submodules UNPOPULATED (no
// exec.Command, no BCFTOOLS_PLUGINS): they exercise the grouping decision, the
// per-allele variant classification (bcf_set_variant_type), the exon cursor
// (bcf_sr_regions_overlap), and the OOF computation directly.

// TestUnitAlleleVariant checks the per-allele classification (signed length n and
// the VCF_INDEL bit) against the htslib bcf_set_variant_type values verified with
// a C probe against the vendored htslib.
func TestUnitAlleleVariant(t *testing.T) {
	cases := []struct {
		ref, alt string
		n        int
		indel    bool
	}{
		{"A", "ATTT", 3, true},       // 3-base insertion
		{"ATTT", "A", -3, true},      // 3-base deletion
		{"A", "AT", 1, true},         // 1-base insertion
		{"AT", "A", -1, true},        // 1-base deletion
		{"A", "ATTTTT", 5, true},     // 5-base insertion
		{"GGGGGG", "G", -5, true},    // 5-base deletion
		{"C", "CGG", 2, true},        // 2-base insertion
		{"ACGT", "A", -3, true},      // 3-base deletion
		{"AT", "ATG", 1, true},       // mid insertion
		{"ATG", "AT", -1, true},      // mid deletion
		{"a", "at", 1, true},         // lower-case prefix still matches (toupper)
		{"A", "G", 1, false},         // SNP
		{"AC", "GT", 2, false},       // MNP (not an indel)
		{"A", "<NON_REF>", 0, false}, // symbolic NON_REF => REF
		{"A", "*", 0, false},         // overlap allele
		{"A", ".", 0, false},         // missing ALT => REF
	}
	for _, c := range cases {
		av := alleleVariant(c.ref, c.alt)
		if av.n != c.n || av.isIndel != c.indel {
			t.Errorf("alleleVariant(%q,%q) = {n=%d indel=%v}, want {n=%d indel=%v}",
				c.ref, c.alt, av.n, av.isIndel, c.n, c.indel)
		}
	}
}

// TestUnitOofForAllele covers BOTH the shipped upstream behaviour (dead-code
// guard => always -1) and the corrected --fix-oof exon-trim + mod-3 computation.
// Exon coordinates are 0-based, inclusive end (exStart..exEnd), matching the
// reg->start/reg->end the cursor exposes.
func TestUnitOofForAllele(t *testing.T) {
	// Exon spanning 0-based 100..199 (BED `chr1 100 200`).
	const exStart, exEnd = 100, 199

	t.Run("upstream-shipped-always-minus1", func(t *testing.T) {
		// Every indel allele overlapping the exon yields -1 in the real binary.
		cases := []struct {
			ref, alt string
			pos0     int
		}{
			{"A", "ATTT", 149}, // insertion, in-frame if fixed
			{"ATTT", "A", 154}, // deletion 3, in-frame if fixed
			{"A", "AT", 159},   // insertion 1, OOF if fixed
			{"AT", "A", 164},   // deletion 1, OOF if fixed
		}
		for _, c := range cases {
			if got := oofForAllele(c.ref, c.alt, c.pos0, exStart, exEnd, false); got != -1 {
				t.Errorf("oofForAllele(%q,%q,pos0=%d,fix=false) = %d, want -1", c.ref, c.alt, c.pos0, got)
			}
		}
	})

	t.Run("fixed-computation", func(t *testing.T) {
		cases := []struct {
			name     string
			ref, alt string
			pos0     int
			want     int
		}{
			{"insertion 3 in exon => in-frame", "A", "ATTT", 149, 0},
			{"insertion 1 in exon => out-of-frame", "A", "AT", 159, 1},
			{"insertion 2 in exon => out-of-frame", "C", "CGG", 149, 1},
			{"deletion 3 fully in exon => in-frame", "ATTT", "A", 154, 0},
			{"deletion 1 fully in exon => out-of-frame", "AT", "A", 164, 1},
			// Deletion partially past the exon end: pos0=197, n=-5 so the deletion
			// spans 197..202 but the exon ends at 199; tlen=5 trimmed by
			// (202-199)=3 => 2 inside the exon => 2%3 => out-of-frame.
			{"deletion 5 trimmed at end", "GGGGGG", "G", 197, 1},
			// SNP allele (not an indel) => not applicable.
			{"snp => na", "A", "G", 150, -1},
			// Insertion whose anchor base is at the exon end boundary: the
			// insertion branch requires exEnd > pos0, and 199 > 199 is false => -1.
			{"insertion at exon end boundary", "A", "AT", 199, -1},
		}
		for _, c := range cases {
			if got := oofForAllele(c.ref, c.alt, c.pos0, exStart, exEnd, true); got != c.want {
				t.Errorf("%s: oofForAllele(%q,%q,pos0=%d,fix=true) = %d, want %d",
					c.name, c.ref, c.alt, c.pos0, got, c.want)
			}
		}
	})
}

// TestUnitExonCursor exercises the bcf_sr_regions_overlap port: forward
// advancement, the start<=end overlap test, sort+merge, re-seek on a backwards
// query, and the absent-sequence/no-overlap branches.
func TestUnitExonCursor(t *testing.T) {
	newCursor := func() *exonCursor {
		c := &exonCursor{byChrom: map[string][]exonRegion{}, iseqOf: map[string]int{}}
		// BED `chr1 100 200`, `chr1 300 400`, `chr2 100 200` => from++ then -1.
		c.add("chr1", 100, 199)
		c.add("chr1", 300, 399)
		c.add("chr2", 100, 199)
		c.sortAndMerge()
		c.reset()
		return c
	}

	t.Run("forward overlaps and gaps", func(t *testing.T) {
		c := newCursor()
		// Inside first exon.
		if !c.overlap("chr1", 149, 149) {
			t.Fatalf("expected overlap at chr1:149")
		}
		if c.start != 100 || c.end != 199 {
			t.Fatalf("cursor at first exon = [%d,%d], want [100,199]", c.start, c.end)
		}
		// Between the two chr1 exons (0-based 200..299): no overlap.
		if c.overlap("chr1", 250, 250) {
			t.Fatalf("did not expect overlap at chr1:250")
		}
		// Inside second exon.
		if !c.overlap("chr1", 350, 350) {
			t.Fatalf("expected overlap at chr1:350")
		}
		if c.start != 300 || c.end != 399 {
			t.Fatalf("cursor at second exon = [%d,%d], want [300,399]", c.start, c.end)
		}
		// New chromosome.
		if !c.overlap("chr2", 150, 150) {
			t.Fatalf("expected overlap at chr2:150")
		}
		if c.start != 100 || c.end != 199 {
			t.Fatalf("cursor at chr2 exon = [%d,%d], want [100,199]", c.start, c.end)
		}
	})

	t.Run("boundaries", func(t *testing.T) {
		c := newCursor()
		if !c.overlap("chr1", 100, 100) { // first base of exon
			t.Fatalf("expected overlap at exon start")
		}
		c = newCursor()
		if !c.overlap("chr1", 199, 199) { // last base of exon
			t.Fatalf("expected overlap at exon end")
		}
		c = newCursor()
		if c.overlap("chr1", 99, 99) { // one before exon start
			t.Fatalf("did not expect overlap one before exon start")
		}
		c = newCursor()
		// A query that spans into the exon overlaps (end >= exon start).
		if !c.overlap("chr1", 95, 105) {
			t.Fatalf("expected overlap for span crossing the exon start")
		}
	})

	t.Run("backwards query re-seeks", func(t *testing.T) {
		c := newCursor()
		if !c.overlap("chr1", 350, 350) { // advance to second exon
			t.Fatalf("expected overlap at chr1:350")
		}
		// A backwards query (prev_start > start) re-seeks to the chromosome start.
		if !c.overlap("chr1", 150, 150) {
			t.Fatalf("expected overlap after backwards re-seek to chr1:150")
		}
		if c.start != 100 || c.end != 199 {
			t.Fatalf("cursor after re-seek = [%d,%d], want [100,199]", c.start, c.end)
		}
	})

	t.Run("absent sequence", func(t *testing.T) {
		c := newCursor()
		if c.overlap("chrX", 100, 100) {
			t.Fatalf("did not expect overlap on an absent sequence")
		}
	})

	t.Run("sort and merge", func(t *testing.T) {
		c := &exonCursor{byChrom: map[string][]exonRegion{}, iseqOf: map[string]int{}}
		// Out-of-order, overlapping/adjacent regions collapse to one.
		c.add("chr1", 300, 399)
		c.add("chr1", 100, 199)
		c.add("chr1", 150, 250) // overlaps [100,199]
		c.sortAndMerge()
		got := c.byChrom["chr1"]
		want := []exonRegion{{100, 250}, {300, 399}}
		if !reflect.DeepEqual(got, want) {
			t.Fatalf("sortAndMerge => %+v, want %+v", got, want)
		}
	})
}

// TestUnitIsGvcfRefBlock checks the gVCF reference-block predicate.
func TestUnitIsGvcfRefBlock(t *testing.T) {
	cases := []struct {
		alt  []string
		want bool
	}{
		{nil, true},
		{[]string{"<NON_REF>"}, true},
		{[]string{"<*>"}, true},
		{[]string{"."}, true},
		{[]string{"A"}, false},
		{[]string{"<NON_REF>", "A"}, false},
		{[]string{"A", "<NON_REF>"}, false},
	}
	for _, c := range cases {
		v := &vcf.Variant{Alt: c.alt}
		if got := isGvcfRefBlock(v); got != c.want {
			t.Errorf("isGvcfRefBlock(alt=%v) = %v, want %v", c.alt, got, c.want)
		}
	}
}

// TestUnitParseGvcfzGroups checks the -g group-by parser: filter labels, the
// catch-all "-" expression, the quote substitution in the FILTER description, and
// PASS contributing no FILTER line.
func TestUnitParseGvcfzGroups(t *testing.T) {
	hdr := &vcf.Header{
		Samples: []string{"S1"},
		MetaInfo: []string{
			`##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">`,
			`##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">`,
		},
	}
	groupBy := `PASS:GQ>60 & DP<20; Flt1:GQ>20; Flt2:-`
	groups, lines, err := parseGvcfzGroups(groupBy, hdr)
	if err != nil {
		t.Fatalf("parseGvcfzGroups: %v", err)
	}
	if len(groups) != 3 {
		t.Fatalf("got %d groups, want 3", len(groups))
	}
	if groups[0].filterLabel != "" {
		t.Errorf("PASS group filterLabel = %q, want empty", groups[0].filterLabel)
	}
	if groups[0].filter == nil {
		t.Errorf("PASS group expression should compile to a filter")
	}
	if groups[1].filterLabel != "Flt1" || groups[1].filter == nil {
		t.Errorf("Flt1 group = %+v", groups[1])
	}
	if groups[2].filterLabel != "Flt2" || groups[2].filter != nil {
		t.Errorf("Flt2 catch-all should have a nil filter; got %+v", groups[2])
	}
	// Two non-PASS FILTER lines, description = whole group-by string.
	if len(lines) != 2 {
		t.Fatalf("got %d FILTER lines, want 2: %v", len(lines), lines)
	}
	wantDescr := `##FILTER=<ID=Flt1,Description="PASS:GQ>60 & DP<20; Flt1:GQ>20; Flt2:-">`
	if lines[0] != wantDescr {
		t.Errorf("FILTER line[0] = %q, want %q", lines[0], wantDescr)
	}

	// Quote substitution: a double quote in the group-by becomes a single quote.
	groupsQ, linesQ, err := parseGvcfzGroups(`Flt1:GT!="alt"; PASS:GQ>1`, hdr)
	if err != nil {
		t.Fatalf("parseGvcfzGroups(quotes): %v", err)
	}
	if len(groupsQ) != 2 {
		t.Fatalf("got %d groups, want 2", len(groupsQ))
	}
	wantQ := `##FILTER=<ID=Flt1,Description="Flt1:GT!='alt'; PASS:GQ>1">`
	if len(linesQ) != 1 || linesQ[0] != wantQ {
		t.Errorf("quoted FILTER line = %v, want [%q]", linesQ, wantQ)
	}
}

// TestUnitGvcfzGrouping drives the whole block state machine over a small stream
// with NO upstream binary and checks the merged representative records: block
// boundaries, END extension/clamping, the min DP/GQ/PL merges, and a real ALT
// flushing the block.
func TestUnitGvcfzGrouping(t *testing.T) {
	hdr := &vcf.Header{
		Samples: []string{"S1"},
		MetaInfo: []string{
			`##INFO=<ID=END,Number=1,Type=Integer,Description="end">`,
			`##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">`,
			`##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">`,
			`##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">`,
			`##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="min_dp">`,
			`##FORMAT=<ID=PL,Number=G,Type=Integer,Description="pl">`,
		},
	}
	mkRef := func(pos, end, dp, gq, minDP, pl0, pl1, pl2 int) *vcf.Variant {
		return &vcf.Variant{
			Chrom: "chr1", Pos: pos, Ref: "A", Alt: []string{"<NON_REF>"}, Qual: -1,
			Info:      map[string]string{"END": strconv.Itoa(end)},
			InfoOrder: []string{"END"},
			Format:    []string{"GT", "DP", "GQ", "MIN_DP", "PL"},
			Samples: []vcf.Sample{{Name: "S1", Data: map[string]string{
				"GT": "0/0", "DP": strconv.Itoa(dp), "GQ": strconv.Itoa(gq), "MIN_DP": strconv.Itoa(minDP),
				"PL": strconv.Itoa(pl0) + "," + strconv.Itoa(pl1) + "," + strconv.Itoa(pl2),
			}}},
		}
	}

	p := &gvcfzPlugin{}
	out, err := p.Init([]string{"-g", "PASS:-"}, hdr)
	if err != nil {
		t.Fatalf("Init: %v", err)
	}
	_ = out

	variants := []*vcf.Variant{
		mkRef(100, 110, 30, 70, 25, 0, 30, 300),
		mkRef(111, 120, 18, 45, 15, 0, 20, 200),
		// A real variant flushes the block.
		{Chrom: "chr1", Pos: 121, Ref: "T", Alt: []string{"A"}, Qual: 50,
			Filter: []string{"PASS"}, Format: []string{"GT", "DP", "GQ"},
			Samples: []vcf.Sample{{Name: "S1", Data: map[string]string{"GT": "0/1", "DP": "40", "GQ": "99"}}}},
		mkRef(122, 130, 35, 65, 30, 0, 33, 330),
	}
	res, err := p.ProcessAll(variants)
	if err != nil {
		t.Fatalf("ProcessAll: %v", err)
	}
	// Expect: merged [100..120] block, the variant, the [122..130] block.
	if len(res) != 3 {
		t.Fatalf("got %d output records, want 3", len(res))
	}
	blk := res[0]
	if blk.Pos != 100 || blk.Info["END"] != "120" {
		t.Errorf("block 1 = pos %d END %q, want pos 100 END 120", blk.Pos, blk.Info["END"])
	}
	if got := blk.Samples[0].Data["DP"]; got != "15" {
		t.Errorf("block 1 DP = %q, want 15 (min of MIN_DP 25,15)", got)
	}
	if got := blk.Samples[0].Data["GQ"]; got != "45" {
		t.Errorf("block 1 GQ = %q, want 45 (min of 70,45)", got)
	}
	if got := blk.Samples[0].Data["PL"]; got != "0,20,200" {
		t.Errorf("block 1 PL = %q, want 0,20,200 (element-wise min)", got)
	}
	if res[1].Pos != 121 || res[1].Ref != "T" {
		t.Errorf("record 2 should be the passthrough variant, got pos %d", res[1].Pos)
	}
	if res[2].Pos != 122 || res[2].Info["END"] != "130" {
		t.Errorf("block 3 = pos %d END %q, want pos 122 END 130", res[2].Pos, res[2].Info["END"])
	}
}
