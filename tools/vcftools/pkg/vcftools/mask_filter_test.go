package vcftools

import (
	"bytes"
	"path/filepath"
	"strings"
	"testing"
)

// ----- unit tests for the mask cursor state machine ------------------------

func TestMaskFilter_LoadValidation(t *testing.T) {
	cases := []struct {
		name string
		min  int
		want bool // true => want error
	}{
		{"zero ok", 0, false},
		{"nine ok", 9, false},
		{"five ok", 5, false},
		{"negative rejected", -1, true},
		{"ten rejected", 10, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			path := filepath.Join(vcftoolsFixtureDir(t), "mask_fasta.txt")
			_, err := loadMaskFilter(path, false, c.min)
			if (err != nil) != c.want {
				t.Errorf("loadMaskFilter(min=%d) err=%v, wantErr=%v", c.min, err, c.want)
			}
		})
	}
}

func TestMaskFilter_LoadMissingFile(t *testing.T) {
	if _, err := loadMaskFilter("/nonexistent/path/to/mask.txt", false, 0); err == nil {
		t.Fatal("loadMaskFilter(nonexistent) returned nil error, want error")
	}
}

// TestMaskFilter_ParseSlabs verifies the FASTA-style parser segments the file
// into per-chromosome (startPos, line) slabs in input order. This is the
// foundation that the cursor relies on.
func TestMaskFilter_ParseSlabs(t *testing.T) {
	src := `>1
0000000000
5555555555
1111111111
>2 some comment
9999999999
9999999999
`
	chroms, err := parseMaskFile(bytes.NewReader([]byte(src)))
	if err != nil {
		t.Fatalf("parseMaskFile: %v", err)
	}
	if len(chroms) != 2 {
		t.Fatalf("got %d chroms, want 2", len(chroms))
	}
	if chroms[0].name != "1" || chroms[1].name != "2" {
		t.Errorf("chrom names: got %q %q, want \"1\" \"2\"", chroms[0].name, chroms[1].name)
	}
	if len(chroms[0].slabs) != 3 {
		t.Fatalf("chrom 1 slabs: got %d, want 3", len(chroms[0].slabs))
	}
	// startPos cumulative offsets: 1, 11, 21
	wantStart := []int{1, 11, 21}
	for i, slab := range chroms[0].slabs {
		if slab.startPos != wantStart[i] {
			t.Errorf("chrom 1 slab[%d] startPos = %d, want %d", i, slab.startPos, wantStart[i])
		}
	}
	if len(chroms[1].slabs) != 2 {
		t.Errorf("chrom 2 slabs: got %d, want 2", len(chroms[1].slabs))
	}
	// Chrom 2's first slab still starts at 1 (per-chromosome reset).
	if chroms[1].slabs[0].startPos != 1 {
		t.Errorf("chrom 2 slab[0] startPos = %d, want 1", chroms[1].slabs[0].startPos)
	}
}

// TestMaskFilter_PassesDefault exercises the cursor: with --mask-min=0 only
// digit '0' passes. Sites must be visited in chromosome order matching the
// mask file (forward-only cursor).
func TestMaskFilter_PassesDefault(t *testing.T) {
	m := mustMaskFilter(t, "mask_fasta.txt", false, 0)

	cases := []struct {
		chrom string
		pos   int
		want  bool
	}{
		{"1", 5, true},   // mask digit '0'
		{"1", 15, false}, // mask digit '5'
		{"1", 25, false}, // mask digit '1'
		{"2", 5, false},  // mask digit '9'
	}
	for _, c := range cases {
		if got := m.passes(c.chrom, c.pos); got != c.want {
			t.Errorf("passes(%q,%d) = %v, want %v", c.chrom, c.pos, got, c.want)
		}
	}
}

func TestMaskFilter_PassesMin5(t *testing.T) {
	m := mustMaskFilter(t, "mask_fasta.txt", false, 5)
	if !m.passes("1", 5) || !m.passes("1", 15) || !m.passes("1", 25) {
		t.Error("chr1 sites (digits 0,5,1) should all pass with min=5")
	}
	if m.passes("2", 5) {
		t.Error("chr2:5 (digit 9) should fail with min=5")
	}
}

func TestMaskFilter_Invert(t *testing.T) {
	m := mustMaskFilter(t, "mask_fasta.txt", true, 5)
	if m.passes("1", 5) {
		t.Error("invert: chr1:5 (digit 0) should fail")
	}
	if !m.passes("2", 5) {
		t.Error("invert: chr2:5 (digit 9) should pass")
	}
}

// TestMaskFilter_OffEndDrops covers the case where a VCF position is past the
// end of the mask sequence for its chromosome: the partial mask only has
// chr1 1-10 but the VCF asks for chr1:15, chr1:25. Both must drop, and the
// cursor must not get stuck.
func TestMaskFilter_OffEndDrops(t *testing.T) {
	m := mustMaskFilter(t, "mask_partial.txt", false, 0)
	if !m.passes("1", 5) {
		t.Error("chr1:5 (digit 0) should pass with partial mask")
	}
	if m.passes("1", 15) {
		t.Error("chr1:15 is past end of partial mask, should drop")
	}
	if m.passes("1", 25) {
		t.Error("chr1:25 is past end of partial mask, should drop")
	}
	// After the cursor walked off the last chromosome, anything else drops.
	if m.passes("2", 5) {
		t.Error("after exhausted, chr2:5 must drop")
	}
}

// TestMaskFilter_OutOfOrderVCFDrops mirrors the upstream forward-only cursor:
// if the VCF presents chr2 before chr1 (relative to mask order), chr1 sites
// after the cursor advances onto chr2 are dropped.
func TestMaskFilter_OutOfOrderVCFDrops(t *testing.T) {
	m := mustMaskFilter(t, "mask_fasta.txt", false, 9)
	// First call: chr2 — cursor walks past chr1 in the mask.
	if !m.passes("2", 5) {
		t.Error("chr2:5 (digit 9, min=9) should pass")
	}
	// Subsequent chr1 call: cursor cannot rewind, must drop.
	if m.passes("1", 5) {
		t.Error("chr1:5 after cursor advanced to chr2 must drop")
	}
}

func TestMaskFilter_NilPasses(t *testing.T) {
	var m *maskFilter
	if !m.passes("1", 5) {
		t.Error("nil maskFilter must pass everything")
	}
}

// TestMaskFilter_DataBeforeHeader covers the parse path where lines appear
// before any `>NAME` header. Upstream's mask_chr is "" until the first header
// is read, so any leading data lines never match a VCF site; we drop them.
func TestMaskFilter_DataBeforeHeader(t *testing.T) {
	src := "1234567890\n>1\n0000\n"
	chroms, err := parseMaskFile(strings.NewReader(src))
	if err != nil {
		t.Fatalf("parseMaskFile: %v", err)
	}
	if len(chroms) != 1 {
		t.Fatalf("got %d chroms, want 1", len(chroms))
	}
	if len(chroms[0].slabs) != 1 || chroms[0].slabs[0].line != "0000" {
		t.Errorf("data-before-header should be dropped; got slabs=%+v", chroms[0].slabs)
	}
}

// ----- parity tests against upstream-generated goldens ---------------------

// TestParity_Mask_Default — `--mask FILE --recode` with default --mask-min=0.
// Upstream parameters.cpp:280 + entry_filters.cpp:674-752. Only sites whose
// mask digit is exactly '0' survive.
func TestParity_Mask_Default(t *testing.T) {
	prefix := runVcftoolsParity(t, "mask_fixture.vcf", &Params{
		Mask:   filepath.Join(vcftoolsFixtureDir(t), "mask_fasta.txt"),
		Recode: true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "mask_default.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Mask_Min5 — `--mask FILE --mask-min 5 --recode`. Sites with mask
// digit <= 5 are kept.
func TestParity_Mask_Min5(t *testing.T) {
	prefix := runVcftoolsParity(t, "mask_fixture.vcf", &Params{
		Mask:    filepath.Join(vcftoolsFixtureDir(t), "mask_fasta.txt"),
		MaskMin: 5,
		Recode:  true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "mask_min5.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_InvertMask_Min5 — `--invert-mask FILE --mask-min 5 --recode`.
// Inverted: sites with mask digit > 5 are kept.
func TestParity_InvertMask_Min5(t *testing.T) {
	prefix := runVcftoolsParity(t, "mask_fixture.vcf", &Params{
		Mask:       filepath.Join(vcftoolsFixtureDir(t), "mask_fasta.txt"),
		InvertMask: true,
		MaskMin:    5,
		Recode:     true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "invmask_min5.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Mask_Partial — partial mask file that covers only chr1:1-10.
// Sites past the end of the mask sequence (chr1:15, chr1:25) and on
// chromosomes absent from the mask (chr2) all drop.
func TestParity_Mask_Partial(t *testing.T) {
	prefix := runVcftoolsParity(t, "mask_fixture.vcf", &Params{
		Mask:   filepath.Join(vcftoolsFixtureDir(t), "mask_partial.txt"),
		Recode: true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "mask_partial.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// mustMaskFilter is a t.Helper that loads a maskFilter from a parity fixture
// or fails the test.
func mustMaskFilter(t *testing.T, fixture string, invert bool, min int) *maskFilter {
	t.Helper()
	path := filepath.Join(vcftoolsFixtureDir(t), fixture)
	m, err := loadMaskFilter(path, invert, min)
	if err != nil {
		t.Fatalf("loadMaskFilter(%s): %v", fixture, err)
	}
	return m
}
