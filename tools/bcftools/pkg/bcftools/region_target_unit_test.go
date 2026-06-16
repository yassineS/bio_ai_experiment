// Binary-free unit tests for the shared -r/-R/-t/-T region/target selection
// (region_target.go). These pin the region-overlap (-r/-R) vs
// target-start-in-region (-t/-T) distinction, the negation marker, the spec
// parsing, and the BED/region-list file loader using small in-temp-dir
// fixtures. No upstream binary or reference_code submodule is required.
package bcftools

import (
	"os"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// TestUnitRegionTargetStartInAny pins the -t/-T positional semantics: a
// variant is selected only when its 1-based POS falls within [beg,end] of some
// window on the same chromosome (the REF span is irrelevant).
func TestUnitRegionTargetStartInAny(t *testing.T) {
	regions := []region{{chrom: "chr1", beg: 100, end: 200}}
	tests := []struct {
		name  string
		chrom string
		pos   int
		ref   string
		want  bool
	}{
		{"start inside", "chr1", 150, "A", true},
		{"start at beg", "chr1", 100, "A", true},
		{"start at end", "chr1", 200, "A", true},
		{"start below window", "chr1", 99, "A", false},
		{"start above window", "chr1", 201, "A", false},
		{"wrong chrom", "chr2", 150, "A", false},
		{
			// An indel starting at 99 spanning 99..103 still fails startInAny: its
			// START is outside the window even though its span overlaps it.
			name: "span overlaps but start outside", chrom: "chr1", pos: 99, ref: "ACGTT", want: false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &vcf.Variant{Chrom: tt.chrom, Pos: tt.pos, Ref: tt.ref}
			if got := startInAny(v, regions); got != tt.want {
				t.Fatalf("startInAny(%s:%d ref=%q) = %v, want %v", tt.chrom, tt.pos, tt.ref, got, tt.want)
			}
		})
	}
}

// TestUnitRegionTargetKeepRegionVsTarget pins the headline -r vs -t difference:
// an indel at POS=100 spanning 100..104 is kept by `-r chr:102-102` (span
// overlap) but dropped by `-t chr:102-102` (start 100 not in 102..102).
func TestUnitRegionTargetKeep(t *testing.T) {
	// REF "ACGTT" => span [100,104].
	indel := &vcf.Variant{Chrom: "chr1", Pos: 100, Ref: "ACGTT"}

	t.Run("region overlap keeps the indel", func(t *testing.T) {
		f := regionTargetFilter{
			regions:    []region{{chrom: "chr1", beg: 102, end: 102}},
			hasRegions: true,
		}
		if !f.keep(indel) {
			t.Fatal("expected -r chr1:102-102 to keep the overlapping indel")
		}
	})
	t.Run("target start excludes the indel", func(t *testing.T) {
		f := regionTargetFilter{
			targets:    []region{{chrom: "chr1", beg: 102, end: 102}},
			hasTargets: true,
		}
		if f.keep(indel) {
			t.Fatal("expected -t chr1:102-102 to drop the indel whose start is 100")
		}
	})
	t.Run("target start within keeps the variant", func(t *testing.T) {
		f := regionTargetFilter{
			targets:    []region{{chrom: "chr1", beg: 100, end: 100}},
			hasTargets: true,
		}
		if !f.keep(indel) {
			t.Fatal("expected -t chr1:100-100 to keep the variant starting at 100")
		}
	})
}

// TestUnitRegionTargetKeepNegation pins the leading-'^' target negation: a
// negated target list EXCLUDES the matches and keeps everything else.
func TestUnitRegionTargetKeepNegation(t *testing.T) {
	f := regionTargetFilter{
		targets:        []region{{chrom: "chr1", beg: 100, end: 100}},
		targetsNegated: true,
		hasTargets:     true,
	}
	inWindow := &vcf.Variant{Chrom: "chr1", Pos: 100, Ref: "A"}
	if f.keep(inWindow) {
		t.Fatal("negated -t should drop a variant whose start is in the window")
	}
	outWindow := &vcf.Variant{Chrom: "chr1", Pos: 101, Ref: "A"}
	if !f.keep(outWindow) {
		t.Fatal("negated -t should keep a variant whose start is outside the window")
	}
}

// TestUnitRegionTargetKeepAnded pins that -r and -t are ANDed: a variant must
// pass both the region overlap and the target start test.
func TestUnitRegionTargetKeepAnded(t *testing.T) {
	f := regionTargetFilter{
		regions:    []region{{chrom: "chr1", beg: 90, end: 110}},
		hasRegions: true,
		targets:    []region{{chrom: "chr1", beg: 100, end: 100}},
		hasTargets: true,
	}
	// Start 100 is in the target and overlaps the region: kept.
	if !f.keep(&vcf.Variant{Chrom: "chr1", Pos: 100, Ref: "A"}) {
		t.Fatal("variant matching both region and target should be kept")
	}
	// Start 105 overlaps the region but is not in the target window: dropped.
	if f.keep(&vcf.Variant{Chrom: "chr1", Pos: 105, Ref: "A"}) {
		t.Fatal("variant outside target window should be dropped even if region overlaps")
	}
}

// TestUnitRegionTargetActive pins that the zero value is a no-op filter.
func TestUnitRegionTargetActive(t *testing.T) {
	var zero regionTargetFilter
	if zero.active() {
		t.Fatal("zero-value filter must be inactive")
	}
	if !(&regionTargetFilter{hasRegions: true}).active() {
		t.Fatal("filter with hasRegions must be active")
	}
	if !(&regionTargetFilter{hasTargets: true}).active() {
		t.Fatal("filter with hasTargets must be active")
	}
}

// TestUnitRegionTargetSplitNegation pins the '^' stripping helper.
func TestUnitRegionTargetSplitNegation(t *testing.T) {
	tests := []struct {
		in       string
		wantNeg  bool
		wantRest string
	}{
		{"^chr1:100", true, "chr1:100"},
		{"chr1:100", false, "chr1:100"},
		{"^", true, ""},
		{"", false, ""},
	}
	for _, tt := range tests {
		neg, rest := splitNegation(tt.in)
		if neg != tt.wantNeg || rest != tt.wantRest {
			t.Errorf("splitNegation(%q) = (%v,%q), want (%v,%q)", tt.in, neg, rest, tt.wantNeg, tt.wantRest)
		}
	}
}

// TestUnitRegionTargetParseRegionSpecs pins the verbatim-label retention: the
// inline token text is kept as the spec label.
func TestUnitRegionTargetParseRegionSpecs(t *testing.T) {
	specs, err := parseRegionSpecs([]string{"chr1", "chr2:200-405"})
	if err != nil {
		t.Fatalf("parseRegionSpecs error: %v", err)
	}
	if len(specs) != 2 {
		t.Fatalf("got %d specs, want 2", len(specs))
	}
	if specs[0].label != "chr1" || specs[0].region.chrom != "chr1" {
		t.Errorf("spec[0] = %+v, want label chr1", specs[0])
	}
	if specs[1].label != "chr2:200-405" {
		t.Errorf("spec[1].label = %q, want chr2:200-405", specs[1].label)
	}
	if specs[1].region.beg != 200 || specs[1].region.end != 405 {
		t.Errorf("spec[1].region = %+v, want beg 200 end 405", specs[1].region)
	}
}

// TestUnitRegionTargetLoadBED pins the BED (.bed, 0-based half-open) loader: a
// "chr1 99 104" line becomes the 1-based-inclusive "chr1:100-104" spec.
func TestUnitRegionTargetLoadBED(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "regions.bed")
	content := "# a comment line\n" +
		"\n" + // blank line skipped
		"chr1\t99\t104\n" +
		"chr2\t0\t10\textra\tcols\n"
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	got, err := loadRegionTargetFile(path)
	if err != nil {
		t.Fatalf("loadRegionTargetFile error: %v", err)
	}
	want := []string{"chr1:100-104", "chr2:1-10"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("loadRegionTargetFile(BED) = %v, want %v", got, want)
	}
}

// TestUnitRegionTargetLoadRegionList pins the non-BED (1-based) loader: a
// two-column line is a single position, a three+-column line is beg..end (NOT
// shifted, unlike BED). It also pins the exported wrapper.
func TestUnitRegionTargetLoadRegionList(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "regions.txt")
	content := "chr1\t150\n" + // single position -> 150-150
		"chr2\t200\t405\n" // beg..end, 1-based inclusive, no shift
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	got, err := loadRegionTargetFile(path)
	if err != nil {
		t.Fatalf("loadRegionTargetFile error: %v", err)
	}
	want := []string{"chr1:150-150", "chr2:200-405"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("loadRegionTargetFile(list) = %v, want %v", got, want)
	}
	// The exported wrapper must agree byte-for-byte.
	gotExported, err := LoadRegionTargetFile(path)
	if err != nil {
		t.Fatalf("LoadRegionTargetFile error: %v", err)
	}
	if !reflect.DeepEqual(gotExported, want) {
		t.Fatalf("LoadRegionTargetFile = %v, want %v", gotExported, want)
	}
}

// TestUnitRegionTargetLoadSingleColumnError pins the synced-reader rule that a
// single-column line (a bare contig or chr:beg-end string) is a parse error,
// distinct from `bcftools view`'s regidx which would accept it.
func TestUnitRegionTargetLoadSingleColumnError(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "single.txt")
	if err := os.WriteFile(path, []byte("chr1:100-200\n"), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	if _, err := loadRegionTargetFile(path); err == nil {
		t.Fatal("expected a parse error for a single-column region-list line")
	}
}

// TestUnitRegionTargetLoadBadNumber pins that a non-numeric position is an
// error in both the two-column and three-column branches.
func TestUnitRegionTargetLoadBadNumber(t *testing.T) {
	dir := t.TempDir()
	cases := map[string]string{
		"bad_pos.txt": "chr1\tnope\n",
		"bad_beg.bed": "chr1\tnope\t10\n",
		"bad_end.bed": "chr1\t5\tnope\n",
	}
	for name, content := range cases {
		path := filepath.Join(dir, name)
		if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
			t.Fatalf("write %s: %v", name, err)
		}
		if _, err := loadRegionTargetFile(path); err == nil {
			t.Errorf("%s: expected a parse error", name)
		}
	}
}
