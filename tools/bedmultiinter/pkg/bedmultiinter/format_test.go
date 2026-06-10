package bedmultiinter

import (
	"strings"
	"testing"
)

// TestDetectFormat covers upstream BedFile::parseLine's detection precedence:
// BED (cols 2,3 integer) wins over VCF, VCF (col 2 integer, >=8 cols) wins over
// GFF, and GFF requires exactly 8/9 cols with cols 4,5 integer.
func TestDetectFormat(t *testing.T) {
	cases := []struct {
		name   string
		line   string
		expect inputFormat
	}{
		{"bed3", "chr1\t10\t20", formatBED},
		{"bed6", "chr1\t10\t20\tname\t0\t+", formatBED},
		{"bed_wins_over_vcf", "chr1\t10\t20\tn\t0\t+\t1\t2", formatBED},
		{"vcf", "chr1\t32\t.\tA\tT\t0\tPASS\tDP=22", formatVCF},
		{"gff8", "chr1\tsrc\texon\t32\t40\t.\t-\t.", formatGFF},
		{"gff9", "chr1\tsrc\texon\t32\t40\t.\t-\t.\tgene_id 1", formatGFF},
		{"too_few", "chr1\t10", formatUnknown},
		{"non_integer_coords", "chr1\tabc\tdef\tx\ty\tz\tp\tq", formatUnknown},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := detectFormat(strings.Split(c.line, "\t"))
			if got != c.expect {
				t.Fatalf("detectFormat(%q) = %d, want %d", c.line, got, c.expect)
			}
		})
	}
}

// TestParseTyped checks coordinate conversion for each format matches upstream:
// BED is already 0-based half-open; VCF is POS-1 .. POS-1+len(REF); GFF is
// start-1 .. end.
func TestParseTyped(t *testing.T) {
	cases := []struct {
		name             string
		line             string
		format           inputFormat
		wantStart, wantE int
	}{
		{"bed", "chr1\t10\t20", formatBED, 10, 20},
		{"vcf_snp", "chr1\t32\t.\tA\tT\t0\tPASS\t.", formatVCF, 31, 32},
		{"vcf_del", "chr1\t10\t.\tACGT\tA\t0\tPASS\t.", formatVCF, 9, 13},
		{"gff", "chr1\tsrc\texon\t12\t18\t.\t+\t.", formatGFF, 11, 18},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			iv, err := parseTyped(strings.Split(c.line, "\t"), c.format)
			if err != nil {
				t.Fatalf("parseTyped: %v", err)
			}
			if iv.start != c.wantStart || iv.end != c.wantE {
				t.Fatalf("parseTyped(%q) = [%d,%d), want [%d,%d)",
					c.line, iv.start, iv.end, c.wantStart, c.wantE)
			}
		})
	}
}

// TestReadIntervals_ForcedVCFHeader confirms that a `##fileformat=VCF` header
// line forces VCF detection even though the first data line, on its own, could
// look ambiguous.
func TestReadIntervals_ForcedVCFHeader(t *testing.T) {
	in := "##fileformat=VCF4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t5\t.\tAC\tG\t0\tPASS\t.\n"
	ivs, err := readIntervals(strings.NewReader(in))
	if err != nil {
		t.Fatalf("readIntervals: %v", err)
	}
	if len(ivs) != 1 {
		t.Fatalf("want 1 interval, got %d", len(ivs))
	}
	if ivs[0].chrom != "chr1" || ivs[0].start != 4 || ivs[0].end != 6 {
		t.Fatalf("got %+v, want {chr1 4 6}", ivs[0])
	}
}

// TestReadIntervals_UnknownFormat surfaces an error for a non-conforming line.
func TestReadIntervals_UnknownFormat(t *testing.T) {
	if _, err := readIntervals(strings.NewReader("just one column\n")); err == nil {
		t.Fatalf("expected error for unparseable input")
	}
}
