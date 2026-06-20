package vcf

import (
	"bufio"
	"bytes"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// TestFormatVCFFloatG covers the float-to-string rule that VCF output uses for
// QUAL and Float-typed INFO/FORMAT values. Upstream htslib's kputd is
// equivalent to C printf("%g", ...): six significant digits, trailing zeros
// stripped, and a switch to scientific notation (e+NN) for magnitudes outside
// the [0.0001, 999999] window. The differential fuzzer (gap A16) found a large
// QUAL such as 4294967296 was printed verbatim where htslib prints
// "4.29497e+09"; these cases lock that behavior down.
func TestFormatVCFFloatG(t *testing.T) {
	cases := []struct {
		in   float64
		want string
	}{
		// Small integers print without a decimal point.
		{60, "60"},
		{0, "0"},
		{1, "1"},
		{100000, "100000"},
		{999999, "999999"},
		// The %g window boundary: >999999 switches to scientific.
		{1000000, "1e+06"},
		{1234567, "1.23457e+06"},
		// The fuzzer's regression case.
		{4294967296, "4.29497e+09"},
		{1234567890, "1.23457e+09"},
		// Fractions keep six significant digits.
		{123.456, "123.456"},
		{12345.678, "12345.7"},
		{29.99, "29.99"},
		{0.5, "0.5"},
		{0.75, "0.75"},
		// Small-magnitude boundary: <0.0001 switches to scientific.
		{0.0001, "0.0001"},
		{0.00009999, "9.999e-05"},
		{2.5e-05, "2.5e-05"},
		{1e-20, "1e-20"},
		{1e20, "1e+20"},
		// Negatives.
		{-4294967296, "-4.29497e+09"},
		{-123.456, "-123.456"},
		// Special spellings.
		{math.Inf(1), "inf"},
		{math.Inf(-1), "-inf"},
		{math.NaN(), "nan"},
	}
	for _, c := range cases {
		if got := FormatVCFFloat64(c.in); got != c.want {
			t.Errorf("FormatVCFFloat64(%v) = %q, want %q", c.in, got, c.want)
		}
	}

	// Negative zero is the one value where Go's %g differs from C's; both
	// FormatVCFFloat helpers must spell it "-0".
	if got := FormatVCFFloat64(math.Copysign(0, -1)); got != "-0" {
		t.Errorf("FormatVCFFloat64(-0) = %q, want %q", got, "-0")
	}
	if got := FormatVCFFloat32(math.Copysign(0, -1)); got != "-0" {
		t.Errorf("FormatVCFFloat32(-0) = %q, want %q", got, "-0")
	}
}

// TestFormatVCFFloat32Narrowing verifies that FormatVCFFloat32 narrows to a
// 32-bit float before formatting, matching htslib's bcf1_t.qual storage. The
// sixth significant figure of 0.0001157624993 rounds differently in float32
// (…763) than in float64 (…762); upstream bcftools emits the float32 form.
func TestFormatVCFFloat32Narrowing(t *testing.T) {
	const in = 0.0001157624993
	if got, want := FormatVCFFloat32(in), "0.000115763"; got != want {
		t.Errorf("FormatVCFFloat32(%v) = %q, want %q (float32-narrowed)", in, got, want)
	}
	if got, want := FormatVCFFloat64(in), "0.000115762"; got != want {
		t.Errorf("FormatVCFFloat64(%v) = %q, want %q (no narrowing)", in, got, want)
	}
}

// TestQUALByteParityWithUpstream pipes a VCF whose QUAL column spans the small,
// large (scientific), fractional, and boundary ranges through the upstream
// `bcftools view` binary and asserts our QUAL formatter reproduces it
// byte-for-byte. It skips gracefully when the upstream binary is absent.
func TestQUALByteParityWithUpstream(t *testing.T) {
	bin, err := filepath.Abs("../../../reference_code/bcftools/bcftools")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(bin); err != nil {
		t.Skipf("upstream bcftools not built at %s; skipping byte-parity check", bin)
	}

	quals := []float64{
		60, 0, 29.99, 999999, 1000000, 4294967296, 1234567, 123.456,
		0.0001, 0.00009999, 100000, 12345.678, 1e-20, 1234567890,
		0.123456789, 3.14159265, 2.5e-05,
	}

	var vcf bytes.Buffer
	vcf.WriteString("##fileformat=VCFv4.2\n")
	vcf.WriteString("##contig=<ID=1,length=300000000>\n")
	vcf.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for i, q := range quals {
		// Use the full-precision decimal so upstream's atof sees the same value.
		vcf.WriteString("1\t")
		vcf.WriteString(strconv.Itoa((i + 1) * 100))
		vcf.WriteString("\t.\tA\tG\t")
		vcf.WriteString(strconv.FormatFloat(q, 'g', -1, 64))
		vcf.WriteString("\t.\t.\n")
	}

	cmd := exec.Command(bin, "view", "--no-version", "-")
	cmd.Stdin = bytes.NewReader(vcf.Bytes())
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("upstream bcftools view failed: %v", err)
	}

	// Collect upstream QUAL column for each data line.
	var upstreamQuals []string
	sc := bufio.NewScanner(bytes.NewReader(out))
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Split(line, "\t")
		if len(fields) < 6 {
			t.Fatalf("unexpected upstream line: %q", line)
		}
		upstreamQuals = append(upstreamQuals, fields[5])
	}
	if len(upstreamQuals) != len(quals) {
		t.Fatalf("got %d upstream records, want %d", len(upstreamQuals), len(quals))
	}

	for i, q := range quals {
		got := formatQual(q)
		want := upstreamQuals[i]
		if got != want {
			t.Errorf("QUAL %v: formatQual = %q, upstream bcftools = %q", q, got, want)
		}
	}
}
