package bedsummary

import (
	"bytes"
	"strings"
	"testing"
)

// runSummary parses the genome, runs Run, and returns the rendered output.
func runSummary(t *testing.T, bedInput, genomeInput string, opts Options) string {
	t.Helper()
	g, err := ParseGenome(strings.NewReader(genomeInput))
	if err != nil {
		t.Fatalf("ParseGenome: %v", err)
	}
	var out bytes.Buffer
	if err := Run(strings.NewReader(bedInput), g, &out, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return out.String()
}

// TestUnit_Header checks the exact upstream 10-column header.
func TestUnit_Header(t *testing.T) {
	got := runSummary(t, "chr1\t0\t10\n", "chr1\t1000\n", Options{})
	wantHeader := "chrom\tchrom_length\tnum_ivls\ttotal_ivl_bp\tchrom_frac_genome\t" +
		"frac_all_ivls\tfrac_all_bp\tmin\tmax\tmean\n"
	if !strings.HasPrefix(got, wantHeader) {
		t.Fatalf("header mismatch.\nwant prefix: %q\ngot: %q", wantHeader, got)
	}
}

// TestUnit_Aggregation verifies the full per-chrom + all aggregation, including
// the per-data-row trailing tab and the literal "1.0" fraction columns on the
// "all" row. Mirrors the byte layout of upstream `bedtools summary`.
func TestUnit_Aggregation(t *testing.T) {
	bedIn := "chr1\t0\t100\nchr1\t50\t200\nchr2\t10\t20\nchr1\t0\t10\nchr3\t5\t500\n"
	genIn := "chr1\t1000\nchr2\t500\nchr3\t2000\nchr4\t300\n"
	got := runSummary(t, bedIn, genIn, Options{})

	want := "chrom\tchrom_length\tnum_ivls\ttotal_ivl_bp\tchrom_frac_genome\t" +
		"frac_all_ivls\tfrac_all_bp\tmin\tmax\tmean\n" +
		"chr1\t1000\t3\t260\t0.263157895\t0.600000000\t0.339869281\t10\t150\t86.666666667\t\n" +
		"chr2\t500\t1\t10\t0.131578947\t0.200000000\t0.013071895\t10\t10\t10.000000000\t\n" +
		"chr3\t2000\t1\t495\t0.526315789\t0.200000000\t0.647058824\t495\t495\t495.000000000\t\n" +
		"chr4\t300\t0\t0\t0.078947368\t0.000000000\t0.000000000\t-1\t-1\t-1\n" +
		"all\t3800\t5\t765\t1.0\t1.0\t1.0\t10\t495\t153.000000000\n"
	if got != want {
		t.Fatalf("aggregation mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestUnit_ChromOrderFollowsGenome verifies chromosomes are reported in the
// order they appear in the genome file, not in input order.
func TestUnit_ChromOrderFollowsGenome(t *testing.T) {
	// BED lists chr2 before chr1; genome lists chr1 first.
	got := runSummary(t, "chr2\t0\t10\nchr1\t0\t5\n", "chr1\t100\nchr2\t100\n", Options{NoHeader: true})
	lines := strings.Split(strings.TrimRight(got, "\n"), "\n")
	if !strings.HasPrefix(lines[0], "chr1\t") || !strings.HasPrefix(lines[1], "chr2\t") {
		t.Fatalf("chrom order should follow genome (chr1 then chr2), got:\n%s", got)
	}
}

// TestUnit_MissingChromError verifies an input chromosome absent from the
// genome file is a hard error (matching upstream).
func TestUnit_MissingChromError(t *testing.T) {
	g, err := ParseGenome(strings.NewReader("chr1\t1000\n"))
	if err != nil {
		t.Fatalf("ParseGenome: %v", err)
	}
	var out bytes.Buffer
	err = Run(strings.NewReader("chrZ\t0\t50\n"), g, &out, Options{})
	if err == nil {
		t.Fatal("expected error for chromosome not in genome file")
	}
	if !strings.Contains(err.Error(), "chrZ") {
		t.Errorf("error should name the missing chromosome, got: %v", err)
	}
}

// TestUnit_NilGenome verifies Run rejects a nil genome.
func TestUnit_NilGenome(t *testing.T) {
	if err := Run(strings.NewReader(""), nil, &bytes.Buffer{}, Options{}); err == nil {
		t.Fatal("expected error for nil genome")
	}
}

// TestUnit_NoHeader suppresses the header line.
func TestUnit_NoHeader(t *testing.T) {
	got := runSummary(t, "chr1\t0\t10\n", "chr1\t1000\n", Options{NoHeader: true})
	if strings.HasPrefix(got, "chrom\t") {
		t.Errorf("expected no header, got: %q", got)
	}
}

// TestUnit_EmptyInput emits every genome chromosome as a default (-1) row plus
// the "all" row.
func TestUnit_EmptyInput(t *testing.T) {
	got := runSummary(t, "", "chr1\t100\nchr2\t200\n", Options{NoHeader: true})
	want := "chr1\t100\t0\t0\t0.333333333\t0.000000000\t0.000000000\t-1\t-1\t-1\n" +
		"chr2\t200\t0\t0\t0.666666667\t0.000000000\t0.000000000\t-1\t-1\t-1\n" +
		"all\t300\t0\t0\t1.0\t1.0\t1.0\t0\t0\t0.000000000\n"
	if got != want {
		t.Fatalf("empty-input mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestUnit_BadBED surfaces a parse error from a malformed coordinate.
func TestUnit_BadBED(t *testing.T) {
	g, _ := ParseGenome(strings.NewReader("chr1\t1000\n"))
	if err := Run(strings.NewReader("chr1\tnotanumber\t100\n"), g, &bytes.Buffer{}, Options{}); err == nil {
		t.Error("expected error on bad start coordinate")
	}
}

// TestUnit_BadGenome surfaces a parse error from a 1-column genome file.
func TestUnit_BadGenome(t *testing.T) {
	if _, err := ParseGenome(strings.NewReader("chr1\n")); err == nil {
		t.Error("expected error on 1-column genome file")
	}
}
