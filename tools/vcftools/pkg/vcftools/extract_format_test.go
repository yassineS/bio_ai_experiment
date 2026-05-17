package vcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestParity_ExtractFormatInfo_DP verifies --extract-FORMAT-info DP
// against an upstream golden built from sample.vcf:
//
//   - Sites whose FORMAT lacks DP (rows 19:111, 19:112, 20:1235237,
//     X:9, X:10, X:12 in sample.vcf) are SKIPPED entirely.
//   - Sites where a sample's value vector is too short (e.g. 20:1234567
//     S1 = "0/1:.:4" has DP at index 2 = "4"; 20:1230237 S1 = "0|0:54:.:..."
//     has DP at index 2 = ".") emit the raw colon-token, or "."
//     when the token slot is absent in the sample string.
//
// Ported from variant_file_format_convert.cpp:1204-1263 +
// vcf_entry.cpp:610-639.
func TestParity_ExtractFormatInfo_DP(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		ExtractFormatInfo: "DP",
	})
	got := readFileBytes(t, prefix+".DP.FORMAT")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "extract_format_dp.expected.DP.FORMAT"))
	if !bytes.Equal(got, want) {
		t.Errorf(".DP.FORMAT mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_ExtractFormatInfo_HQ verifies --extract-FORMAT-info HQ.
// Notable: HQ is a comma-separated pair so the per-sample column
// contains "10,10" / "." / ".,." raw.
func TestParity_ExtractFormatInfo_HQ(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		ExtractFormatInfo: "HQ",
	})
	got := readFileBytes(t, prefix+".HQ.FORMAT")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "extract_format_hq.expected.HQ.FORMAT"))
	if !bytes.Equal(got, want) {
		t.Errorf(".HQ.FORMAT mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_ExtractFormatInfo_GQ verifies --extract-FORMAT-info GQ.
// Covers the case where a sample's value vector is truncated and the
// runner emits "." (e.g. X:11 S2 = "./." with no further fields).
func TestParity_ExtractFormatInfo_GQ(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		ExtractFormatInfo: "GQ",
	})
	got := readFileBytes(t, prefix+".GQ.FORMAT")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "extract_format_gq.expected.GQ.FORMAT"))
	if !bytes.Equal(got, want) {
		t.Errorf(".GQ.FORMAT mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_ExtractFormatInfo_EdgeCases pins behaviour on a dedicated
// fixture that exercises every value-vector layout: (a) full vector
// (1:100 GT:DP:GQ), (b) missing trailing field (1:300 GT only — FORMAT
// lacks DP, row skipped), (c) partial vector (1:400 S2 = "0/1" with
// DP/GQ absent → both emit ".").
func TestParity_ExtractFormatInfo_EdgeCases(t *testing.T) {
	prefix := runVcftoolsParity(t, "extract_format_fixture.vcf", &Params{
		ExtractFormatInfo: "DP",
	})
	got := readFileBytes(t, prefix+".DP.FORMAT")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "extract_format_fixture_dp.expected.DP.FORMAT"))
	if !bytes.Equal(got, want) {
		t.Errorf(".DP.FORMAT mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_ExtractFormatInfo_EdgeCases_GQ — companion to the DP test
// above; verifies sites where FORMAT contains GQ but the per-sample
// value is missing/truncated.
func TestParity_ExtractFormatInfo_EdgeCases_GQ(t *testing.T) {
	prefix := runVcftoolsParity(t, "extract_format_fixture.vcf", &Params{
		ExtractFormatInfo: "GQ",
	})
	got := readFileBytes(t, prefix+".GQ.FORMAT")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "extract_format_fixture_gq.expected.GQ.FORMAT"))
	if !bytes.Equal(got, want) {
		t.Errorf(".GQ.FORMAT mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestExtractFormatInfo_UnknownTagIsEmpty pins that asking for a tag
// that no site lists in FORMAT produces only the header (no data
// rows). Mirrors upstream's `continue` at
// variant_file_format_convert.cpp:1247-1248.
func TestExtractFormatInfo_UnknownTagIsEmpty(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		ExtractFormatInfo: "ZZZ",
	})
	lines := readFileLines(t, prefix+".ZZZ.FORMAT")
	if len(lines) != 1 {
		t.Errorf("unknown FORMAT tag should yield header-only file; got %d lines: %v", len(lines), lines)
	}
	if len(lines) > 0 && !strings.HasPrefix(lines[0], "CHROM\tPOS") {
		t.Errorf("expected header CHROM\\tPOS..., got %q", lines[0])
	}
}

// TestExtractFormatInfo_EmptyNameRejected pins that an explicitly empty
// FORMAT name surfaces an error rather than producing
// `<prefix>..FORMAT`.
func TestExtractFormatInfo_EmptyNameRejected(t *testing.T) {
	// We bypass Run() and construct the runner directly so the test
	// pins the error message contract of newExtractFormatRunner.
	tmp := t.TempDir()
	if _, err := newExtractFormatRunner(filepath.Join(tmp, "out"), "", []string{"S1"}); err == nil {
		t.Errorf("expected error for empty FORMAT name; got nil")
	}
}

// TestExtractFormatInfo_NoSamplesHeader pins that the header still
// emits at least `CHROM\tPOS` when the VCF has no samples (e.g. a
// site-only VCF). The runner is informational in that case: no data
// rows can be emitted, but the file exists.
func TestExtractFormatInfo_NoSamplesHeader(t *testing.T) {
	tmp := t.TempDir()
	r, err := newExtractFormatRunner(filepath.Join(tmp, "out"), "DP", nil)
	if err != nil {
		t.Fatalf("newExtractFormatRunner: %v", err)
	}
	if err := r.close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	got, err := os.ReadFile(filepath.Join(tmp, "out.DP.FORMAT"))
	if err != nil {
		t.Fatalf("ReadFile: %v", err)
	}
	if string(got) != "CHROM\tPOS\n" {
		t.Errorf("header-only file mismatch; got %q", string(got))
	}
}
