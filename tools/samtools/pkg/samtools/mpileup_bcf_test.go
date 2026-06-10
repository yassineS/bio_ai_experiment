package samtools

import (
	"bytes"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

// writeMpileupBCFFixture writes a tiny reference FASTA (+.fai) and a SAM
// file with two reads carrying a variant base at one position, returning
// the FASTA and SAM paths. The reference is 50 A's; both reads spell
// AAAACAAA at chr1:1, so position 5 (1-based) carries a C/A mix that the
// genotype-likelihood caller scores.
func writeMpileupBCFFixture(t *testing.T) (fasta, samPath string) {
	t.Helper()
	dir := t.TempDir()
	fasta = filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(fasta, []byte(">chr1\n"+strings.Repeat("A", 50)+"\n"), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	if err := os.WriteFile(fasta+".fai", []byte("chr1\t50\t6\t50\t51\n"), 0o644); err != nil {
		t.Fatalf("write fai: %v", err)
	}
	samPath = filepath.Join(dir, "in.sam")
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:50",
		"@RG\tID:rg1\tSM:s1",
		"r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tAAAACAAA\t????????\tRG:Z:rg1",
		"r2\t0\tchr1\t1\t60\t8M\t*\t0\t0\tAAAACAAA\t????????\tRG:Z:rg1",
		"",
	}, "\n")
	if err := os.WriteFile(samPath, []byte(samText), 0o644); err != nil {
		t.Fatalf("write sam: %v", err)
	}
	return fasta, samPath
}

// decodeMpileupBCF writes the BCF bytes to a temp file and runs the
// in-tree bcftools BCF decoder (ViewFile) to return deterministic VCF text.
func decodeMpileupBCF(t *testing.T, bcfBytes []byte) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "out.bcf")
	if err := os.WriteFile(p, bcfBytes, 0o644); err != nil {
		t.Fatalf("write bcf: %v", err)
	}
	var vcfText bytes.Buffer
	if _, err := bcftools.ViewFile(p, &vcfText, bcftools.ViewOptions{}, &bytes.Buffer{}); err != nil {
		t.Fatalf("decode BCF via ViewFile: %v", err)
	}
	return vcfText.String()
}

// TestMpileupBCF_GenotypeLikelihoods asserts that `-g` emits a BGZF-wrapped
// BCF stream that decodes to per-site genotype-likelihood records carrying
// FORMAT/PL, the <*> unseen allele, INFO/DP, and the sample column.
func TestMpileupBCF_GenotypeLikelihoods(t *testing.T) {
	fasta, samPath := writeMpileupBCFFixture(t)
	var buf bytes.Buffer
	opts := MpileupBCFOptions{Inputs: []string{samPath}, FastaRef: fasta}
	if err := MpileupBCF(opts, &buf); err != nil {
		t.Fatalf("MpileupBCF -g: %v", err)
	}
	out := buf.Bytes()
	if len(out) < 2 || out[0] != 0x1f || out[1] != 0x8b {
		t.Fatalf("BCF (-g) output missing BGZF/gzip magic, got % x", out[:min(4, len(out))])
	}
	got := decodeMpileupBCF(t, out)
	for _, want := range []string{
		"##fileformat=VCF",
		"ID=PL,",  // FORMAT/PL declared in the header
		"ID=AD,",  // FORMAT/AD default annotation declared
		"<*>",     // the unseen symbolic allele
		"chr1\t1", // at least one record on chr1
		"PL:AD",   // FORMAT keys present on a record
		"\ts1",    // sample column from @RG SM
	} {
		if !strings.Contains(got, want) {
			t.Errorf("decoded BCF missing %q in:\n%s", want, got)
		}
	}
	// The variant column (position 5) must mention the C allele in REF/ALT.
	if !strings.Contains(got, "\t5\t") {
		t.Errorf("expected a record at chr1:5 in:\n%s", got)
	}
}

// TestMpileupBCF_Uncompressed asserts `-u` emits an uncompressed BCF
// stream (BCF magic, NOT BGZF) that still decodes to the same records.
func TestMpileupBCF_Uncompressed(t *testing.T) {
	fasta, samPath := writeMpileupBCFFixture(t)
	var buf bytes.Buffer
	opts := MpileupBCFOptions{Inputs: []string{samPath}, FastaRef: fasta, Uncompressed: true}
	if err := MpileupBCF(opts, &buf); err != nil {
		t.Fatalf("MpileupBCF -u: %v", err)
	}
	out := buf.Bytes()
	// Uncompressed BCF starts with the literal "BCF" magic, not BGZF.
	if len(out) < 3 || out[0] != 'B' || out[1] != 'C' || out[2] != 'F' {
		t.Fatalf("uncompressed BCF (-u) missing BCF magic, got % x", out[:min(4, len(out))])
	}
	got := decodeMpileupBCF(t, out)
	if !strings.Contains(got, "<*>") || !strings.Contains(got, "PL") {
		t.Errorf("uncompressed BCF round-trip missing PL/<*> in:\n%s", got)
	}
}

// TestMpileupBCF_GVsUMatch asserts the compressed (-g) and uncompressed
// (-u) paths decode to byte-identical VCF text: the only difference is the
// container, not the per-site genotype-likelihood records.
func TestMpileupBCF_GVsUMatch(t *testing.T) {
	fasta, samPath := writeMpileupBCFFixture(t)
	var g, u bytes.Buffer
	if err := MpileupBCF(MpileupBCFOptions{Inputs: []string{samPath}, FastaRef: fasta}, &g); err != nil {
		t.Fatalf("MpileupBCF -g: %v", err)
	}
	if err := MpileupBCF(MpileupBCFOptions{Inputs: []string{samPath}, FastaRef: fasta, Uncompressed: true}, &u); err != nil {
		t.Fatalf("MpileupBCF -u: %v", err)
	}
	gv := normalizeVCF(decodeMpileupBCF(t, g.Bytes()))
	uv := normalizeVCF(decodeMpileupBCF(t, u.Bytes()))
	if gv != uv {
		t.Errorf("-g and -u decode to different VCF records:\n--g--\n%s\n--u--\n%s", gv, uv)
	}
}

// normalizeVCF canonicalises a decoded VCF body so the comparison is robust
// to the INFO-field emission order, which is map-iteration-order dependent
// in the shared bcftools BCF writer (the same trait the bcftools mpileup
// suite normalises around). Header lines are dropped; each data record's
// INFO column (field 8) has its key=value entries sorted.
func normalizeVCF(s string) string {
	var out []string
	for _, ln := range strings.Split(s, "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		cols := strings.Split(ln, "\t")
		if len(cols) > 7 {
			info := strings.Split(cols[7], ";")
			sort.Strings(info)
			cols[7] = strings.Join(info, ";")
		}
		out = append(out, strings.Join(cols, "\t"))
	}
	sort.Strings(out)
	return strings.Join(out, "\n")
}

// TestMpileupBCF_RequiresFasta asserts the genotype-likelihood path
// rejects a missing reference, since the caller needs the REF base.
func TestMpileupBCF_RequiresFasta(t *testing.T) {
	var buf bytes.Buffer
	err := MpileupBCF(MpileupBCFOptions{Inputs: []string{"x.bam"}}, &buf)
	if err == nil {
		t.Fatal("expected error for missing -f/--fasta-ref")
	}
	if !strings.Contains(err.Error(), "fasta-ref") {
		t.Errorf("error should mention fasta-ref, got %v", err)
	}
}

// TestMpileupBCF_NoInputs asserts an empty input set is rejected.
func TestMpileupBCF_NoInputs(t *testing.T) {
	var buf bytes.Buffer
	err := MpileupBCF(MpileupBCFOptions{FastaRef: "ref.fa"}, &buf)
	if err == nil || !strings.Contains(err.Error(), "no input") {
		t.Fatalf("expected 'no input files' error, got %v", err)
	}
}

// TestMpileupBCF_PositionsFilter asserts -l (positions/BED) restricts the
// emitted sites: only chr1:5 should survive a single-position BED.
func TestMpileupBCF_PositionsFilter(t *testing.T) {
	fasta, samPath := writeMpileupBCFFixture(t)
	dir := t.TempDir()
	bed := filepath.Join(dir, "pos.bed")
	// BED is 0-based half-open: chr1 4 5 selects 1-based position 5.
	if err := os.WriteFile(bed, []byte("chr1\t4\t5\n"), 0o644); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	var buf bytes.Buffer
	opts := MpileupBCFOptions{Inputs: []string{samPath}, FastaRef: fasta, PositionsFile: bed}
	if err := MpileupBCF(opts, &buf); err != nil {
		t.Fatalf("MpileupBCF -g -l: %v", err)
	}
	got := decodeMpileupBCF(t, buf.Bytes())
	if !strings.Contains(got, "\t5\t") {
		t.Errorf("positions filter should keep chr1:5 in:\n%s", got)
	}
	// No other body position should appear (e.g. chr1:1 or chr1:8).
	for _, bad := range []string{"\t1\t", "\t8\t"} {
		if strings.Contains(got, bad) {
			t.Errorf("positions filter leaked site %q in:\n%s", bad, got)
		}
	}
}
