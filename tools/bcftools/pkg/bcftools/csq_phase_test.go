package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// The GFF phase-consistency parity tests below mirror upstream gff.c
// tscript_init_cds. They were verified byte-for-byte against
// bcftools 1.23.1 (bin/upstream/bcftools) on the same crafted fixtures.

// writeCSQFixture writes a FASTA + GFF3 + VCF trio into dir and returns
// their paths. The FASTA is a fixed 600 bp chr1 so the crafted CDS
// coordinates below always resolve.
func writeCSQFixture(t *testing.T, dir, gff, vcf string) (fa, g, v string) {
	t.Helper()
	// Deterministic 600 bp reference; content is irrelevant to the phase
	// check (which is purely coordinate/phase arithmetic) but must be long
	// enough for the crafted CDS/variant coordinates.
	var b strings.Builder
	b.WriteString(">chr1\n")
	const row = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" // 60 bp
	for i := 0; i < 10; i++ {
		b.WriteString(row)
		b.WriteByte('\n')
	}
	fa = filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(fa, []byte(b.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	g = filepath.Join(dir, "anno.gff3")
	if err := os.WriteFile(g, []byte(gff), 0o644); err != nil {
		t.Fatal(err)
	}
	v = filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(v, []byte(vcf), 0o644); err != nil {
		t.Fatal(err)
	}
	return fa, g, v
}

// forwardBadGFF is a forward-strand transcript whose second CDS exon has
// a phase (2) that disagrees with the cumulative CDS length (30, i.e.
// phase should be 0). Upstream errors on this without --force.
const forwardBadGFF = `##gff-version 3
chr1	test	gene	11	80	.	+	.	ID=gene:G1;Name=GENE1;biotype=protein_coding
chr1	test	mRNA	11	80	.	+	.	ID=transcript:T1;Parent=gene:G1;biotype=protein_coding
chr1	test	CDS	11	40	.	+	0	Parent=transcript:T1
chr1	test	CDS	51	80	.	+	2	Parent=transcript:T1
`

// forwardBadVCF places a het SNP inside the (bad) second CDS exon.
const forwardBadVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	55	.	G	A	50	PASS	.	GT	0/1
`

// TestCSQPhaseInconsistentErrorsByDefault pins the upstream default:
// an inconsistent GFF3 CDS phase is a hard error (no --force). The error
// text (after the wrapper's "bcftools csq: " prefix) is byte-identical to
// upstream 1.23.1's gff.c message.
func TestCSQPhaseInconsistentErrorsByDefault(t *testing.T) {
	dir := t.TempDir()
	fa, g, _ := writeCSQFixture(t, dir, forwardBadGFF, forwardBadVCF)

	_, err := loadCSQIndexUnified(fa, g, "", "", "", csqPhaseCheck{Force: false, Verbosity: 1})
	if err == nil {
		t.Fatal("expected an error for inconsistent GFF phase, got nil")
	}
	const want = "Error: GFF3 assumption failed for transcript T1, CDS=51: phase!=len%3 (phase=1, len=30). Use the --force option to proceed anyway (at your own risk)."
	if got := err.Error(); got != want {
		t.Errorf("error mismatch\n got: %q\nwant: %q", got, want)
	}
}

// TestCSQPhaseInconsistentErrorsWithVerbose confirms that -v/--verbose
// alone does NOT downgrade the error: upstream only warns-and-skips when
// --force is set (verbosity merely gates the warning text). A higher
// verbosity without --force still errors.
func TestCSQPhaseInconsistentErrorsWithVerbose(t *testing.T) {
	dir := t.TempDir()
	fa, g, _ := writeCSQFixture(t, dir, forwardBadGFF, forwardBadVCF)

	_, err := loadCSQIndexUnified(fa, g, "", "", "", csqPhaseCheck{Force: false, Verbosity: 2})
	if err == nil {
		t.Fatal("expected an error under -v without --force, got nil")
	}
	if !strings.Contains(err.Error(), "GFF3 assumption failed for transcript T1") {
		t.Errorf("unexpected error under -v: %v", err)
	}
}

// TestCSQPhaseInconsistentForceSkips confirms that --force downgrades the
// error to a warn-and-skip: the offending transcript is dropped (its CDS
// removed) so the overlapping variant degrades to an intron consequence,
// exactly as upstream 1.23.1 emits. The warning text is also asserted.
func TestCSQPhaseInconsistentForceSkips(t *testing.T) {
	dir := t.TempDir()
	fa, g, v := writeCSQFixture(t, dir, forwardBadGFF, forwardBadVCF)

	// Capture the warning written to stderr by the package.
	origStderr := os.Stderr
	r, w, _ := os.Pipe()
	os.Stderr = w
	var out bytes.Buffer
	_, err := CSQFile(v, &out, CSQOptions{
		FastaRef:     fa,
		GFFAnnot:     g,
		Phase:        'a',
		CustomTag:    "BCSQ",
		Force:        true,
		Verbosity:    1,
		OutputFormat: OutputVCF,
	})
	w.Close()
	os.Stderr = origStderr
	var errbuf bytes.Buffer
	_, _ = errbuf.ReadFrom(r)

	if err != nil {
		t.Fatalf("CSQFile under --force returned error: %v", err)
	}
	const wantWarn = "Warning: The GFF has inconsistent phase column in transcript T1, skipping. CDS pos=51: phase!=len%3 (phase=1, len=30)"
	if !strings.Contains(errbuf.String(), wantWarn) {
		t.Errorf("warning mismatch\n got: %q\nwant substring: %q", errbuf.String(), wantWarn)
	}
	// The skipped transcript degrades to an intron container, so the SNP
	// is annotated as intron|GENE1||protein_coding (matching upstream).
	if !strings.Contains(out.String(), "BCSQ=intron|GENE1||protein_coding") {
		var recs []string
		for _, line := range strings.Split(out.String(), "\n") {
			if strings.HasPrefix(line, "chr1") {
				recs = append(recs, line)
			}
		}
		t.Errorf("expected intron consequence under --force, got records: %v", recs)
	}
}

// reverseWalkbackGFF is a reverse-strand transcript with a 5' incomplete
// CDS whose leading phase (2) exceeds the 1 bp 5' CDS exon, forcing the
// multi-exon walk-back (upstream gff.c STRAND_REV branch). The 5' end on
// the reverse strand is the highest-coordinate CDS (101..101), so the
// residual phase after consuming that 1 bp exon spills into 40..99.
const reverseWalkbackGFF = `##gff-version 3
chr1	test	gene	40	101	.	-	.	ID=gene:GR;Name=REVG;biotype=protein_coding
chr1	test	mRNA	40	101	.	-	.	ID=transcript:TR;Parent=gene:GR;biotype=protein_coding
chr1	test	CDS	40	99	.	-	1	Parent=transcript:TR
chr1	test	CDS	101	101	.	-	2	Parent=transcript:TR
`

// reverseWalkbackVCF places het SNPs inside the reverse CDS. The REF
// bases match the ACGT-repeat reference used by the test (base at 1-based
// pos p is "ACGT"[(p-1)%4]): pos 50=C, 60=T, 70=C.
const reverseWalkbackVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=600>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	50	.	C	A	50	PASS	.	GT	0/1
chr1	60	.	T	A	50	PASS	.	GT	0/1
chr1	70	.	C	A	50	PASS	.	GT	0/1
`

// TestCSQReverseStrandMultiExonWalkback pins FIX B: the reverse-strand 5'
// incomplete-CDS phase walk-back across multiple exons. The expected
// consequences were verified byte-for-byte against upstream bcftools
// 1.23.1 (`bcftools csq -p a --force`) on this exact fixture; a
// single-exon-only trim would put the whole CDS out of frame and produce
// different amino-acid positions.
func TestCSQReverseStrandMultiExonWalkback(t *testing.T) {
	dir := t.TempDir()
	// Use a fixed reference whose bases pin the amino-acid changes below.
	fa := filepath.Join(dir, "ref.fa")
	var b strings.Builder
	b.WriteString(">chr1\n")
	const row = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	for i := 0; i < 10; i++ {
		b.WriteString(row)
		b.WriteByte('\n')
	}
	if err := os.WriteFile(fa, []byte(b.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	g := filepath.Join(dir, "rev.gff3")
	if err := os.WriteFile(g, []byte(reverseWalkbackGFF), 0o644); err != nil {
		t.Fatal(err)
	}
	v := filepath.Join(dir, "rev.vcf")
	if err := os.WriteFile(v, []byte(reverseWalkbackVCF), 0o644); err != nil {
		t.Fatal(err)
	}

	var out bytes.Buffer
	if _, err := CSQFile(v, &out, CSQOptions{
		FastaRef:     fa,
		GFFAnnot:     g,
		Phase:        'a',
		CustomTag:    "BCSQ",
		Force:        true,
		Verbosity:    1,
		OutputFormat: OutputVCF,
	}); err != nil {
		t.Fatalf("CSQFile: %v", err)
	}

	// Expected BCSQ strings (upstream 1.23.1, -p a --force). The amino-acid
	// positions (17, 13, 10) are only correct if the walk-back trimmed the
	// leading phase across both exons.
	wantSub := []string{
		"BCSQ=missense|REVG|TR|protein_coding|-|17V>17L|50C>A",
		"BCSQ=synonymous|REVG|TR|protein_coding|-|13V|60T>A",
		"BCSQ=missense|REVG|TR|protein_coding|-|10R>10L|70C>A",
	}
	got := out.String()
	for _, w := range wantSub {
		if !strings.Contains(got, w) {
			t.Errorf("walk-back consequence missing: %q\nfull output:\n%s", w, got)
		}
	}
}

// reverseWalkbackNonMaskingGFF is a reverse-strand 5'-incomplete CDS whose
// RESIDUAL exon phase DIFFERS from the loop-residual phase, so it exercises
// the walk-back's final trim without the masking coincidence present in
// reverseWalkbackGFF. Here the 1 bp 5' exon (101..101) carries phase 2, but
// so does the residual exon (40..99): after consuming the 1 bp exon the loop
// residual is 1, whereas the residual exon's OWN phase is 2. Upstream
// gff.c:780 subtracts the exon's own phase (2), so a walk-back that instead
// subtracts the loop residual (1) mis-frames the whole CDS. This fixture
// therefore fails when the bug (End -= ph) is present and passes once the
// fix (End -= exon.Phase) matches upstream.
const reverseWalkbackNonMaskingGFF = `##gff-version 3
chr1	test	gene	40	101	.	-	.	ID=gene:GR;Name=REVG;biotype=protein_coding
chr1	test	mRNA	40	101	.	-	.	ID=transcript:TR;Parent=gene:GR;biotype=protein_coding
chr1	test	CDS	40	99	.	-	2	Parent=transcript:TR
chr1	test	CDS	101	101	.	-	2	Parent=transcript:TR
`

// TestCSQReverseStrandWalkbackNonMasking is the NON-MASKING regression for
// the reverse-strand multi-exon 5' phase walk-back. The expected BCSQ
// strings were regenerated byte-for-byte from upstream bcftools 1.23.1
// (`bcftools csq -p a --force`) on reverseWalkbackNonMaskingGFF +
// reverseWalkbackVCF. With the pre-fix code (subtracting the loop residual
// phase 1 instead of the residual exon's own phase 2) ours produced the
// wrong frame (17V>17L / 13V / 10R>10L); after the fix it matches upstream
// (16T / 13Y>13F / 10V>10L).
func TestCSQReverseStrandWalkbackNonMasking(t *testing.T) {
	dir := t.TempDir()
	fa := filepath.Join(dir, "ref.fa")
	var b strings.Builder
	b.WriteString(">chr1\n")
	const row = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	for i := 0; i < 10; i++ {
		b.WriteString(row)
		b.WriteByte('\n')
	}
	if err := os.WriteFile(fa, []byte(b.String()), 0o644); err != nil {
		t.Fatal(err)
	}
	g := filepath.Join(dir, "rev.gff3")
	if err := os.WriteFile(g, []byte(reverseWalkbackNonMaskingGFF), 0o644); err != nil {
		t.Fatal(err)
	}
	v := filepath.Join(dir, "rev.vcf")
	if err := os.WriteFile(v, []byte(reverseWalkbackVCF), 0o644); err != nil {
		t.Fatal(err)
	}

	var out bytes.Buffer
	if _, err := CSQFile(v, &out, CSQOptions{
		FastaRef:     fa,
		GFFAnnot:     g,
		Phase:        'a',
		CustomTag:    "BCSQ",
		Force:        true,
		Verbosity:    1,
		OutputFormat: OutputVCF,
	}); err != nil {
		t.Fatalf("CSQFile: %v", err)
	}

	// Expected BCSQ strings (upstream 1.23.1, -p a --force). These are only
	// correct if the final walk-back trim subtracts the residual exon's OWN
	// phase (2), matching gff.c:780; the buggy loop-residual (1) yields
	// 17V>17L / 13V / 10R>10L instead.
	wantSub := []string{
		"BCSQ=synonymous|REVG|TR|protein_coding|-|16T|50C>A",
		"BCSQ=missense|REVG|TR|protein_coding|-|13Y>13F|60T>A",
		"BCSQ=missense|REVG|TR|protein_coding|-|10V>10L|70C>A",
	}
	got := out.String()
	for _, w := range wantSub {
		if !strings.Contains(got, w) {
			t.Errorf("non-masking walk-back consequence missing: %q\nfull output:\n%s", w, got)
		}
	}
}

// TestCSQPhaseConsistentUnaffected guards the common case: a fully
// phase-consistent GFF (all phases agree with cumulative CDS length) must
// pass the check without error under the default (no --force) settings.
func TestCSQPhaseConsistentUnaffected(t *testing.T) {
	const okGFF = `##gff-version 3
chr1	test	gene	11	80	.	+	.	ID=gene:G1;Name=GENE1;biotype=protein_coding
chr1	test	mRNA	11	80	.	+	.	ID=transcript:T1;Parent=gene:G1;biotype=protein_coding
chr1	test	CDS	11	40	.	+	0	Parent=transcript:T1
chr1	test	CDS	51	80	.	+	0	Parent=transcript:T1
`
	dir := t.TempDir()
	fa, g, _ := writeCSQFixture(t, dir, okGFF, forwardBadVCF)
	idx, err := loadCSQIndexUnified(fa, g, "", "", "", csqPhaseCheck{Force: false, Verbosity: 1})
	if err != nil {
		t.Fatalf("consistent GFF should not error: %v", err)
	}
	if len(idx.Transcripts) != 1 {
		t.Fatalf("expected exactly one transcript, got %d", len(idx.Transcripts))
	}
	for _, tx := range idx.Transcripts {
		if tx.ID != "T1" || len(tx.CDSExons) != 2 {
			t.Fatalf("consistent transcript should keep both CDS exons, got %#v", tx)
		}
	}
}
