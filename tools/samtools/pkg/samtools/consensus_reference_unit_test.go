package samtools

import (
	"bufio"
	"bytes"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// writeTestFASTA writes a plain FASTA with its sibling .fai index into a temp
// dir and returns the FASTA path. It needs no upstream binary.
func writeTestFASTA(t *testing.T, body string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(path, []byte(body), 0o600); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	idx, err := fasta.BuildIndex(path)
	if err != nil {
		t.Fatalf("build fai: %v", err)
	}
	if err := idx.Save(path + ".fai"); err != nil {
		t.Fatalf("save fai: %v", err)
	}
	return path
}

// TestUnitConsensusRefBase verifies the per-position reference lookup
// (consensusRef.base): correct 0-based indexing, original byte-case
// preservation for a plain FASTA, and the 'N' fallback for out-of-range
// coordinates and contigs absent from the reference. Binary-free.
func TestUnitConsensusRefBase(t *testing.T) {
	path := writeTestFASTA(t, ">chr1\nACGTacgtNN\n>chr2\nTTTT\n")
	ref, err := loadConsensusRef(path)
	if err != nil {
		t.Fatalf("loadConsensusRef: %v", err)
	}
	defer ref.close()

	want := "ACGTacgtNN"
	for i := 0; i < len(want); i++ {
		if got := ref.base("chr1", i); got != want[i] {
			t.Errorf("chr1[%d] = %q, want %q (case must be preserved)", i, got, want[i])
		}
	}
	if got := ref.base("chr2", 0); got != 'T' {
		t.Errorf("chr2[0] = %q, want 'T'", got)
	}
	// Out-of-range and unknown-contig coordinates fall back to 'N'
	// (upstream update_ref returns <0 -> caller emits 'N').
	if got := ref.base("chr1", 10); got != 'N' {
		t.Errorf("chr1[10] (past end) = %q, want 'N'", got)
	}
	if got := ref.base("chr1", -1); got != 'N' {
		t.Errorf("chr1[-1] = %q, want 'N'", got)
	}
	if got := ref.base("chrX", 0); got != 'N' {
		t.Errorf("chrX[0] (unknown) = %q, want 'N'", got)
	}
}

// TestUnitWriteEmptyPileupRowsRefSubstitution pins the empty-row reference
// substitution in pileup output (writeEmptyPileupRows): with no reference the
// no-coverage call is 'N'; with a reference it is the reference base at that
// position; a reference position beyond the loaded contig falls back to 'N'.
// The depth (0), quality (0) and seq/qual ('*') columns are unchanged either
// way. Binary-free.
func TestUnitWriteEmptyPileupRowsRefSubstitution(t *testing.T) {
	// Reference base mapping (1-based): 1=A 2=C 3=G 4=T 5=A.
	path := writeTestFASTA(t, ">chr1\nACGTA\n")
	ref, err := loadConsensusRef(path)
	if err != nil {
		t.Fatalf("loadConsensusRef: %v", err)
	}
	defer ref.close()

	render := func(r *consensusRef, start, end int) string {
		var buf bytes.Buffer
		bw := bufio.NewWriter(&buf)
		if err := writeEmptyPileupRows(bw, "chr1", start, end, nil, r); err != nil {
			t.Fatalf("writeEmptyPileupRows: %v", err)
		}
		bw.Flush()
		return buf.String()
	}

	// No reference: every no-coverage call is 'N'.
	gotN := render(nil, 1, 3)
	wantN := "chr1\t1\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t2\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t3\t0\t0\tN\t0\t*\t*\n"
	if gotN != wantN {
		t.Errorf("no-ref empty rows:\n got %q\nwant %q", gotN, wantN)
	}

	// With reference: the call is the reference base (ACG at 1..3).
	gotR := render(ref, 1, 3)
	wantR := "chr1\t1\t0\t0\tA\t0\t*\t*\n" +
		"chr1\t2\t0\t0\tC\t0\t*\t*\n" +
		"chr1\t3\t0\t0\tG\t0\t*\t*\n"
	if gotR != wantR {
		t.Errorf("ref empty rows:\n got %q\nwant %q", gotR, wantR)
	}

	// Position 6 is past the 5bp contig: fall back to 'N'.
	gotPast := render(ref, 6, 6)
	wantPast := "chr1\t6\t0\t0\tN\t0\t*\t*\n"
	if gotPast != wantPast {
		t.Errorf("ref past-end row:\n got %q\nwant %q", gotPast, wantPast)
	}
}
