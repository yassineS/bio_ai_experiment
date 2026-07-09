package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// writeFastaWithFai writes a plain FASTA plus its sibling .fai to dir and
// returns the FASTA path. The line width is fixed at 6 bases per line so the
// index geometry is exercised (multi-line contigs, a short final line).
func writeFastaWithFai(t *testing.T, dir string, recs []*fasta.Record) string {
	t.Helper()
	faPath := filepath.Join(dir, "ref.fa")
	var b bytes.Buffer
	const wrap = 6
	for _, r := range recs {
		b.WriteString(">" + r.ID + "\n")
		seq := r.Sequence
		for i := 0; i < len(seq); i += wrap {
			j := i + wrap
			if j > len(seq) {
				j = len(seq)
			}
			b.Write(seq[i:j])
			b.WriteByte('\n')
		}
	}
	if err := os.WriteFile(faPath, b.Bytes(), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	idx, err := fasta.BuildIndex(faPath)
	if err != nil {
		t.Fatalf("build index: %v", err)
	}
	if err := idx.Save(faPath + ".fai"); err != nil {
		t.Fatalf("save fai: %v", err)
	}
	return faPath
}

// TestConsensusFileStreamMatchesInMemory proves the faidx-streaming
// ConsensusFile path is byte-identical to the in-memory Consensus path across
// multiple contigs, soft-masked (lowercase) reference bases, a SNP, an
// insertion and a deletion. The streaming path holds only one contig in memory
// at a time; this guards against any divergence introduced by that refactor.
func TestConsensusFileStreamMatchesInMemory(t *testing.T) {
	recs := []*fasta.Record{
		// Mixed-case to verify FetchRaw preserves soft-masking exactly.
		{ID: "chr1", Sequence: []byte("AAAAccccGGGGttttAAAA")},
		{ID: "chr2", Sequence: []byte("TTTTggggCCCCaaaa")},
	}
	vcfText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=20>\n" +
		"##contig=<ID=chr2,length=16>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t1\t.\tA\tG\t.\tPASS\t.\n" + // SNP at start
		"chr1\t9\t.\tG\tGTT\t.\tPASS\t.\n" + // insertion
		"chr1\t13\t.\ttttt\tt\t.\tPASS\t.\n" + // deletion (lowercase ref)
		"chr2\t5\t.\tg\tC\t.\tPASS\t.\n" // SNP into lowercase span

	// In-memory reference path.
	var memOut bytes.Buffer
	memN, err := Consensus(strings.NewReader(vcfText), &memOut, ConsensusOptions{
		Reference: recs,
		LineWidth: 60,
	})
	if err != nil {
		t.Fatalf("Consensus (in-memory): %v", err)
	}

	// Streaming file path.
	dir := t.TempDir()
	faPath := writeFastaWithFai(t, dir, recs)
	vcfPath := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(vcfPath, []byte(vcfText), 0o644); err != nil {
		t.Fatalf("write vcf: %v", err)
	}
	var fileOut bytes.Buffer
	fileN, err := ConsensusFile(vcfPath, faPath, &fileOut, ConsensusOptions{
		LineWidth: 60,
	})
	if err != nil {
		t.Fatalf("ConsensusFile (streaming): %v", err)
	}

	if memN != fileN {
		t.Errorf("applied-count mismatch: in-memory %d, streaming %d", memN, fileN)
	}
	if memOut.String() != fileOut.String() {
		t.Errorf("streaming output diverged from in-memory:\n in-memory: %q\n streaming: %q",
			memOut.String(), fileOut.String())
	}
}

// TestConsensusFileStreamMask exercises the streaming path with a BED mask and
// a prefix, again requiring byte-identity with the in-memory path. The mask is
// applied in place on the single fetched contig buffer.
func TestConsensusFileStreamMask(t *testing.T) {
	recs := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("ACGTACGTACGTACGT")},
	}
	vcfText := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=16>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t2\t.\tC\tT\t.\tPASS\t.\n"
	mask := []MaskRegion{{Chrom: "chr1", Beg: 5, End: 8}} // 1-based inclusive

	opts := ConsensusOptions{
		Mask:      mask,
		MaskWith:  MarkSpec{Mode: MarkNone}, // default 'N'
		Prefix:    "px_",
		LineWidth: 60,
	}

	var memOut bytes.Buffer
	memOpts := opts
	memOpts.Reference = recs
	if _, err := Consensus(strings.NewReader(vcfText), &memOut, memOpts); err != nil {
		t.Fatalf("Consensus (in-memory): %v", err)
	}

	dir := t.TempDir()
	faPath := writeFastaWithFai(t, dir, recs)
	vcfPath := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(vcfPath, []byte(vcfText), 0o644); err != nil {
		t.Fatalf("write vcf: %v", err)
	}
	var fileOut bytes.Buffer
	if _, err := ConsensusFile(vcfPath, faPath, &fileOut, opts); err != nil {
		t.Fatalf("ConsensusFile (streaming): %v", err)
	}

	if memOut.String() != fileOut.String() {
		t.Errorf("streaming(mask) diverged:\n in-memory: %q\n streaming: %q",
			memOut.String(), fileOut.String())
	}
}
