package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bam"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// writeIndexedTestBAM builds a coordinate-sorted BAM on disk from SAM text
// and writes a .bai sibling using the in-tree index builder, so MpileupFile's
// indexed -r fast path (openBAMRegionReader) engages without depending on the
// upstream samtools binary. It returns the BAM path.
func writeIndexedTestBAM(t *testing.T, dir, samText string) string {
	t.Helper()
	// Parse the SAM text into a header + records.
	rd, err := sam.NewReader(strings.NewReader(samText))
	if err != nil {
		t.Fatalf("sam.NewReader: %v", err)
	}
	hdr := rd.Header()
	bamPath := filepath.Join(dir, "in.bam")
	bf, err := os.Create(bamPath)
	if err != nil {
		t.Fatalf("create bam: %v", err)
	}
	bw := sam.NewBAMWriter(bf)
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for {
		rec, rerr := rd.Read()
		if rerr != nil {
			break
		}
		if err := bw.Write(rec); err != nil {
			t.Fatalf("BAM write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("BAM close: %v", err)
	}
	if err := bf.Close(); err != nil {
		t.Fatalf("bam file close: %v", err)
	}

	// Build the .bai sibling from the on-disk BAM.
	rbf, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("reopen bam: %v", err)
	}
	defer rbf.Close()
	br, err := sam.NewBAMReader(rbf)
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	idx, err := bam.BuildBAI(br, len(br.Header().Refs))
	if err != nil {
		t.Fatalf("BuildBAI: %v", err)
	}
	var idxBuf bytes.Buffer
	if err := bam.WriteBAI(&idxBuf, idx); err != nil {
		t.Fatalf("WriteBAI: %v", err)
	}
	if err := os.WriteFile(bamPath+".bai", idxBuf.Bytes(), 0o644); err != nil {
		t.Fatalf("write .bai: %v", err)
	}
	return bamPath
}

// mpileupDataRecords runs MpileupFile and returns the non-header VCF lines.
func mpileupDataRecords(t *testing.T, opts MpileupOptions) []string {
	t.Helper()
	var buf bytes.Buffer
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}
	var data []string
	for _, ln := range strings.Split(buf.String(), "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		data = append(data, ln)
	}
	return data
}

// TestMpileupIndexedRegionMatchesLinear is the regression guard for task #47:
// a `-r` query that takes the indexed BAI/CSI seek path (openBAMRegionReader)
// must yield byte-identical VCF records to the linear whole-file scan, on the
// same records and the same region. The indexed reader keeps left-overlapping
// reads (reads that start before the region but span into it), so the two
// paths agree column-for-column inside the requested window. This locks the
// seek fast path to the linear semantics that the upstream-parity tests pin.
func TestMpileupIndexedRegionMatchesLinear(t *testing.T) {
	dir := t.TempDir()

	// 200bp of A on chr1; .fai sidecar for OpenRandomAccess.
	famPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(famPath, []byte(">chr1\n"+strings.Repeat("A", 200)+"\n"), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	if err := os.WriteFile(famPath+".fai", []byte("chr1\t200\t6\t200\t201\n"), 0o644); err != nil {
		t.Fatalf("write fai: %v", err)
	}

	// Reads scattered across chr1, coordinate-sorted. Two reads (r3, r4)
	// carry a C at column 55 (0-based 54) so there is a variant-ish site
	// inside the target window. r_left starts at pos 48 and spans into the
	// window [50,60], exercising the left-overlap invariant.
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:200",
		"@RG\tID:rg1\tSM:sample1",
		"r_left\t0\tchr1\t48\t60\t20M\t*\t0\t0\t" + strings.Repeat("A", 20) + "\t" + strings.Repeat("?", 20) + "\tRG:Z:rg1",
		"r1\t0\tchr1\t50\t60\t20M\t*\t0\t0\t" + strings.Repeat("A", 20) + "\t" + strings.Repeat("?", 20) + "\tRG:Z:rg1",
		"r2\t0\tchr1\t50\t60\t20M\t*\t0\t0\t" + strings.Repeat("A", 20) + "\t" + strings.Repeat("?", 20) + "\tRG:Z:rg1",
		"r3\t0\tchr1\t50\t60\t20M\t*\t0\t0\tAAAACAAAAAAAAAAAAAAA\t" + strings.Repeat("?", 20) + "\tRG:Z:rg1",
		"r4\t0\tchr1\t50\t60\t20M\t*\t0\t0\tAAAACAAAAAAAAAAAAAAA\t" + strings.Repeat("?", 20) + "\tRG:Z:rg1",
		"r_far\t0\tchr1\t150\t60\t20M\t*\t0\t0\t" + strings.Repeat("A", 20) + "\t" + strings.Repeat("?", 20) + "\tRG:Z:rg1",
		"",
	}, "\n")

	bamPath := writeIndexedTestBAM(t, dir, samText)
	region := "chr1:50-60"

	// Indexed path: the .bai sibling exists, so openBAMRegionReader seeks.
	indexed := mpileupDataRecords(t, MpileupOptions{
		Inputs:   []string{bamPath},
		FastaRef: famPath,
		Regions:  []string{region},
		NoBAQ:    true,
	})

	// Linear path: same records, no index sidecar, so the whole-file scan
	// runs and the region is applied as a post-filter.
	linDir := t.TempDir()
	linBAM := filepath.Join(linDir, "lin.bam")
	if err := copyFileForTest(bamPath, linBAM); err != nil {
		t.Fatalf("copy bam: %v", err)
	}
	linear := mpileupDataRecords(t, MpileupOptions{
		Inputs:   []string{linBAM},
		FastaRef: famPath,
		Regions:  []string{region},
		NoBAQ:    true,
	})

	if len(indexed) == 0 {
		t.Fatalf("indexed -r produced no records for %s", region)
	}
	if len(indexed) != len(linear) {
		t.Fatalf("record count differs: indexed=%d linear=%d", len(indexed), len(linear))
	}
	for i := range indexed {
		if indexed[i] != linear[i] {
			t.Fatalf("record %d differs:\n indexed: %s\n linear : %s", i, indexed[i], linear[i])
		}
	}
}

// copyFileForTest copies src to dst byte-for-byte.
func copyFileForTest(src, dst string) error {
	b, err := os.ReadFile(src)
	if err != nil {
		return err
	}
	return os.WriteFile(dst, b, 0o644)
}
