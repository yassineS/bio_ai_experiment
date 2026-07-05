package runner

import (
	"bytes"
	"compress/gzip"
	"crypto/md5"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/matrix"
)

// TestStripProvenanceSAM checks @PG/@CO removal without touching data lines.
func TestStripProvenanceSAM(t *testing.T) {
	in := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@PG\tID:samtools\tPN:samtools\tVN:1.22\n" +
		"read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\tIIII\n"
	want := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n" +
		"read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\tIIII\n"
	if got := string(stripProvenance([]byte(in))); got != want {
		t.Errorf("stripProvenance SAM:\n got=%q\nwant=%q", got, want)
	}
}

// TestStripProvenanceDictUR checks that the machine-specific "UR:" reference-URI
// field is stripped from @SQ lines (as "samtools dict" emits) while every other
// @SQ field is preserved, so two dict outputs that differ ONLY by their local
// file:// path compare equal but a genuine field difference still diverges.
func TestStripProvenanceDictUR(t *testing.T) {
	ours := "@HD\tVN:1.0\tSO:unsorted\n" +
		"@SQ\tSN:chr1\tLN:100\tM5:abc\tUR:file:///home/ours/ref.fa\tAS:GRCh38\n"
	up := "@HD\tVN:1.0\tSO:unsorted\n" +
		"@SQ\tSN:chr1\tLN:100\tM5:abc\tUR:file:///opt/upstream/ref.fa\tAS:GRCh38\n"
	if a, b := string(stripProvenance([]byte(ours))), string(stripProvenance([]byte(up))); a != b {
		t.Errorf("UR-only diff should compare equal after strip:\n ours=%q\n  up =%q", a, b)
	}
	want := "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\tM5:abc\tAS:GRCh38\n"
	if got := string(stripProvenance([]byte(ours))); got != want {
		t.Errorf("dict UR strip:\n got=%q\nwant=%q", got, want)
	}
	// A genuine field difference (different M5) must still diverge.
	bad := "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:100\tM5:DIFFERENT\tUR:file:///x\tAS:GRCh38\n"
	if string(stripProvenance([]byte(ours))) == string(stripProvenance([]byte(bad))) {
		t.Error("genuine M5 difference must not be masked by UR stripping")
	}
	// The streaming filter must match the batch path byte-for-byte.
	var buf bytes.Buffer
	pf := newProvenanceFilter(&buf)
	if _, err := io.Copy(pf, bytes.NewReader([]byte(ours))); err != nil {
		t.Fatal(err)
	}
	if err := pf.Close(); err != nil {
		t.Fatal(err)
	}
	if got := buf.String(); got != string(stripProvenance([]byte(ours))) {
		t.Errorf("streaming vs batch drift:\n stream=%q\n  batch=%q", got, string(stripProvenance([]byte(ours))))
	}
}

// TestStripProvenanceVCF checks command/version header removal.
func TestStripProvenanceVCF(t *testing.T) {
	in := "##fileformat=VCFv4.2\n##bcftools_viewCommand=view a.vcf; Date=...\n##contig=<ID=chr1>\nchr1\t1\t.\tA\tG\t60\tPASS\t.\n"
	got := string(stripProvenance([]byte(in)))
	if want := "##fileformat=VCFv4.2\n##contig=<ID=chr1>\nchr1\t1\t.\tA\tG\t60\tPASS\t.\n"; got != want {
		t.Errorf("stripProvenance VCF:\n got=%q\nwant=%q", got, want)
	}
}

// TestStripProvenanceStatsBlock checks the samtools/bcftools stats-style
// provenance comment block ("# This file was produced by ...", the command-line
// echo, the working-directory block, the bare "#" separator, and the gtcheck
// timing line) is removed while the data-describing comment rows survive.
func TestStripProvenanceStatsBlock(t *testing.T) {
	in := "# This file was produced by bcftools stats (1.23+htslib-1.23)\n" +
		"# The command line was:\tbcftools stats a.vcf\n" +
		"#\n" +
		"# ID\t[2]id\t[3]file names\n" +
		"SN\t0\tnumber of records:\t400\n" +
		"INFO\tTime required to process one record .. 0.000003 seconds\n"
	want := "# ID\t[2]id\t[3]file names\n" +
		"SN\t0\tnumber of records:\t400\n"
	if got := string(stripProvenance([]byte(in))); got != want {
		t.Errorf("stripProvenance stats block:\n got=%q\nwant=%q", got, want)
	}
}

// TestStripProvenanceFilterPass checks the auto-inserted ##FILTER=PASS
// boilerplate is dropped (its position differs between ours/upstream) while a
// real ##FILTER definition is preserved.
func TestStripProvenanceFilterPass(t *testing.T) {
	in := "##fileformat=VCFv4.2\n" +
		"##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
		"##FILTER=<ID=q10,Description=\"Quality below 10\">\n" +
		"chr1\t1\t.\tA\tG\t60\tPASS\t.\n"
	want := "##fileformat=VCFv4.2\n" +
		"##FILTER=<ID=q10,Description=\"Quality below 10\">\n" +
		"chr1\t1\t.\tA\tG\t60\tPASS\t.\n"
	if got := string(stripProvenance([]byte(in))); got != want {
		t.Errorf("stripProvenance FILTER=PASS:\n got=%q\nwant=%q", got, want)
	}
}

// provenanceTestInputs is the shared corpus exercising every provenance pattern
// plus the trailing-newline / chunk-boundary edge cases. It is the table the
// streaming-equivalence test iterates.
func provenanceTestInputs() map[string][]byte {
	return map[string][]byte{
		"empty":                  []byte(""),
		"single_newline":         []byte("\n"),
		"only_data_no_trailing":  []byte("chr1\t1\tA"),
		"only_data_trailing":     []byte("chr1\t1\tA\n"),
		"only_provenance_pg":     []byte("@PG\tID:x\n"),
		"only_provenance_no_nl":  []byte("@PG\tID:x"),
		"two_provenance":         []byte("@PG\tID:x\n@CO\tcomment\n"),
		"prov_then_data":         []byte("@PG\tID:x\nchr1\t1\tA\n"),
		"data_then_prov":         []byte("chr1\t1\tA\n@PG\tID:x\n"),
		"interleaved":            []byte("@HD\tVN:1.6\n@PG\tID:x\n@SQ\tSN:chr1\nr1\t0\tchr1\n@CO\tc\nr2\t0\tchr1\n"),
		"vcf_command":            []byte("##fileformat=VCFv4.2\n##bcftools_viewCommand=view a.vcf; Date=x\n##contig=<ID=chr1>\nchr1\t1\t.\tA\tG\n"),
		"vcf_source_filedate":    []byte("##source=myTool\n##fileDate=20200101\n##reference=hs37\n##contig=<ID=1>\n1\t1\t.\tA\tG\n"),
		"filter_pass":            []byte("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##FILTER=<ID=q10,Description=\"q\">\nchr1\t1\n"),
		"stats_banner":           []byte("# This file was produced by bcftools stats (1.23)\n# The command line was:\tbcftools stats a.vcf\n# \t/work\n#\n# ID\t[2]id\nSN\t0\t400\n"),
		"stats_bare_hash":        []byte("#\n# CHK, Checksum\ndata\n"),
		"stats_tab_echo":         []byte("# \tbcftools stats a.vcf\n# \t/some/dir\nSN\t0\t1\n"),
		"gtcheck_timing":         []byte("INFO\tTime required to process one record .. 0.000003 seconds\nDC\tsample\t0.1\n"),
		"working_dir_banner":     []byte("# and the working directory was:\t/work\nSN\t0\t1\n"),
		"contains_stats_banner2": []byte("# This file contains statistics for all reads.\nSN\t0\t1\n"),
		"trailing_blank_lines":   []byte("chr1\t1\n\n\n@PG\tID:x\n\n"),
		"prov_only_then_nl":      []byte("@PG\ta\n@CO\tb\n\n"),
		"crlf_like":              []byte("chr1\t1\r\n@PG\tx\r\nchr2\t2\r\n"),
		"long_line":              append(append([]byte("@PG\t"), bytes.Repeat([]byte("X"), 200000)...), '\n'),
	}
}

// TestStreamDigestEquivalence is the GATE: for every input in the corpus, the
// streaming provenance filter (StreamDigest) must produce the SAME md5 as
// md5(stripProvenance(b)) — i.e. the streaming normalization is byte-for-byte
// identical to what CompareByteExact compares. It also asserts that feeding the
// bytes at arbitrary chunk boundaries (1-byte, 3-byte, and whole) yields the
// same digest, proving the partial-line buffering is correct.
func TestStreamDigestEquivalence(t *testing.T) {
	for name, b := range provenanceTestInputs() {
		want := md5.Sum(stripProvenance(b))

		// Whole-slice path.
		gotSum, _, err := StreamDigest(bytes.NewReader(b))
		if err != nil {
			t.Fatalf("%s: StreamDigest error: %v", name, err)
		}
		if gotSum != want {
			t.Errorf("%s: whole StreamDigest md5=%x want md5(stripProvenance)=%x\n  stripped=%q",
				name, gotSum, want, stripProvenance(b))
		}

		// Arbitrary chunk boundaries: 1-byte, 3-byte, and whole reader.
		for _, chunk := range []int{1, 3, 7, len(b) + 1} {
			if chunk < 1 {
				chunk = 1
			}
			sum, _, err := StreamDigest(&chunkReader{b: b, chunk: chunk})
			if err != nil {
				t.Fatalf("%s chunk=%d: StreamDigest error: %v", name, chunk, err)
			}
			if sum != want {
				t.Errorf("%s chunk=%d: md5=%x want %x", name, chunk, sum, want)
			}
		}

		// Also feed via repeated direct Write calls at byte-by-byte boundaries
		// straight into the filter, mirroring how a pipe delivers bytes.
		for _, chunk := range []int{1, 2, 5} {
			h := md5.New()
			pf := newProvenanceFilter(h)
			for off := 0; off < len(b); off += chunk {
				end := off + chunk
				if end > len(b) {
					end = len(b)
				}
				if _, err := pf.Write(b[off:end]); err != nil {
					t.Fatalf("%s write chunk=%d: %v", name, chunk, err)
				}
			}
			if err := pf.Close(); err != nil {
				t.Fatalf("%s close chunk=%d: %v", name, chunk, err)
			}
			var got [md5.Size]byte
			copy(got[:], h.Sum(nil))
			if got != want {
				t.Errorf("%s direct-write chunk=%d: md5=%x want %x", name, chunk, got, want)
			}
		}
	}
}

// TestStreamDigestHead checks the captured head equals the first 64 KiB of the
// provenance-stripped output and is itself capped.
func TestStreamDigestHead(t *testing.T) {
	// Data that exceeds the head window after stripping.
	var sb bytes.Buffer
	sb.WriteString("@PG\tID:x\n") // dropped
	for i := 0; i < 5000; i++ {
		sb.WriteString("chr1\t1\tA\tdata\n")
	}
	b := sb.Bytes()
	_, head, err := StreamDigest(bytes.NewReader(b))
	if err != nil {
		t.Fatal(err)
	}
	stripped := stripProvenance(b)
	wantHeadLen := streamHeadCap
	if len(stripped) < wantHeadLen {
		wantHeadLen = len(stripped)
	}
	if len(head) != wantHeadLen {
		t.Fatalf("head len=%d want %d (cap %d, stripped %d)", len(head), wantHeadLen, streamHeadCap, len(stripped))
	}
	if !bytes.Equal(head, stripped[:wantHeadLen]) {
		t.Errorf("head does not match prefix of stripProvenance output")
	}
}

// chunkReader yields b in fixed-size chunks to exercise partial-line buffering
// across Read boundaries.
type chunkReader struct {
	b     []byte
	chunk int
	off   int
}

func (c *chunkReader) Read(p []byte) (int, error) {
	if c.off >= len(c.b) {
		return 0, io.EOF
	}
	n := c.chunk
	if n > len(p) {
		n = len(p)
	}
	if c.off+n > len(c.b) {
		n = len(c.b) - c.off
	}
	copy(p, c.b[c.off:c.off+n])
	c.off += n
	return n, nil
}

// TestCompareByteExact covers match and mismatch.
func TestCompareByteExact(t *testing.T) {
	a := []byte("@PG\tID:x\nchr1\t1\n")
	b := []byte("@PG\tID:y\nchr1\t1\n") // differs only in stripped @PG
	if r := CompareByteExact(a, b); !r.Equal {
		t.Errorf("expected equal after provenance strip, got %+v", r)
	}
	c := []byte("chr1\t2\n")
	if r := CompareByteExact(a, c); r.Equal {
		t.Errorf("expected mismatch, got equal")
	}
}

// TestCompareSimilarity covers numeric tolerance and structural mismatch.
func TestCompareSimilarity(t *testing.T) {
	a := []byte("chr1\t0.1000000\t2\n")
	b := []byte("chr1\t0.1000001\t2\n")
	r := CompareSimilarity(a, b, similarityEpsilon)
	if !r.Equal {
		t.Errorf("expected within-tolerance equal, got %+v", r)
	}
	if r.MaxDeviation == 0 {
		t.Errorf("expected non-zero recorded deviation")
	}

	d := CompareSimilarity([]byte("chr1\t1.0\n"), []byte("chr2\t1.0\n"), similarityEpsilon)
	if d.Equal {
		t.Errorf("expected non-numeric field mismatch to diverge")
	}

	e := CompareSimilarity([]byte("chr1\t1.0\n"), []byte("chr1\t2.0\n"), similarityEpsilon)
	if e.Equal {
		t.Errorf("expected out-of-tolerance numeric to diverge")
	}

	// A per-entry tolerance widens acceptance: a deviation that fails at the
	// default epsilon passes when the entry opts into a looser bound (the
	// bcftools call QUAL libm-last-ULP case).
	f := CompareSimilarity([]byte("chr1\t15.6999\n"), []byte("chr1\t15.6998\n"), resolveEpsilon(2e-5))
	if !f.Equal {
		t.Errorf("expected within widened tolerance, got %+v", f)
	}
	g := CompareSimilarity([]byte("chr1\t15.6999\n"), []byte("chr1\t15.6998\n"), resolveEpsilon(0))
	if g.Equal {
		t.Errorf("expected default tolerance to reject the QUAL last-ULP deviation")
	}
}

// TestCompareOutputFiles_ByteExact covers matching, mismatch, and gzip handling
// of the output-file comparison path used by vcftools and mosdepth.
func TestCompareOutputFiles_ByteExact(t *testing.T) {
	dir := t.TempDir()
	write := func(name, content string) {
		if err := os.WriteFile(filepath.Join(dir, name), []byte(content), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	writeGz := func(name, content string) {
		f, err := os.Create(filepath.Join(dir, name))
		if err != nil {
			t.Fatal(err)
		}
		gw := gzip.NewWriter(f)
		if _, err := gw.Write([]byte(content)); err != nil {
			t.Fatal(err)
		}
		gw.Close()
		f.Close()
	}
	// Plain text, identical.
	write("a.frq", "chr1\t1\tA\n")
	write("b.frq", "chr1\t1\tA\n")
	if r := CompareOutputFiles(filepath.Join(dir, "a"), filepath.Join(dir, "b"), []string{".frq"}, matrix.ByteExact, similarityEpsilon); !r.Equal {
		t.Errorf("identical .frq should match: %+v", r)
	}
	// Gzipped, identical payload (different framing is irrelevant after decode).
	writeGz("a.bed.gz", "chr1\t0\t100\t5\n")
	writeGz("b.bed.gz", "chr1\t0\t100\t5\n")
	if r := CompareOutputFiles(filepath.Join(dir, "a"), filepath.Join(dir, "b"), []string{".bed.gz"}, matrix.ByteExact, similarityEpsilon); !r.Equal {
		t.Errorf("identical gzip payload should match: %+v", r)
	}
	// Mismatch.
	write("a.diff", "x\n")
	write("b.diff", "y\n")
	if r := CompareOutputFiles(filepath.Join(dir, "a"), filepath.Join(dir, "b"), []string{".diff"}, matrix.ByteExact, similarityEpsilon); r.Equal {
		t.Errorf("differing files should diverge")
	}
	// Presence mismatch (one side missing).
	write("a.only", "x\n")
	if r := CompareOutputFiles(filepath.Join(dir, "a"), filepath.Join(dir, "b"), []string{".only"}, matrix.ByteExact, similarityEpsilon); r.Equal {
		t.Errorf("missing-on-one-side should diverge")
	}
}
