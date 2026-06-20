package cram

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// CRAM v4.0 parity. CRAM v4 keeps the v3 container/block/slice structure
// but encodes every variable-length integer as a uint7 varint and uses a
// new data-series codec set (VARINT_UNSIGNED/SIGNED, CONST_INT/BYTE,
// EXTERNAL-for-bytes-only). These tests prove our v4 decode reproduces
// `samtools view` byte-for-byte.
//
// The fixtures are produced at test time by the vendored upstream samtools
// (`samtools view -C --output-fmt-option version=4.0`). The upstream in
// this tree CAN write CRAM 4.0 (it prints a "still a draft" warning but
// emits a valid file), so the oracle is live, not canned. As with the
// other parity tests, a missing/un-buildable upstream is a hard failure,
// never a skip.

// writeV4CRAM runs `samtools view -C --output-fmt-option version=4.0` to
// produce a CRAM v4.0 file at dst from src, optionally with an external
// reference FASTA. It fails the test if the produced file is not actually
// CRAM major version 4 — proving the upstream really emitted v4.0 rather
// than silently falling back.
func writeV4CRAM(t *testing.T, samtools, src, dst, refFA string) {
	t.Helper()
	args := []string{"view", "-C", "--output-fmt-option", "version=4.0"}
	if refFA != "" {
		args = append(args, "-T", refFA)
	}
	args = append(args, "-o", dst, src)
	cmd := exec.Command(samtools, args...)
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("samtools %s: %v\n%s", strings.Join(args, " "), err, out)
	}
	rd, err := Open(dst)
	if err != nil {
		t.Fatalf("Open %s: %v", dst, err)
	}
	defer rd.Close()
	if got := rd.FileDefinition().Major; got != 4 {
		t.Fatalf("upstream wrote CRAM major version %d, not 4 — this tree's samtools cannot write v4.0; %s",
			got, strings.TrimSpace(string(out)))
	}
}

// TestV4UpstreamWritesV4 documents and locks in that the vendored upstream
// samtools+htslib can write CRAM v4.0, which the rest of the v4 parity
// suite depends on. If a future submodule bump drops v4 write support this
// fails loudly here rather than mysteriously elsewhere.
func TestV4UpstreamWritesV4(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	srcSAM := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.sam")
	cramPath := filepath.Join(t.TempDir(), "probe.v4.cram")
	writeV4CRAM(t, samtools, srcSAM, cramPath, "")

	rd, err := Open(cramPath)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer rd.Close()
	if v := rd.FileDefinition().VersionString(); v != "4.0" {
		t.Fatalf("fixture version = %q, want 4.0", v)
	}
}

// TestV4EmbeddedReferenceParity is the headline v4 decode test: it encodes
// a known SAM to CRAM v4.0 (samtools enables embed_ref because the @SQ
// lines carry no M5 and no external reference is supplied) and asserts our
// decode matches `samtools view` field-by-field (QNAME/FLAG/RNAME/POS/
// MAPQ/CIGAR/RNEXT/PNEXT/TLEN/SEQ/QUAL/tags). The dat/test_input_1_a.sam
// fixture exercises the v4 deltas end-to-end: uint7 container/block/slice
// headers, the VARINT/CONST codecs, the v4 RG* tag-dictionary placeholder,
// the v4 read-name dedup, and the QO quality-orientation flag.
func TestV4EmbeddedReferenceParity(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	srcSAM := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.sam")
	if _, err := os.Stat(srcSAM); err != nil {
		t.Skipf("source SAM fixture missing: %v", err)
	}
	cramPath := filepath.Join(t.TempDir(), "embed.v4.cram")
	writeV4CRAM(t, samtools, srcSAM, cramPath, "")
	assertVersionAndEmbeddedRef(t, cramPath, "4.0")

	want := samtoolsViewRecords(t, samtools, cramPath)
	got := ourViewRecords(t, cramPath)
	if len(got) != len(want) {
		t.Fatalf("decoded %d records, samtools decoded %d", len(got), len(want))
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("v4.0 record %d mismatch:\n got=%q\nwant=%q", i, got[i], want[i])
		}
	}
}

// TestV4UsesV4Codecs confirms the v4 fixture's compression header really
// declares the v4-only codecs (VARINT_UNSIGNED/SIGNED, CONST_INT), so the
// parity above exercises the new decode paths rather than passing on a
// file that happened to avoid them.
func TestV4UsesV4Codecs(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	srcSAM := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.sam")
	cramPath := filepath.Join(t.TempDir(), "codecs.v4.cram")
	writeV4CRAM(t, samtools, srcSAM, cramPath, "")

	rd, err := Open(cramPath)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer rd.Close()
	conts, err := rd.Containers()
	if err != nil {
		t.Fatalf("Containers: %v", err)
	}
	seen := map[EncodingID]bool{}
	for _, c := range conts {
		dc, derr := ParseDataContainer(c)
		if derr != nil {
			continue
		}
		for _, enc := range dc.Compression.DataSeries {
			seen[enc.ID] = true
		}
	}
	for _, want := range []EncodingID{EncodingVarintUnsigned, EncodingConstInt} {
		if !seen[want] {
			t.Errorf("v4 fixture declared no %s data series — that v4 codec path is untested", want)
		}
	}
}

// TestV4ExternalReferenceParity exercises the v4 decode against an
// external reference: a self-consistent reference and SAM are written,
// encoded to CRAM v4.0 with -T, and decoded with SetReferenceFASTA. The
// result is compared to the v3.0 decode of the same source, so MD/NM
// regeneration, substitution-feature reconstruction and the
// reference-resolve path are all asserted identical across versions. (A
// raw `samtools view` of either CRAM omits MD/NM because the file does not
// store them; both our decodes regenerate them from the reference, which
// is the behaviour under test.)
func TestV4ExternalReferenceParity(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	dir := t.TempDir()

	faPath := filepath.Join(dir, "ref.fa")
	// A 104bp reference (two 52bp lines) so a record can substitute,
	// insert and span it cleanly.
	if err := os.WriteFile(faPath,
		[]byte(">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"),
		0o644); err != nil {
		t.Fatalf("writing reference: %v", err)
	}
	// Compute the contig MD5 with samtools so the @SQ M5 matches (samtools
	// validates it before embedding/encoding).
	dictOut, err := exec.Command(samtools, "dict", faPath).Output()
	if err != nil {
		t.Fatalf("samtools dict: %v", err)
	}
	m5 := ""
	for _, line := range strings.Split(string(dictOut), "\n") {
		if strings.HasPrefix(line, "@SQ") {
			for _, f := range strings.Split(line, "\t") {
				if strings.HasPrefix(f, "M5:") {
					m5 = strings.TrimPrefix(f, "M5:")
				}
			}
		}
	}
	if m5 == "" {
		t.Fatalf("samtools dict produced no M5 for the reference")
	}

	samPath := filepath.Join(dir, "reads.sam")
	sam := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:104\tM5:" + m5 + "\n" +
		"@RG\tID:rg1\tSM:sample1\n" +
		// read3 substitutes a T where the reference has A at position 5
		// (1-based 8), so MD/NM regeneration is non-trivial.
		"read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1\n" +
		"read2\t16\tchr1\t11\t60\t10M\t*\t0\t0\tGTACGTACGT\tIIIIIIIIII\tRG:Z:rg1\n" +
		"read3\t0\tchr1\t5\t60\t8M\t*\t0\t0\tACGTTCGT\tIIIIIIII\tRG:Z:rg1\n"
	if err := os.WriteFile(samPath, []byte(sam), 0o644); err != nil {
		t.Fatalf("writing SAM: %v", err)
	}

	v4Path := filepath.Join(dir, "ext.v4.cram")
	writeV4CRAM(t, samtools, samPath, v4Path, faPath)

	v3Path := filepath.Join(dir, "ext.v3.cram")
	cmd := exec.Command(samtools, "view", "-C", "-T", faPath,
		"--output-fmt-option", "version=3.0", "-o", v3Path, samPath)
	if out, cerr := cmd.CombinedOutput(); cerr != nil {
		t.Fatalf("samtools view -C version=3.0: %v\n%s", cerr, out)
	}

	v4Recs := decodeWithReferenceFASTA(t, v4Path, faPath)
	v3Recs := decodeWithReferenceFASTA(t, v3Path, faPath)
	if len(v4Recs) != len(v3Recs) {
		t.Fatalf("v4 decoded %d records, v3 decoded %d", len(v4Recs), len(v3Recs))
	}
	for i := range v3Recs {
		if v4Recs[i] != v3Recs[i] {
			t.Fatalf("ext-ref record %d differs between v4 and v3:\nv4=%q\nv3=%q", i, v4Recs[i], v3Recs[i])
		}
	}
	// Sanity: the substitution must have produced a non-trivial MD/NM, so
	// the comparison is not vacuously matching two empty tag sets.
	sawSubst := false
	for _, r := range v4Recs {
		if strings.Contains(r, "MD:Z:4A3") && strings.Contains(r, "NM:i:1") {
			sawSubst = true
		}
	}
	if !sawSubst {
		t.Fatalf("v4 ext-ref decode did not regenerate the expected MD:Z:4A3 / NM:i:1 for the substituted read; records=%v", v4Recs)
	}
}

// decodeWithReferenceFASTA decodes path with the given reference FASTA
// attached and returns the alignment record lines (no header).
func decodeWithReferenceFASTA(t *testing.T, path, faPath string) []string {
	t.Helper()
	rr, err := OpenRecords(path)
	if err != nil {
		t.Fatalf("OpenRecords %s: %v", path, err)
	}
	defer rr.Close()
	if err := rr.SetReferenceFASTA(faPath); err != nil {
		t.Fatalf("SetReferenceFASTA %s: %v", faPath, err)
	}
	var recs []string
	for _, line := range ourSAMBody(t, rr) {
		recs = append(recs, line)
	}
	return recs
}

// ourSAMBody renders rr's records to SAM and returns the non-header lines.
func ourSAMBody(t *testing.T, rr *RecordReader) []string {
	t.Helper()
	var buf strings.Builder
	if err := rr.WriteSAM(&buf); err != nil {
		t.Fatalf("WriteSAM: %v", err)
	}
	var out []string
	for _, line := range strings.Split(buf.String(), "\n") {
		if line != "" && !strings.HasPrefix(line, "@") {
			out = append(out, line)
		}
	}
	return out
}
