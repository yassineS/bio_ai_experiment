package bcftools

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// streamVCF is a small multi-record VCF used by the streaming reheader tests.
const streamVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	rs1	A	T	30	PASS	DP=10	GT	0/1	0/0
chr1	200	rs2	C	G	30	PASS	DP=20	GT	1/1	0/1
chr1	300	rs3	G	A	.	.	DP=5	GT	./.	0/0
`

// streamBodyOf returns the body bytes (everything after the #CHROM header line)
// of a VCF text block.
func streamBodyOf(t *testing.T, vcfText string) []byte {
	t.Helper()
	_, body := splitVCFHeader([]byte(vcfText))
	if len(body) == 0 {
		t.Fatalf("no body extracted from test VCF")
	}
	return body
}

// TestReheaderStreamBodyVerbatim confirms the streaming path copies the body
// records byte-for-byte and swaps the header (via -h) without decoding records.
func TestReheaderStreamBodyVerbatim(t *testing.T) {
	dir := t.TempDir()
	hpath := filepath.Join(dir, "newhdr.vcf")
	const newHdr = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##source=StreamTest
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`
	if err := os.WriteFile(hpath, []byte(newHdr), 0644); err != nil {
		t.Fatal(err)
	}

	wantBody := streamBodyOf(t, streamVCF)

	var out bytes.Buffer
	n, err := Reheader(strings.NewReader(streamVCF), &out, ReheaderOptions{HeaderFile: hpath})
	if err != nil {
		t.Fatal(err)
	}
	if n != 0 {
		t.Errorf("streaming path record count = %d, want 0 (records not materialised)", n)
	}

	gotHeader, gotBody := splitVCFHeader(out.Bytes())
	// (a) Body byte-identical to the input body.
	if !bytes.Equal(gotBody, wantBody) {
		t.Errorf("body not byte-identical:\n got %q\nwant %q", gotBody, wantBody)
	}
	// (b) Header was replaced.
	if !bytes.Contains(gotHeader, []byte("##source=StreamTest")) {
		t.Errorf("header not replaced:\n%s", gotHeader)
	}
	if bytes.Equal(gotHeader, []byte("##fileformat=VCFv4.2\n")) {
		t.Errorf("header unexpectedly empty")
	}
}

// TestReheaderStreamVCFGz confirms both the plain-VCF and VCF.gz streaming paths
// produce a byte-identical body. The .vcf.gz input is built via pkg/htsgo/bgzf.
func TestReheaderStreamVCFGz(t *testing.T) {
	dir := t.TempDir()
	hpath := filepath.Join(dir, "newhdr.vcf")
	const newHdr = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##source=StreamTestGz
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`
	if err := os.WriteFile(hpath, []byte(newHdr), 0644); err != nil {
		t.Fatal(err)
	}
	wantBody := streamBodyOf(t, streamVCF)

	// Build a BGZF-compressed .vcf.gz input.
	gzPath := filepath.Join(dir, "in.vcf.gz")
	var gzBuf bytes.Buffer
	bw := bgzip.NewWriter(&gzBuf)
	if _, err := bw.Write([]byte(streamVCF)); err != nil {
		t.Fatal(err)
	}
	if err := bw.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(gzPath, gzBuf.Bytes(), 0644); err != nil {
		t.Fatal(err)
	}

	// Reheader the .vcf.gz -> VCF.gz output, then decompress and check the body.
	var out bytes.Buffer
	if _, err := ReheaderFile(gzPath, &out, ReheaderOptions{HeaderFile: hpath}); err != nil {
		t.Fatal(err)
	}
	if len(out.Bytes()) < 2 || out.Bytes()[0] != 0x1f || out.Bytes()[1] != 0x8b {
		t.Fatalf("VCF.gz input did not produce compressed output")
	}
	dec, err := iohelper.OpenReader(gzPath) // sanity: reader works on our .gz
	if err != nil {
		t.Fatal(err)
	}
	dec.Close()

	br, err := bgzip.NewReader(bytes.NewReader(out.Bytes()))
	if err != nil {
		t.Fatal(err)
	}
	plain, err := io.ReadAll(br)
	if err != nil {
		t.Fatal(err)
	}
	_, gotBody := splitVCFHeader(plain)
	if !bytes.Equal(gotBody, wantBody) {
		t.Errorf("VCF.gz body not byte-identical:\n got %q\nwant %q", gotBody, wantBody)
	}
	if !bytes.Contains(plain, []byte("##source=StreamTestGz")) {
		t.Errorf("VCF.gz header not replaced:\n%s", plain)
	}
}

// makeMultiBlockVCF returns a text VCF whose body exceeds one 64 KiB BGZF block,
// forcing the raw passthrough path to raw-copy at least one whole compressed
// block after the header/straddle block. It returns the full VCF text.
func makeMultiBlockVCF(nRecords int) string {
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	b.WriteString("##contig=<ID=chr1,length=100000000>\n")
	b.WriteString("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n")
	b.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
	for i := 0; i < nRecords; i++ {
		fmt.Fprintf(&b, "chr1\t%d\trs%d\tA\tT\t30\tPASS\tDP=%d\tGT\t0/1\t0/0\n", i+1, i+1, i%50)
	}
	return b.String()
}

// TestReheaderRawHeaderOnlyGz confirms the raw passthrough path handles a
// header-only VCF.gz (no body records): the header is replaced and the
// (empty) body stays empty.
func TestReheaderRawHeaderOnlyGz(t *testing.T) {
	dir := t.TempDir()
	const hdrOnly = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`
	hpath := filepath.Join(dir, "newhdr.vcf")
	const newHdr = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##source=HeaderOnlyGz
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`
	if err := os.WriteFile(hpath, []byte(newHdr), 0644); err != nil {
		t.Fatal(err)
	}

	gzPath := filepath.Join(dir, "hdronly.vcf.gz")
	writeBGZF(t, gzPath, hdrOnly)

	var out bytes.Buffer
	if _, err := ReheaderFile(gzPath, &out, ReheaderOptions{HeaderFile: hpath}); err != nil {
		t.Fatal(err)
	}
	plain := decodeBGZF(t, out.Bytes())
	_, gotBody := splitVCFHeader(plain)
	if len(gotBody) != 0 {
		t.Errorf("header-only input produced a non-empty body: %q", gotBody)
	}
	if !bytes.Contains(plain, []byte("##source=HeaderOnlyGz")) {
		t.Errorf("header not replaced:\n%s", plain)
	}
}

// TestReheaderRawHeaderOnlyPlain confirms the plain-VCF path handles a
// header-only VCF (no records) too.
func TestReheaderRawHeaderOnlyPlain(t *testing.T) {
	dir := t.TempDir()
	const hdrOnly = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`
	hpath := filepath.Join(dir, "newhdr.vcf")
	if err := os.WriteFile(hpath, []byte("##fileformat=VCFv4.2\n##source=HeaderOnlyPlain\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"), 0644); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := Reheader(strings.NewReader(hdrOnly), &out, ReheaderOptions{HeaderFile: hpath}); err != nil {
		t.Fatal(err)
	}
	_, gotBody := splitVCFHeader(out.Bytes())
	if len(gotBody) != 0 {
		t.Errorf("header-only input produced a non-empty body: %q", gotBody)
	}
	if !bytes.Contains(out.Bytes(), []byte("##source=HeaderOnlyPlain")) {
		t.Errorf("header not replaced:\n%s", out.Bytes())
	}
}

// TestReheaderRawMultiBlockGz exercises the raw compressed-block passthrough on
// a VCF.gz whose body spans several 64 KiB BGZF blocks. It asserts the body is
// byte-identical after reheader and the header was replaced — i.e. the raw
// block copy and the header/straddle-block re-compression compose correctly.
func TestReheaderRawMultiBlockGz(t *testing.T) {
	dir := t.TempDir()
	vcfText := makeMultiBlockVCF(4000) // ~250 KiB body -> multiple BGZF blocks
	wantBody := streamBodyOf(t, vcfText)
	if len(wantBody) < 64<<10 {
		t.Fatalf("test body too small (%d bytes) to span multiple BGZF blocks", len(wantBody))
	}

	gzPath := filepath.Join(dir, "multi.vcf.gz")
	writeBGZF(t, gzPath, vcfText)

	hpath := filepath.Join(dir, "newhdr.vcf")
	const newHdr = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000000>
##source=MultiBlockGz
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`
	if err := os.WriteFile(hpath, []byte(newHdr), 0644); err != nil {
		t.Fatal(err)
	}

	var out bytes.Buffer
	if _, err := ReheaderFile(gzPath, &out, ReheaderOptions{HeaderFile: hpath}); err != nil {
		t.Fatal(err)
	}
	if out.Len() < 2 || out.Bytes()[0] != 0x1f || out.Bytes()[1] != 0x8b {
		t.Fatalf("VCF.gz input did not produce compressed output")
	}
	plain := decodeBGZF(t, out.Bytes())
	gotHeader, gotBody := splitVCFHeader(plain)
	if !bytes.Equal(gotBody, wantBody) {
		t.Errorf("multi-block body not byte-identical (got %d bytes, want %d)", len(gotBody), len(wantBody))
	}
	if !bytes.Contains(gotHeader, []byte("##source=MultiBlockGz")) {
		t.Errorf("header not replaced:\n%s", gotHeader)
	}
	if bytes.Contains(gotHeader, []byte("chr1\t1\trs1")) {
		t.Errorf("a body record leaked into the header")
	}
}

// writeBGZF BGZF-compresses text into a file at path.
func writeBGZF(t *testing.T, path, text string) {
	t.Helper()
	var buf bytes.Buffer
	bw := bgzip.NewWriter(&buf)
	if _, err := bw.Write([]byte(text)); err != nil {
		t.Fatal(err)
	}
	if err := bw.Close(); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(path, buf.Bytes(), 0644); err != nil {
		t.Fatal(err)
	}
}

// decodeBGZF fully inflates a BGZF byte stream.
func decodeBGZF(t *testing.T, gz []byte) []byte {
	t.Helper()
	br, err := bgzip.NewReader(bytes.NewReader(gz))
	if err != nil {
		t.Fatal(err)
	}
	plain, err := io.ReadAll(br)
	if err != nil {
		t.Fatal(err)
	}
	return plain
}

// TestReheaderStreamSampleRenameChromOnly asserts that a -s sample rename
// rewrites only the #CHROM line and leaves every body record byte-identical.
func TestReheaderStreamSampleRenameChromOnly(t *testing.T) {
	dir := t.TempDir()
	srename := filepath.Join(dir, "names.txt")
	if err := os.WriteFile(srename, []byte("NEW1\nNEW2\n"), 0644); err != nil {
		t.Fatal(err)
	}
	wantBody := streamBodyOf(t, streamVCF)

	var out bytes.Buffer
	if _, err := Reheader(strings.NewReader(streamVCF), &out, ReheaderOptions{SamplesFile: srename}); err != nil {
		t.Fatal(err)
	}
	gotHeader, gotBody := splitVCFHeader(out.Bytes())

	// Body untouched.
	if !bytes.Equal(gotBody, wantBody) {
		t.Errorf("sample rename changed the body:\n got %q\nwant %q", gotBody, wantBody)
	}
	// Only the #CHROM line changed: every ## meta line is preserved verbatim.
	inHdr, _ := splitVCFHeader([]byte(streamVCF))
	for _, line := range bytes.Split(inHdr, []byte("\n")) {
		if bytes.HasPrefix(line, []byte("##")) && len(line) > 0 {
			if !bytes.Contains(gotHeader, line) {
				t.Errorf("meta line dropped/changed: %q", line)
			}
		}
	}
	if !bytes.Contains(gotHeader, []byte("\tNEW1\tNEW2\n")) {
		t.Errorf("#CHROM line not renamed:\n%s", gotHeader)
	}
	if bytes.Contains(gotHeader, []byte("\tS1\tS2\n")) {
		t.Errorf("old sample names still present in #CHROM line:\n%s", gotHeader)
	}
}
