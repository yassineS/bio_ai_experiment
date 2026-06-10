package bcftools

// Live upstream-binary parity for the bcftools `-@/--threads` threading tail:
// the subcommands wired to route their bgzipped output through the shared
// bgzf.MultiWriter chokepoint when Threads > 1 (isec, convert, reheader,
// mendelian, mendelian2, csq, filter/vcffilter and the plugin output path).
//
// Two invariants are asserted per subcommand:
//
//  1. Thread-independence: -@ N output decodes byte-identically to -@ 1
//     output. Because every BGZF block is an independent gzip member, the
//     framed bytes may legitimately differ at block boundaries between thread
//     counts — so the comparison is always on DECODED plaintext.
//
//  2. Upstream parity (where an upstream surface exists): the decoded records
//     match the live upstream C binary's `-O z` output on the same fixture.
//
// A separate test pins finding #1: `view -Ob` / `-Oz` now flush the header
// into its OWN BGZF block (matching htslib's bgzf_flush after vcf_hdr_write /
// bcf_hdr_write), keeping tabix/.csi offsets clean.
//
// Per project rule the tests never t.Skip: when the upstream binary cannot be
// built they fail loudly via t.Fatalf through the shared builder.

import (
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

// upstreamBcftoolsThreadTail is the uniquely-named, sync.Once-guarded entry
// point the threading-tail tests use to obtain the freshly built upstream
// bcftools binary. It delegates to the shared builder (itself memoised) and
// calls t.Fatalf — never t.Skip — when the binary cannot be produced.
var threadTailBinOnce sync.Once
var threadTailBinPath string

func upstreamBcftoolsThreadTail(t *testing.T) string {
	t.Helper()
	threadTailBinOnce.Do(func() {
		threadTailBinPath = upstreamBcftools(t)
	})
	if threadTailBinPath == "" {
		t.Fatalf("upstream bcftools binary unavailable")
	}
	return threadTailBinPath
}

// splitBGZFBlocks walks a BGZF stream and returns the decompressed payload of
// each non-empty block in order. Every BGZF block is a self-contained gzip
// member, so it can be inflated independently with compress/gzip. The trailing
// 28-byte EOF block (ISIZE 0) is skipped. It t.Fatalf's on any malformed frame.
func splitBGZFBlocks(t *testing.T, b []byte) [][]byte {
	t.Helper()
	var blocks [][]byte
	for off := 0; off < len(b); {
		if off+18 > len(b) {
			t.Fatalf("bgzf: truncated header at offset %d", off)
		}
		if b[off] != 0x1f || b[off+1] != 0x8b || b[off+2] != 0x08 || b[off+3] != 0x04 {
			t.Fatalf("bgzf: bad block magic at offset %d", off)
		}
		xlen := int(binary.LittleEndian.Uint16(b[off+10 : off+12]))
		// Locate the BC subfield (SI1='B', SI2='C', SLEN=2) to read BSIZE.
		bsize := -1
		extra := b[off+12 : off+12+xlen]
		for i := 0; i+4 <= len(extra); {
			si1, si2 := extra[i], extra[i+1]
			slen := int(binary.LittleEndian.Uint16(extra[i+2 : i+4]))
			if si1 == 'B' && si2 == 'C' && slen == 2 {
				bsize = int(binary.LittleEndian.Uint16(extra[i+4 : i+6]))
				break
			}
			i += 4 + slen
		}
		if bsize < 0 {
			t.Fatalf("bgzf: BC subfield not found at offset %d", off)
		}
		blockLen := bsize + 1
		if off+blockLen > len(b) {
			t.Fatalf("bgzf: block overruns stream at offset %d", off)
		}
		// ISIZE (uncompressed size) is the last 4 bytes of the block.
		isize := binary.LittleEndian.Uint32(b[off+blockLen-4 : off+blockLen])
		if isize != 0 {
			zr, err := gzip.NewReader(bytes.NewReader(b[off : off+blockLen]))
			if err != nil {
				t.Fatalf("bgzf: gzip.NewReader at offset %d: %v", off, err)
			}
			payload, err := io.ReadAll(zr)
			if err != nil {
				t.Fatalf("bgzf: inflate block at offset %d: %v", off, err)
			}
			_ = zr.Close()
			blocks = append(blocks, payload)
		}
		off += blockLen
	}
	return blocks
}

// TestViewHeaderOwnBGZFBlock_Oz pins finding #1 for the -O z (VCF.gz) path:
// after the fix, openOutput flushes the header into its own BGZF block, so the
// first decompressed block ends exactly at the header/record boundary (its last
// byte is the newline terminating the #CHROM line) and carries no data record.
func TestViewHeaderOwnBGZFBlock_Oz(t *testing.T) {
	in := []byte(makeLargeVCF(2000))
	out := runParityView(t, in, ViewOptions{OutputFormat: OutputVCFGz, Threads: 1})

	blocks := splitBGZFBlocks(t, out)
	if len(blocks) < 2 {
		t.Fatalf("expected the header to occupy its own block plus >=1 record block, got %d blocks", len(blocks))
	}
	first := string(blocks[0])
	// The header block must end with the #CHROM column line and contain no
	// data records (data lines do not start with '#').
	if len(first) == 0 || first[len(first)-1] != '\n' {
		t.Fatalf("first block does not end on a line boundary")
	}
	// Every line in the first block must be a header line.
	for _, line := range splitLines(first) {
		if line == "" {
			continue
		}
		if line[0] != '#' {
			t.Fatalf("first BGZF block contains a non-header line %q — header was not flushed into its own block", line)
		}
	}
	// And the last header line must be the #CHROM line (so the WHOLE header,
	// not a prefix of it, lives in the first block).
	lines := splitLines(first)
	last := ""
	for i := len(lines) - 1; i >= 0; i-- {
		if lines[i] != "" {
			last = lines[i]
			break
		}
	}
	if len(last) < 6 || last[:6] != "#CHROM" {
		t.Fatalf("header block does not end on the #CHROM line, got %q", last)
	}
}

// TestViewHeaderOwnBGZFBlock_Ob pins finding #1 for the -O b (BCF) path. BCF is
// binary, but the same structural invariant holds: the BCF header section
// (magic 'BCF\2\2' + l_text + header text) must be flushed into the first BGZF
// block, separate from the record blocks. We assert the first decompressed
// block starts with the BCF magic and that there is at least one further block.
func TestViewHeaderOwnBGZFBlock_Ob(t *testing.T) {
	in := []byte(makeLargeVCF(2000))
	out := runParityView(t, in, ViewOptions{OutputFormat: OutputBCF, Threads: 1})

	blocks := splitBGZFBlocks(t, out)
	if len(blocks) < 2 {
		t.Fatalf("expected the BCF header in its own block plus >=1 record block, got %d blocks", len(blocks))
	}
	first := blocks[0]
	if len(first) < 5 || string(first[:3]) != "BCF" {
		t.Fatalf("first BGZF block does not start with the BCF magic: %q", first[:minInt(5, len(first))])
	}
	// The header block must contain exactly the magic + l_text + header text:
	// l_text (uint32 LE at offset 5) plus 9 prefix bytes equals the block len.
	if len(first) < 9 {
		t.Fatalf("BCF header block too short: %d bytes", len(first))
	}
	lText := binary.LittleEndian.Uint32(first[5:9])
	if int(lText)+9 != len(first) {
		t.Fatalf("BCF header block length %d != 9 + l_text %d — header not flushed into its own block",
			len(first), lText)
	}
}

func splitLines(s string) []string {
	var out []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == '\n' {
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	if start < len(s) {
		out = append(out, s[start:])
	}
	return out
}

// --- threading-tail helpers -------------------------------------------------

// writeFixture writes b to a temp file and returns its path.
func writeFixture(t *testing.T, name string, b []byte) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatalf("write fixture %s: %v", name, err)
	}
	return p
}

// runUpstream runs the upstream bcftools binary with the given args and returns
// raw stdout (which may be BGZF/BCF binary). It t.Fatalf's on failure.
func runUpstreamRaw(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools %v failed: %v\n%s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// convertThreadFixture is a small, valid VCF used by the threading-tail tests
// that re-emit input through openOutput unchanged. It is intentionally large
// enough (via makeLargeVCF) to span multiple BGZF blocks at -O z.

// TestConvertThreadsByteIdentical asserts `convert -O z` (which here is an
// identity VCF->VCF.gz conversion) decodes byte-identically across thread
// counts, and matches the live upstream binary.
func TestConvertThreadsByteIdentical(t *testing.T) {
	bin := upstreamBcftoolsThreadTail(t)
	in := []byte(makeLargeVCF(4000))
	path := writeFixture(t, "in.vcf", in)

	base := convertGz(t, path, 1)
	wantPlain := gunzipBytes(t, base)
	if len(wantPlain) < 2*65280 {
		t.Fatalf("fixture too small to span multiple BGZF blocks: %d bytes", len(wantPlain))
	}
	for _, n := range []int{2, 4, 8} {
		got := convertGz(t, path, n)
		if !bytes.Equal(gunzipBytes(t, got), wantPlain) {
			t.Fatalf("convert -@ %d decoded plaintext differs from -@ 1", n)
		}
	}

	// Upstream parity on decoded records.
	upstream := runUpstreamRaw(t, bin, "convert", "--no-version", "-O", "z", "--threads", "4", path)
	upRecs := dataLines(string(gunzipBytes(t, upstream)))
	gotRecs := dataLines(string(wantPlain))
	if !equalStringSlices(gotRecs, upRecs) {
		t.Fatalf("convert: record content mismatch vs upstream (got %d, upstream %d)", len(gotRecs), len(upRecs))
	}
}

func convertGz(t *testing.T, path string, threads int) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := ConvertFile(path, &out, ConvertOptions{OutputFormat: OutputVCFGz, CompressLevel: -1, Threads: threads}); err != nil {
		t.Fatalf("ConvertFile(threads=%d): %v", threads, err)
	}
	return out.Bytes()
}

// TestReheaderThreadsByteIdentical asserts `reheader -O z` output is
// thread-count independent on decoded plaintext. Upstream reheader preserves
// the input container (it has no -O), only multithreads BCF output, and on a
// no-op edit of a plain uncompressed VCF declines to run — so a direct
// upstream `-O z` comparison is not meaningful here. Record-level parity of
// reheader's pass-through is already covered by reheader_test.go; this test
// pins the threading invariant that the -@ wiring must not perturb output.
func TestReheaderThreadsByteIdentical(t *testing.T) {
	in := []byte(makeLargeVCF(4000))
	path := writeFixture(t, "in.vcf", in)

	base := reheaderGz(t, path, 1)
	wantPlain := gunzipBytes(t, base)
	if len(wantPlain) < 2*65280 {
		t.Fatalf("fixture too small to span multiple BGZF blocks: %d bytes", len(wantPlain))
	}
	for _, n := range []int{2, 4, 8} {
		got := reheaderGz(t, path, n)
		if !bytes.Equal(gunzipBytes(t, got), wantPlain) {
			t.Fatalf("reheader -@ %d decoded plaintext differs from -@ 1", n)
		}
	}
}

func reheaderGz(t *testing.T, path string, threads int) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := ReheaderFile(path, &out, ReheaderOptions{OutputFormat: OutputVCFGz, CompressLevel: -1, Threads: threads}); err != nil {
		t.Fatalf("ReheaderFile(threads=%d): %v", threads, err)
	}
	return out.Bytes()
}

// TestVCFFilterThreadsByteIdentical asserts `filter -O z` output is
// thread-count independent on decoded plaintext. No soft-filter is set so the
// records pass through unchanged; upstream parity is checked on records.
func TestVCFFilterThreadsByteIdentical(t *testing.T) {
	bin := upstreamBcftoolsThreadTail(t)
	in := []byte(makeLargeVCF(4000))
	path := writeFixture(t, "in.vcf", in)

	base := filterGz(t, in, 1)
	wantPlain := gunzipBytes(t, base)
	if len(wantPlain) < 2*65280 {
		t.Fatalf("fixture too small to span multiple BGZF blocks: %d bytes", len(wantPlain))
	}
	for _, n := range []int{2, 4, 8} {
		got := filterGz(t, in, n)
		if !bytes.Equal(gunzipBytes(t, got), wantPlain) {
			t.Fatalf("filter -@ %d decoded plaintext differs from -@ 1", n)
		}
	}

	// Upstream `filter -i 'QUAL>=0'` (always true) emits every record; compare
	// records only.
	upstream := runUpstreamRaw(t, bin, "filter", "--no-version", "-i", "QUAL>=0", "-O", "z", "--threads", "4", path)
	upRecs := dataLines(string(gunzipBytes(t, upstream)))
	gotFull := filterGzInclude(t, in, 4)
	gotRecs := dataLines(string(gunzipBytes(t, gotFull)))
	if !equalStringSlices(gotRecs, upRecs) {
		t.Fatalf("filter: record content mismatch vs upstream (got %d, upstream %d)", len(gotRecs), len(upRecs))
	}
}

func filterGz(t *testing.T, in []byte, threads int) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := VCFFilter(bytes.NewReader(in), &out, VCFFilterOptions{OutputFormat: OutputVCFGz, CompressLevel: -1, Threads: threads, NoVersion: true}); err != nil {
		t.Fatalf("VCFFilter(threads=%d): %v", threads, err)
	}
	return out.Bytes()
}

func filterGzInclude(t *testing.T, in []byte, threads int) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := VCFFilter(bytes.NewReader(in), &out, VCFFilterOptions{OutputFormat: OutputVCFGz, CompressLevel: -1, Threads: threads, IncludeExpr: "QUAL>=0", NoVersion: true}); err != nil {
		t.Fatalf("VCFFilter include(threads=%d): %v", threads, err)
	}
	return out.Bytes()
}

// TestIsecThreadsByteIdentical asserts `isec -w1 -O z` stdout output is
// thread-count independent on decoded plaintext.
func TestIsecThreadsByteIdentical(t *testing.T) {
	a := []byte(makeLargeVCF(4000))
	// Build a second input that shares every site with the first so the
	// intersection (default -w1 dump) is the full record set.
	pathA := writeFixture(t, "a.vcf", a)
	pathB := writeFixture(t, "b.vcf", a)

	base := isecGz(t, []string{pathA, pathB}, 1)
	wantPlain := gunzipBytes(t, base)
	if len(wantPlain) < 2*65280 {
		t.Fatalf("fixture too small to span multiple BGZF blocks: %d bytes", len(wantPlain))
	}
	for _, n := range []int{2, 4, 8} {
		got := isecGz(t, []string{pathA, pathB}, n)
		if !bytes.Equal(gunzipBytes(t, got), wantPlain) {
			t.Fatalf("isec -@ %d decoded plaintext differs from -@ 1", n)
		}
	}
}

func isecGz(t *testing.T, paths []string, threads int) []byte {
	t.Helper()
	var out bytes.Buffer
	opts := IsecOptions{
		Nfiles:       NfilesSpec{},
		OutputFormat: OutputVCFGz,
		Threads:      threads,
		Write:        []int{1},
	}
	if _, err := IsecFiles(paths, &out, opts); err != nil {
		t.Fatalf("IsecFiles(threads=%d): %v", threads, err)
	}
	return out.Bytes()
}

// makeLargeTrioVCF builds a valid trio VCF (samples CHILD/FATHER/MOTHER) large
// enough to span multiple BGZF blocks, for the mendelian2 threading test.
func makeLargeTrioVCF(nRecords int) string {
	var b bytes.Buffer
	b.WriteString("##fileformat=VCFv4.2\n")
	b.WriteString(`##FILTER=<ID=PASS,Description="All filters passed">` + "\n")
	b.WriteString("##contig=<ID=chr1,length=300000000>\n")
	b.WriteString(`##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">` + "\n")
	b.WriteString(`##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">` + "\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHILD\tFATHER\tMOTHER\n")
	bases := []string{"A", "C", "G", "T"}
	gts := []string{"0/0", "0/1", "1/1"}
	for i := 0; i < nRecords; i++ {
		pos := (i + 1) * 7
		ref := bases[i%4]
		alt := bases[(i+1)%4]
		c := gts[i%3]
		f := gts[(i+1)%3]
		m := gts[(i+2)%3]
		b.WriteString("chr1\t")
		b.WriteString(itoa(pos))
		b.WriteString("\t.\t" + ref + "\t" + alt + "\t.\tPASS\tDP=" + itoa(10+(i%90)) + "\tGT\t" + c + "\t" + f + "\t" + m + "\n")
	}
	return b.String()
}

// TestMendelian2ThreadsByteIdentical asserts mendelian2's VCF.gz output (here
// in annotate mode, which streams every record back out through openOutput) is
// thread-count independent on decoded plaintext.
func TestMendelian2ThreadsByteIdentical(t *testing.T) {
	in := []byte(makeLargeTrioVCF(6000))

	base := mendelian2Gz(t, in, 1)
	wantPlain := gunzipBytes(t, base)
	if len(wantPlain) < 2*65280 {
		t.Fatalf("fixture too small to span multiple BGZF blocks: %d bytes", len(wantPlain))
	}
	for _, n := range []int{2, 4, 8} {
		got := mendelian2Gz(t, in, n)
		if !bytes.Equal(gunzipBytes(t, got), wantPlain) {
			t.Fatalf("mendelian2 -@ %d decoded plaintext differs from -@ 1", n)
		}
	}
}

func mendelian2Gz(t *testing.T, in []byte, threads int) []byte {
	t.Helper()
	var out bytes.Buffer
	_, err := Mendelian2(bytes.NewReader(in), &out, Mendelian2Options{
		Trios:         []Mendelian2Trio{{Child: "CHILD", Father: "FATHER", Mother: "MOTHER"}},
		Mode:          Mendelian2Annotate,
		OutputFormat:  OutputVCFGz,
		CompressLevel: -1,
		Threads:       threads,
	})
	if err != nil {
		t.Fatalf("Mendelian2(threads=%d): %v", threads, err)
	}
	return out.Bytes()
}
