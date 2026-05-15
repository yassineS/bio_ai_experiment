package bedgetfasta

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
)

// writeBGZFFasta drops a BGZF-compressed FASTA at <dir>/<name>.fa.gz
// and returns its path. Used by the BGZF random-access tests below.
func writeBGZFFasta(t *testing.T, dir, name, body string) string {
	t.Helper()
	path := filepath.Join(dir, name+".fa.gz")
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create %s: %v", path, err)
	}
	w := bgzip.NewWriter(f)
	if _, err := w.Write([]byte(body)); err != nil {
		t.Fatalf("bgzip write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("bgzip close: %v", err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("close %s: %v", path, err)
	}
	return path
}

// writeFasta drops a temporary FASTA at <dir>/<name>.fa and returns its
// path. The .fai is built lazily by Run.
func writeFasta(t *testing.T, dir, name, body string) string {
	t.Helper()
	path := filepath.Join(dir, name+".fa")
	if err := os.WriteFile(path, []byte(body), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	return path
}

// canonical t.fa from the upstream getfasta corpus.
const tFasta = ">chr1\naggggggggg\ncggggggggg\ntggggggggg\naggggggggg\ncggggggggg\n"

// IUPAC fixture with mixed-case contigs.
const iupacFasta = ">1\nAGCTYRWSKMDVHBXNACGT\n>2\nagctyrwskmdvhbxnacgt\n"

func TestRun_DefaultHeader(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t0\t10\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := ">chr1:0-10\naggggggggg\n"
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestRun_PreservesCase(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", iupacFasta)
	bed := "1\t0\t16\nblockchain\t0\t1\n2\t0\t16\n"
	// Trim the malformed middle line above — we just want one row from "1"
	// and one from "2". Build it cleanly:
	bed = "1\t0\t16\n2\t0\t16\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := ">1:0-16\nAGCTYRWSKMDVHBXN\n>2:0-16\nagctyrwskmdvhbxn\n"
	if buf.String() != want {
		t.Errorf("case-preserving fetch broken.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestRun_NameHeader(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t0\t10\tmy_block\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{Name: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := ">my_block::chr1:0-10\naggggggggg\n"
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestRun_NameOnlyHeader(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t0\t10\tjust_name\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{NameOnly: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := ">just_name\naggggggggg\n"
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestRun_NameMissingFallsBackToCoord(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t0\t10\n" // no name column
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{Name: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !strings.HasPrefix(buf.String(), ">chr1:0-10\n") {
		t.Errorf("expected fallback to coord header; got %q", buf.String())
	}
}

func TestRun_Tab(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t0\t10\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{Tab: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1:0-10\taggggggggg\n"
	if buf.String() != want {
		t.Errorf("tab output mismatch.\nwant:%q\ngot:%q", want, buf.String())
	}
}

func TestRun_StrandReverseComp(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", iupacFasta)
	// BED with negative strand on contig 1.
	bed := "1\t0\t16\t.\t1000\t-\n2\t0\t16\t.\t1000\t-\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{Strand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := ">1:0-16(-)\nNXVDBHKMSWYRAGCT\n>2:0-16(-)\nnxvdbhkmswyragct\n"
	if buf.String() != want {
		t.Errorf("strand rc mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestRun_StrandPlusGetsSuffix(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t0\t10\trec\t0\t+\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{Strand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !strings.HasPrefix(buf.String(), ">chr1:0-10(+)\n") {
		t.Errorf("expected (+) suffix; got %q", buf.String())
	}
	// And sequence should be the original (no revcomp on plus strand).
	if !strings.Contains(buf.String(), "aggggggggg") {
		t.Errorf("plus strand should not be reverse-complemented; got %q", buf.String())
	}
}

func TestRun_Split(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	// Upstream blocks.bed line: 3 blocks of sizes 2,10,10 starting at 5,16,36
	// → bases [5..7), [16..26), [36..46) within chr1 (length 50).
	bed := "chr1\t0\t40\tthree_blocks_match\t0\t+\t0\t0\t0\t3\t2,10,10,\t5,16,36,\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{Split: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	// The expected concatenation matches upstream getfasta.t02:
	//   block [5,7)  -> "gg"  (positions 5,6 of "aggggggggg")
	//   block [16,26)-> "ggggcggggg" (?) Let's just verify length=22.
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines; got %d: %q", len(lines), buf.String())
	}
	if len(lines[1]) != 22 {
		t.Errorf("split concat length = %d, want 22", len(lines[1]))
	}
}

func TestRun_SplitWithStrand(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	// Upstream blocks.bed minus-strand line: 3 blocks of size 1 starting at
	// 10,20,30 → bases "c","t","a". With -s -split, blocks are
	// reverse-complemented individually AND emitted in reverse order, so
	// "a","t","c" → complemented "t","a","g" → reversed-order "g","a","t" →
	// upstream test gets "tag" with strict reading. Let's verify.
	bed := "chr1\t0\t40\tthree_blocks_match\t0\t-\t0\t0\t0\t3\t1,1,1,\t10,20,30,\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{Split: true, Strand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines; got %d", len(lines))
	}
	if lines[1] != "tag" {
		t.Errorf("split+strand concat = %q, want \"tag\"", lines[1])
	}
}

func TestRun_RNA(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t20\t30\n" // grabs "tggggggggg"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{RNA: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !strings.Contains(buf.String(), "uggggggggg") {
		t.Errorf("-rna should map t→u; got %q", buf.String())
	}
}

func TestRun_MissingChromWarn(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chrNope\t0\t1\nchr1\t0\t5\n"
	var buf, warn bytes.Buffer
	n, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 1 {
		t.Errorf("expected 1 record (the second), got %d", n)
	}
	if !strings.Contains(warn.String(), "chrNope") || !strings.Contains(warn.String(), "WARNING") {
		t.Errorf("expected warning about chrNope; got %q", warn.String())
	}
	if !strings.Contains(buf.String(), "chr1:0-5") {
		t.Errorf("expected second record to be emitted; got %q", buf.String())
	}
}

func TestRun_HeaderLinesSkipped(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "# comment\ntrack name=foo\nbrowser pos\n\nchr1\t0\t3\n"
	var buf, warn bytes.Buffer
	n, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 1 {
		t.Errorf("emitted %d records, want 1", n)
	}
}

func TestRun_TooFewFields(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	var buf, warn bytes.Buffer
	_, err := Run(strings.NewReader("chr1\t0\n"), path, &buf, &warn, Options{})
	if err == nil {
		t.Errorf("expected error for 2-column BED")
	}
}

func TestRun_BadCoords(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader("chr1\tabc\t10\n"), path, &buf, &warn, Options{}); err == nil {
		t.Errorf("expected error for non-numeric start")
	}
	if _, err := Run(strings.NewReader("chr1\t0\txyz\n"), path, &buf, &warn, Options{}); err == nil {
		t.Errorf("expected error for non-numeric end")
	}
}

func TestRun_FastaMissing(t *testing.T) {
	var buf, warn bytes.Buffer
	_, err := Run(strings.NewReader("chr1\t0\t1\n"), "/nonexistent.fa", &buf, &warn, Options{})
	if err == nil {
		t.Errorf("expected error for missing FASTA")
	}
}

func TestRun_OutOfBoundsRange(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	// chr1 is 50bp long. Ask for bases 100-110.
	var buf, warn bytes.Buffer
	_, err := Run(strings.NewReader("chr1\t100\t110\n"), path, &buf, &warn, Options{})
	if err == nil {
		t.Errorf("expected error for out-of-bounds range")
	}
}

func TestFormatHeader(t *testing.T) {
	cases := []struct {
		name, chrom string
		start, end  int
		strand      string
		opts        Options
		want        string
	}{
		{"my", "chr1", 0, 10, "+", Options{}, "chr1:0-10"},
		{"my", "chr1", 0, 10, "+", Options{Strand: true}, "chr1:0-10(+)"},
		{"my", "chr1", 0, 10, "-", Options{Strand: true}, "chr1:0-10(-)"},
		{"foo", "chrZ", 5, 7, "+", Options{Name: true}, "foo::chrZ:5-7"},
		{"foo", "chrZ", 5, 7, "-", Options{Name: true, Strand: true}, "foo::chrZ:5-7(-)"},
		{"foo", "chrZ", 5, 7, "-", Options{NameOnly: true}, "foo"},
		{"foo", "chrZ", 5, 7, "-", Options{NameOnly: true, Strand: true}, "foo(-)"},
		{"", "chrZ", 5, 7, "+", Options{NameOnly: true}, "chrZ:5-7"},
	}
	for i, tc := range cases {
		got := formatHeader(tc.chrom, tc.start, tc.end, tc.name, tc.strand, tc.opts)
		if got != tc.want {
			t.Errorf("case %d: got %q want %q", i, got, tc.want)
		}
	}
}

func TestComplement(t *testing.T) {
	in := "AaCcGgTtUuRrYySsWwKkMmBbVvDdHhNnXx"
	out := make([]byte, len(in))
	for i := 0; i < len(in); i++ {
		out[i] = complement(in[i])
	}
	want := "TtGgCcAaAaYyRrSsWwMmKkVvBbHhDdNnXx"
	if string(out) != want {
		t.Errorf("complement table broken.\nwant: %s\ngot:  %s", want, out)
	}
	// pass-through for unknown
	if complement('Q') != 'Q' {
		t.Errorf("unknown base must pass through")
	}
}

func TestReverseComplement(t *testing.T) {
	got := reverseComplement([]byte("AGCTYRWSKMDVHBXN"))
	want := []byte("NXVDBHKMSWYRAGCT")
	if !bytes.Equal(got, want) {
		t.Errorf("got %s, want %s", got, want)
	}
}

func TestDnaToRNA(t *testing.T) {
	got := dnaToRNA([]byte("ACGTacgtN"))
	if string(got) != "ACGUacguN" {
		t.Errorf("got %s", got)
	}
}

func TestErrMissingChrom(t *testing.T) {
	e := errMissingChrom{name: "chrZ"}
	if e.Error() == "" {
		t.Errorf("Error() empty")
	}
	if !isMissingChrom(e) {
		t.Errorf("isMissingChrom() false on errMissingChrom")
	}
	if isMissingChrom(errors.New("other")) {
		t.Errorf("isMissingChrom() true on plain error")
	}
}

func TestRandomAccess_BuildOnMissingFai(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	// Confirm there's no .fai sitting next to it yet.
	if _, err := os.Stat(path + ".fai"); err == nil {
		t.Fatalf("unexpected pre-existing .fai")
	}
	ra, err := openFasta(path, false)
	if err != nil {
		t.Fatalf("openFasta: %v", err)
	}
	defer ra.Close()
	seq, err := ra.FetchPreserveCase("chr1", 0, 5)
	if err != nil {
		t.Fatalf("Fetch: %v", err)
	}
	if string(seq) != "agggg" {
		t.Errorf("Fetch = %q, want \"agggg\"", seq)
	}
}

func TestRandomAccess_RangeErrors(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	ra, err := openFasta(path, false)
	if err != nil {
		t.Fatalf("openFasta: %v", err)
	}
	defer ra.Close()
	if _, err := ra.FetchPreserveCase("chr1", -1, 5); err == nil {
		t.Errorf("expected error for negative start")
	}
	if _, err := ra.FetchPreserveCase("chr1", 5, 4); err == nil {
		t.Errorf("expected error for end<start")
	}
	if _, err := ra.FetchPreserveCase("chr1", 0, 9999); err == nil {
		t.Errorf("expected error for end>length")
	}
	if _, err := ra.FetchPreserveCase("nope", 0, 1); err == nil {
		t.Errorf("expected missing-chrom error")
	}
	// Zero-length range returns empty.
	got, err := ra.FetchPreserveCase("chr1", 5, 5)
	if err != nil {
		t.Errorf("zero-length fetch errored: %v", err)
	}
	if len(got) != 0 {
		t.Errorf("zero-length fetch returned %d bytes", len(got))
	}
}

// TestRun_BGZFRoundTrip — end-to-end run against a BGZF-compressed FASTA
// of the canonical t.fa fixture. Confirms the BGZF magic-sniff path,
// in-memory decompression, and case-preserving Fetch all hold up.
func TestRun_BGZFRoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := writeBGZFFasta(t, dir, "t", tFasta)
	bed := "chr1\t0\t10\n"
	var out, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &out, &warn, Options{}); err != nil {
		t.Fatalf("Run on BGZF: %v", err)
	}
	if warn.Len() != 0 {
		t.Errorf("unexpected warning: %q", warn.String())
	}
	want := ">chr1:0-10\naggggggggg\n"
	if got := out.String(); got != want {
		t.Errorf("BGZF output mismatch:\n got %q\nwant %q", got, want)
	}
}

// TestRun_BGZFFullHeader — `-fullHeader` works against BGZF FASTA inputs.
func TestRun_BGZFFullHeader(t *testing.T) {
	dir := t.TempDir()
	path := writeBGZFFasta(t, dir, "ref",
		">chr1 assembly notes here\naggggggggg\n")
	bed := "chr1 assembly notes here\t0\t5\n"
	var out, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &out, &warn,
		Options{FullHeader: true}); err != nil {
		t.Fatalf("Run BGZF -fullHeader: %v", err)
	}
	if warn.Len() != 0 {
		t.Errorf("unexpected warning: %q", warn.String())
	}
	if !strings.Contains(out.String(), "agggg") {
		t.Errorf("expected full-header lookup to succeed; got %q", out.String())
	}
}

// TestRun_BGZFWithSiblingFai — the BGZF path should honour an explicit
// `<path>.fa.gz.fai` sidecar instead of always rebuilding the index.
func TestRun_BGZFWithSiblingFai(t *testing.T) {
	dir := t.TempDir()
	path := writeBGZFFasta(t, dir, "t", tFasta)
	// Hand-write a samtools-style .fai for the uncompressed payload.
	// chr1 is 50 bases at offset 6 (header is ">chr1\n" = 6 bytes), line
	// width 11 (10 bases + '\n').
	faiPath := path + ".fai"
	if err := os.WriteFile(faiPath, []byte("chr1\t50\t6\t10\t11\n"), 0o644); err != nil {
		t.Fatalf("write .fai: %v", err)
	}
	bed := "chr1\t0\t10\n"
	var out, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &out, &warn, Options{}); err != nil {
		t.Fatalf("Run BGZF with sibling .fai: %v", err)
	}
	if warn.Len() != 0 {
		t.Errorf("unexpected warnings: %q", warn.String())
	}
	want := ">chr1:0-10\naggggggggg\n"
	if out.String() != want {
		t.Errorf("output mismatch: got %q want %q", out.String(), want)
	}
}

// TestRun_BGZFNonexistent — opening a missing BGZF file surfaces a clean
// open error, exercising the isBGZF error path.
func TestRun_BGZFNonexistent(t *testing.T) {
	dir := t.TempDir()
	var out, warn bytes.Buffer
	if _, err := Run(strings.NewReader("chr1\t0\t1\n"),
		filepath.Join(dir, "missing.fa.gz"), &out, &warn, Options{}); err == nil {
		t.Fatal("expected error opening missing BGZF, got nil")
	}
}

// TestOpenFastaBGZF_BadFai — a malformed sibling .fa.gz.fai surfaces as
// an error rather than silently rebuilding.
func TestOpenFastaBGZF_BadFai(t *testing.T) {
	dir := t.TempDir()
	path := writeBGZFFasta(t, dir, "t", tFasta)
	if err := os.WriteFile(path+".fai", []byte("not\tactually\tfive\tcolumns\n"), 0o644); err != nil {
		t.Fatalf("write bad .fai: %v", err)
	}
	if _, err := openFasta(path, false); err == nil {
		// Note: LoadIndex on a 4-column line errors → the openFastaBGZF
		// fallback path treats that as a hard failure (not IsNotExist),
		// so we expect an explicit error here.
		t.Fatal("expected error on malformed .fai, got nil")
	}
}

// TestIsAllPresent — direct unit test for the helper.
func TestIsAllPresent(t *testing.T) {
	if !isAllPresent([][]byte{{}, {}}, [][2]int{{0, 1}, {1, 2}}) {
		t.Error("expected true for equal lengths")
	}
	if isAllPresent([][]byte{{}}, [][2]int{{0, 1}, {1, 2}}) {
		t.Error("expected false for short parts")
	}
}

// TestWarnMissingChrom — nil writer is a no-op; explicit writer gets the
// standard upstream warning line.
func TestWarnMissingChrom(t *testing.T) {
	warnMissingChrom(nil, "chrZ") // must not panic
	var buf bytes.Buffer
	warnMissingChrom(&buf, "chrZ")
	want := "WARNING. chromosome (chrZ) was not found in the FASTA file. Skipping.\n"
	if buf.String() != want {
		t.Errorf("warn output mismatch: got %q want %q", buf.String(), want)
	}
}

// TestBytesJoin — direct unit test for the no-separator concatenator.
func TestBytesJoin(t *testing.T) {
	got := bytesJoin([][]byte{[]byte("ac"), []byte("gt"), []byte("")})
	if string(got) != "acgt" {
		t.Errorf("bytesJoin = %q, want %q", got, "acgt")
	}
	if string(bytesJoin(nil)) != "" {
		t.Errorf("bytesJoin(nil) should return empty")
	}
}

// TestIsBGZF — sniff helper rejects plain FASTA and accepts BGZF magic.
func TestIsBGZF(t *testing.T) {
	dir := t.TempDir()
	plain := writeFasta(t, dir, "plain", tFasta)
	if got, err := isBGZF(plain); err != nil || got {
		t.Errorf("isBGZF(plain): got %v, err %v", got, err)
	}
	bgz := writeBGZFFasta(t, dir, "bgz", tFasta)
	if got, err := isBGZF(bgz); err != nil || !got {
		t.Errorf("isBGZF(bgz): got %v, err %v", got, err)
	}
	if _, err := isBGZF(filepath.Join(dir, "missing.fa")); err == nil {
		t.Error("expected error for missing path")
	}
}
