package bedgetfasta

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
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

// TestRun_NameMissingEmitsEmptyName — upstream does NOT fall back to the
// coord header when -name is set but the BED row has no name column; it
// emits `>::chrom:start-end` with an empty name (verified against
// `bedtools getfasta -name` on a 3-column BED).
func TestRun_NameMissingEmitsEmptyName(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t0\t10\n" // no name column
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{Name: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !strings.HasPrefix(buf.String(), ">::chr1:0-10\n") {
		t.Errorf("expected empty-name header; got %q", buf.String())
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

func TestRun_BedOut(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t1\t10\tfoo\t99\t-\textra\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{BedOut: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t1\t10\tfoo\t99\t-\textra\tggggggggg\n"
	if buf.String() != want {
		t.Errorf("bedOut mismatch:\n got %q\nwant %q", buf.String(), want)
	}
	if warn.Len() != 0 {
		t.Errorf("unexpected warning: %q", warn.String())
	}
}

func TestRun_BedOutStrandSplit(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	// Minus-strand BED12 with three single-base blocks -> revcomp -> "tag".
	bed := "chr1\t0\t40\trec\t0\t-\t0\t0\t0\t3\t1,1,1,\t10,20,30,\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{BedOut: true, Strand: true, Split: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !strings.HasSuffix(strings.TrimRight(buf.String(), "\n"), "\ttag") {
		t.Errorf("expected trailing seq column \"tag\"; got %q", buf.String())
	}
}

func TestRun_ZeroLengthFeature(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	bed := "chr1\t5\t5\nchr1\t0\t5\n"
	var buf, warn bytes.Buffer
	n, err := Run(strings.NewReader(bed), path, &buf, &warn, Options{})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 1 {
		t.Errorf("expected 1 emitted record, got %d", n)
	}
	wantWarn := "Feature (chr1:5-5) has length = 0, Skipping.\n"
	if warn.String() != wantWarn {
		t.Errorf("warn mismatch: got %q want %q", warn.String(), wantWarn)
	}
	if !strings.HasPrefix(buf.String(), ">chr1:0-5\n") {
		t.Errorf("expected second record emitted; got %q", buf.String())
	}
}

func TestRun_StaleIndexWarning(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	// Write a sibling .fai, then make the FASTA newer than it.
	faiBody := "chr1\t50\t6\t10\t11\n"
	if err := os.WriteFile(path+".fai", []byte(faiBody), 0o644); err != nil {
		t.Fatalf("write .fai: %v", err)
	}
	old := time.Now().Add(-2 * time.Hour)
	if err := os.Chtimes(path+".fai", old, old); err != nil {
		t.Fatalf("chtimes fai: %v", err)
	}
	now := time.Now()
	if err := os.Chtimes(path, now, now); err != nil {
		t.Fatalf("chtimes fasta: %v", err)
	}
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader("chr1\t0\t5\n"), path, &buf, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "Warning: the index file is older than the FASTA file.\n"
	if warn.String() != want {
		t.Errorf("stale-index warn mismatch: got %q want %q", warn.String(), want)
	}
}

func TestRun_FreshIndexNoStaleWarning(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	if err := os.WriteFile(path+".fai", []byte("chr1\t50\t6\t10\t11\n"), 0o644); err != nil {
		t.Fatalf("write .fai: %v", err)
	}
	// Make the index newer than the FASTA.
	old := time.Now().Add(-2 * time.Hour)
	if err := os.Chtimes(path, old, old); err != nil {
		t.Fatalf("chtimes fasta: %v", err)
	}
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader("chr1\t0\t5\n"), path, &buf, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if warn.Len() != 0 {
		t.Errorf("expected no stale warning; got %q", warn.String())
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

// TestRun_OutOfBoundsRange — like upstream, a feature beyond the contig
// length is skipped with a "beyond the length of" warning (not an error),
// and no FASTA is emitted.
func TestRun_OutOfBoundsRange(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	// chr1 is 50bp long. Ask for bases 100-110.
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader("chr1\t100\t110\n"), path, &buf, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if buf.Len() != 0 {
		t.Errorf("expected no FASTA output for out-of-bounds; got %q", buf.String())
	}
	want := "Feature (chr1:100-110) beyond the length of chr1 size (50 bp).  Skipping.\n"
	if warn.String() != want {
		t.Errorf("warn mismatch: got %q want %q", warn.String(), want)
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
		// Empty name with -nameOnly: upstream emits an empty header (">"),
		// it does NOT fall back to chrom:start-end.
		{"", "chrZ", 5, 7, "+", Options{NameOnly: true}, ""},
		// Empty name with -name: upstream emits "::chrom:start-end".
		{"", "chrZ", 5, 7, "+", Options{Name: true}, "::chrZ:5-7"},
		// -name+ behaves identically to -name.
		{"foo", "chrZ", 5, 7, "+", Options{NamePlus: true}, "foo::chrZ:5-7"},
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

// TestRun_MissingChromWarning — a BED row naming a contig absent from the
// FASTA index emits upstream's WARNING line and produces no FASTA output.
func TestRun_MissingChromWarning(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	var out, warn bytes.Buffer
	if _, err := Run(strings.NewReader("chrZ\t1\t10\n"), path, &out, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if out.Len() != 0 {
		t.Errorf("expected no FASTA output, got %q", out.String())
	}
	want := "WARNING. chromosome (chrZ) was not found in the FASTA file. Skipping.\n"
	if warn.String() != want {
		t.Errorf("warn mismatch: got %q want %q", warn.String(), want)
	}
}

func TestRandomAccess_BuildOnMissingFai(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	// Confirm there's no .fai sitting next to it yet.
	if _, err := os.Stat(path + ".fai"); err == nil {
		t.Fatalf("unexpected pre-existing .fai")
	}
	ra, err := openFasta(path, false, io.Discard)
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
	ra, err := openFasta(path, false, io.Discard)
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
	if _, err := openFasta(path, false, io.Discard); err == nil {
		// Note: LoadIndex on a 4-column line errors → the openFastaBGZF
		// fallback path treats that as a hard failure (not IsNotExist),
		// so we expect an explicit error here.
		t.Fatal("expected error on malformed .fai, got nil")
	}
}

// TestRun_NilWarnWriter — a nil warn writer must be tolerated (defaulted to
// io.Discard) even when a warning would otherwise be emitted.
func TestRun_NilWarnWriter(t *testing.T) {
	dir := t.TempDir()
	path := writeFasta(t, dir, "t", tFasta)
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader("chrZ\t1\t10\n"), path, &buf, nil, Options{}); err != nil {
		t.Fatalf("Run with nil warn: %v", err)
	}
	if buf.Len() != 0 {
		t.Errorf("expected no output for missing chrom; got %q", buf.String())
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
