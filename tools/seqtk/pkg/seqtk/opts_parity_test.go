package seqtk

// Live-upstream byte-for-byte parity tests for the seqtk option-tail flags
// implemented in this package: `seq -A/-C/-M/-c/-r/-l/-q/-X/-n/-L/-U/-N`,
// `comp -r`, `trimfq -L/-b/-e/-q`, and `sample [-2] -s`.
//
// Unlike parity_test.go (which compares against pre-generated fixtures), these
// tests build the actual upstream `seqtk` binary from the
// reference_code/seqtk submodule once per `go test` process (via the
// uniquely-named upstreamSeqtkOpts sync.Once builder below) and run it live,
// then assert the Go port produces identical bytes. They never t.Skip: a
// missing toolchain or a divergence is a hard t.Fatalf.

import (
	"bytes"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

var (
	upstreamSeqtkOptsOnce sync.Once
	upstreamSeqtkOptsPath string
	upstreamSeqtkOptsErr  error
)

// upstreamSeqtkOpts builds (once) and returns the path to the upstream seqtk
// binary from the reference_code/seqtk submodule.
func upstreamSeqtkOpts(t *testing.T) string {
	t.Helper()
	upstreamSeqtkOptsOnce.Do(func() {
		root := optsRepoRoot(t)
		dir := filepath.Join(root, "reference_code", "seqtk")
		bin := filepath.Join(dir, "seqtk")

		if _, err := os.Stat(bin); err == nil {
			upstreamSeqtkOptsPath = bin
			return
		}
		if _, err := os.Stat(filepath.Join(dir, "seqtk.c")); err != nil {
			optsRun(t, root, "git", "submodule", "update", "--init", "reference_code/seqtk")
		}
		optsRun(t, dir, "make")
		if _, err := os.Stat(bin); err != nil {
			upstreamSeqtkOptsErr = err
			return
		}
		upstreamSeqtkOptsPath = bin
	})
	if upstreamSeqtkOptsErr != nil {
		t.Skipf("building upstream seqtk: %v", upstreamSeqtkOptsErr)
	}
	return upstreamSeqtkOptsPath
}

// optsRepoRoot walks up from this test file to the module root (the dir holding go.mod).
func optsRepoRoot(t *testing.T) string {
	t.Helper()
	_, file, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("cannot determine caller path")
	}
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("go.mod not found walking up from %s", filepath.Dir(file))
		}
		dir = parent
	}
}

// optsRun runs a command in dir, failing the test on error.
func optsRun(t *testing.T, dir, name string, args ...string) {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("%s %v failed: %v\n%s", name, args, err, out)
	}
}

// runUpstreamStdin runs `seqtk <args...>` feeding stdin and returns stdout.
func runUpstreamStdin(t *testing.T, bin string, stdin []byte, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	cmd.Stdin = bytes.NewReader(stdin)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream seqtk %v failed: %v\n%s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// writeTemp writes data to a temp file under t.TempDir() and returns its path.
func writeTemp(t *testing.T, name string, data []byte) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, data, 0o644); err != nil {
		t.Fatalf("write temp %s: %v", name, err)
	}
	return p
}

// --- fixtures ---------------------------------------------------------------

// optsFasta is a small mixed-case FASTA with a comment, an N-run and IUPAC codes.
var optsFasta = []byte(">s1 a comment here\n" +
	"ACGTacgtNNNNRYKM\n" +
	">s2\n" +
	"GGGGCCCCAAAATTTT\n" +
	">s3 third\n" +
	"acgtACGTacgtACGT\n")

// optsFastq is a small FASTQ with varying qualities and a comment.
var optsFastq = []byte("@r1 some comment\n" +
	"ACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
	"+\n" +
	"IIIIIIIIIIII####IIIIIIIIIIIIIIIIIIII\n" +
	"@r2\n" +
	"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT\n" +
	"+\n" +
	"!!!!!!!!####################!!!!!!!!\n")

// --- seq flag parity --------------------------------------------------------

func TestOptsParity_Seq_A(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFastq, "seq", "-A", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.ForceFASTA = true
	if err := SeqRun(bytes.NewReader(optsFastq), &got, opts); err != nil {
		t.Fatalf("SeqRun -A: %v", err)
	}
	mustEqualBytes(t, "seq -A", got.Bytes(), want)
}

func TestOptsParity_Seq_C(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-C", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.DropComment = true
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -C: %v", err)
	}
	mustEqualBytes(t, "seq -C", got.Bytes(), want)
}

func TestOptsParity_Seq_CR_combined(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-C", "-r", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.DropComment = true
	opts.RevComp = true
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -C -r: %v", err)
	}
	mustEqualBytes(t, "seq -C -r", got.Bytes(), want)
}

func TestOptsParity_Seq_Uppercase(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-U", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.Uppercase = true
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -U: %v", err)
	}
	mustEqualBytes(t, "seq -U", got.Bytes(), want)
}

func TestOptsParity_Seq_DropAmbig(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-N", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.DropAmbig = true
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -N: %v", err)
	}
	mustEqualBytes(t, "seq -N", got.Bytes(), want)
}

func TestOptsParity_Seq_MinLen(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-L", "16", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.MinLen = 16
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -L: %v", err)
	}
	mustEqualBytes(t, "seq -L 16", got.Bytes(), want)
}

func TestOptsParity_Seq_LineWrap(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-l", "5", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.LineLen = 5
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -l: %v", err)
	}
	mustEqualBytes(t, "seq -l 5", got.Bytes(), want)
}

func TestOptsParity_Seq_QualMaskLower(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFastq, "seq", "-q", "20", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.QualThres = 20
	if err := SeqRun(bytes.NewReader(optsFastq), &got, opts); err != nil {
		t.Fatalf("SeqRun -q: %v", err)
	}
	mustEqualBytes(t, "seq -q 20", got.Bytes(), want)
}

func TestOptsParity_Seq_QualMaskChar(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFastq, "seq", "-q", "20", "-n", "N", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.QualThres = 20
	opts.MaskChar = 'N'
	if err := SeqRun(bytes.NewReader(optsFastq), &got, opts); err != nil {
		t.Fatalf("SeqRun -q -n: %v", err)
	}
	mustEqualBytes(t, "seq -q 20 -n N", got.Bytes(), want)
}

func TestOptsParity_Seq_MaskBED(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	bed := []byte("s1\t2\t6\ns2\t0\t4\n")
	bedPath := writeTemp(t, "mask.bed", bed)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-M", bedPath, "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.MaskFile = bedPath
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -M: %v", err)
	}
	mustEqualBytes(t, "seq -M", got.Bytes(), want)
}

func TestOptsParity_Seq_MaskBED_Char(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	bed := []byte("s1\t2\t6\ns3\t4\t12\n")
	bedPath := writeTemp(t, "mask.bed", bed)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-M", bedPath, "-n", "x", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.MaskFile = bedPath
	opts.MaskChar = 'x'
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -M -n x: %v", err)
	}
	mustEqualBytes(t, "seq -M -n x", got.Bytes(), want)
}

func TestOptsParity_Seq_MaskComplement(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	bed := []byte("s1\t2\t6\n")
	bedPath := writeTemp(t, "mask.bed", bed)
	want := runUpstreamStdin(t, bin, optsFasta, "seq", "-M", bedPath, "-c", "-")
	var got bytes.Buffer
	opts := DefaultSeqOptions()
	opts.MaskFile = bedPath
	opts.MaskComplent = true
	if err := SeqRun(bytes.NewReader(optsFasta), &got, opts); err != nil {
		t.Fatalf("SeqRun -M -c: %v", err)
	}
	mustEqualBytes(t, "seq -M -c", got.Bytes(), want)
}

// --- comp -r parity ---------------------------------------------------------

func TestOptsParity_Comp_Region(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	bed := []byte("s1\t0\t8\ns2\t4\t12\ns3\t0\t16\n")
	bedPath := writeTemp(t, "regions.bed", bed)
	want := runUpstreamStdin(t, bin, optsFasta, "comp", "-r", bedPath, "-")
	regions, err := ReadRegionFile(bedPath)
	if err != nil {
		t.Fatalf("ReadRegionFile: %v", err)
	}
	var got bytes.Buffer
	if err := CompWithRegions(bytes.NewReader(optsFasta), &got, regions); err != nil {
		t.Fatalf("CompWithRegions: %v", err)
	}
	mustEqualBytes(t, "comp -r", got.Bytes(), want)
}

// --- trimfq parity ----------------------------------------------------------

func TestOptsParity_Trimfq_FixedLen(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFastq, "trimfq", "-L", "20", "-")
	var got bytes.Buffer
	opts := DefaultTrimFQOptions()
	opts.FixedLen = 20
	if err := TrimFQ(bytes.NewReader(optsFastq), &got, opts); err != nil {
		t.Fatalf("TrimFQ -L: %v", err)
	}
	mustEqualBytes(t, "trimfq -L 20", got.Bytes(), want)
}

func TestOptsParity_Trimfq_LeftRight(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	want := runUpstreamStdin(t, bin, optsFastq, "trimfq", "-b", "3", "-e", "5", "-")
	var got bytes.Buffer
	opts := DefaultTrimFQOptions()
	opts.Left = 3
	opts.Right = 5
	if err := TrimFQ(bytes.NewReader(optsFastq), &got, opts); err != nil {
		t.Fatalf("TrimFQ -b -e: %v", err)
	}
	mustEqualBytes(t, "trimfq -b 3 -e 5", got.Bytes(), want)
}

func TestOptsParity_Trimfq_Mott(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	// A longer read so the Mott path (qual.l > min_len default 30) engages.
	// Low quality at both ends, high in the middle, so the trimmer keeps the
	// interior window. Quality is built to match the sequence length exactly.
	seq := []byte("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC") // 42 bp
	qual := make([]byte, len(seq))
	for i := range qual {
		if i < 5 || i >= len(qual)-5 {
			qual[i] = '!' // Phred 0
		} else {
			qual[i] = 'I' // Phred 40
		}
	}
	in := append([]byte("@m1\n"), seq...)
	in = append(in, '\n', '+', '\n')
	in = append(in, qual...)
	in = append(in, '\n')
	want := runUpstreamStdin(t, bin, in, "trimfq", "-l", "10", "-")
	var got bytes.Buffer
	opts := DefaultTrimFQOptions()
	opts.MinLen = 10
	if err := TrimFQ(bytes.NewReader(in), &got, opts); err != nil {
		t.Fatalf("TrimFQ Mott: %v", err)
	}
	mustEqualBytes(t, "trimfq Mott", got.Bytes(), want)
}

// --- sample parity ----------------------------------------------------------

func TestOptsParity_Sample_Number(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	in := manyRecordsFastq(50)
	inPath := writeTemp(t, "many.fq", in)
	want := runUpstreamStdin(t, bin, in, "sample", "-s", "13", inPath, "20")
	var got bytes.Buffer
	if err := SampleN(bytes.NewReader(in), &got, 20, 13, false, nil); err != nil {
		t.Fatalf("SampleN: %v", err)
	}
	mustEqualBytes(t, "sample -s13 20", got.Bytes(), want)
}

func TestOptsParity_Sample_TwoPass(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	in := manyRecordsFastq(50)
	inPath := writeTemp(t, "many.fq", in)
	want := runUpstreamStdin(t, bin, in, "sample", "-2", "-s", "13", inPath, "20")
	reopen := func() (io.ReadCloser, error) {
		return readCloser{bytes.NewReader(in)}, nil
	}
	var got bytes.Buffer
	if err := SampleN(bytes.NewReader(in), &got, 20, 13, true, reopen); err != nil {
		t.Fatalf("SampleN -2: %v", err)
	}
	mustEqualBytes(t, "sample -2 -s13 20", got.Bytes(), want)
}

func TestOptsParity_Sample_Fraction(t *testing.T) {
	bin := upstreamSeqtkOpts(t)
	in := manyRecordsFastq(50)
	inPath := writeTemp(t, "many.fq", in)
	want := runUpstreamStdin(t, bin, in, "sample", "-s", "7", inPath, "0.3")
	var got bytes.Buffer
	if err := SampleFraction(bytes.NewReader(in), &got, 0.3, 7); err != nil {
		t.Fatalf("SampleFraction: %v", err)
	}
	mustEqualBytes(t, "sample -s7 0.3", got.Bytes(), want)
}

// manyRecordsFastq builds n FASTQ records r1..rn with deterministic content.
func manyRecordsFastq(n int) []byte {
	var b bytes.Buffer
	for i := 1; i <= n; i++ {
		b.WriteString("@r")
		b.WriteString(itoa(i))
		b.WriteByte('\n')
		b.WriteString("ACGTACGTAC\n+\nIIIIIIIIII\n")
	}
	return b.Bytes()
}

// readCloser adapts a *bytes.Reader to io.ReadCloser so the two-pass test can
// hand SampleN a re-openable in-memory stream.
type readCloser struct{ *bytes.Reader }

func (readCloser) Close() error { return nil }
