package samtools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// smallFasta is a tiny multi-contig FASTA with mixed case and a 10-base wrap,
// used by the self-contained faidx tests. chr1 is 25 bases (uneven last line),
// chr2 is exactly 20 bases.
const smallFasta = ">chr1 desc one\n" +
	"ACGTacgtAC\n" +
	"GTACGTacgt\n" +
	"ACGTa\n" +
	">chr2\n" +
	"TTTTGGGGCC\n" +
	"CCAAAANNNN\n"

// writeTempFasta writes content to a temp FASTA and returns its path.
func writeTempFasta(t *testing.T, content string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	return p
}

// runFaidx is a small in-process driver: it builds DefaultFaidxOptions, applies
// the mutators, and captures stdout/stderr/exit.
func runFaidxInProc(t *testing.T, path string, regions []string, mutate func(*FaidxOptions)) (stdout, stderr string, exit int) {
	t.Helper()
	opts := DefaultFaidxOptions(FaidxFASTA)
	if mutate != nil {
		mutate(&opts)
	}
	var out, errBuf bytes.Buffer
	exit = Faidx(path, regions, opts, &out, &errBuf)
	return out.String(), errBuf.String(), exit
}

// TestFaidxBuildIndex verifies the .fai build matches the expected 5-column
// layout for a known fixture.
func TestFaidxBuildIndex(t *testing.T) {
	path := writeTempFasta(t, smallFasta)
	if err := FaidxBuild(path, DefaultFaidxOptions(FaidxFASTA)); err != nil {
		t.Fatalf("FaidxBuild: %v", err)
	}
	got, err := os.ReadFile(path + ".fai")
	if err != nil {
		t.Fatalf("read .fai: %v", err)
	}
	// chr1: 25 bases, offset 15 (">chr1 desc one\n"), 10 bases/line, 11 bytes.
	// chr2: 20 bases. chr1 occupies 15 + (11+11+6) = 15 + 28 = 43; chr2 header
	// ">chr2\n" is 6 bytes → chr2 seq offset 43 + 6 = 49.
	want := "chr1\t25\t15\t10\t11\n" +
		"chr2\t20\t49\t10\t11\n"
	if string(got) != want {
		t.Errorf("fai mismatch:\n got: %q\nwant: %q", got, want)
	}
}

// TestFaidxExtractTable drives a matrix of extraction cases through the
// in-process API, checking stdout, stderr and the exit code.
func TestFaidxExtractTable(t *testing.T) {
	path := writeTempFasta(t, smallFasta)
	if err := FaidxBuild(path, DefaultFaidxOptions(FaidxFASTA)); err != nil {
		t.Fatalf("FaidxBuild: %v", err)
	}

	cases := []struct {
		name    string
		regions []string
		mutate  func(*FaidxOptions)
		stdout  string
		stderr  string
		exit    int
	}{
		{
			name:    "single base preserves case",
			regions: []string{"chr1:4-4"},
			stdout:  ">chr1:4-4\nT\n",
		},
		{
			name:    "range wraps at contig width (10)",
			regions: []string{"chr1:1-25"},
			stdout:  ">chr1:1-25\nACGTacgtAC\nGTACGTacgt\nACGTa\n",
		},
		{
			name:    "whole contig",
			regions: []string{"chr2"},
			stdout:  ">chr2\nTTTTGGGGCC\nCCAAAANNNN\n",
		},
		{
			name:    "chr:0 is whole contig",
			regions: []string{"chr2:0"},
			stdout:  ">chr2:0\nTTTTGGGGCC\nCCAAAANNNN\n",
		},
		{
			name:    "open-ended start",
			regions: []string{"chr1:21"},
			stdout:  ">chr1:21\nACGTa\n",
		},
		{
			name:    "leading-dash end is 1..M",
			regions: []string{"chr1:-5"},
			stdout:  ">chr1:-5\nACGTa\n",
		},
		{
			name:    "end clamp emits truncated warning, exit 0",
			regions: []string{"chr1:20-40"},
			stdout:  ">chr1:20-40\ntACGTa\n",
			stderr:  "[faidx] Truncated sequence: chr1:20-40\n",
		},
		{
			name:    "single line with -n 0",
			regions: []string{"chr1:1-25"},
			mutate:  func(o *FaidxOptions) { o.LineLen = 0 },
			stdout:  ">chr1:1-25\nACGTacgtACGTACGTacgtACGTa\n",
		},
		{
			name:    "custom wrap -n 5",
			regions: []string{"chr2:1-20"},
			mutate:  func(o *FaidxOptions) { o.LineLen = 5 },
			stdout:  ">chr2:1-20\nTTTTG\nGGGCC\nCCAAA\nANNNN\n",
		},
		{
			name:    "reverse complement adds /rc, preserves case",
			regions: []string{"chr1:1-10"},
			mutate:  func(o *FaidxOptions) { o.ReverseComplement = true },
			// ACGTacgtAC reverse-complemented = GTacgtACGT
			stdout: ">chr1:1-10/rc\nGTacgtACGT\n",
		},
		{
			name:    "mark-strand sign on + strand",
			regions: []string{"chr1:1-5"},
			mutate:  func(o *FaidxOptions) { _ = ParseMarkStrand(o, "sign") },
			stdout:  ">chr1:1-5(+)\nACGTa\n",
		},
		{
			name:    "mark-strand sign on - strand",
			regions: []string{"chr1:1-5"},
			mutate: func(o *FaidxOptions) {
				_ = ParseMarkStrand(o, "sign")
				o.ReverseComplement = true
			},
			// ACGTa rc = tACGT
			stdout: ">chr1:1-5(-)\ntACGT\n",
		},
		{
			name:    "multi-region in one call",
			regions: []string{"chr1:1-5", "chr2:1-4"},
			stdout:  ">chr1:1-5\nACGTa\n>chr2:1-4\nTTTT\n",
		},
		{
			name:    "unknown contig aborts with exit 1",
			regions: []string{"chrZ:1-5"},
			stdout:  ">chrZ:1-5\n",
			stderr: "[W::fai_get_val] Reference chrZ:1-5 not found in FASTA file, returning empty sequence\n" +
				"[W::fai_get_val] Reference chrZ:1-5 not found in FASTA file, returning empty sequence\n" +
				"[faidx] Failed to fetch sequence in chrZ:1-5\n",
			exit: 1,
		},
		{
			name:    "unknown contig with -c continues",
			regions: []string{"chrZ:1-5", "chr1:1-3"},
			mutate:  func(o *FaidxOptions) { o.Continue = true },
			stdout:  ">chrZ:1-5\n>chr1:1-3\nACG\n",
			stderr: "[W::fai_get_val] Reference chrZ:1-5 not found in FASTA file, returning empty sequence\n" +
				"[W::fai_get_val] Reference chrZ:1-5 not found in FASTA file, returning empty sequence\n" +
				"[faidx] Failed to fetch sequence in chrZ:1-5\n",
			exit: 0,
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			stdout, stderr, exit := runFaidxInProc(t, path, c.regions, c.mutate)
			if stdout != c.stdout {
				t.Errorf("stdout:\n got: %q\nwant: %q", stdout, c.stdout)
			}
			if stderr != c.stderr {
				t.Errorf("stderr:\n got: %q\nwant: %q", stderr, c.stderr)
			}
			if exit != c.exit {
				t.Errorf("exit: got %d, want %d", exit, c.exit)
			}
		})
	}
}

// TestParseFaidxRegion checks the region parser against the index for the
// coordinate-defaulting branches that mirror htslib's hts_parse_region.
func TestParseFaidxRegion(t *testing.T) {
	idx, err := fasta.BuildIndexBytes([]byte(smallFasta))
	if err != nil {
		t.Fatalf("build index: %v", err)
	}
	cases := []struct {
		region   string
		wantOK   bool
		wantName string
		wantBeg  int64
		wantEnd  int64 // clamped end (0-based exclusive)
	}{
		{"chr1", true, "chr1", 0, 25},
		{"chr1:", true, "chr1", 0, 25},
		{"chr1:0", true, "chr1", 0, 25},
		{"chr1:5", true, "chr1", 4, 25},
		{"chr1:5-", true, "chr1", 4, 25},
		{"chr1:-5", true, "chr1", 0, 5},
		{"chr1:3-8", true, "chr1", 2, 8},
		{"chr1:20-40", true, "chr1", 19, 25}, // clamped to length
		{"chr1:5-4", false, "", 0, 0},        // end < start → parse fails
		{"chrZ", false, "", 0, 0},
		{"chrZ:1-5", false, "", 0, 0},
		{"chr1:1,0-1,5", true, "chr1", 9, 15}, // commas: 10-15
	}
	for _, c := range cases {
		t.Run(c.region, func(t *testing.T) {
			pr, ok := parseFaidxRegion(idx, c.region)
			if ok != c.wantOK {
				t.Fatalf("ok: got %v want %v", ok, c.wantOK)
			}
			if !ok {
				return
			}
			if pr.name != c.wantName || pr.beg0 != c.wantBeg || pr.end0 != c.wantEnd {
				t.Errorf("got {%s %d %d}, want {%s %d %d}", pr.name, pr.beg0, pr.end0, c.wantName, c.wantBeg, c.wantEnd)
			}
		})
	}
}

// TestParseMarkStrand checks the strand-marker configuration variants.
func TestParseMarkStrand(t *testing.T) {
	cases := []struct {
		typ     string
		pos     string
		neg     string
		wantErr bool
	}{
		{"rc", "", "/rc", false},
		{"no", "", "", false},
		{"sign", "(+)", "(-)", false},
		{"custom,_f,_r", "_f", "_r", false},
		{"custom,_only", "_only", "", false},
		{"bogus", "", "", true},
	}
	for _, c := range cases {
		t.Run(c.typ, func(t *testing.T) {
			opts := DefaultFaidxOptions(FaidxFASTA)
			err := ParseMarkStrand(&opts, c.typ)
			if (err != nil) != c.wantErr {
				t.Fatalf("err: got %v, wantErr %v", err, c.wantErr)
			}
			if c.wantErr {
				return
			}
			if opts.PosStrandName != c.pos || opts.NegStrandName != c.neg {
				t.Errorf("markers: got (%q,%q), want (%q,%q)", opts.PosStrandName, opts.NegStrandName, c.pos, c.neg)
			}
		})
	}
}

// TestReverseComplementInPlace checks the in-place RC against htslib's
// comp_base semantics (including the centre base on odd lengths).
func TestReverseComplementInPlace(t *testing.T) {
	cases := []struct{ in, want string }{
		{"ACGT", "ACGT"},
		{"AAAA", "TTTT"},
		{"acgt", "acgt"},
		{"ACGTN", "NACGT"},
		{"A", "T"},
		{"", ""},
	}
	for _, c := range cases {
		b := []byte(c.in)
		reverseComplementInPlace(b)
		if string(b) != c.want {
			t.Errorf("rc(%q) = %q, want %q", c.in, b, c.want)
		}
	}
}

// TestFaidxFastqRoundTrip builds a FASTQ index and extracts, checking the
// 6-column index and the per-region output (sequence + quality) in process.
func TestFaidxFastqRoundTrip(t *testing.T) {
	const fq = "@read1\nACGTACGTAC\nGTACGTACGT\n+\nIIIIIIIIII\nJJJJJJJJJJ\n" +
		"@read2 d\nTTTTGGGG\n+\nFFFFFFFF\n"
	dir := t.TempDir()
	p := filepath.Join(dir, "reads.fq")
	if err := os.WriteFile(p, []byte(fq), 0o644); err != nil {
		t.Fatalf("write fastq: %v", err)
	}
	if err := FaidxBuild(p, DefaultFaidxOptions(FaidxFASTQ)); err != nil {
		t.Fatalf("FaidxBuild fastq: %v", err)
	}
	fai, err := os.ReadFile(p + ".fai")
	if err != nil {
		t.Fatalf("read fastq fai: %v", err)
	}
	// read1: 20 bases, seq offset 7, 10 bases/line, 11 bytes/line, qual offset
	// = 7 + 22 (seq incl newlines) + 2 ("+\n") = 31. read2 header "@read2 d\n"
	// (9 bytes) begins at 53, so read2 seq offset = 62, qual offset = 73.
	wantFai := "read1\t20\t7\t10\t11\t31\n" +
		"read2\t8\t62\t8\t9\t73\n"
	if string(fai) != wantFai {
		t.Errorf("fastq fai:\n got: %q\nwant: %q", fai, wantFai)
	}

	opts := DefaultFaidxOptions(FaidxFASTQ)
	var out, errBuf bytes.Buffer
	if exit := Faidx(p, []string{"read1:3-12"}, opts, &out, &errBuf); exit != 0 {
		t.Fatalf("extract exit=%d stderr=%s", exit, errBuf.String())
	}
	want := "@read1:3-12\nGTACGTACGT\n+\nIIIIIIIIJJ\n"
	if out.String() != want {
		t.Errorf("fastq extract:\n got: %q\nwant: %q", out.String(), want)
	}
}

// ---- live upstream parity -------------------------------------------------

// TestFaidxUpstreamParity diffs the Go port against the live upstream
// `samtools faidx` binary on a temp FASTA, covering index build and a region
// extraction matrix. It builds the upstream binary on demand (skipping only
// in -short mode).
func TestFaidxUpstreamParity(t *testing.T) {
	bin := upstreamSamtools(t)

	dir := t.TempDir()
	usPath := filepath.Join(dir, "us.fa")
	ourPath := filepath.Join(dir, "our.fa")
	if err := os.WriteFile(usPath, []byte(smallFasta), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(ourPath, []byte(smallFasta), 0o644); err != nil {
		t.Fatal(err)
	}

	// Index build parity.
	if out, err := exec.Command(bin, "faidx", usPath).CombinedOutput(); err != nil {
		t.Fatalf("upstream faidx build: %v\n%s", err, out)
	}
	if err := FaidxBuild(ourPath, DefaultFaidxOptions(FaidxFASTA)); err != nil {
		t.Fatalf("our FaidxBuild: %v", err)
	}
	usFai, _ := os.ReadFile(usPath + ".fai")
	ourFai, _ := os.ReadFile(ourPath + ".fai")
	if !bytes.Equal(usFai, ourFai) {
		t.Errorf(".fai differs:\n upstream: %q\n ours:     %q", usFai, ourFai)
	}

	// Extraction parity matrix (stdout only — the upstream binary mixes htslib
	// log lines onto stderr that we compare separately in the CLI tests).
	matrix := [][]string{
		{"chr1:4-4"},
		{"chr1:1-25"},
		{"chr2"},
		{"chr1:0"},
		{"chr1:-5"},
		{"chr1:20-40"},
		{"chr1:1-5", "chr2:1-4"},
		{"-i", "chr1:1-10"},
		{"-n", "0", "chr1:1-25"},
		{"-n", "5", "chr2:1-20"},
		{"--mark-strand", "sign", "-i", "chr1:1-10"},
	}
	for _, args := range matrix {
		t.Run(strings.Join(args, "_"), func(t *testing.T) {
			usArgs := append([]string{"faidx", usPath}, args...)
			usOut, _ := exec.Command(bin, usArgs...).Output()

			// Drive our port through the same flag surface by parsing the
			// leading flags out of args.
			ourOut := ourFaidxOut(t, ourPath, args)
			if !bytes.Equal(usOut, []byte(ourOut)) {
				t.Errorf("args %v:\n upstream: %q\n ours:     %q", args, usOut, ourOut)
			}
		})
	}
}

// ourFaidxOut runs our in-process Faidx for a flag+region argument vector
// (the same form passed to the upstream CLI) and returns stdout.
func ourFaidxOut(t *testing.T, path string, args []string) string {
	t.Helper()
	opts := DefaultFaidxOptions(FaidxFASTA)
	var regions []string
	for i := 0; i < len(args); i++ {
		switch args[i] {
		case "-i":
			opts.ReverseComplement = true
		case "-n":
			i++
			n := 0
			neg := false
			s := args[i]
			j := 0
			if j < len(s) && s[j] == '-' {
				neg = true
				j++
			}
			for ; j < len(s); j++ {
				n = n*10 + int(s[j]-'0')
			}
			if neg {
				n = -n
			}
			opts.LineLen = n
		case "--mark-strand":
			i++
			if err := ParseMarkStrand(&opts, args[i]); err != nil {
				t.Fatalf("mark-strand: %v", err)
			}
		default:
			regions = append(regions, args[i])
		}
	}
	var out, errBuf bytes.Buffer
	Faidx(path, regions, opts, &out, &errBuf)
	return out.String()
}
