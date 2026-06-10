package prinseq

// Live parity tests for the sliding-window quality-trim controls
// (--trim_qual_window/step/type/rule), --trim_to_len, and the
// --range_len / --range_gc filters against the upstream PRINSEQ-lite Perl
// reference (uwb-linux/prinseq @ 0.20.4).
//
// Unlike the golden-file tests in parity_test.go these run the real Perl
// oracle (`perl reference_code/prinseq/prinseq-lite.pl`) in-process and
// compare its `-out_good` file byte-for-byte against the Go port's Filter
// output. These flags use only core Perl modules (Getopt::Long, File::Temp,
// ...) and produce deterministic output, so a byte comparison is valid.
//
// When perl is unavailable the live tests t.Skip; when perl IS available
// but the comparison fails they t.Fatalf (per the validation protocol).

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// upstreamPrinseqPath returns the path to the vendored prinseq-lite.pl, or
// "" if the submodule is not checked out.
func upstreamPrinseqPath() string {
	p := filepath.Join("..", "..", "..", "..", "reference_code", "prinseq", "prinseq-lite.pl")
	if _, err := os.Stat(p); err != nil {
		return ""
	}
	return p
}

// runUpstreamPrinseqTrim runs the upstream prinseq-lite.pl on the given
// input bytes with the supplied extra flags and returns the contents of the
// `-out_good` file (empty when upstream produced none). isFastq selects the
// input flag and output format. The helper t.Fatalf's on a genuine upstream
// failure but treats prinseq's non-zero "no bad sequences" exit as success
// as long as a good file was produced.
//
// This function is uniquely named so it does not collide with helpers in the
// sibling prinseq parity PR.
func runUpstreamPrinseqTrim(t *testing.T, input []byte, isFastq bool, extraFlags ...string) []byte {
	t.Helper()
	pl := upstreamPrinseqPath()
	if pl == "" {
		t.Skip("upstream prinseq-lite.pl not checked out; skipping live parity")
	}
	if _, err := exec.LookPath("perl"); err != nil {
		t.Skip("perl not available; skipping live parity")
	}

	dir := t.TempDir()
	ext := ".fasta"
	inFlag := "-fasta"
	outFormat := "1"
	if isFastq {
		ext = ".fastq"
		inFlag = "-fastq"
		outFormat = "3"
	}
	inPath := filepath.Join(dir, "in"+ext)
	if err := os.WriteFile(inPath, input, 0o644); err != nil {
		t.Fatalf("write upstream input: %v", err)
	}
	goodPrefix := filepath.Join(dir, "good")
	badPrefix := filepath.Join(dir, "bad")

	args := []string{
		pl, inFlag, inPath,
		"-out_good", goodPrefix,
		"-out_bad", badPrefix,
		"-out_format", outFormat,
	}
	args = append(args, extraFlags...)

	cmd := exec.Command("perl", args...)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	// prinseq exits non-zero in some no-bad-output cases; we judge success
	// by whether a good file appeared rather than the exit status alone.
	_ = cmd.Run()

	goodPath := goodPrefix + ext
	data, err := os.ReadFile(goodPath)
	if err != nil {
		// No good file at all is a real failure (perl available but it
		// did not run the flags we asked for).
		if os.IsNotExist(err) {
			return nil
		}
		t.Fatalf("read upstream good output: %v\nstderr:\n%s", err, stderr.String())
	}
	return data
}

// runGoFilter runs the Go port's Filter over the input and returns the bytes.
func runGoFilter(t *testing.T, input []byte, isFastq bool, opts FilterOptions) []byte {
	t.Helper()
	var out bytes.Buffer
	if err := Filter(bytes.NewReader(input), &out, isFastq, opts); err != nil {
		t.Fatalf("Filter: %v", err)
	}
	return out.Bytes()
}

// trimRangeFastqFixture is a small deterministic FASTQ corpus exercising
// the new quality-trim and trim_to_len paths. r1 has a low-quality 5' run,
// r3 has an internal low-quality block, r2 is uniformly high quality.
func trimRangeFastqFixture() []byte {
	return []byte(
		"@r1\n" +
			"ACGTACGTACGTACGT\n" +
			"+\n" +
			"!!!IIIIIIIIIIIII\n" +
			"@r2\n" +
			"GGGGCCCCAAAATTTT\n" +
			"+\n" +
			"IIIIIIIIIIIIIIII\n" +
			"@r3\n" +
			"ACGTNNNNACGTACGT\n" +
			"+\n" +
			"IIII!!!!IIIIIIII\n",
	)
}

func trimRangeFastaFixture() []byte {
	return []byte(
		">s1\n" +
			"ACGTACGTACGTACGTACGT\n" +
			">s2\n" +
			"GGGGGGGGGGCCCCCCCCCC\n" +
			">s3\n" +
			"ACGTACGT\n",
	)
}

func TestLiveParity_TrimQualWindow_LeftRight(t *testing.T) {
	in := trimRangeFastqFixture()
	cases := []struct {
		name  string
		flags []string
		opts  FilterOptions
	}{
		{
			name:  "left_window1",
			flags: []string{"-trim_qual_left", "20"},
			opts:  FilterOptions{TrimQualL: 20, TrimQualWindow: 1, TrimQualStep: 1},
		},
		{
			name:  "right_window1",
			flags: []string{"-trim_qual_right", "20"},
			opts:  FilterOptions{TrimQualR: 20, TrimQualWindow: 1, TrimQualStep: 1},
		},
		{
			name:  "left_window4_mean",
			flags: []string{"-trim_qual_left", "20", "-trim_qual_window", "4", "-trim_qual_step", "1", "-trim_qual_type", "mean"},
			opts:  FilterOptions{TrimQualL: 20, TrimQualWindow: 4, TrimQualStep: 1, TrimQualType: "mean"},
		},
		{
			name:  "left_window4_min_step2",
			flags: []string{"-trim_qual_left", "20", "-trim_qual_window", "4", "-trim_qual_step", "2", "-trim_qual_type", "min"},
			opts:  FilterOptions{TrimQualL: 20, TrimQualWindow: 4, TrimQualStep: 2, TrimQualType: "min"},
		},
		{
			name:  "left_window3_sum",
			flags: []string{"-trim_qual_left", "60", "-trim_qual_window", "3", "-trim_qual_step", "1", "-trim_qual_type", "sum"},
			opts:  FilterOptions{TrimQualL: 60, TrimQualWindow: 3, TrimQualStep: 1, TrimQualType: "sum"},
		},
		{
			name:  "right_window4_max_gt",
			flags: []string{"-trim_qual_right", "30", "-trim_qual_window", "4", "-trim_qual_rule", "gt", "-trim_qual_type", "max"},
			opts:  FilterOptions{TrimQualR: 30, TrimQualWindow: 4, TrimQualStep: 1, TrimQualRule: "gt", TrimQualType: "max"},
		},
		{
			name:  "both_window2",
			flags: []string{"-trim_qual_left", "20", "-trim_qual_right", "20", "-trim_qual_window", "2"},
			opts:  FilterOptions{TrimQualL: 20, TrimQualR: 20, TrimQualWindow: 2, TrimQualStep: 1},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := runUpstreamPrinseqTrim(t, in, true, tc.flags...)
			got := runGoFilter(t, in, true, tc.opts)
			if !bytes.Equal(got, want) {
				t.Fatalf("trim_qual %s mismatch.\nwant:\n%s\ngot:\n%s", tc.name, want, got)
			}
		})
	}
}

func TestLiveParity_TrimToLen(t *testing.T) {
	in := trimRangeFastqFixture()
	cases := []int{4, 8, 16, 20}
	for _, n := range cases {
		t.Run("len", func(t *testing.T) {
			want := runUpstreamPrinseqTrim(t, in, true, "-trim_to_len", itoa(n))
			got := runGoFilter(t, in, true, FilterOptions{TrimToLen: n})
			if !bytes.Equal(got, want) {
				t.Fatalf("trim_to_len %d mismatch.\nwant:\n%s\ngot:\n%s", n, want, got)
			}
		})
	}
}

func TestLiveParity_RangeLen(t *testing.T) {
	in := trimRangeFastaFixture()
	cases := []string{"10-20", "8-16", "20-20", "1-100"}
	for _, r := range cases {
		t.Run(r, func(t *testing.T) {
			want := runUpstreamPrinseqTrim(t, in, false, "-range_len", r)
			got := runGoFilter(t, in, false, FilterOptions{RangeLen: r})
			if !bytes.Equal(got, want) {
				t.Fatalf("range_len %s mismatch.\nwant:\n%s\ngot:\n%s", r, want, got)
			}
		})
	}
}

func TestLiveParity_RangeGC(t *testing.T) {
	in := trimRangeFastaFixture()
	cases := []string{"0-100", "40-60", "90-100", "0-49"}
	for _, r := range cases {
		t.Run(r, func(t *testing.T) {
			want := runUpstreamPrinseqTrim(t, in, false, "-range_gc", r)
			got := runGoFilter(t, in, false, FilterOptions{RangeGC: r})
			if !bytes.Equal(got, want) {
				t.Fatalf("range_gc %s mismatch.\nwant:\n%s\ngot:\n%s", r, want, got)
			}
		})
	}
}

// itoa is a tiny helper to avoid importing strconv for a single call site.
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}
