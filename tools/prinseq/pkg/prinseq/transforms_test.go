package prinseq

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// --- Unit tests for the sequence/header transforms and misc knobs ---

func TestApplySeqCase(t *testing.T) {
	tests := []struct {
		mode, in, want string
	}{
		{"upper", "acGTn", "ACGTN"},
		{"lower", "ACgtN", "acgtn"},
		{"", "AcGt", "AcGt"},
		{"bogus", "AcGt", "AcGt"},
	}
	for _, tt := range tests {
		if got := applySeqCase(tt.in, tt.mode); got != tt.want {
			t.Errorf("applySeqCase(%q,%q)=%q want %q", tt.in, tt.mode, got, tt.want)
		}
	}
}

func TestApplyDNARNA(t *testing.T) {
	tests := []struct {
		mode, in, want string
	}{
		{"rna", "ACGTacgt", "ACGUacgu"},
		{"dna", "ACGUacgu", "ACGTacgt"},
		{"rna", "NNNN", "NNNN"},
		{"", "ACGT", "ACGT"},
	}
	for _, tt := range tests {
		if got := applyDNARNA(tt.in, tt.mode); got != tt.want {
			t.Errorf("applyDNARNA(%q,%q)=%q want %q", tt.in, tt.mode, got, tt.want)
		}
	}
}

func TestSplitJoinHeader(t *testing.T) {
	tests := []struct {
		desc, id, comment string
	}{
		{"read1 sample=A more", "read1", "sample=A more"},
		{"read1", "read1", ""},
		{"read1   trailing", "read1", "trailing"},
	}
	for _, tt := range tests {
		id, comment := splitHeader(tt.desc)
		if id != tt.id || comment != tt.comment {
			t.Errorf("splitHeader(%q)=(%q,%q) want (%q,%q)", tt.desc, id, comment, tt.id, tt.comment)
		}
		// joinHeader(splitHeader) collapses internal whitespace to a single
		// space, which matches upstream's `$sid.' '.$header`.
	}
	if got := joinHeader("id", ""); got != "id" {
		t.Errorf("joinHeader empty comment = %q want id", got)
	}
	if got := joinHeader("id", "c d"); got != "id c d" {
		t.Errorf("joinHeader = %q want 'id c d'", got)
	}
}

func TestParseCustomParams(t *testing.T) {
	rules := ParseCustomParams("AT 10 ; A 70% ; bad ; GG xyz ; ACX 5")
	if len(rules) != 2 {
		t.Fatalf("got %d rules, want 2: %#v", len(rules), rules)
	}
	if rules[0] != (CustomParam{Bases: "AT", Count: 10, Percent: false}) {
		t.Errorf("rule0=%#v", rules[0])
	}
	if rules[1] != (CustomParam{Bases: "A", Count: 70, Percent: true}) {
		t.Errorf("rule1=%#v", rules[1])
	}
}

func TestParseCustomParamsLowercaseSkipped(t *testing.T) {
	// Upstream's tr/ACGTN// counts only upper-case bases, so a lower-case
	// (or mixed-case) bases token fails the length check and is skipped.
	for _, s := range []string{"at 5", "aT 5", "Acgt 5"} {
		if rules := ParseCustomParams(s); len(rules) != 0 {
			t.Errorf("ParseCustomParams(%q) = %#v, want no rules (lowercase bases must be skipped)", s, rules)
		}
	}
	// Upper-case still parses.
	if rules := ParseCustomParams("AT 5"); len(rules) != 1 {
		t.Errorf("ParseCustomParams(\"AT 5\") = %#v, want 1 rule", rules)
	}
}

func TestShouldFilterCustomParams(t *testing.T) {
	// Repeat rule: reject when "AT" x3 = "ATATAT" appears.
	repeat := []CustomParam{{Bases: "AT", Count: 3, Percent: false}}
	if !shouldFilterCustomParams("GGATATATGG", repeat) {
		t.Error("expected repeat rule to reject ATATAT")
	}
	if shouldFilterCustomParams("GGATATGG", repeat) {
		t.Error("did not expect rejection without 3 repeats")
	}
	// Percentage rule: reject when >50% A.
	pct := []CustomParam{{Bases: "A", Count: 50, Percent: true}}
	if !shouldFilterCustomParams("AAAAAG", pct) { // 5/6 ~= 83%
		t.Error("expected percentage rule to reject A-rich seq")
	}
	if shouldFilterCustomParams("ACGTACGT", pct) { // 2/8 = 25%
		t.Error("did not expect rejection at 25% A")
	}
}

func TestFilterSeqNum(t *testing.T) {
	in := ">a\nACGT\n>b\nACGT\n>c\nACGT\n"
	var out bytes.Buffer
	if err := Filter(strings.NewReader(in), &out, false, FilterOptions{SeqNum: 2}); err != nil {
		t.Fatal(err)
	}
	got := out.String()
	if strings.Count(got, ">") != 2 {
		t.Errorf("seq_num=2 kept %d records: %q", strings.Count(got, ">"), got)
	}
	if !strings.Contains(got, ">a") || !strings.Contains(got, ">b") || strings.Contains(got, ">c") {
		t.Errorf("seq_num kept wrong records: %q", got)
	}
}

func TestFilterRmHeaderAndSeqID(t *testing.T) {
	in := ">read1 comment here\nACGT\n"
	// seq_id preserves comment.
	var out bytes.Buffer
	if err := Filter(strings.NewReader(in), &out, false, FilterOptions{SeqID: "S_"}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), ">S_1 comment here\n") {
		t.Errorf("seq_id should preserve comment: %q", out.String())
	}
	// rm_header drops comment.
	out.Reset()
	if err := Filter(strings.NewReader(in), &out, false, FilterOptions{RmHeader: true}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), ">read1\n") || strings.Contains(out.String(), "comment") {
		t.Errorf("rm_header should drop comment: %q", out.String())
	}
}

func TestFilterLineWidthFasta(t *testing.T) {
	in := ">a\nAAAACCCCGGGG\n"
	var out bytes.Buffer
	if err := Filter(strings.NewReader(in), &out, false, FilterOptions{QualLineWidth: 4}); err != nil {
		t.Fatal(err)
	}
	if out.String() != ">a\nAAAA\nCCCC\nGGGG\n" {
		t.Errorf("line_width wrap = %q", out.String())
	}
}

func TestFilterNoQualHeader(t *testing.T) {
	in := "@r1 c\nACGT\n+r1 c\nIIII\n"
	var out bytes.Buffer
	if err := Filter(strings.NewReader(in), &out, true, FilterOptions{NoQualHeader: true}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), "\n+\nIIII\n") {
		t.Errorf("no_qual_header should emit bare +: %q", out.String())
	}
}

func TestFilterSeqCaseDNARNAFastq(t *testing.T) {
	in := "@r1\nacgt\n+r1\nIIII\n"
	var out bytes.Buffer
	if err := Filter(strings.NewReader(in), &out, true, FilterOptions{SeqCase: "upper", DNARNA: "rna"}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(out.String(), "\nACGU\n") {
		t.Errorf("seq_case upper + dna_rna rna = %q", out.String())
	}
}

// --- Live parity tests against the upstream Perl prinseq-lite.pl ---

// upstreamPrinseqPath returns the path to the vendored prinseq-lite.pl, or ""
// if the submodule is not checked out.
func upstreamPrinseqPath(t *testing.T) string {
	t.Helper()
	// Walk up to the repo root looking for reference_code/prinseq.
	dir, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 8; i++ {
		p := filepath.Join(dir, "reference_code", "prinseq", "prinseq-lite.pl")
		if _, err := os.Stat(p); err == nil {
			return p
		}
		dir = filepath.Dir(dir)
	}
	return ""
}

// runUpstreamPrinseqXform runs the upstream Perl prinseq-lite.pl on the given
// FASTA or FASTQ fixture with extra flags and returns the contents of the
// good-output file. It writes the fixture into a temp dir and uses an explicit
// -out_good prefix so the output is deterministic and easy to read back. The
// helper is uniquely named to avoid colliding with the sibling prinseq PR's
// oracle helper.
func runUpstreamPrinseqXform(t *testing.T, perl, script, fixture string, isFastq bool, extra ...string) string {
	t.Helper()
	dir := t.TempDir()
	var inName, outPrefix string
	if isFastq {
		inName = filepath.Join(dir, "in.fastq")
	} else {
		inName = filepath.Join(dir, "in.fasta")
	}
	if err := os.WriteFile(inName, []byte(fixture), 0o644); err != nil {
		t.Fatal(err)
	}
	outPrefix = filepath.Join(dir, "good")

	args := []string{script}
	if isFastq {
		args = append(args, "-fastq", inName)
	} else {
		args = append(args, "-fasta", inName)
	}
	args = append(args, "-out_good", outPrefix, "-out_bad", "null")
	args = append(args, extra...)

	cmd := exec.Command(perl, args...)
	cmd.Dir = dir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream prinseq failed: %v\noutput: %s", err, out)
	}

	// Upstream names the good file <prefix>.fasta or <prefix>.fastq.
	candidates := []string{outPrefix + ".fastq", outPrefix + ".fasta"}
	for _, c := range candidates {
		if data, err := os.ReadFile(c); err == nil {
			return string(data)
		}
	}
	t.Fatalf("upstream produced no good output file under %s", outPrefix)
	return ""
}

// perlAvailable reports whether perl can run the upstream script (core modules
// like Getopt::Long, Digest::MD5, Cwd are part of every standard Perl).
func perlAvailable(t *testing.T) (perl, script string, ok bool) {
	t.Helper()
	perl, err := exec.LookPath("perl")
	if err != nil {
		return "", "", false
	}
	script = upstreamPrinseqPath(t)
	if script == "" {
		return "", "", false
	}
	return perl, script, true
}

func TestParityTransforms(t *testing.T) {
	perl, script, ok := perlAvailable(t)
	if !ok {
		t.Skip("perl or upstream prinseq submodule unavailable")
	}

	fasta := ">read1 sampleA\nacgtACGTttttaaaa\n>read2 sampleB\nGGGGccccTTTT\n"
	// A sequence longer than the 60-column default wrap so the parity check
	// exercises upstream's implicit FASTA wrap (no -line_width given).
	longFasta := ">long desc\n" + strings.Repeat("ACGT", 40) + "\n"
	fastq := "@read1 sampleA\nacgtACGT\n+read1 sampleA\nIIIIIIII\n@read2 sampleB\nGGGGcccc\n+read2 sampleB\nIIIIIIII\n"

	cases := []struct {
		name    string
		isFastq bool
		fixture string
		flags   []string
	}{
		{"seq_case_upper", false, fasta, []string{"-seq_case", "upper"}},
		{"seq_case_lower", false, fasta, []string{"-seq_case", "lower"}},
		{"dna_rna_rna", false, fasta, []string{"-dna_rna", "rna"}},
		{"dna_rna_dna", false, fasta, []string{"-dna_rna", "dna"}},
		{"rm_header_fasta", false, fasta, []string{"-rm_header"}},
		{"seq_id_fasta", false, fasta, []string{"-seq_id", "S_"}},
		{"line_width_fasta", false, fasta, []string{"-line_width", "4"}},
		{"seq_num_fasta", false, fasta, []string{"-seq_num", "1"}},
		{"seq_case_upper_fastq", true, fastq, []string{"-seq_case", "upper"}},
		{"dna_rna_rna_fastq", true, fastq, []string{"-dna_rna", "rna"}},
		{"rm_header_fastq", true, fastq, []string{"-rm_header"}},
		{"no_qual_header_fastq", true, fastq, []string{"-no_qual_header"}},
		{"seq_id_fastq", true, fastq, []string{"-seq_id", "S_"}},
		{"seq_num_fastq", true, fastq, []string{"-seq_num", "1"}},
		{"custom_params_fasta", false, ">a\nAAAAAAAAAAAA\n>b\nACGTACGTACGT\n", []string{"-custom_params", "A 50%"}},
		{"default_wrap_fasta", false, longFasta, []string{"-min_len", "1"}},
		{"line_width_zero_fasta", false, longFasta, []string{"-line_width", "0"}},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			want := runUpstreamPrinseqXform(t, perl, script, tc.fixture, tc.isFastq, tc.flags...)
			got := runGoFilter(t, tc.fixture, tc.isFastq, tc.flags)
			if got != want {
				t.Fatalf("parity diverged for %s\n--- go ---\n%q\n--- perl ---\n%q", tc.name, got, want)
			}
		})
	}
}

// runGoFilter runs the Go Filter with options derived from upstream-style
// flags, returning the good output as a string. It replicates the CLI's
// effective-line-width resolution (default 60 unless -line_width is given)
// so the parity comparison genuinely exercises upstream's default FASTA wrap.
func runGoFilter(t *testing.T, fixture string, isFastq bool, flags []string) string {
	t.Helper()
	opts := FilterOptions{}
	lineWidthSet := false
	for i := 0; i < len(flags); i++ {
		switch flags[i] {
		case "-seq_case":
			i++
			opts.SeqCase = flags[i]
		case "-dna_rna":
			i++
			opts.DNARNA = flags[i]
		case "-rm_header":
			opts.RmHeader = true
		case "-no_qual_header":
			opts.NoQualHeader = true
		case "-line_width":
			i++
			n, _ := parseUint(flags[i])
			opts.QualLineWidth = n
			lineWidthSet = true
		case "-seq_num":
			i++
			n, _ := parseUint(flags[i])
			opts.SeqNum = n
		case "-seq_id":
			i++
			opts.SeqID = flags[i]
		case "-custom_params":
			i++
			opts.CustomParams = ParseCustomParams(flags[i])
		case "-min_len":
			i++
			n, _ := parseUint(flags[i])
			opts.MinLen = n
		default:
			t.Fatalf("runGoFilter: unhandled flag %q", flags[i])
		}
	}
	// Mirror upstream $linelen default of 60 for FASTA output when
	// -line_width is not given (prinseq-lite.pl:931-939). FASTQ-only output
	// (no FASTA/QUAL here) is unaffected because the FASTQ writer never wraps.
	if !lineWidthSet && !isFastq {
		opts.QualLineWidth = 60
	}
	var out bytes.Buffer
	if err := Filter(strings.NewReader(fixture), &out, isFastq, opts); err != nil {
		t.Fatal(err)
	}
	return out.String()
}
