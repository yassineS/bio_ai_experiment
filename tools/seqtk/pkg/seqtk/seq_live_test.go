package seqtk

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// seqLiveCase describes one `seqtk seq` invocation: a key (used to name
// the committed golden file), the input fixture, and the SeqOptions that
// reproduce the flags. maskFile, when set, is loaded into the options.
type seqLiveCase struct {
	key     string
	input   string // "input.fa" or "input.fq"
	build   func(o *SeqOptions)
	mask    string // mask fixture filename, "" if none
	cliArgs []string
}

// seqLiveCases enumerates the matrix exercised by both the golden-file
// comparison and (when the genuine binary is present) the live oracle.
// The cliArgs mirror the flags exactly; build() configures the Go port.
func seqLiveCases() []seqLiveCase {
	return []seqLiveCase{
		{"fa_none", "input.fa", func(o *SeqOptions) {}, "", nil},
		{"fa_A", "input.fa", func(o *SeqOptions) { o.ForceFasta = true }, "", []string{"-A"}},
		{"fa_l4", "input.fa", func(o *SeqOptions) { o.LineLen = 4 }, "", []string{"-l", "4"}},
		{"fa_U", "input.fa", func(o *SeqOptions) { o.Uppercase = true }, "", []string{"-U"}},
		{"fa_C", "input.fa", func(o *SeqOptions) { o.DropComment = true }, "", []string{"-C"}},
		{"fa_S", "input.fa", func(o *SeqOptions) { o.StripSpace = true }, "", []string{"-S"}},
		{"fa_r", "input.fa", func(o *SeqOptions) { o.RevComp = true }, "", []string{"-r"}},
		{"fa_R", "input.fa", func(o *SeqOptions) { o.BothStrands = true }, "", []string{"-R"}},
		{"fa_N", "input.fa", func(o *SeqOptions) { o.DropAmbig = true }, "", []string{"-N"}},
		{"fa_r_l4", "input.fa", func(o *SeqOptions) { o.RevComp = true; o.LineLen = 4 }, "", []string{"-r", "-l", "4"}},
		{"fa_C_S_U", "input.fa", func(o *SeqOptions) { o.DropComment = true; o.StripSpace = true; o.Uppercase = true }, "", []string{"-C", "-S", "-U"}},
		{"fa_M_bed", "input.fa", func(o *SeqOptions) {}, "mask.bed", []string{"-M"}},
		{"fa_M_bed_nN", "input.fa", func(o *SeqOptions) { o.MaskChar = 'N' }, "mask.bed", []string{"-M", "_M_", "-n", "N"}},
		{"fa_M_bed_c", "input.fa", func(o *SeqOptions) { o.MaskComp = true }, "mask.bed", []string{"-M", "_M_", "-c"}},
		{"fa_M_names", "input.fa", func(o *SeqOptions) {}, "mask.names", []string{"-M"}},
		{"fa_x_nN", "input.fa", func(o *SeqOptions) { o.LowerToMask = true; o.MaskChar = 'N' }, "", []string{"-x", "-n", "N"}},

		{"fq_none", "input.fq", func(o *SeqOptions) {}, "", nil},
		{"fq_A", "input.fq", func(o *SeqOptions) { o.ForceFasta = true }, "", []string{"-A"}},
		{"fq_l4", "input.fq", func(o *SeqOptions) { o.LineLen = 4 }, "", []string{"-l", "4"}},
		{"fq_U", "input.fq", func(o *SeqOptions) { o.Uppercase = true }, "", []string{"-U"}},
		{"fq_C", "input.fq", func(o *SeqOptions) { o.DropComment = true }, "", []string{"-C"}},
		{"fq_S", "input.fq", func(o *SeqOptions) { o.StripSpace = true }, "", []string{"-S"}},
		{"fq_r", "input.fq", func(o *SeqOptions) { o.RevComp = true }, "", []string{"-r"}},
		{"fq_R", "input.fq", func(o *SeqOptions) { o.BothStrands = true }, "", []string{"-R"}},
		{"fq_q20", "input.fq", func(o *SeqOptions) { o.QualThres = 20 }, "", []string{"-q", "20"}},
		{"fq_q20_nN", "input.fq", func(o *SeqOptions) { o.QualThres = 20; o.MaskChar = 'N' }, "", []string{"-q", "20", "-n", "N"}},
		{"fq_X40_nN", "input.fq", func(o *SeqOptions) { o.MaxQual = 40; o.MaskChar = 'N' }, "", []string{"-X", "40", "-n", "N"}},
		{"fq_Q64", "input.fq", func(o *SeqOptions) { o.QualShift = 64 }, "", []string{"-Q", "64"}},
		{"fq_V_Q64", "input.fq", func(o *SeqOptions) { o.ShiftQual = true; o.QualShift = 64 }, "", []string{"-V", "-Q", "64"}},
		{"fq_FI", "input.fq", func(o *SeqOptions) { o.FakeQual = int('I') }, "", []string{"-F", "I"}},
		{"fq_1", "input.fq", func(o *SeqOptions) { o.Odd = true }, "", []string{"-1"}},
		{"fq_2", "input.fq", func(o *SeqOptions) { o.Even = true }, "", []string{"-2"}},
		{"fq_L10", "input.fq", func(o *SeqOptions) { o.MinLen = 10 }, "", []string{"-L", "10"}},
		{"fq_r_l4", "input.fq", func(o *SeqOptions) { o.RevComp = true; o.LineLen = 4 }, "", []string{"-r", "-l", "4"}},
	}
}

// runGoSeq applies a case's options to the Go port and returns the output.
func runGoSeq(t *testing.T, tc seqLiveCase, tdDir string) []byte {
	t.Helper()
	in, err := os.ReadFile(filepath.Join(tdDir, tc.input))
	if err != nil {
		t.Fatalf("read input %s: %v", tc.input, err)
	}
	opts := NewSeqOptions()
	tc.build(&opts)
	if tc.mask != "" {
		mf, err := os.Open(filepath.Join(tdDir, tc.mask))
		if err != nil {
			t.Fatalf("open mask %s: %v", tc.mask, err)
		}
		defer mf.Close()
		if err := opts.LoadMaskFile(mf); err != nil {
			t.Fatalf("load mask %s: %v", tc.mask, err)
		}
	}
	var out bytes.Buffer
	if err := Seq(bytes.NewReader(in), &out, opts); err != nil {
		t.Fatalf("Seq(%s): %v", tc.key, err)
	}
	return out.Bytes()
}

// TestSeqtkSeqGolden compares the Go port against committed golden outputs
// captured from the genuine seqtk 1.5-r133 binary. It always runs.
func TestSeqtkSeqGolden(t *testing.T) {
	tdDir := filepath.Join("testdata", "seq_live")
	for _, tc := range seqLiveCases() {
		t.Run(tc.key, func(t *testing.T) {
			got := runGoSeq(t, tc, tdDir)
			golden := filepath.Join(tdDir, tc.key+".golden")
			want, err := os.ReadFile(golden)
			if err != nil {
				t.Fatalf("read golden %s: %v", golden, err)
			}
			if !bytes.Equal(got, want) {
				t.Errorf("seq %s: output differs from golden\n--- got ---\n%s\n--- want ---\n%s",
					tc.key, got, want)
			}
		})
	}
}

// TestSeqtkSeqLiveOracle runs the same matrix through the genuine seqtk
// binary (if it can be located and is executable) and compares its output
// to the Go port byte-for-byte. It is skipped when the binary is absent so
// CI without submodules still passes.
func TestSeqtkSeqLiveOracle(t *testing.T) {
	bin := findGenuineSeqtk()
	if bin == "" {
		t.Skip("genuine seqtk binary not found; skipping live oracle")
	}
	tdDir := filepath.Join("testdata", "seq_live")
	for _, tc := range seqLiveCases() {
		t.Run(tc.key, func(t *testing.T) {
			args := []string{"seq"}
			for _, a := range tc.cliArgs {
				if a == "_M_" {
					args = append(args, filepath.Join(tdDir, tc.mask))
				} else {
					args = append(args, a)
				}
			}
			// For cases that pass "-M" alone, append the mask path.
			if tc.mask != "" && !containsToken(tc.cliArgs, "_M_") {
				args = append(args, filepath.Join(tdDir, tc.mask))
			}
			args = append(args, filepath.Join(tdDir, tc.input))

			cmd := exec.Command(bin, args...)
			ref, err := cmd.Output()
			if err != nil {
				t.Fatalf("genuine seqtk %v: %v", args, err)
			}
			got := runGoSeq(t, tc, tdDir)
			if !bytes.Equal(got, ref) {
				t.Errorf("seq %s: Go output differs from genuine binary\n--- go ---\n%s\n--- ref ---\n%s",
					tc.key, got, ref)
			}
		})
	}
}

func containsToken(s []string, tok string) bool {
	for _, x := range s {
		if x == tok {
			return true
		}
	}
	return false
}

// findGenuineSeqtk locates the vendored genuine seqtk binary, walking up
// from the package directory to the repo root's reference_code submodule.
func findGenuineSeqtk() string {
	candidates := []string{
		filepath.Join("..", "..", "..", "..", "reference_code", "seqtk", "seqtk"),
	}
	for _, c := range candidates {
		if fi, err := os.Stat(c); err == nil && !fi.IsDir() && fi.Mode()&0111 != 0 {
			return c
		}
	}
	if p, err := exec.LookPath("seqtk"); err == nil {
		return p
	}
	return ""
}
