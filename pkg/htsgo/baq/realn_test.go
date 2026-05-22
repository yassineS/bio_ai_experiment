package baq

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// loadRefs reads every contig of a FASTA file into a name->sequence map. The
// realn test FASTAs are tiny, so reading them whole is fine.
func loadRefs(t *testing.T, path string) map[string][]byte {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	recs, err := fasta.NewReader(f).ReadAll()
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	m := make(map[string][]byte, len(recs))
	for _, r := range recs {
		m[r.ID] = r.Sequence
	}
	return m
}

// runRealn replays an input SAM file through SamProbRealn with the given flag
// and returns the resulting SAM text, mirroring htslib's test_realn driver
// (which only calls sam_prob_realn — it does not touch MD/NM).
func runRealn(t *testing.T, samPath, faPath string, flag int) string {
	t.Helper()
	refs := loadRefs(t, faPath)

	in, err := os.Open(samPath)
	if err != nil {
		t.Fatalf("open %s: %v", samPath, err)
	}
	defer in.Close()

	r, err := sam.NewReader(in)
	if err != nil {
		t.Fatalf("sam reader: %v", err)
	}
	var buf bytes.Buffer
	w := sam.NewSAMWriter(&buf)
	if err := w.WriteHeader(r.Header()); err != nil {
		t.Fatalf("write header: %v", err)
	}
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		if rec.RName != "" && rec.RName != "*" {
			ref := refs[rec.RName]
			res := SamProbRealn(rec, ref, flag)
			if res < -3 {
				t.Fatalf("SamProbRealn failed (%d) for read %q", res, rec.QName)
			}
		}
		if err := w.Write(rec); err != nil {
			t.Fatalf("write record: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("close writer: %v", err)
	}
	return buf.String()
}

// TestSamProbRealnGoldens diffs our SamProbRealn output byte-for-byte against
// htslib's realn0* golden fixtures. The flag combinations match test.pl's
// test_realn() invocations: plain -r (flag 0, computes BQ), -a (FlagApply),
// and -e (FlagExtend).
func TestSamProbRealnGoldens(t *testing.T) {
	cases := []struct {
		name string
		in   string // input SAM (relative to testdata/)
		fa   string // reference FASTA
		flag int
		exp  string // expected golden SAM
	}{
		{"realn01_r", "realn01.sam", "realn01.fa", 0, "realn01_exp.sam"},
		{"realn01_a", "realn01.sam", "realn01.fa", FlagApply, "realn01_exp-a.sam"},
		{"realn01_e", "realn01.sam", "realn01.fa", FlagExtend, "realn01_exp-e.sam"},
		{"realn02_r", "realn02.sam", "realn02.fa", 0, "realn02_exp.sam"},
		{"realn02_a", "realn02.sam", "realn02.fa", FlagApply, "realn02_exp-a.sam"},
		{"realn02_e", "realn02.sam", "realn02.fa", FlagExtend, "realn02_exp-e.sam"},
		{"realn03_e", "realn03.sam", "realn03.fa", FlagExtend, "realn03_exp.sam"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := runRealn(t,
				filepath.Join("testdata", tc.in),
				filepath.Join("testdata", tc.fa),
				tc.flag)
			wantBytes, err := os.ReadFile(filepath.Join("testdata", tc.exp))
			if err != nil {
				t.Fatalf("read golden %s: %v", tc.exp, err)
			}
			want := string(wantBytes)
			if got != want {
				t.Errorf("output mismatch for %s\n%s", tc.exp, firstDiff(got, want))
			}
		})
	}
}

// TestSamProbRealnRedo checks the FlagRedo path: realn02-r.sam carries a
// pre-existing (now stale) BQ tag; a plain recompute (FlagRedo) must discard
// it and reproduce realn02_exp.sam.
func TestSamProbRealnRedo(t *testing.T) {
	got := runRealn(t,
		filepath.Join("testdata", "realn02-r.sam"),
		filepath.Join("testdata", "realn02.fa"),
		FlagRedo)
	wantBytes, err := os.ReadFile(filepath.Join("testdata", "realn02_exp.sam"))
	if err != nil {
		t.Fatalf("read golden: %v", err)
	}
	if got != string(wantBytes) {
		t.Errorf("redo output mismatch\n%s", firstDiff(got, string(wantBytes)))
	}
}

// firstDiff returns a human-readable description of the first line where two
// SAM texts differ.
func firstDiff(got, want string) string {
	gl := strings.Split(got, "\n")
	wl := strings.Split(want, "\n")
	for i := 0; i < len(gl) || i < len(wl); i++ {
		var g, w string
		if i < len(gl) {
			g = gl[i]
		}
		if i < len(wl) {
			w = wl[i]
		}
		if g != w {
			return "line " + itoa(i+1) + ":\n  got:  " + g + "\n  want: " + w
		}
	}
	return "(no line-level difference found)"
}

// itoa is a tiny strconv.Itoa stand-in to keep the import list minimal.
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var b [20]byte
	i := len(b)
	for n > 0 {
		i--
		b[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		b[i] = '-'
	}
	return string(b[i:])
}
