package bcf

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// upstreamBcftoolsFmtKey locates (building if necessary) the bcftools binary
// vendored under reference_code/bcftools. The build runs at most once per test
// process. These tests assert FORMAT-key decode parity against live upstream
// bcftools output, so an unavailable binary is a hard failure, never a skip.
var (
	upstreamBcftoolsOnce sync.Once
	upstreamBcftoolsPath string
	upstreamBcftoolsErr  error
)

func upstreamBcftoolsFmtKey(t *testing.T) string {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamBcftoolsOnce.Do(func() {
		bcfDir, err := filepath.Abs("../../../reference_code/bcftools")
		if err != nil {
			upstreamBcftoolsErr = err
			return
		}
		bin := filepath.Join(bcfDir, "bcftools")
		if _, statErr := os.Stat(bin); statErr == nil {
			upstreamBcftoolsPath = bin
			return
		}
		htslibDir, err := filepath.Abs("../../../reference_code/htslib")
		if err != nil {
			upstreamBcftoolsErr = err
			return
		}
		if _, statErr := os.Stat(filepath.Join(bcfDir, "config.mk")); statErr != nil {
			for _, args := range [][]string{
				{"autoheader"},
				{"autoconf"},
				{"./configure", "--with-htslib=" + htslibDir},
			} {
				cmd := exec.Command(args[0], args[1:]...)
				cmd.Dir = bcfDir
				if out, runErr := cmd.CombinedOutput(); runErr != nil {
					upstreamBcftoolsErr = fmt.Errorf("%v: %v\n%s", args, runErr, out)
					return
				}
			}
		}
		cmd := exec.Command("make", "-j4", "bcftools")
		cmd.Dir = bcfDir
		if out, runErr := cmd.CombinedOutput(); runErr != nil {
			upstreamBcftoolsErr = fmt.Errorf("make bcftools: %v\n%s", runErr, out)
			return
		}
		upstreamBcftoolsPath = bin
	})
	if upstreamBcftoolsErr != nil {
		t.Skipf("locating/building upstream bcftools: %v", upstreamBcftoolsErr)
	}
	if upstreamBcftoolsPath == "" {
		t.Skipf("upstream bcftools not available")
	}
	return upstreamBcftoolsPath
}

// fmtKeyCase is one VCF body exercising a FORMAT-key edge case. The header is
// shared and supplied per-case so the dictionary indices for the FORMAT keys
// can be pushed high (forcing int16-encoded keys) where relevant.
type fmtKeyCase struct {
	name   string
	header string
	body   string
}

// sharedHeader is a small multi-FORMAT, multi-sample header used by most cases.
const sharedHeader = `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="x">
##FILTER=<ID=q10,Description="x">
##FILTER=<ID=LowQual,Description="x">
##contig=<ID=1,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="d">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="pl">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="ad">
##FORMAT=<ID=HQ,Number=.,Type=Float,Description="hq">
##FORMAT=<ID=FT,Number=1,Type=String,Description="ft">
##FORMAT=<ID=BC,Number=1,Type=String,Description="bc">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
`

// bigDictHeader pads the dictionary with 130 INFO tags before the FORMAT keys,
// forcing the BCF encoder to emit the FORMAT keys as int16-typed dictionary
// indices (size class 2) rather than int8 — the exact path the historical
// "only the first FORMAT key resolves" bug touched.
func bigDictHeader() string {
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	b.WriteString("##FILTER=<ID=PASS,Description=\"x\">\n")
	b.WriteString("##contig=<ID=1,length=100000>\n")
	for i := 0; i < 130; i++ {
		fmt.Fprintf(&b, "##INFO=<ID=I%d,Number=1,Type=Integer,Description=\"x\">\n", i)
	}
	b.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n")
	b.WriteString("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"dp\">\n")
	b.WriteString("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"gq\">\n")
	b.WriteString("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"pl\">\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
	return b.String()
}

func fmtKeyCases() []fmtKeyCase {
	return []fmtKeyCase{
		{
			name:   "many_numeric_keys",
			header: sharedHeader,
			body: "1\t100\t.\tA\tG\t50\tPASS\tDP=20\tGT:DP:GQ:PL\t0/1:10:30:50,0,60\t1/1:12:40:80,10,0\n" +
				"1\t200\t.\tC\tT\t60\tPASS\tDP=22\tGT:DP:GQ:PL\t0/0:15:35:0,20,90\t0/1:18:45:30,0,70\n",
		},
		{
			name:   "string_fields_mixed_ploidy",
			header: sharedHeader,
			body: "1\t300\t.\tA\tG\t50\tPASS\t.\tGT:FT:DP\t0/1:PASS:10\t1/1:LowQual:12\n" +
				"1\t400\t.\tC\tT\t60\tPASS\t.\tGT:FT:DP\t0:q10:15\t0/1/1:.:18\n",
		},
		{
			name:   "ragged_vectors_and_floats",
			header: sharedHeader,
			body:   "1\t500\t.\tA\tG\t50\tPASS\t.\tGT:AD:HQ\t0/1:10,5:1.5\t1/1:3,4,9:2.0,3.5,4.5\n",
		},
		{
			name:   "phased_and_all_missing_sample",
			header: sharedHeader,
			body:   "1\t600\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0|1:10\t.:.\n",
		},
		{
			name:   "ragged_per_sample_strings",
			header: sharedHeader,
			body:   "1\t700\t.\tA\tG\t50\tPASS\t.\tGT:BC\t0/1:AAA\t1/1:B\n",
		},
		{
			name:   "int16_encoded_keys",
			header: bigDictHeader(),
			body: "1\t100\t.\tA\tG\t50\tPASS\tI0=5\tGT:DP:GQ:PL\t0/1:10:30:1,2,3\t1/1:12:40:4,5,6\n" +
				"1\t200\t.\tC\tT\t60\tPASS\tI1=7\tGT:DP:GQ:PL\t0/0:15:35:7,8,9\t0/1:18:45:10,11,12\n",
		},
	}
}

// TestBCF_FormatKeyParity verifies, for each edge case, that our BCF reader
// reconstructs every FORMAT key and per-sample value such that the resulting
// VCF text matches what upstream bcftools emits from the same BCF — pinning the
// historical "only the first FORMAT key resolves" regression as fixed.
func TestBCF_FormatKeyParity(t *testing.T) {
	bcftools := upstreamBcftoolsFmtKey(t)
	for _, tc := range fmtKeyCases() {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			vcfPath := filepath.Join(dir, "in.vcf")
			if err := os.WriteFile(vcfPath, []byte(tc.header+tc.body), 0o644); err != nil {
				t.Fatalf("write vcf: %v", err)
			}
			bcfPath := filepath.Join(dir, "in.bcf")
			if out, err := exec.Command(bcftools, "view", "-O", "b", "-o", bcfPath, vcfPath).CombinedOutput(); err != nil {
				t.Fatalf("bcftools view -O b: %v\n%s", err, out)
			}

			// Upstream VCF body (excluding ## header lines) is the oracle.
			wantOut, err := exec.Command(bcftools, "view", bcfPath).Output()
			if err != nil {
				t.Fatalf("bcftools view: %v", err)
			}
			wantBody := bodyLines(string(wantOut))

			// Our reader → VCF body.
			f, err := iohelper.OpenReader(bcfPath)
			if err != nil {
				t.Fatalf("open bcf: %v", err)
			}
			r, err := NewReader(f)
			if err != nil {
				t.Fatalf("bcf NewReader: %v", err)
			}
			recs, err := r.ReadAll()
			if err != nil {
				t.Fatalf("ReadAll: %v", err)
			}
			var gotBody []string
			for _, rec := range recs {
				v := rec.ToVariant(r.Header())
				gotBody = append(gotBody, variantToVCF(v, r.Header().Samples))
				// Hard assertion: every FORMAT key must resolve to a header tag.
				for i, k := range rec.FmtKeys {
					if r.Header().FmtTag(k) == nil {
						t.Fatalf("FORMAT key %d (index %d) did not resolve to a header tag", i, k)
					}
				}
			}
			if len(gotBody) != len(wantBody) {
				t.Fatalf("record count mismatch: got %d, want %d\ngot=%v\nwant=%v", len(gotBody), len(wantBody), gotBody, wantBody)
			}
			for i := range wantBody {
				if gotBody[i] != wantBody[i] {
					t.Fatalf("record %d mismatch:\n got=%q\nwant=%q", i, gotBody[i], wantBody[i])
				}
			}
		})
	}
}

// bodyLines returns the non-header (non-"##", non-"#CHROM") lines.
func bodyLines(s string) []string {
	var out []string
	for _, ln := range strings.Split(strings.TrimRight(s, "\n"), "\n") {
		if strings.HasPrefix(ln, "#") {
			continue
		}
		if ln == "" {
			continue
		}
		out = append(out, ln)
	}
	return out
}

// variantToVCF renders a vcf.Variant body line (CHROM..FORMAT + samples) in the
// canonical bcftools layout so it can be compared against upstream output.
func variantToVCF(v *vcf.Variant, samples []string) string {
	cols := []string{
		v.Chrom,
		fmt.Sprint(v.Pos),
		v.ID,
		v.Ref,
		joinOrDot(v.Alt, ","),
		qualStr(v.Qual),
		strings.Join(v.Filter, ";"),
		infoStr(v),
	}
	if len(v.Format) > 0 {
		cols = append(cols, strings.Join(v.Format, ":"))
		for _, name := range samples {
			cols = append(cols, sampleStr(v, name))
		}
	}
	return strings.Join(cols, "\t")
}

func joinOrDot(xs []string, sep string) string {
	if len(xs) == 0 {
		return "."
	}
	return strings.Join(xs, sep)
}

func qualStr(q float64) string {
	if q < 0 {
		return "."
	}
	if q == float64(int64(q)) {
		return fmt.Sprint(int64(q))
	}
	return fmt.Sprintf("%g", q)
}

func infoStr(v *vcf.Variant) string {
	if len(v.InfoOrder) == 0 {
		return "."
	}
	var parts []string
	for _, k := range v.InfoOrder {
		val := v.Info[k]
		if val == "" {
			parts = append(parts, k)
		} else {
			parts = append(parts, k+"="+val)
		}
	}
	if len(parts) == 0 {
		return "."
	}
	return strings.Join(parts, ";")
}

func sampleStr(v *vcf.Variant, name string) string {
	for _, s := range v.Samples {
		if s.Name == name {
			var vals []string
			for _, key := range v.Format {
				val, ok := s.Data[key]
				if !ok || val == "" {
					val = "."
				}
				vals = append(vals, val)
			}
			return strings.Join(vals, ":")
		}
	}
	return "."
}
