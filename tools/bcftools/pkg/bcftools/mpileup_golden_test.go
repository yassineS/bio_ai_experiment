package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// mpileupGoldenDir is where the upstream `bcftools mpileup` fixtures are
// vendored: the mpileup.{1,2,3}.bam read sets, the mpileup.ref.fa
// reference (with its .fai), and the mpileup.{1,11}.out goldens.
const mpileupGoldenDir = "../../testdata/mpileup"

// mpileupFixture returns the absolute path to a vendored mpileup fixture
// and skips the test if it is missing.
func mpileupFixture(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join(mpileupGoldenDir, name))
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	if _, err := os.Stat(abs); err != nil {
		t.Skipf("vendored fixture %s missing: %v", name, err)
	}
	return abs
}

// splitMpileupVCF splits VCF text into its header lines (the leading
// "##"/"#" block) and its data records. The upstream goldens were
// produced with `--no-version`, so the header is compared verbatim.
func splitMpileupVCF(s string) (header, data []string) {
	for _, ln := range strings.Split(s, "\n") {
		if ln == "" {
			continue
		}
		if strings.HasPrefix(ln, "#") {
			header = append(header, ln)
		} else {
			data = append(data, ln)
		}
	}
	return header, data
}

// TestMpileupSNPGoldens is the slice-4 byte-for-byte parity check for
// the SNP MAQ path. With the per-site bias annotations (VDB, SGB, RPBZ,
// MQBZ, BQBZ, MQSBZ, SCBZ), the MQ0F fraction, the INFO/QS float32
// rounding and the MPLP_SMART_OVERLAPS read-pair quality merge all
// ported, `bcftools mpileup` must reproduce the upstream goldens
// exactly — header and every data record, INFO bias tags included.
//
// Two upstream invocations from reference_code/bcftools/test/test.pl
// are replayed:
//
//   - mpileup.11.out: `mpileup -a -AD mpileup.3.bam`, the full 4200 bp
//     contig. 4001 covered positions, 87 SNP ALT records, two
//     overlapping mate pairs (17:1118-1142 and 17:3785-3836) that
//     exercise smart-overlaps.
//   - mpileup.1.out: `mpileup -r17:100-150 -a -AD mpileup.{1,2,3}.bam`,
//     the three-sample multi-BAM path over a 51 bp window.
//
// `-a -AD` removes FORMAT/AD from the default tag set, so the goldens
// carry FORMAT=PL only — exactly what the SNP path emits today.
//
// The single upstream INDEL record in mpileup.11.out (17:302 TA) is the
// only golden line not reproduced: indel calling (bam2bcf_indel.c) is
// the remaining deferred mpileup work. The comparison aligns records by
// CHROM:POS and skips that INDEL position.
func TestMpileupSNPGoldens(t *testing.T) {
	ref := mpileupFixture(t, "mpileup.ref.fa")
	mpileupFixture(t, "mpileup.ref.fa.fai") // sidecar required by the FASTA reader

	cases := []struct {
		name     string
		golden   string
		inputs   []string
		regions  []string
		annotate string
	}{
		{
			name:     "single-bam-full-contig",
			golden:   "mpileup.11.out",
			inputs:   []string{mpileupFixture(t, "mpileup.3.bam")},
			annotate: "-AD", // upstream: -a -AD removes FORMAT/AD from defaults
		},
		{
			name:     "multi-bam-region",
			golden:   "mpileup.1.out",
			inputs:   []string{mpileupFixture(t, "mpileup.1.bam"), mpileupFixture(t, "mpileup.2.bam"), mpileupFixture(t, "mpileup.3.bam")},
			regions:  []string{"17:100-150"},
			annotate: "-AD",
		},
		{
			// mpileup.12.out exercises the default `-a` set, which
			// includes FORMAT/AD. The upstream test.pl invocation
			// (mpileup.1+mpileup.2+mpileup.3, region list 17:100-105)
			// covers six positions including the multi-sample SNP rows
			// at 17:103 and 17:104.
			name:     "multi-bam-region-format-AD",
			golden:   "mpileup.12.out",
			inputs:   []string{mpileupFixture(t, "mpileup.1.bam"), mpileupFixture(t, "mpileup.2.bam"), mpileupFixture(t, "mpileup.3.bam")},
			regions:  []string{"17:100-102", "17:102-103", "17:103-104", "17:104-105", "17:100-105"},
			annotate: "",
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			goldenBytes, err := os.ReadFile(mpileupFixture(t, tc.golden))
			if err != nil {
				t.Fatalf("read golden %s: %v", tc.golden, err)
			}
			wantHeader, wantData := splitMpileupVCF(string(goldenBytes))

			var buf bytes.Buffer
			opts := MpileupOptions{
				Inputs:    tc.inputs,
				FastaRef:  ref,
				Regions:   tc.regions,
				Annotate:  tc.annotate,
				NoVersion: true, // upstream goldens were made with --no-version
			}
			if err := MpileupFile(opts, &buf); err != nil {
				t.Fatalf("MpileupFile: %v", err)
			}
			gotHeader, gotData := splitMpileupVCF(buf.String())

			// Header must match byte-for-byte.
			if len(gotHeader) != len(wantHeader) {
				t.Fatalf("header line count: got %d, want %d\n got:  %v\n want: %v",
					len(gotHeader), len(wantHeader), gotHeader, wantHeader)
			}
			for i := range gotHeader {
				if gotHeader[i] != wantHeader[i] {
					t.Errorf("header line %d:\n got:  %s\n want: %s", i, gotHeader[i], wantHeader[i])
				}
			}

			// Index the golden data records by CHROM:POS, splitting off
			// the deferred INDEL records so the SNP comparison is clean.
			wantByPos := make(map[string]string, len(wantData))
			indelSkipped := 0
			for _, ln := range wantData {
				f := strings.Split(ln, "\t")
				if len(f) < 8 {
					t.Fatalf("malformed golden record %q", ln)
				}
				if strings.HasPrefix(f[7], "INDEL") || strings.Contains(f[7], ";INDEL;") {
					indelSkipped++
					continue
				}
				wantByPos[f[0]+":"+f[1]] = ln
			}

			checked, diffs := 0, 0
			for _, ln := range gotData {
				f := strings.Split(ln, "\t")
				key := f[0] + ":" + f[1]
				want, ok := wantByPos[key]
				if !ok {
					t.Errorf("emitted an extra record at %s:\n %s", key, ln)
					continue
				}
				checked++
				delete(wantByPos, key)
				if ln != want {
					diffs++
					if diffs <= 10 {
						t.Errorf("record %s mismatch:\n got:  %s\n want: %s", key, ln, want)
					}
				}
			}
			// Every non-INDEL golden record must have been produced.
			for key := range wantByPos {
				t.Errorf("missing record at %s (golden has it, we did not emit it)", key)
			}
			if diffs > 10 {
				t.Errorf("... and %d more record mismatches", diffs-10)
			}
			t.Logf("%s: %d SNP records byte-for-byte identical, %d upstream INDEL record(s) skipped (indel calling deferred)",
				tc.golden, checked, indelSkipped)
		})
	}
}

// TestMpileupGoldensDeferred documents the upstream mpileup goldens that
// still do not byte-match and the precise reason, so the remaining work
// stays visible in the test output.
func TestMpileupGoldensDeferred(t *testing.T) {
	deferred := []struct{ golden, reason string }{
		{
			"mpileup/mpileup.11.out (17:302 TA record)",
			"the single INDEL record; indel calling (bam2bcf_indel.c) is " +
				"the only remaining deferred mpileup path. Every SNP record " +
				"of this golden, bias tags included, matches byte-for-byte.",
		},
		{
			"mpileup/mpileup.2.out, mpileup.4.out, mpileup.5.out, mpileup.6.out",
			"produced with FORMAT tags beyond PL (DP, DV, DP4, SP, AD, ADF, " +
				"ADR and the gVCF mode). The per-sample FORMAT tag emission " +
				"is a separate slice; the SNP INFO columns already match.",
		},
		{
			"mpileup/mpileup.{3,7,8,9,10}.out",
			"produced with `--ff` FLAG filtering, `-s/-S/-G` sample/read-" +
				"group selection or `-t` targets — accepted flags whose " +
				"read/sample subsetting is a separate parity item.",
		},
		{
			"mpileup/iupac.1.out",
			"reference FASTA carries IUPAC ambiguity codes; rendering an " +
				"ambiguous REF base is a separate parity gap.",
		},
		{
			"mpileup/indel-AD.*.out, mpileup-SCR.out",
			"exercise indel calling and INFO/FMT SCR; both depend on the " +
				"deferred bam2bcf_indel.c port.",
		},
	}
	for _, d := range deferred {
		t.Logf("DEFERRED golden %s: %s", d.golden, d.reason)
	}
	t.Skip("documented above; these goldens are deferred to follow-up slices")
}
