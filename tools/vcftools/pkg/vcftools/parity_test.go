package vcftools

// Parity tests for vcftools against the upstream test corpus
// (reference_code/vcftools/examples/) plus hand-computed cases pinned to
// the same VCF fixture.
//
// Methodology:
//
//   - Vendor a small VCF fixture (tools/vcftools/testdata/parity/sample.vcf)
//     copied byte-for-byte from reference_code/vcftools/examples/valid-4.0.vcf,
//     which has 12 sites across 3 chromosomes and 3 samples, covering
//     biallelic SNPs, multi-allelic SNPs, indels, monomorphic sites,
//     filtered sites, and partially-missing genotypes.
//
//   - For every option in our port that produces a textual output file
//     (.frq, .ldepth*, .lqual, .lmiss, .imiss, .idepth, .hwe, .sites.pi,
//     .singletons, .TsTv.summary, .TsTv.count, .TsTv.qual, .het,
//     .recode.vcf, .012*, .ped/.map, .tped/.tfam, .geno.ld, .hap.ld), this
//     file runs the port via Run and diffs the output against a
//     hand-computed `<case>.expected.<ext>` golden file (or, where
//     upstream and our port have a known format / formula gap, asserts
//     the file shape and skips byte parity).
//
// Two notable deviations from upstream documented elsewhere:
//
//   1. `--site-pi` formula. Upstream emits a per-genotype pairwise
//      quantity; our port uses the textbook
//      `(n^2 - sum(c_a^2)) / (n*(n-1))`. See
//      docs/UPSTREAM_BUGS.md#vcftools-site-pi. Affected tests skip byte
//      parity against upstream and instead diff against the textbook
//      golden file.
//
//   2. Many small format gaps where our port emits a subset of upstream's
//      columns (e.g. .ldepth missing SUMSQ_DEPTH, .hwe missing the
//      directional P-values). These are documented in
//      docs/PARITY_ROADMAP.md#vcftools and tracked with t.Skip rather
//      than masked.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func vcftoolsFixtureDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity"))
	if err != nil {
		t.Fatalf("Abs: %v", err)
	}
	return abs
}

// runVcftoolsParity opens the named VCF fixture under
// tools/vcftools/testdata/parity/, runs the port with params, and returns
// the OutPrefix used (so the caller can read individual output files).
func runVcftoolsParity(t *testing.T, vcfName string, params *Params) string {
	t.Helper()
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "out")
	params.OutPrefix = prefix
	in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), vcfName))
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	if err := Run(in, params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return prefix
}

// readFileBytes reads a parity output file and returns its bytes.
func readFileBytes(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return b
}

// readFileLines reads a parity output file and returns its lines with the
// trailing blank line trimmed.
func readFileLines(t *testing.T, path string) []string {
	return strings.Split(strings.TrimRight(string(readFileBytes(t, path)), "\n"), "\n")
}

// assertLinesEqual fails the test with a side-by-side diff if got/want
// differ.
func assertLinesEqual(t *testing.T, got, want []string) {
	t.Helper()
	if len(got) != len(want) {
		t.Errorf("line count mismatch: got %d, want %d", len(got), len(want))
	}
	maxN := len(got)
	if len(want) > maxN {
		maxN = len(want)
	}
	for i := 0; i < maxN; i++ {
		var g, w string
		if i < len(got) {
			g = got[i]
		}
		if i < len(want) {
			w = want[i]
		}
		if g != w {
			t.Errorf("line %d:\n  got:  %q\n  want: %q", i+1, g, w)
		}
	}
}

// -----------------------------------------------------------------------------
// Filtering: --chr, --from-bp/--to-bp, --maf, --mac, --minQ, --max-missing,
// --remove-indels, --keep-only-indels
// -----------------------------------------------------------------------------

// TestParity_Chr_19 — `--chr 19` keeps only the two chr19 records.
func TestParity_Chr_19(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Chr:    "19",
		Recode: true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	var data int
	for _, ln := range got {
		if ln != "" && !strings.HasPrefix(ln, "#") {
			data++
		}
	}
	if data != 2 {
		t.Errorf("want 2 chr19 rows, got %d", data)
	}
	for _, ln := range got {
		if !strings.HasPrefix(ln, "#") && ln != "" && !strings.HasPrefix(ln, "19\t") {
			t.Errorf("non-chr19 row leaked: %q", ln)
		}
	}
}

// TestParity_FromTo_BP — `--chr 20 --from-bp 14000 --to-bp 20000` selects
// the two chr20 sites in that range (14370 and 17330).
func TestParity_FromTo_BP(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Chr:    "20",
		FromBp: 14000,
		ToBp:   20000,
		Recode: true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	want := map[string]bool{"20\t14370": true, "20\t17330": true}
	seen := map[string]bool{}
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		key := strings.Join(strings.SplitN(ln, "\t", 3)[:2], "\t")
		seen[key] = true
	}
	for k := range want {
		if !seen[k] {
			t.Errorf("missing site %q", k)
		}
	}
	for k := range seen {
		if !want[k] {
			t.Errorf("unexpected site %q", k)
		}
	}
}

// TestParity_Maf — `--maf 0.4` keeps sites whose minor allele frequency
// >= 0.4. In sample.vcf the qualifying sites are 20:14370 (G/A 3/3 — MAF
// 0.5), X:9 (A/T 3/2 — MAF 0.4), X:11 (T/<DEL:ME:ALU> 1/1 — MAF 0.5,
// only one sample has non-missing data), and X:12 (T/A 2/3 — MAF 0.4).
// We also confirm low-MAF sites like 19:111 (MAF 0.167) are dropped.
func TestParity_Maf(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Maf:    0.4,
		Recode: true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	keepers := map[string]bool{
		"20\t14370": true,
		"X\t9":      true,
		"X\t11":     true,
		"X\t12":     true,
	}
	seen := 0
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		key := strings.Join(strings.SplitN(ln, "\t", 3)[:2], "\t")
		if !keepers[key] {
			t.Errorf("unexpected --maf site: %q", key)
		}
		seen++
	}
	if seen != len(keepers) {
		t.Errorf("want %d --maf sites, got %d", len(keepers), seen)
	}
	// Sanity: a low-MAF site must be dropped.
	for _, ln := range got {
		if strings.HasPrefix(ln, "19\t111\t") {
			t.Errorf("--maf 0.4 should have dropped 19:111 (MAF 0.167)")
		}
	}
}

// TestParity_Mac — `--mac 3` keeps biallelic sites whose minor allele
// count is >= 3.
func TestParity_Mac(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Mac:    3,
		Recode: true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	seen := 0
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		seen++
	}
	if seen < 1 {
		t.Errorf("expected at least one --mac 3 record, got %d", seen)
	}
}

// TestParity_NonRefAF_03 — `--non-ref-af 0.3 --recode` byte-for-byte
// against upstream. Upstream registration: parameters.cpp:303. Filter
// logic ported from entry_filters.cpp:801-815: every ALT's frequency
// (count/non-missing-chr) must be >= 0.3, and a monomorphic site
// (ALT == "." or no ALT) is also dropped because
// N_failed == N_alleles-1 == 0 fires the line-814 fallback.
//
// For sample.vcf the expected kept sites are 20:14370 (A freq 0.5),
// 20:1110696 (G freq 0.333, T freq 0.5 — both >= 0.3), X:9 (T 0.4),
// X:12 (A 0.5). All other sites have at least one ALT below 0.3.
func TestParity_NonRefAF_03(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MinNonRefAF: 0.3,
		Recode:      true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "non_ref_af_03.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_NonRefAF_05 — `--non-ref-af 0.5 --recode` byte-for-byte.
// Pin: at threshold 0.5 only 20:14370 (A freq 0.5) and X:12 (A freq 0.5)
// survive — both rare biallelic SNPs whose single ALT exactly meets the
// threshold. 20:1110696 drops because G's freq is 0.333.
func TestParity_NonRefAF_05(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MinNonRefAF: 0.5,
		Recode:      true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "non_ref_af_05.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_NonRefAC_2_ChrX — `--chr X --non-ref-ac 2 --recode`
// byte-for-byte. Upstream registration: parameters.cpp:302. Filter
// ported from entry_filters.cpp:902-907: every ALT count must be >= 2.
// We restrict to chr X to avoid two pre-existing baseline-recode
// discrepancies on chr 20 (the port currently drops 20:1235237 and
// X:11 even without filters); those gaps are documented in
// docs/PARITY_ROADMAP.md and tracked separately. Within chr X, X:9
// (T count 2) and X:12 (A count 3) pass while X:10 and X:11 fail
// because some ALT has count < 2.
func TestParity_NonRefAC_2_ChrX(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Chr:         "X",
		MinNonRefAC: 2,
		Recode:      true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "non_ref_ac_2_chrX.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_NonRefAC_3_ChrX — `--chr X --non-ref-ac 3 --recode`
// byte-for-byte. At threshold 3 only X:12 (A count 3) survives in chr X.
func TestParity_NonRefAC_3_ChrX(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Chr:         "X",
		MinNonRefAC: 3,
		Recode:      true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "non_ref_ac_3_chrX.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_NonRefAF_DropsMonomorphic — pin the upstream-quirk
// behaviour that `--non-ref-af > 0` drops monomorphic (no-ALT) sites
// via the N_failed == N_alleles-1 fallback. Counterpart: --non-ref-ac
// does NOT drop them because the analogous fallback is gated on
// `_any`. See entry_filters.cpp:814 vs 912.
func TestParity_NonRefAF_DropsMonomorphic(t *testing.T) {
	// With --non-ref-af set very low, every real ALT passes; only the
	// upstream monomorphic fallback can drop sites here.
	prefixAF := runVcftoolsParity(t, "sample.vcf", &Params{
		MinNonRefAF: 0.01,
		Recode:      true,
	})
	gotAF := readFileLines(t, prefixAF+".recode.vcf")
	for _, ln := range gotAF {
		if strings.HasPrefix(ln, "20\t1230237") || strings.HasPrefix(ln, "20\t1235237") {
			t.Errorf("--non-ref-af 0.01 should drop monomorphic site, got: %q", ln)
		}
	}
	// Counterpart: --non-ref-ac 1 should NOT drop monomorphic sites,
	// because upstream's count branch keys the fallback on `_any` only.
	// We can only check the non-baseline-buggy mono site here:
	// 20:1230237 (the X:11 / 20:1235237 baseline gap is separate).
	prefixAC := runVcftoolsParity(t, "sample.vcf", &Params{
		MinNonRefAC: 1,
		Recode:      true,
	})
	gotAC := readFileLines(t, prefixAC+".recode.vcf")
	saw1230237 := false
	for _, ln := range gotAC {
		if strings.HasPrefix(ln, "20\t1230237") {
			saw1230237 = true
		}
	}
	if !saw1230237 {
		t.Errorf("--non-ref-ac 1 should keep monomorphic site 20:1230237 (upstream behaviour)")
	}
}

// TestParity_NonRefACAny_2 — `--non-ref-ac-any 2 --recode` byte-for-byte.
// Upstream registration: parameters.cpp:304. Filter ported from
// entry_filters.cpp:902-913: each ALT's count gets compared against the
// `_any` threshold and contributes to N_failed; the site is dropped only
// when N_failed == N_alleles-1 (every ALT failed). Monomorphic sites
// (N_alleles == 1) satisfy 0 == 0 and so the fallback at line 912 drops
// them too — keyed on the `_any` thresholds.
//
// Expected kept sites for sample.vcf: 20:14370 (A count 3), 20:1110696
// (G=2,T=4 — T passes), 20:1234567 (GA=3,GAC=1 — GA passes), X:9 (T=2),
// X:12 (A=3). Sites dropped: 19:111 (C=1), 19:112 (G=1), 20:17330 (A=1),
// 20:1230237 (mono), 20:1235237 (mono), X:10 (A=1,ATG=1 both fail), X:11
// (A=1,DEL=1 both fail).
func TestParity_NonRefACAny_2(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MinNonRefACAny: 2,
		Recode:         true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "non_ref_ac_any_2.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_NonRefACAny_1_Chr20 — `--chr 20 --non-ref-ac-any 1 --recode`
// byte-for-byte. Pins the upstream-quirk monomorphic drop for the AC
// `_any` branch (entry_filters.cpp:912) — counterpart of plain
// `--non-ref-ac` which does NOT drop monomorphic sites. We restrict to
// chr 20 because the chr-X record X:11 (T -> A,<DEL:ME:ALU>) is dropped
// by our port's `MinAlleles=2` default and would mismatch upstream; the
// only chr-20 sites kept by the filter (the four ALT-bearing records)
// pass our baseline as well. The two monomorphic chr-20 sites
// (20:1230237 and 20:1235237) are dropped by the `_any` fallback,
// matching upstream.
func TestParity_NonRefACAny_1_Chr20(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Chr:            "20",
		MinNonRefACAny: 1,
		Recode:         true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "non_ref_ac_any_1_chr20.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_MaxNonRefAF_03_Chr20 — `--chr 20 --max-non-ref-af 0.3
// --recode` byte-for-byte. Upstream registration: parameters.cpp:288.
// Filter ported from entry_filters.cpp:807 (`freq > max_non_ref_af`).
// At threshold 0.3 the only chr-20 site whose every ALT has freq ≤ 0.3
// is 20:17330 (A freq 0.167). Monomorphic chr-20 sites are dropped via
// the line-814 fallback (gated on `max_non_ref_af < 1.0`, which our
// internal nonzero MaxNonRefAF mirrors).
func TestParity_MaxNonRefAF_03_Chr20(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Chr:         "20",
		MaxNonRefAF: 0.3,
		Recode:      true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "max_non_ref_af_03_chr20.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_MaxNonRefAC_2_Chr19 — `--chr 19 --max-non-ref-ac 2 --recode`
// byte-for-byte. Upstream registration: parameters.cpp:287. Filter ported
// from entry_filters.cpp:905 (`count > max_non_ref_ac`). Both chr-19 sites
// have a single ALT with count 1 (≤ 2) so both are kept. Counterpart pin:
// plain `--max-non-ref-ac` does NOT drop monomorphic sites (the
// line-912 fallback is keyed on `_any`); chr 19 has no monomorphic
// sites so this test focuses on the plain per-ALT max-check.
func TestParity_MaxNonRefAC_2_Chr19(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Chr:         "19",
		MaxNonRefAC: 2,
		Recode:      true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "max_non_ref_ac_2_chr19.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_MaxNonRefACAny_2_Chr20 — `--chr 20 --max-non-ref-ac-any 2
// --recode` byte-for-byte. Upstream registration: parameters.cpp:289.
// Filter ported from entry_filters.cpp:908 (`count > max_non_ref_ac_any`
// increments N_failed) and line 912 (drop if N_failed == N_alleles-1
// when the `_any` thresholds are active). At threshold 2 a chr-20 site
// is dropped only when EVERY ALT has count > 2:
//
//   - 20:14370 (A=3) → all ALTs fail max → drop
//   - 20:17330 (A=1) → A passes max → keep
//   - 20:1110696 (G=2,T=4) → G passes max → keep
//   - 20:1230237 (mono, N_alleles=1, N_failed=0, fallback 0==0) → drop
//   - 20:1234567 (GA=3,GAC=1) → GAC passes max → keep
//   - 20:1235237 (mono) → fallback drops
func TestParity_MaxNonRefACAny_2_Chr20(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Chr:            "20",
		MaxNonRefACAny: 2,
		Recode:         true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "max_non_ref_ac_any_2_chr20.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_NonRefAF_03_Any_06 — `--non-ref-af 0.3 --non-ref-af-any 0.6
// --recode` byte-for-byte. Documents that the AF `_any` flag has an
// observable effect ONLY when the plain `--non-ref-af` flag is also set
// (upstream parameters.cpp:303,305 + entry_filters.cpp:814 fallback
// gated on `min_non_ref_af > 0`, NOT on `min_non_ref_af_any > 0`).
//
// With plain ≥ 0.3 alone the kept sites are 20:14370, 20:1110696, X:9,
// X:12 (see TestParity_NonRefAF_03). Adding `--non-ref-af-any 0.6`
// fires the post-loop fallback for biallelic sites whose single ALT
// has freq < 0.6 (N_failed=1, N_alleles-1=1):
//
//   - 20:14370 (A freq 0.5): plain passes; N_failed=1; fallback drops.
//   - 20:1110696 (G=0.333,T=0.667): plain passes; N_failed=1 of 2
//     (only G < 0.6); fallback NOT fired (1 != 2) → KEEP.
//   - X:9 (T freq 0.4): N_failed=1; fallback drops.
//   - X:12 (A freq exactly 0.6): `freq < 0.6` is false → N_failed=0;
//     not == N_alleles-1=1 → KEEP.
func TestParity_NonRefAF_03_Any_06(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MinNonRefAF:    0.3,
		MinNonRefAFAny: 0.6,
		Recode:         true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "non_ref_af_03_any_06.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_NonRefAFAny_NoOp — `--non-ref-af-any` and
// `--max-non-ref-af-any` are observably NO-OPS when used alone in
// upstream vcftools because the line-814 fallback gate is keyed on the
// PLAIN thresholds (`min_non_ref_af > 0.0` / `max_non_ref_af < 1.0`),
// not on the `_any` thresholds. Pinning this upstream quirk prevents
// future refactors from "fixing" the bug and breaking parity.
//
// Strategy: run the port with extreme `_any` thresholds and assert the
// result is bit-identical to a no-filter recode. We do NOT diff against
// upstream because of pre-existing baseline-recode gaps (20:1235237 and
// X:11 in the port's `MinAlleles=2` default). Comparing port-to-port is
// the right invariant here.
func TestParity_NonRefAFAny_NoOp(t *testing.T) {
	// Baseline: just --recode.
	basePrefix := runVcftoolsParity(t, "sample.vcf", &Params{Recode: true})
	base := readFileBytes(t, basePrefix+".recode.vcf")

	// With --non-ref-af-any at an extreme threshold — should match baseline.
	anyMinPrefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MinNonRefAFAny: 0.99,
		Recode:         true,
	})
	if got := readFileBytes(t, anyMinPrefix+".recode.vcf"); !bytes.Equal(got, base) {
		t.Errorf("--non-ref-af-any alone should be a no-op (upstream quirk); diff vs baseline")
	}

	// With --max-non-ref-af-any at a very low threshold — should match baseline.
	anyMaxPrefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MaxNonRefAFAny: 0.001,
		Recode:         true,
	})
	if got := readFileBytes(t, anyMaxPrefix+".recode.vcf"); !bytes.Equal(got, base) {
		t.Errorf("--max-non-ref-af-any alone should be a no-op (upstream quirk); diff vs baseline")
	}
}

// TestParity_MinQ — `--minQ 20` drops sites with QUAL < 20.
func TestParity_MinQ(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MinQ:   20.0,
		Recode: true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	wantKeys := map[string]bool{
		"20\t14370":   true,
		"20\t1110696": true,
		"20\t1230237": true,
		"20\t1234567": true,
	}
	seen := map[string]bool{}
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		key := strings.Join(strings.SplitN(ln, "\t", 3)[:2], "\t")
		seen[key] = true
	}
	for k := range wantKeys {
		if !seen[k] {
			t.Errorf("missing site %q after --minQ 20", k)
		}
	}
	for k := range seen {
		if !wantKeys[k] {
			t.Errorf("unexpected site %q after --minQ 20", k)
		}
	}
}

// TestParity_RemoveIndels — drops indel sites.
func TestParity_RemoveIndels(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		RemoveIndels: true,
		Recode:       true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		key := strings.Join(strings.SplitN(ln, "\t", 3)[:2], "\t")
		if key == "20\t1234567" || key == "X\t10" {
			t.Errorf("--remove-indels left indel %q in output", key)
		}
	}
}

// TestParity_KeepOnlyIndels — inverse. The accepted indels are
// 20:1234567 (G -> GA,GAC), X:10 (AC -> A,ATG), and X:11 (T ->
// A,<DEL:ME:ALU>) because the symbolic <DEL:ME:ALU> ALT also tests as
// "different length from REF" — matching upstream behaviour, since
// upstream also treats symbolic ALTs as indels for this filter (it
// rejects only when every ALT is a single nucleotide).
func TestParity_KeepOnlyIndels(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		KeepOnlyIndels: true,
		Recode:         true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	wantKeys := map[string]bool{
		"20\t1234567": true,
		"X\t10":       true,
		"X\t11":       true,
	}
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		key := strings.Join(strings.SplitN(ln, "\t", 3)[:2], "\t")
		if !wantKeys[key] {
			t.Errorf("--keep-only-indels left non-indel %q", key)
		}
	}
}

// TestParity_MaxMissing — drops X:11 (partially-missing).
func TestParity_MaxMissing(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MaxMissing: 1.0,
		Recode:     true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		if strings.HasPrefix(ln, "X\t11\t") {
			t.Errorf("--max-missing 1.0 kept partially-missing X:11: %q", ln)
		}
	}
}

// -----------------------------------------------------------------------------
// Sample management: --indv, --remove-indv, --keep
// -----------------------------------------------------------------------------

// TestParity_IndvKeep — `--indv NA00001` keeps only that sample column.
func TestParity_IndvKeep(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		IndvList: []string{"NA00001"},
		Recode:   true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	for _, ln := range got {
		if !strings.HasPrefix(ln, "#CHROM") {
			continue
		}
		fields := strings.Split(ln, "\t")
		if len(fields) != 10 {
			t.Errorf("expected 10 columns (1 sample), got %d: %q", len(fields), ln)
		}
		if len(fields) >= 10 && fields[9] != "NA00001" {
			t.Errorf("expected only NA00001, got %q", fields[9])
		}
	}
}

// TestParity_RemoveIndv — drops a sample column.
func TestParity_RemoveIndv(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		RemoveIndvList: []string{"NA00003"},
		Recode:         true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	for _, ln := range got {
		if !strings.HasPrefix(ln, "#CHROM") {
			continue
		}
		fields := strings.Split(ln, "\t")
		for _, f := range fields {
			if f == "NA00003" {
				t.Errorf("--remove-indv NA00003 left it in #CHROM line: %q", ln)
			}
		}
	}
}

// TestParity_Keep — `--keep FILE` keeps sample names listed in FILE.
func TestParity_Keep(t *testing.T) {
	tmp := t.TempDir()
	kf := filepath.Join(tmp, "keep.txt")
	if err := os.WriteFile(kf, []byte("NA00001\nNA00002\n"), 0644); err != nil {
		t.Fatalf("write keep: %v", err)
	}
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		KeepFile: kf,
		Recode:   true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	for _, ln := range got {
		if !strings.HasPrefix(ln, "#CHROM") {
			continue
		}
		f := strings.Split(ln, "\t")
		// 9 fixed + 2 samples.
		if len(f) != 11 {
			t.Errorf("want 11 #CHROM cols (2 samples), got %d: %q", len(f), ln)
		}
	}
}

// -----------------------------------------------------------------------------
// Site stats: --freq, --counts, --site-pi, --hardy, --missing-site,
// --missing-indv, --depth, --site-depth, --site-mean-depth, --het,
// --singletons, --TsTv-summary, --TsTv-by-count, --TsTv-by-qual
// -----------------------------------------------------------------------------

// TestParity_Freq — `--freq` byte-for-byte against an upstream-format
// golden file. Only biallelic SNPs appear because our port restricts
// --freq to biallelic loci (PARITY_ROADMAP.md#vcftools tracks the
// multi-allelic gap).
func TestParity_Freq(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Freq: true})
	got := readFileBytes(t, prefix+".frq")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "freq.expected.frq"))
	if !bytes.Equal(got, want) {
		t.Errorf(".frq mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Counts — `--counts` byte-for-byte.
func TestParity_Counts(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Counts: true})
	got := readFileBytes(t, prefix+".frq.count")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "counts.expected.frq.count"))
	if !bytes.Equal(got, want) {
		t.Errorf(".frq.count mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_SitePi — `--site-pi` byte-for-byte against an upstream
// golden generated by reference_code/vcftools. The two implementations
// compute mathematically identical values (mismatches /
// total*(total-1) is the same quantity as the textbook
// (n^2 - sum(c_a^2)) / (n*(n-1))); after this PR we also match
// upstream's C++-default float formatting and skip non-diploid sites,
// so byte parity holds.
func TestParity_SitePi(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{SitePi: true})
	got := readFileBytes(t, prefix+".sites.pi")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "site_pi.expected.sites.pi"))
	if !bytes.Equal(got, want) {
		t.Errorf(".sites.pi mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_SitePi_TextbookFormula — sanity-check three hand-traced
// values against the C++-default formatted output. For 20:14370
// (G:3/A:3, n=6), pi = (36 - 9 - 9) / (6*5) = 18/30 = 0.6, which
// upstream emits as "0.6" (six-significant-digit shortest representation).
// X:9 is non-diploid in this fixture and is filtered out by the
// is_diploid gate, matching upstream.
func TestParity_SitePi_TextbookFormula(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{SitePi: true})
	lines := readFileLines(t, prefix+".sites.pi")
	wanted := map[string]string{
		"20\t14370": "0.6",
		"19\t111":   "0.333333",
	}
	for _, ln := range lines[1:] {
		fields := strings.SplitN(ln, "\t", 3)
		if len(fields) != 3 {
			continue
		}
		key := fields[0] + "\t" + fields[1]
		if want, ok := wanted[key]; ok && fields[2] != want {
			t.Errorf("pi at %s: got %s, want %s", key, fields[2], want)
		}
		if fields[0] == "X" {
			t.Errorf("non-diploid X chrom site at %s should have been skipped by is_diploid filter", key)
		}
	}
}

// TestParity_Hardy — byte-for-byte against an upstream-format golden file.
// The directional P-values are placeholders (see PARITY_ROADMAP.md).
func TestParity_Hardy(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Hardy: true})
	got := readFileBytes(t, prefix+".hwe")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hardy.expected.hwe"))
	if !bytes.Equal(got, want) {
		t.Errorf(".hwe mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_MissingSite — byte-for-byte.
func TestParity_MissingSite(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{MissingSite: true})
	got := readFileBytes(t, prefix+".lmiss")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "missing_site.expected.lmiss"))
	if !bytes.Equal(got, want) {
		t.Errorf(".lmiss mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_MissingIndv — byte-for-byte.
func TestParity_MissingIndv(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{MissingIndv: true})
	got := readFileBytes(t, prefix+".imiss")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "missing_indv.expected.imiss"))
	if !bytes.Equal(got, want) {
		t.Errorf(".imiss mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Depth — byte-for-byte.
func TestParity_Depth(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Depth: true})
	got := readFileBytes(t, prefix+".idepth")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "depth.expected.idepth"))
	if !bytes.Equal(got, want) {
		t.Errorf(".idepth mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_SiteDepth — byte-for-byte. SUMSQ_DEPTH is a literal 0 — see
// docs/PARITY_ROADMAP.md#vcftools.
func TestParity_SiteDepth(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{SiteDepth: true})
	got := readFileBytes(t, prefix+".ldepth")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "site_depth.expected.ldepth"))
	if !bytes.Equal(got, want) {
		t.Errorf(".ldepth mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_SiteMeanDepth — byte-for-byte. VAR_DEPTH is a literal 0.
func TestParity_SiteMeanDepth(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{SiteMeanDepth: true})
	got := readFileBytes(t, prefix+".ldepth.mean")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "site_mean_depth.expected.ldepth.mean"))
	if !bytes.Equal(got, want) {
		t.Errorf(".ldepth.mean mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Het — byte-for-byte.
func TestParity_Het(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Het: true})
	got := readFileBytes(t, prefix+".het")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "het.expected.het"))
	if !bytes.Equal(got, want) {
		t.Errorf(".het mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_Singletons — byte-for-byte.
func TestParity_Singletons(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Singletons: true})
	got := readFileBytes(t, prefix+".singletons")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "singletons.expected.singletons"))
	if !bytes.Equal(got, want) {
		t.Errorf(".singletons mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_TsTvSummary — byte-for-byte.
func TestParity_TsTvSummary(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{TsTvSummary: true})
	got := readFileBytes(t, prefix+".TsTv.summary")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "tstv_summary.expected.TsTv.summary"))
	if !bytes.Equal(got, want) {
		t.Errorf(".TsTv.summary mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestVcftoolsLiveOracle runs --freq, --counts, --het, --missing-site,
// --site-pi and --TsTv-summary on testdata/live/in.vcf and compares each
// output file byte-for-byte against expected.* — golden files generated
// by a freshly-built genuine VCFtools 0.1.18 binary
// (reference_code/vcftools/src/cpp/vcftools). All six must match exactly.
func TestVcftoolsLiveOracle(t *testing.T) {
	liveDir, err := filepath.Abs(filepath.Join("testdata", "live"))
	if err != nil {
		t.Fatalf("Abs: %v", err)
	}
	in := filepath.Join(liveDir, "in.vcf")

	cases := []struct {
		name     string
		params   *Params
		out      string // output suffix
		expected string // expected.* file name
	}{
		{"freq", &Params{Freq: true}, ".frq", "expected.frq"},
		{"counts", &Params{Counts: true}, ".frq.count", "expected.frq.count"},
		{"het", &Params{Het: true}, ".het", "expected.het"},
		{"missing-site", &Params{MissingSite: true}, ".lmiss", "expected.lmiss"},
		{"site-pi", &Params{SitePi: true}, ".sites.pi", "expected.sites.pi"},
		{"TsTv-summary", &Params{TsTvSummary: true}, ".TsTv.summary", "expected.TsTv.summary"},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			tmp := t.TempDir()
			prefix := filepath.Join(tmp, "out")
			tc.params.OutPrefix = prefix
			f, err := os.Open(in)
			if err != nil {
				t.Fatalf("open in.vcf: %v", err)
			}
			defer f.Close()
			if err := Run(f, tc.params); err != nil {
				t.Fatalf("Run: %v", err)
			}
			got := readFileBytes(t, prefix+tc.out)
			want := readFileBytes(t, filepath.Join(liveDir, tc.expected))
			if !bytes.Equal(got, want) {
				t.Errorf("%s mismatch\nwant:\n%s\ngot:\n%s", tc.out, want, got)
			}
		})
	}
}

// genuineVcftoolsBinary returns the absolute path of the upstream
// vcftools binary built from the vendored submodule. Tests that need to
// drive it should call this and t.Skip when it is missing — the parity
// suite is designed to run without it, but the live-oracle tests below
// fail closed when both the binary and the port are available.
func genuineVcftoolsBinary(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join(
		"..", "..", "..", "..",
		"reference_code", "vcftools", "src", "cpp", "vcftools",
	))
	if err != nil {
		t.Fatalf("Abs: %v", err)
	}
	if _, err := os.Stat(abs); err != nil {
		t.Skipf("genuine vcftools binary not built at %s: %v", abs, err)
	}
	return abs
}

// TestVcftoolsLiveOracleFull is the broad expansion of TestVcftoolsLiveOracle.
// For every vcftools option wired into our port that the genuine C binary
// also supports on testdata/live/in.vcf, this test:
//
//  1. invokes the genuine binary
//     (reference_code/vcftools/src/cpp/vcftools) at test time;
//  2. invokes our port on the same input with the same Params;
//  3. diffs each non-".log" output file byte-for-byte.
//
// Goldens are NOT vendored — they are regenerated every run from the live
// binary. This is the only way to keep us honest as the port evolves; the
// pre-existing vendored expected.* files (which TestVcftoolsLiveOracle
// still relies on) drifted out of date and hid bugs in --freq/--het/
// --missing-site that landed in commit dea8695.
//
// Subtests t.Skip when the genuine binary refuses the option (e.g.
// --BEAGLE-GL requires GL/PL FORMAT tags absent from in.vcf, --LROH
// requires --chr explicitly, etc.). They t.Errorf — not t.Skip — on
// mismatch: the goal is to surface gaps.
func TestVcftoolsLiveOracleFull(t *testing.T) {
	bin := genuineVcftoolsBinary(t)

	liveDir, err := filepath.Abs(filepath.Join("testdata", "live"))
	if err != nil {
		t.Fatalf("Abs: %v", err)
	}
	in := filepath.Join(liveDir, "in.vcf")

	// Each case lists (subtest name, the upstream CLI args, the matching
	// Params, and the output-file suffixes to diff). The args slice does
	// NOT include --vcf / --out; runOracle adds those. The suffixes do
	// NOT include the ".log" suffix; that file embeds timestamps and
	// upstream version banners and is never compared.
	type oracle struct {
		name     string
		args     []string
		params   *Params
		suffixes []string
	}

	cases := []oracle{
		// ------ Site / allele frequency family ------
		{"freq", []string{"--freq"}, &Params{Freq: true}, []string{".frq"}},
		{"counts", []string{"--counts"}, &Params{Counts: true}, []string{".frq.count"}},
		{"freq2", []string{"--freq2"}, &Params{Freq2: true}, []string{".frq"}},
		{"counts2", []string{"--counts2"}, &Params{Counts2: true}, []string{".frq.count"}},

		// ------ Site stats ------
		{"site-pi", []string{"--site-pi"}, &Params{SitePi: true}, []string{".sites.pi"}},
		{"site-quality", []string{"--site-quality"}, &Params{SiteQuality: true}, []string{".lqual"}},
		{"site-depth", []string{"--site-depth"}, &Params{SiteDepth: true}, []string{".ldepth"}},
		{"site-mean-depth", []string{"--site-mean-depth"}, &Params{SiteMeanDepth: true}, []string{".ldepth.mean"}},
		{"missing-site", []string{"--missing-site"}, &Params{MissingSite: true}, []string{".lmiss"}},
		{"hardy", []string{"--hardy"}, &Params{Hardy: true}, []string{".hwe"}},
		{"singletons", []string{"--singletons"}, &Params{Singletons: true}, []string{".singletons"}},
		{"hist-indel-len", []string{"--hist-indel-len"}, &Params{HistIndelLen: true}, []string{".indel.hist"}},
		{"FILTER-summary", []string{"--FILTER-summary"}, &Params{FilterSummary: true}, []string{".FILTER.summary"}},
		{"SNPdensity-100", []string{"--SNPdensity", "100"}, &Params{SNPDensity: 100}, []string{".snpden"}},
		{"kept-sites", []string{"--kept-sites"}, &Params{KeptSites: true}, []string{".kept.sites"}},
		{"removed-sites-with-filter", []string{"--minQ", "100", "--removed-sites"}, &Params{MinQ: 100, RemovedSites: true}, []string{".removed.sites"}},

		// ------ Per-individual stats ------
		{"depth", []string{"--depth"}, &Params{Depth: true}, []string{".idepth"}},
		{"missing-indv", []string{"--missing-indv"}, &Params{MissingIndv: true}, []string{".imiss"}},
		{"het", []string{"--het"}, &Params{Het: true}, []string{".het"}},
		{"geno-depth", []string{"--geno-depth"}, &Params{GenoDepth: true}, []string{".gdepth"}},
		{"indv-burden", []string{"--indv-burden"}, &Params{IndvBurden: true}, []string{".iburden"}},
		{"indv-freq-burden", []string{"--indv-freq-burden"}, &Params{IndvFreqBurden: true}, []string{".ifreqburden"}},
		{"indv-freq-burden2", []string{"--indv-freq-burden2"}, &Params{IndvFreqBurden2: true}, []string{".ifreqburden"}},

		// ------ TsTv family ------
		{"TsTv-summary", []string{"--TsTv-summary"}, &Params{TsTvSummary: true}, []string{".TsTv.summary"}},
		{"TsTv-by-count", []string{"--TsTv-by-count"}, &Params{TsTvByCount: true}, []string{".TsTv.count"}},
		{"TsTv-by-qual", []string{"--TsTv-by-qual"}, &Params{TsTvByQual: true}, []string{".TsTv.qual"}},
		{"TsTv-binned-100", []string{"--TsTv", "100"}, &Params{TsTvBinSize: 100}, []string{".TsTv"}},

		// ------ Windowed stats ------
		{"TajimaD-100", []string{"--TajimaD", "100"}, &Params{TajimaD: 100}, []string{".Tajima.D"}},
		{"window-pi-1000", []string{"--window-pi", "1000"}, &Params{WindowPi: 1000}, []string{".windowed.pi"}},

		// ------ LD ------
		{"geno-r2", []string{"--geno-r2"}, &Params{GenoR2: true}, []string{".geno.ld"}},
		{"hap-r2", []string{"--hap-r2"}, &Params{HapR2: true}, []string{".hap.ld"}},
		{"interchrom-geno-r2", []string{"--interchrom-geno-r2"}, &Params{InterchromGenoR2: true}, []string{".interchrom.geno.ld"}},
		{"interchrom-hap-r2", []string{"--interchrom-hap-r2"}, &Params{InterchromHapR2: true}, []string{".interchrom.hap.ld"}},
		{"geno-chisq", []string{"--geno-chisq"}, &Params{GenoChiSq: true}, []string{".geno.chisq"}},

		// ------ Relatedness ------
		{"relatedness", []string{"--relatedness"}, &Params{Relatedness: true}, []string{".relatedness"}},
		{"relatedness2", []string{"--relatedness2"}, &Params{Relatedness2: true}, []string{".relatedness2"}},

		// ------ LROH (needs --chr) ------
		{"LROH-chr1", []string{"--chr", "1", "--LROH"}, &Params{Chr: "1", LROH: true}, []string{".LROH"}},

		// ------ INFO extraction ------
		{"get-INFO-DP", []string{"--get-INFO", "DP"}, &Params{GetINFO: "DP"}, []string{".INFO"}},

		// ------ Format conversions ------
		{"012", []string{"--012"}, &Params{Output012: true}, []string{".012", ".012.indv", ".012.pos"}},
		{"plink", []string{"--plink"}, &Params{OutputPlink: true}, []string{".ped", ".map"}},
		{"plink-tped", []string{"--plink-tped"}, &Params{OutputPlinkTped: true}, []string{".tped", ".tfam"}},

		// ------ Recode passthroughs ------
		{"recode", []string{"--recode"}, &Params{Recode: true}, []string{".recode.vcf"}},
		{"recode-INFO-all", []string{"--recode", "--recode-INFO-all"}, &Params{Recode: true, RecodeInfoAll: true}, []string{".recode.vcf"}},

		// ------ Site-set filters (paired with --recode for an observable
		// output) ------
		{"chr-filter", []string{"--chr", "1", "--recode"}, &Params{Chr: "1", Recode: true}, []string{".recode.vcf"}},
		{"from-to-bp", []string{"--chr", "1", "--from-bp", "150", "--to-bp", "350", "--recode"}, &Params{Chr: "1", FromBp: 150, ToBp: 350, Recode: true}, []string{".recode.vcf"}},
		{"maf", []string{"--maf", "0.3", "--recode"}, &Params{Maf: 0.3, Recode: true}, []string{".recode.vcf"}},
		{"max-maf", []string{"--max-maf", "0.4", "--recode"}, &Params{MaxMaf: 0.4, Recode: true}, []string{".recode.vcf"}},
		{"mac", []string{"--mac", "2", "--recode"}, &Params{Mac: 2, Recode: true}, []string{".recode.vcf"}},
		{"max-mac", []string{"--max-mac", "4", "--recode"}, &Params{MaxMac: 4, Recode: true}, []string{".recode.vcf"}},
		{"minQ", []string{"--minQ", "20", "--recode"}, &Params{MinQ: 20, Recode: true}, []string{".recode.vcf"}},
		{"remove-indels", []string{"--remove-indels", "--recode"}, &Params{RemoveIndels: true, Recode: true}, []string{".recode.vcf"}},
		{"keep-only-indels", []string{"--keep-only-indels", "--recode"}, &Params{KeepOnlyIndels: true, Recode: true}, []string{".recode.vcf"}},
		{"max-missing-1", []string{"--max-missing", "1.0", "--recode"}, &Params{MaxMissing: 1.0, Recode: true}, []string{".recode.vcf"}},
		{"max-missing-count-0", []string{"--max-missing-count", "0", "--recode"}, &Params{MaxMissingCount: 0, MaxMissingCountSet: true, Recode: true}, []string{".recode.vcf"}},
		{"remove-filtered-all", []string{"--remove-filtered-all", "--recode"}, &Params{RemoveFilteredAll: true, Recode: true}, []string{".recode.vcf"}},
		{"non-ref-af", []string{"--non-ref-af", "0.3", "--recode"}, &Params{MinNonRefAF: 0.3, Recode: true}, []string{".recode.vcf"}},
		{"non-ref-ac", []string{"--non-ref-ac", "1", "--recode"}, &Params{MinNonRefAC: 1, Recode: true}, []string{".recode.vcf"}},
		{"max-non-ref-af", []string{"--max-non-ref-af", "0.9", "--recode"}, &Params{MaxNonRefAF: 0.9, Recode: true}, []string{".recode.vcf"}},
		{"max-non-ref-ac", []string{"--max-non-ref-ac", "5", "--recode"}, &Params{MaxNonRefAC: 5, Recode: true}, []string{".recode.vcf"}},
		{"minDP", []string{"--minDP", "6", "--recode"}, &Params{MinDP: 6, Recode: true}, []string{".recode.vcf"}},
		{"maxDP", []string{"--maxDP", "50", "--recode"}, &Params{MaxDP: 50, Recode: true}, []string{".recode.vcf"}},
		{"min-meanDP", []string{"--min-meanDP", "5", "--recode"}, &Params{MinMeanDP: 5, Recode: true}, []string{".recode.vcf"}},
		{"max-meanDP", []string{"--max-meanDP", "50", "--recode"}, &Params{MaxMeanDP: 50, Recode: true}, []string{".recode.vcf"}},
		{"minGQ", []string{"--minGQ", "40", "--recode"}, &Params{MinGQ: 40, Recode: true}, []string{".recode.vcf"}},
		{"thin-50", []string{"--thin", "50", "--recode"}, &Params{Thin: 50, Recode: true}, []string{".recode.vcf"}},
		{"max-alleles-2", []string{"--max-alleles", "2", "--recode"}, &Params{MaxAlleles: 2, Recode: true}, []string{".recode.vcf"}},
		{"phased", []string{"--phased", "--recode"}, &Params{Phased: true, Recode: true}, []string{".recode.vcf"}},
		{"hwe", []string{"--hwe", "0.001", "--recode"}, &Params{MinHWEPvalue: 0.001, MaxAlleles: 2, Recode: true}, []string{".recode.vcf"}},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			runOracleCase(t, bin, in, tc.args, tc.params, tc.suffixes)
		})
	}
}

// runOracleCase drives one live-oracle subtest. It runs the genuine
// vcftools binary in a fresh temp dir; if that succeeds, it runs our
// port on the same input with the supplied Params; then it diffs each
// requested output suffix byte-for-byte. If the genuine binary fails
// (non-zero exit) the subtest is skipped — the option is not exercisable
// on this fixture and there is nothing to compare. Port failures are
// surfaced via t.Fatalf so they cannot hide behind a missing golden.
func runOracleCase(t *testing.T, bin, in string, args []string, params *Params, suffixes []string) {
	t.Helper()
	tmp := t.TempDir()

	// Drive the genuine binary in `goldDir` so its log/aux files cannot
	// collide with the port's output.
	goldDir := filepath.Join(tmp, "gold")
	if err := os.Mkdir(goldDir, 0o755); err != nil {
		t.Fatalf("mkdir gold: %v", err)
	}
	goldPrefix := filepath.Join(goldDir, "out")
	cmdArgs := append([]string{"--vcf", in, "--out", goldPrefix}, args...)
	cmd := exec.Command(bin, cmdArgs...)
	cmd.Dir = goldDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Skipf("genuine vcftools failed for %v: %v\n%s", args, err, out)
	}

	// Drive the port.
	portPrefix := filepath.Join(tmp, "port_out")
	params.OutPrefix = portPrefix
	f, err := os.Open(in)
	if err != nil {
		t.Fatalf("open in.vcf: %v", err)
	}
	defer f.Close()
	if err := Run(f, params); err != nil {
		t.Fatalf("port Run failed: %v", err)
	}

	// Diff each requested suffix byte-for-byte.
	for _, suf := range suffixes {
		gold, err := os.ReadFile(goldPrefix + suf)
		if err != nil {
			t.Skipf("genuine binary did not produce %s%s: %v", "out", suf, err)
			return
		}
		got, err := os.ReadFile(portPrefix + suf)
		if err != nil {
			t.Errorf("port did not produce %s: %v", suf, err)
			continue
		}
		if !bytes.Equal(got, gold) {
			t.Errorf("%s mismatch\n--- want (genuine) ---\n%s\n--- got (port) ---\n%s",
				suf, gold, got)
		}
	}
}

// TestVcftoolsLiveOracleLROH drives the --LROH detector against the genuine
// vcftools binary on a runs-of-homozygosity-rich fixture
// (testdata/live/lroh_in.vcf) and asserts the .LROH output is byte-identical.
//
// Unlike the LROH-chr1 subtest of TestVcftoolsLiveOracleFull — whose fixture
// (in.vcf) carries no autozygous runs, so BOTH binaries emit only the header —
// this fixture contains long homozygous stretches that force the forward-
// backward HMM to call non-empty runs. The test therefore proves the ported
// algorithm (states / emission + transition probabilities / decode / run
// extraction in relatedness.go) actually matches upstream, rather than passing
// on an empty-rows tautology. It explicitly fails if the genuine binary's
// output has no data rows, so the fixture can never silently regress to empty.
func TestVcftoolsLiveOracleLROH(t *testing.T) {
	bin := genuineVcftoolsBinary(t)

	liveDir, err := filepath.Abs(filepath.Join("testdata", "live"))
	if err != nil {
		t.Fatalf("Abs: %v", err)
	}
	in := filepath.Join(liveDir, "lroh_in.vcf")

	tmp := t.TempDir()
	goldDir := filepath.Join(tmp, "gold")
	if err := os.Mkdir(goldDir, 0o755); err != nil {
		t.Fatalf("mkdir gold: %v", err)
	}
	goldPrefix := filepath.Join(goldDir, "out")
	cmd := exec.Command(bin, "--vcf", in, "--chr", "1", "--LROH", "--out", goldPrefix)
	cmd.Dir = goldDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Skipf("genuine vcftools --LROH failed: %v\n%s", err, out)
	}

	gold, err := os.ReadFile(goldPrefix + ".LROH")
	if err != nil {
		t.Fatalf("read genuine .LROH: %v", err)
	}
	// Guard against fixture rot: the genuine binary must emit data rows, not
	// just the header, or this test degrades into the very tautology it exists
	// to replace.
	if dataRows := bytes.Count(gold, []byte("\n")) - 1; dataRows < 1 {
		t.Fatalf("fixture produced no LROH data rows in genuine output:\n%s", gold)
	}

	portPrefix := filepath.Join(tmp, "port_out")
	f, err := os.Open(in)
	if err != nil {
		t.Fatalf("open lroh_in.vcf: %v", err)
	}
	defer f.Close()
	if err := Run(f, &Params{Chr: "1", LROH: true, OutPrefix: portPrefix}); err != nil {
		t.Fatalf("port Run failed: %v", err)
	}
	got, err := os.ReadFile(portPrefix + ".LROH")
	if err != nil {
		t.Fatalf("read port .LROH: %v", err)
	}
	if !bytes.Equal(got, gold) {
		t.Errorf(".LROH mismatch\n--- want (genuine) ---\n%s\n--- got (port) ---\n%s", gold, got)
	}
}

// TestParity_TsTvByCount_Header — header byte-for-byte.
func TestParity_TsTvByCount_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{TsTvByCount: true})
	lines := readFileLines(t, prefix+".TsTv.count")
	if len(lines) == 0 {
		t.Fatalf("empty TsTv.count")
	}
	want := "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv"
	if lines[0] != want {
		t.Errorf("header mismatch.\nwant: %q\ngot:  %q", want, lines[0])
	}
}

// TestParity_TsTvByCount_FullRows — `--TsTv-by-count` byte-for-byte
// against an upstream golden. Exercises the dense 0..2*N_kept_indv-1
// enumeration with empty bins rendered as "0\t0\t-nan" to match
// upstream's `double(0)/0` glibc literal. Mirrors
// variant_file_output.cpp:3220-3225.
func TestParity_TsTvByCount_FullRows(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{TsTvByCount: true})
	got := readFileBytes(t, prefix+".TsTv.count")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "tstv_by_count.expected.TsTv.count"))
	if !bytes.Equal(got, want) {
		t.Errorf(".TsTv.count mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_TsTvByQual_Header — `--TsTv-by-qual` header byte-for-byte.
func TestParity_TsTvByQual_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{TsTvByQual: true})
	lines := readFileLines(t, prefix+".TsTv.qual")
	if len(lines) == 0 {
		t.Fatalf("empty TsTv.qual")
	}
	want := "QUAL_THRESHOLD\tN_Ts_LT_QUAL_THRESHOLD\tN_Tv_LT_QUAL_THRESHOLD\tTs/Tv_LT_QUAL_THRESHOLD\tN_Ts_GT_QUAL_THRESHOLD\tN_Tv_GT_QUAL_THRESHOLD\tTs/Tv_GT_QUAL_THRESHOLD"
	if lines[0] != want {
		t.Errorf("header mismatch.\nwant: %q\ngot:  %q", want, lines[0])
	}
}

// TestParity_TsTv_Binned — `--TsTv N` (binned) byte-for-byte against
// an upstream golden generated with --TsTv 1000000. Exercises the
// per-chromosome dense bin layout and CHROM/BinStart/SNP_count/Ts/Tv
// column order from variant_file_output.cpp:3057-3068.
func TestParity_TsTv_Binned(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{TsTvBinSize: 1000000})
	got := readFileBytes(t, prefix+".TsTv")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "tstv_binned_1m.expected.TsTv"))
	if !bytes.Equal(got, want) {
		t.Errorf(".TsTv mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// -----------------------------------------------------------------------------
// Population genetics: --weir-fst-pop, --fst-window-size, --fst-window-step
// -----------------------------------------------------------------------------

// TestParity_WeirFstPop_HeaderOnly — exercises --weir-fst-pop family and
// asserts the per-site Fst file has the expected upstream header. Byte
// parity on values is skipped because three samples / one variant is too
// sparse to produce a stable expected value without invoking upstream.
func TestParity_WeirFstPop_HeaderOnly(t *testing.T) {
	tmp := t.TempDir()
	pop1 := filepath.Join(tmp, "pop1.txt")
	pop2 := filepath.Join(tmp, "pop2.txt")
	if err := os.WriteFile(pop1, []byte("NA00001\n"), 0644); err != nil {
		t.Fatalf("write pop1: %v", err)
	}
	if err := os.WriteFile(pop2, []byte("NA00002\nNA00003\n"), 0644); err != nil {
		t.Fatalf("write pop2: %v", err)
	}
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		WeirFstPop: []string{pop1, pop2},
	})
	lines := readFileLines(t, prefix+".weir.fst")
	if len(lines) == 0 {
		t.Fatalf("empty .weir.fst output")
	}
	if lines[0] != "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST" {
		t.Errorf("header mismatch: %q", lines[0])
	}
}

// TestParity_WeirFst_Windowed_HeaderOnly — header check for
// --fst-window-size.
func TestParity_WeirFst_Windowed_HeaderOnly(t *testing.T) {
	tmp := t.TempDir()
	pop1 := filepath.Join(tmp, "pop1.txt")
	pop2 := filepath.Join(tmp, "pop2.txt")
	if err := os.WriteFile(pop1, []byte("NA00001\n"), 0644); err != nil {
		t.Fatalf("write pop1: %v", err)
	}
	if err := os.WriteFile(pop2, []byte("NA00002\nNA00003\n"), 0644); err != nil {
		t.Fatalf("write pop2: %v", err)
	}
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		WeirFstPop:    []string{pop1, pop2},
		FstWindowSize: 1000,
	})
	lines := readFileLines(t, prefix+".windowed.weir.fst")
	if len(lines) == 0 {
		t.Fatalf("empty windowed Fst output")
	}
	want := "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST"
	if lines[0] != want {
		t.Errorf("header mismatch:\nwant: %q\ngot:  %q", want, lines[0])
	}
}

// -----------------------------------------------------------------------------
// LD: --geno-r2, --hap-r2
// -----------------------------------------------------------------------------

// TestParity_GenoR2_Header — header byte-for-byte.
func TestParity_GenoR2_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{GenoR2: true})
	lines := readFileLines(t, prefix+".geno.ld")
	if len(lines) == 0 {
		t.Fatalf("empty .geno.ld")
	}
	want := "CHR\tPOS1\tPOS2\tN_INDV\tR^2"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// TestParity_HapR2_Header — header byte-for-byte.
func TestParity_HapR2_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{HapR2: true})
	lines := readFileLines(t, prefix+".hap.ld")
	if len(lines) == 0 {
		t.Fatalf("empty .hap.ld")
	}
	want := "CHR\tPOS1\tPOS2\tN_CHR\tR^2\tD\tDprime"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// -----------------------------------------------------------------------------
// VCF recoding: --recode, --recode-INFO-all
// -----------------------------------------------------------------------------

// TestParity_Recode_AllSites — `--recode` with no filter emits 12 rows.
func TestParity_Recode_AllSites(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Recode: true})
	got := readFileLines(t, prefix+".recode.vcf")
	var data int
	for _, ln := range got {
		if !strings.HasPrefix(ln, "#") && ln != "" {
			data++
		}
	}
	if data != 12 {
		t.Errorf("--recode no-filter: want 12 data rows, got %d", data)
	}
}

// TestParity_Recode_InfoAll — `--recode-INFO-all` preserves the original
// INFO column. We sample one site and check that every INFO key from the
// source line is present in the output. We can't pin field order because
// our VCF writer iterates a map — see PARITY_ROADMAP.md#vcftools.
func TestParity_Recode_InfoAll(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		Recode:        true,
		RecodeInfoAll: true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	wantParts := []string{"NS=3", "DP=14", "AF=0.5", "DB", "H2"}
	for _, ln := range got {
		if !strings.HasPrefix(ln, "20\t14370\t") {
			continue
		}
		for _, p := range wantParts {
			if !strings.Contains(ln, p) {
				t.Errorf("20:14370 INFO should contain %q, got:\n%q", p, ln)
			}
		}
		return
	}
	t.Errorf("20:14370 not found in --recode output")
}

// -----------------------------------------------------------------------------
// Format conversions: --012, --plink, --plink-tped
// -----------------------------------------------------------------------------

// TestParity_012_Indv — `--012` emits one row per sample in .012.indv.
func TestParity_012_Indv(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Output012: true})
	got := readFileLines(t, prefix+".012.indv")
	want := []string{"NA00001", "NA00002", "NA00003"}
	assertLinesEqual(t, got, want)
}

// TestParity_012_RowPrefix — data row first column is the 0-based sample
// index, matching upstream.
func TestParity_012_RowPrefix(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Output012: true})
	got := readFileLines(t, prefix+".012")
	wantPrefix := []string{"0", "1", "2"}
	if len(got) != 3 {
		t.Fatalf("want 3 .012 rows, got %d", len(got))
	}
	for i, ln := range got {
		fields := strings.SplitN(ln, "\t", 2)
		if fields[0] != wantPrefix[i] {
			t.Errorf("row %d prefix: got %q, want %q", i, fields[0], wantPrefix[i])
		}
	}
}

// TestParity_012_Biallelic — `.012` only emits biallelic sites. The 12
// sites in sample.vcf reduce to 8 biallelic sites after dropping the four
// multi-allelic ones (20:1110696, 20:1234567, X:10, X:11).
func TestParity_012_Biallelic(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Output012: true})
	got := readFileLines(t, prefix+".012.pos")
	if len(got) != 8 {
		t.Errorf("want 8 biallelic rows in .012.pos, got %d:\n%s",
			len(got), strings.Join(got, "\n"))
	}
}

// TestParity_Plink_FilesExist — `--plink` emits a .ped and a .map.
func TestParity_Plink_FilesExist(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{OutputPlink: true})
	if _, err := os.Stat(prefix + ".ped"); err != nil {
		t.Errorf("missing .ped: %v", err)
	}
	if _, err := os.Stat(prefix + ".map"); err != nil {
		t.Errorf("missing .map: %v", err)
	}
}

// TestParity_PlinkTped_FilesExist — `--plink-tped` emits a .tped and
// a .tfam.
func TestParity_PlinkTped_FilesExist(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{OutputPlinkTped: true})
	if _, err := os.Stat(prefix + ".tped"); err != nil {
		t.Errorf("missing .tped: %v", err)
	}
	if _, err := os.Stat(prefix + ".tfam"); err != nil {
		t.Errorf("missing .tfam: %v", err)
	}
}

// -----------------------------------------------------------------------------
// Tajima's D, windowed pi, SNP density headers.
// -----------------------------------------------------------------------------

// TestParity_TajimaD_Header — `--TajimaD` header byte-for-byte.
func TestParity_TajimaD_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{TajimaD: 1000000})
	lines := readFileLines(t, prefix+".Tajima.D")
	if len(lines) == 0 {
		t.Fatalf("empty Tajima.D output")
	}
	want := "CHROM\tBIN_START\tN_SNPS\tTajimaD"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// TestParity_WindowPi_Header — `--window-pi` header byte-for-byte.
func TestParity_WindowPi_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{WindowPi: 1000000})
	lines := readFileLines(t, prefix+".windowed.pi")
	if len(lines) == 0 {
		t.Fatalf("empty windowed.pi output")
	}
	want := "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// -----------------------------------------------------------------------------
// Wave 1 long-tail closures: interchrom LD, chi-square, relatedness,
// relatedness2, LROH, phased-blocks, get-INFO, remove-filtered.
// -----------------------------------------------------------------------------

// TestParity_InterchromGenoR2_Header — header byte-for-byte against upstream
// `<prefix>.interchrom.geno.ld` layout: CHR1 POS1 CHR2 POS2 N_INDV R^2.
func TestParity_InterchromGenoR2_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{InterchromGenoR2: true})
	b, err := os.ReadFile(prefix + ".interchrom.geno.ld")
	if err != nil {
		t.Fatalf("read .interchrom.geno.ld: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	want := "CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// TestParity_InterchromHapR2_Header — header for `--interchrom-hap-r2`.
func TestParity_InterchromHapR2_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{InterchromHapR2: true})
	b, err := os.ReadFile(prefix + ".interchrom.hap.ld")
	if err != nil {
		t.Fatalf("read .interchrom.hap.ld: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	want := "CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2\tD\tDprime"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// TestParity_GenoChiSq_Header — header for `--geno-chisq`.
func TestParity_GenoChiSq_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{GenoChiSq: true})
	b, err := os.ReadFile(prefix + ".geno.chisq")
	if err != nil {
		t.Fatalf("read .geno.chisq: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(b), "\n"), "\n")
	want := "CHR\tPOS1\tPOS2\tN_INDV\tCHI^2\tDOF\tPVAL"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// TestParity_Relatedness_Header — header for `--relatedness`.
func TestParity_Relatedness_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Relatedness: true})
	lines := readFileLines(t, prefix+".relatedness")
	if len(lines) == 0 {
		t.Fatalf("empty .relatedness output")
	}
	want := "INDV1\tINDV2\tRELATEDNESS_AJK"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
	// 3 samples -> 3 self + 3 cross = 6 rows.
	if got := len(lines) - 1; got != 6 {
		t.Errorf("got %d rows, want 6", got)
	}
}

// TestParity_Relatedness2_Header — header for `--relatedness2`.
func TestParity_Relatedness2_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Relatedness2: true})
	lines := readFileLines(t, prefix+".relatedness2")
	if len(lines) == 0 {
		t.Fatalf("empty .relatedness2 output")
	}
	want := "INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// TestParity_LROH_Header — header for `--LROH`. Upstream requires a single
// chromosome, so the fixture is restricted to contig 20.
func TestParity_LROH_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{LROH: true, Chr: "20"})
	lines := readFileLines(t, prefix+".LROH")
	if len(lines) == 0 {
		t.Fatalf("empty .LROH output")
	}
	want := "CHROM\tAUTO_START\tAUTO_END\tMIN_START\tMAX_END\tN_VARIANTS_BETWEEN_MAX_BOUNDARIES\tN_MISMATCHES\tINDV"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// TestParity_PhasedBlocks_Header — header for `--phased-blocks`.
func TestParity_PhasedBlocks_Header(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{PhasedBlocks: true})
	lines := readFileLines(t, prefix+".blocks")
	if len(lines) == 0 {
		t.Fatalf("empty .blocks output")
	}
	want := "CHROM\tBLOCK_START\tBLOCK_END\tN_VARIANTS\tINDV"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
}

// TestParity_GetINFO — extract DP and AF tags from the fixture.
func TestParity_GetINFO(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{GetINFO: "DP,AF"})
	lines := readFileLines(t, prefix+".INFO")
	want := "CHROM\tPOS\tREF\tALT\tDP\tAF"
	if lines[0] != want {
		t.Errorf("header mismatch: %q vs %q", lines[0], want)
	}
	// Spot check: 20:14370 has DP=14, AF=0.5.
	for _, ln := range lines[1:] {
		fields := strings.Split(ln, "\t")
		if len(fields) >= 6 && fields[0] == "20" && fields[1] == "14370" {
			if fields[4] != "14" || fields[5] != "0.5" {
				t.Errorf("20:14370 got DP=%s AF=%s; want 14 / 0.5", fields[4], fields[5])
			}
			return
		}
	}
	t.Errorf("20:14370 not found in:\n%s", strings.Join(lines, "\n"))
}

// TestParity_RemoveFiltered — --remove-filtered q10 drops sites listing q10.
// In sample.vcf the q10-filtered sites are 20:17330 and X:11; both should
// drop. The remaining 10 records (PASS, ".", and "q10;s50" is also dropped
// via overlap) should survive less the q10 ones.
func TestParity_RemoveFiltered(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		RemoveFiltered: "q10",
		Recode:         true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		fields := strings.Split(ln, "\t")
		if len(fields) < 7 {
			continue
		}
		// FILTER is column 7 (1-indexed) = index 6.
		filter := fields[6]
		filters := strings.Split(filter, ";")
		for _, f := range filters {
			if f == "q10" {
				t.Errorf("q10-filtered row leaked: %q", ln)
			}
		}
	}
}

// TestParity_KeepFiltered — --keep-filtered q10 keeps only sites listing q10.
func TestParity_KeepFiltered(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		KeepFiltered: "q10",
		Recode:       true,
	})
	got := readFileLines(t, prefix+".recode.vcf")
	rowsSeen := 0
	for _, ln := range got {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		rowsSeen++
		fields := strings.Split(ln, "\t")
		filter := fields[6]
		filters := strings.Split(filter, ";")
		ok := false
		for _, f := range filters {
			if f == "q10" {
				ok = true
				break
			}
		}
		if !ok {
			t.Errorf("non-q10 row leaked with --keep-filtered q10: %q", ln)
		}
	}
	// Expect 2 rows: 20:17330 (q10) and X:11 (q10;s50).
	if rowsSeen != 2 {
		t.Errorf("expected 2 q10 rows, got %d", rowsSeen)
	}
}

// -----------------------------------------------------------------------------
// --hwe FLOAT (parameters.cpp:254): biallelic exact-test HWE p-value filter
// (Wigginton 2005). Upstream couples --hwe with `max_alleles = 2`, so even on
// fixtures with multi-allelic sites the filter looks biallelic-only.
// -----------------------------------------------------------------------------

// TestParity_HWE_005_sample — `--hwe 0.05 --recode` on sample.vcf
// (3 individuals). Sample.vcf has too few samples for any biallelic site to
// fail the exact HWE test (every site is in the "in-equilibrium" regime),
// so the only filtering here is the implicit `max_alleles=2` that upstream
// applies — the four multi-allelic / non-biallelic sites (1110696, 1234567,
// X:10, X:11) get dropped and the remaining 8 biallelic sites survive.
//
// We still pin this case because it exercises the
// `--hwe → max_alleles=2` coupling in the CLI->Params adapter, which is
// the path that is most likely to regress.
func TestParity_HWE_005_sample(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MinHWEPvalue: 0.05,
		MaxAlleles:   2, // upstream's parameters.cpp:254 coupling
		Recode:       true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hwe_005.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_HWE_005_fixture — `--hwe 0.05 --recode` on
// hwe_fixture.vcf (20 individuals, 4 sites). This fixture is designed
// to make the exact HWE test actually fire: site 1:200 has counts
// (HOM1=10, HET=0, HOM2=10), an extreme deficit of heterozygotes with
// exact p ~= 1.34e-6 (verified via `--hardy` against upstream). The
// other two biallelic sites are in equilibrium (p == 1.0). Site 1:300
// is multi-allelic and gets dropped by the upstream max_alleles=2
// coupling.
//
// Expected: sites 1:100 ("in_hwe") and 1:400 ("good") survive; 1:200
// ("out_hwe") drops on the HWE test; 1:300 ("multi") drops on
// max_alleles=2.
func TestParity_HWE_005_fixture(t *testing.T) {
	prefix := runVcftoolsParity(t, "hwe_fixture.vcf", &Params{
		MinHWEPvalue: 0.05,
		MaxAlleles:   2, // upstream's parameters.cpp:254 coupling
		Recode:       true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hwe_fixture_005.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_MaxMissingCount_1 — `--max-missing-count 1 --recode` on
// sample.vcf. Upstream check: `(N_chr - N_non_missing_chr) > 1` drops the
// site. The two sites with missing data:
//   - 20:1235237 (`0/0  0|0  ./.`)        → 2 missing chr  → drops (2 > 1)
//   - X:11      (`.:3:10  ./.  0|2:3`)    → 3 missing chr  → drops (3 > 1)
//
// All other 10 sites pass. Pinned byte-for-byte vs upstream golden.
func TestParity_MaxMissingCount_1(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MaxMissingCount:    1,
		MaxMissingCountSet: true,
		Recode:             true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "max_missing_count_1.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_MaxMissingCount_2 — `--max-missing-count 2 --recode` on
// sample.vcf. The boundary case: 20:1235237 has exactly 2 missing chrs.
// Upstream uses strict `>` so it KEEPS the site (2 > 2 is false). Our
// implementation matches. X:11 still drops (3 > 2). Pinned to catch any
// off-by-one in the comparison.
func TestParity_MaxMissingCount_2(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{
		MaxMissingCount:    2,
		MaxMissingCountSet: true,
		Recode:             true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "max_missing_count_2.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// PCA parity tests live in pca_test.go alongside the algorithm-level
// unit tests. The wave-19 implementation replaces the prior deferral.

// TestSNPHWE_Boundaries — unit test for the SNPHWE port against the
// hand-computable cases.
func TestSNPHWE_Boundaries(t *testing.T) {
	cases := []struct {
		name             string
		hom1, het, hom2  int
		wantMin, wantMax float64
	}{
		// Empty site: p_hwe is vacuously 1.0.
		{"empty", 0, 0, 0, 1.0, 1.0},
		// All homozygous-ref: p_hwe is 1.0 (no rare allele observed).
		{"all-hom1", 10, 0, 0, 1.0, 1.0},
		// Exact equilibrium (5 hom-ref, 10 het, 5 hom-alt; p=q=0.5):
		// observed het is the mode of the distribution, p_hwe == 1.0.
		{"in-hwe-balanced", 5, 10, 5, 0.99, 1.0},
		// Extreme deviation (10 hom-ref, 0 het, 10 hom-alt): exact
		// p_hwe << 1e-5. Verified against upstream's --hardy output
		// in TestParity_HWE_005_fixture's docstring.
		{"out-hwe-extreme", 10, 0, 10, 0.0, 1e-5},
		// One het, otherwise homozygous: a rare het count is well
		// within the typical regime (no extreme excess or deficit).
		{"single-het", 4, 1, 0, 0.5, 1.0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			p := snpHWE(tc.het, tc.hom1, tc.hom2)
			if p < tc.wantMin || p > tc.wantMax {
				t.Errorf("snpHWE(%d,%d,%d) = %g, want in [%g, %g]",
					tc.het, tc.hom1, tc.hom2, p, tc.wantMin, tc.wantMax)
			}
		})
	}
}

// -----------------------------------------------------------------------------
// --kept-sites / --removed-sites: 2-column (CHROM, POS) trace of which sites
// pass / fail filtering. Upstream registration: parameters.cpp:268, 330.
// Implementation: variant_file_output.cpp:4285-4373 (output_kept_sites and
// output_removed_sites). Each writer emits a `CHROM\tPOS` header and one
// row per site, in input file order.
// -----------------------------------------------------------------------------

// TestParity_KeptSites_NoFilter — `--kept-sites` against the 4-site
// hwe_fixture.vcf with no filtering. Every site survives, so the file is
// the header plus all four rows in input order.
func TestParity_KeptSites_NoFilter(t *testing.T) {
	prefix := runVcftoolsParity(t, "hwe_fixture.vcf", &Params{KeptSites: true})
	got := readFileBytes(t, prefix+".kept.sites")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "kept_sites_nofilter.expected.kept.sites"))
	if !bytes.Equal(got, want) {
		t.Errorf(".kept.sites mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_KeptSites_HWE — `--hwe 0.05 --kept-sites` against
// hwe_fixture.vcf. `--hwe` implies max_alleles = 2 (parameters.cpp:254).
// Sites 1:100 and 1:400 are in HWE; 1:200 and 1:300 fail the exact test.
// Golden generated from upstream binary (no LAPACK required).
func TestParity_KeptSites_HWE(t *testing.T) {
	prefix := runVcftoolsParity(t, "hwe_fixture.vcf", &Params{
		MinHWEPvalue: 0.05,
		MaxAlleles:   2, // upstream parameters.cpp:254 forces this when --hwe is set
		KeptSites:    true,
	})
	got := readFileBytes(t, prefix+".kept.sites")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "kept_sites_hwe.expected.kept.sites"))
	if !bytes.Equal(got, want) {
		t.Errorf(".kept.sites mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_RemovedSites_HWE — counterpart to TestParity_KeptSites_HWE.
// The two HWE-failing sites (1:200 and 1:300) appear in
// .removed.sites; the two HWE-passing sites do not.
func TestParity_RemovedSites_HWE(t *testing.T) {
	prefix := runVcftoolsParity(t, "hwe_fixture.vcf", &Params{
		MinHWEPvalue: 0.05,
		MaxAlleles:   2,
		RemovedSites: true,
	})
	got := readFileBytes(t, prefix+".removed.sites")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "removed_sites_hwe.expected.removed.sites"))
	if !bytes.Equal(got, want) {
		t.Errorf(".removed.sites mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_KeptSites_PosFilter — `--from-bp 150 --to-bp 350
// --kept-sites` on hwe_fixture.vcf. Position filter keeps 1:200 and 1:300;
// the file lists only those two rows in input order.
func TestParity_KeptSites_PosFilter(t *testing.T) {
	prefix := runVcftoolsParity(t, "hwe_fixture.vcf", &Params{
		Chr:       "1",
		FromBp:    150,
		ToBp:      350,
		KeptSites: true,
	})
	got := readFileBytes(t, prefix+".kept.sites")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "kept_sites_pos.expected.kept.sites"))
	if !bytes.Equal(got, want) {
		t.Errorf(".kept.sites mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_RemovedSites_PosFilter — counterpart of
// TestParity_KeptSites_PosFilter. The two out-of-range sites (1:100 and
// 1:400) appear in .removed.sites in input order.
func TestParity_RemovedSites_PosFilter(t *testing.T) {
	prefix := runVcftoolsParity(t, "hwe_fixture.vcf", &Params{
		Chr:          "1",
		FromBp:       150,
		ToBp:         350,
		RemovedSites: true,
	})
	got := readFileBytes(t, prefix+".removed.sites")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "removed_sites_pos.expected.removed.sites"))
	if !bytes.Equal(got, want) {
		t.Errorf(".removed.sites mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestKeptRemoved_Disjoint_And_Complete — sanity-check that --kept-sites
// and --removed-sites partition the input perfectly: every input site
// appears in exactly one of the two files, no duplicates, and the union
// equals the input.
//
// Upstream forbids both flags in one invocation (parameters.cpp:685
// "Only one output function may be called"), so this is a port-only
// invariant — we deliberately do not replicate that constraint per
// CLAUDE.md's "don't replicate upstream bugs" rule. The combined
// invocation is strictly more useful than either alone.
func TestKeptRemoved_Disjoint_And_Complete(t *testing.T) {
	prefix := runVcftoolsParity(t, "hwe_fixture.vcf", &Params{
		MinHWEPvalue: 0.05,
		MaxAlleles:   2,
		KeptSites:    true,
		RemovedSites: true,
	})
	kept := readFileLines(t, prefix+".kept.sites")
	removed := readFileLines(t, prefix+".removed.sites")
	// Strip headers.
	if len(kept) == 0 || kept[0] != "CHROM\tPOS" {
		t.Fatalf(".kept.sites: bad header %q", kept[0])
	}
	if len(removed) == 0 || removed[0] != "CHROM\tPOS" {
		t.Fatalf(".removed.sites: bad header %q", removed[0])
	}
	kept = kept[1:]
	removed = removed[1:]
	seen := make(map[string]int)
	for _, s := range kept {
		seen[s]++
	}
	for _, s := range removed {
		seen[s]++
	}
	// 4 input sites, each must appear exactly once across the two files.
	want := []string{"1\t100", "1\t200", "1\t300", "1\t400"}
	for _, w := range want {
		if seen[w] != 1 {
			t.Errorf("site %q appears %d times (want 1)", w, seen[w])
		}
	}
	if len(seen) != len(want) {
		t.Errorf("got %d unique sites, want %d", len(seen), len(want))
	}
}

// TestKeptRemoved_Disabled_NoFiles — with neither flag set, neither file
// is created (we don't leak empty `.kept.sites` / `.removed.sites` files).
func TestKeptRemoved_Disabled_NoFiles(t *testing.T) {
	prefix := runVcftoolsParity(t, "hwe_fixture.vcf", &Params{
		MinHWEPvalue: 0.05,
		MaxAlleles:   2,
	})
	if _, err := os.Stat(prefix + ".kept.sites"); !os.IsNotExist(err) {
		t.Errorf(".kept.sites was created when --kept-sites was not requested")
	}
	if _, err := os.Stat(prefix + ".removed.sites"); !os.IsNotExist(err) {
		t.Errorf(".removed.sites was created when --removed-sites was not requested")
	}
}

// -----------------------------------------------------------------------------
// Wave 10: --remove-filtered-geno, --remove-filtered-geno-all, --max-indv,
// --keep-INFO-all, --version
// -----------------------------------------------------------------------------

// TestParity_RemoveFilteredGenoAll — `--remove-filtered-geno-all --recode`
// rewrites GT to ./. for any kept genotype whose FORMAT FT field is not
// "PASS" or ".". Ported from upstream parameters.cpp:323 +
// vcf_entry.cpp:580-608 (filter_genotypes_by_filter_status with
// remove_all=true). The fixture covers all four FT shapes used by upstream
// (PASS / explicit flag / "." / multi-flag) — see ft_geno.vcf header.
func TestParity_RemoveFilteredGenoAll(t *testing.T) {
	prefix := runVcftoolsParity(t, "ft_geno.vcf", &Params{
		RemoveFilteredGenoAll: true,
		Recode:                true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "ft_geno_all.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_RemoveFilteredGenoQ10 — `--remove-filtered-geno q10 --recode`
// drops only genotypes whose FT lists `q10`; others (including the
// `lowDP`-tagged genotype on site 2) are kept. Mirrors upstream's
// vcf_entry.cpp:601-605 (loop over FT entries, set include_genotype=false
// on match).
func TestParity_RemoveFilteredGenoQ10(t *testing.T) {
	prefix := runVcftoolsParity(t, "ft_geno.vcf", &Params{
		RemoveFilteredGenoList: []string{"q10"},
		Recode:                 true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "ft_geno_q10.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_RemoveFilteredGenoMulti — `--remove-filtered-geno q10
// --remove-filtered-geno lowDP --recode` accepts the flag twice and ORs
// the two named flags into the drop set. Pins the upstream behaviour of
// parameters.cpp:324, which inserts each value into the same set.
func TestParity_RemoveFilteredGenoMulti(t *testing.T) {
	prefix := runVcftoolsParity(t, "ft_geno.vcf", &Params{
		RemoveFilteredGenoList: []string{"q10", "lowDP"},
		Recode:                 true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "ft_geno_multi.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestRemoveFilteredGeno_NoFT_NoOp — sites without a FORMAT FT column are
// left untouched, matching upstream's filter_genotypes_by_filter_status
// in entry_filters.cpp:94-108 (returns early when FT_idx == -1).
func TestRemoveFilteredGeno_NoFT_NoOp(t *testing.T) {
	basePrefix := runVcftoolsParity(t, "sample.vcf", &Params{Recode: true})
	base := readFileBytes(t, basePrefix+".recode.vcf")

	allPrefix := runVcftoolsParity(t, "sample.vcf", &Params{
		RemoveFilteredGenoAll: true,
		Recode:                true,
	})
	if got := readFileBytes(t, allPrefix+".recode.vcf"); !bytes.Equal(got, base) {
		t.Errorf("--remove-filtered-geno-all should be a no-op when no FT column is present; diff vs baseline")
	}

	namedPrefix := runVcftoolsParity(t, "sample.vcf", &Params{
		RemoveFilteredGenoList: []string{"q10"},
		Recode:                 true,
	})
	if got := readFileBytes(t, namedPrefix+".recode.vcf"); !bytes.Equal(got, base) {
		t.Errorf("--remove-filtered-geno NAME should be a no-op when no FT column is present; diff vs baseline")
	}
}

// TestMaxIndv_Count — `--max-indv N` caps the number of kept individuals
// at N. Upstream's filter_individuals_randomly uses srand(time(NULL))
// + random_shuffle so the kept identity is non-deterministic; this port
// keeps the first N in input order. We pin the COUNT only (the strongest
// claim that's portable across upstream's randomness).
func TestMaxIndv_Count(t *testing.T) {
	tests := []struct {
		name string
		n    int
		want int
	}{
		{"cap-1", 1, 1},
		{"cap-2", 2, 2},
		{"cap-3-exact", 3, 3},
		{"cap-above-noop", 5, 3},
		{"cap-0", 0, 0},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			prefix := runVcftoolsParity(t, "sample.vcf", &Params{
				MaxIndv:    tc.n,
				MaxIndvSet: true,
				Recode:     true,
			})
			lines := readFileLines(t, prefix+".recode.vcf")
			var header string
			for _, ln := range lines {
				if strings.HasPrefix(ln, "#CHROM") {
					header = ln
					break
				}
			}
			if header == "" {
				t.Fatalf("no #CHROM header in recode output")
			}
			cols := strings.Split(header, "\t")
			// columns 0..8 are fixed (CHROM..FORMAT); samples follow.
			gotSamples := 0
			if len(cols) > 9 {
				gotSamples = len(cols) - 9
			}
			if gotSamples != tc.want {
				t.Errorf("--max-indv %d kept %d samples, want %d", tc.n, gotSamples, tc.want)
			}
		})
	}
}

// TestMaxIndv_Unset_NoOp — when --max-indv is not supplied, all samples
// survive (MaxIndvSet must gate the cap, since Go's zero value for the
// underlying int would otherwise look like an explicit `--max-indv 0`).
func TestMaxIndv_Unset_NoOp(t *testing.T) {
	prefix := runVcftoolsParity(t, "sample.vcf", &Params{Recode: true})
	lines := readFileLines(t, prefix+".recode.vcf")
	var header string
	for _, ln := range lines {
		if strings.HasPrefix(ln, "#CHROM") {
			header = ln
			break
		}
	}
	cols := strings.Split(header, "\t")
	gotSamples := 0
	if len(cols) > 9 {
		gotSamples = len(cols) - 9
	}
	if gotSamples != 3 {
		t.Errorf("default (no --max-indv) kept %d samples, want 3", gotSamples)
	}
}

// `--keep-INFO-all` is wired as a CLI-only synonym for `--recode-INFO-all`
// in tools/vcftools/cmd/vcftools/main.go (an OR into the same
// Params.RecodeInfoAll bit). There is no package-level test for the
// synonym because both flags collapse to the same struct field before
// reaching `Run` — a unit test of `Run` cannot distinguish them.
// Reviewer-flagged on PR #134: the previous TestKeepINFOAll_Synonym
// was tautological (ran the same params twice). The synonym wiring is
// covered by manual smoke-test (PR #134 body) and by the explicit OR
// at main.go:540. A future CLI-exec test could pin it; tracked as a
// next-PR follow-up.

// -----------------------------------------------------------------------------
// Wave 12: --positions-overlap / --exclude-positions-overlap
// -----------------------------------------------------------------------------
//
// Fixture (pos_overlap_fixture.vcf) has six sites engineered to exercise
// overlap semantics: 1:200 has REF "ACGT" (covers 200..203), 1:400 has REF
// "TTC" (covers 400..402), and 2:100 has REF "GCGC" (covers 100..103). The
// other three sites (1:100, 1:300, 1:500) are 1-base REFs that coincide
// with plain --positions behaviour.
//
// pos_overlap_keep.txt lists 1:202 (mid-REF of 1:200), 1:300 (exact),
// 1:402 (last base of 1:400), 1:600 (no site), and 2:101 (mid-REF of
// 2:100). pos_overlap_exclude.txt lists 1:202, 1:300, and 2:103 (last
// base of 2:100).
//
// Goldens (pos_overlap_keep.expected.recode.vcf and the exclude
// counterpart) were produced by the upstream binary at
// reference_code/vcftools/src/cpp/vcftools (built with
// CXXFLAGS='-O0 -g -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0').

// TestParity_PositionsOverlap_Keep — `--positions-overlap FILE --recode`
// keeps a record when ANY base in [POS, POS+len(REF)-1] matches a position
// in the file. Ported from upstream parameters.cpp:315 +
// entry_filters.cpp:408-531.
func TestParity_PositionsOverlap_Keep(t *testing.T) {
	prefix := runVcftoolsParity(t, "pos_overlap_fixture.vcf", &Params{
		PositionsOverlapFile: filepath.Join(vcftoolsFixtureDir(t), "pos_overlap_keep.txt"),
		Recode:               true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "pos_overlap_keep.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_PositionsOverlap_Exclude — `--exclude-positions-overlap FILE
// --recode` drops a record when ANY base in [POS, POS+len(REF)-1] matches
// a position in the file. Ported from upstream parameters.cpp:221 +
// entry_filters.cpp:533-547.
func TestParity_PositionsOverlap_Exclude(t *testing.T) {
	prefix := runVcftoolsParity(t, "pos_overlap_fixture.vcf", &Params{
		ExcludePositionsOverlapFile: filepath.Join(vcftoolsFixtureDir(t), "pos_overlap_exclude.txt"),
		Recode:                      true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "pos_overlap_exclude.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestPositionsOverlap_VsPlain_DivergesOnMultiBaseRef — pins the
// behavioural difference between --positions (POS-only) and
// --positions-overlap (POS..POS+len(REF)-1) on multi-base REF records.
// The overlap form should keep 1:200 (REF=ACGT) given a positions file
// listing 1:202; the plain form should NOT. This is the whole reason
// upstream ships the overlap variant.
func TestPositionsOverlap_VsPlain_DivergesOnMultiBaseRef(t *testing.T) {
	dir := vcftoolsFixtureDir(t)

	// Build a single-line positions file that targets the interior of
	// the 1:200 REF "ACGT". Tmpfile so the test stays hermetic.
	pos := filepath.Join(t.TempDir(), "interior.txt")
	if err := os.WriteFile(pos, []byte("1\t202\n"), 0o644); err != nil {
		t.Fatalf("write positions: %v", err)
	}

	// Plain --positions: should drop 1:200 (POS != 202).
	plain := runVcftoolsParity(t, "pos_overlap_fixture.vcf", &Params{
		PositionsFile: pos,
		Recode:        true,
	})
	plainLines := readFileLines(t, plain+".recode.vcf")
	for _, ln := range plainLines {
		if strings.HasPrefix(ln, "#") || ln == "" {
			continue
		}
		if strings.HasPrefix(ln, "1\t200\t") {
			t.Errorf("--positions unexpectedly kept 1:200 with pos=1:202")
		}
	}

	// --positions-overlap: should keep 1:200 (interior overlap).
	overlap := runVcftoolsParity(t, "pos_overlap_fixture.vcf", &Params{
		PositionsOverlapFile: pos,
		Recode:               true,
	})
	overlapLines := readFileLines(t, overlap+".recode.vcf")
	kept200 := false
	for _, ln := range overlapLines {
		if strings.HasPrefix(ln, "1\t200\t") {
			kept200 = true
		}
	}
	if !kept200 {
		t.Errorf("--positions-overlap dropped 1:200 even though pos=1:202 is interior to REF=ACGT")
		_ = dir
	}
}

// TestPositionsOverlap_BoundaryHits — pin overlap behaviour at both ends
// of the REF window for both flags. For REF="ACGT" at POS=200, the
// half-open upstream loop is `for ui=POS; ui<POS+REF.size(); ui++` so
// matches are inclusive at both ends: 200 (first base) and 203 (last)
// both qualify, but 199 and 204 do not.
func TestPositionsOverlap_BoundaryHits(t *testing.T) {
	tests := []struct {
		name     string
		posLine  string
		wantKept bool
	}{
		{"first_base_200", "1\t200\n", true},
		{"last_base_203", "1\t203\n", true},
		{"one_before_199", "1\t199\n", false},
		{"one_after_204", "1\t204\n", false},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			pos := filepath.Join(t.TempDir(), "pos.txt")
			if err := os.WriteFile(pos, []byte(tc.posLine), 0o644); err != nil {
				t.Fatalf("write: %v", err)
			}
			prefix := runVcftoolsParity(t, "pos_overlap_fixture.vcf", &Params{
				PositionsOverlapFile: pos,
				Recode:               true,
			})
			lines := readFileLines(t, prefix+".recode.vcf")
			kept200 := false
			for _, ln := range lines {
				if strings.HasPrefix(ln, "1\t200\t") {
					kept200 = true
				}
			}
			if kept200 != tc.wantKept {
				t.Errorf("--positions-overlap pos=%q: kept 1:200 = %v, want %v",
					strings.TrimSpace(tc.posLine), kept200, tc.wantKept)
			}
		})
	}
}

// TestPositionsOverlap_UnknownChromDropped — when a positions-overlap
// file is supplied, sites on chromosomes NOT mentioned in the file are
// dropped. Mirrors upstream entry_filters.cpp:515-516 ("if
// chr_to_idx.find(CHROM) == chr_to_idx.end() passed_filters = false").
func TestPositionsOverlap_UnknownChromDropped(t *testing.T) {
	pos := filepath.Join(t.TempDir(), "pos.txt")
	// Only chromosome "1" is named.
	if err := os.WriteFile(pos, []byte("1\t100\n"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	prefix := runVcftoolsParity(t, "pos_overlap_fixture.vcf", &Params{
		PositionsOverlapFile: pos,
		Recode:               true,
	})
	for _, ln := range readFileLines(t, prefix+".recode.vcf") {
		if strings.HasPrefix(ln, "2\t") {
			t.Errorf("chr 2 site leaked through --positions-overlap: %q", ln)
		}
	}
}

// TestExcludePositionsOverlap_UnknownChromKept — converse: when an
// exclude-overlap file is supplied, sites on chromosomes NOT named are
// passed through unchanged. Mirrors entry_filters.cpp:535 ("if
// chr_to_idx.find(CHROM) != chr_to_idx.end()" — no match means no drop).
func TestExcludePositionsOverlap_UnknownChromKept(t *testing.T) {
	pos := filepath.Join(t.TempDir(), "pos.txt")
	if err := os.WriteFile(pos, []byte("1\t100\n"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	prefix := runVcftoolsParity(t, "pos_overlap_fixture.vcf", &Params{
		ExcludePositionsOverlapFile: pos,
		Recode:                      true,
	})
	saw2 := false
	for _, ln := range readFileLines(t, prefix+".recode.vcf") {
		if strings.HasPrefix(ln, "2\t") {
			saw2 = true
		}
	}
	if !saw2 {
		t.Errorf("chr 2 sites dropped even though chr 2 isn't in --exclude-positions-overlap file")
	}
}

// TestPositionsOverlap_MissingFile — surface an error when the file
// doesn't exist (don't silently produce empty output).
func TestPositionsOverlap_MissingFile(t *testing.T) {
	tmp := t.TempDir()
	in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), "pos_overlap_fixture.vcf"))
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer in.Close()
	params := &Params{
		PositionsOverlapFile: filepath.Join(tmp, "does-not-exist.txt"),
		OutPrefix:            filepath.Join(tmp, "out"),
		Recode:               true,
	}
	if err := Run(in, params); err == nil {
		t.Errorf("expected error for missing --positions-overlap file, got nil")
	}
}

// -----------------------------------------------------------------------------
// Wave 17: --keep-INFO TAG (SITE FILTER)
// -----------------------------------------------------------------------------
//
// Fixture (keep_info_flags.vcf) declares two Flag-type INFO tags
// (FLAG_A, FLAG_B) and one Integer tag (DP). Sites engineered to exercise
// every combination:
//
//   1:100 — FLAG_A only
//   1:200 — FLAG_B only
//   1:300 — both FLAG_A and FLAG_B
//   1:400 — neither (DP only)
//
// Goldens were produced by the upstream binary
// /tmp/vcftools_install/bin/vcftools (VCFtools 0.1.18, built locally
// with the FORTIFY_SOURCE workaround). The fixture intentionally does
// NOT use --recode-INFO-all so that the recoded INFO column is
// upstream's canonical "." (strip-all-INFO) form; this side-steps the
// pre-existing port-level INFO-ordering note documented in
// tools/vcftools/cmd/vcftools/aliases_cli_test.go:138-144 and produces
// byte-for-byte matching .recode.vcf files.

// TestParity_KeepINFO_SingleFlag — `--keep-INFO FLAG_A --recode` keeps
// only sites where FLAG_A is present (sites 100 and 300). Site-set
// parity is byte-for-byte against upstream; the INFO column is "." on
// both sides because --recode-INFO-all is not set. Mirrors upstream's
// entry_filters.cpp:1033-1063 (filter_sites_by_INFO).
func TestParity_KeepINFO_SingleFlag(t *testing.T) {
	prefix := runVcftoolsParity(t, "keep_info_flags.vcf", &Params{
		KeepINFO: "FLAG_A",
		Recode:   true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "keep_info_flagA.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_KeepINFO_OR — multiple --keep-INFO tags compose via OR.
// Upstream entry_filters.cpp:1049-1062 iterates over flags_to_keep
// and sets `keep=true` on the first present tag, then writes
// `passed_filters = keep`. With FLAG_A and FLAG_B as the keep set,
// sites 100, 200, and 300 pass; only site 400 (neither flag) is
// dropped. Byte-for-byte parity against the upstream binary.
func TestParity_KeepINFO_OR(t *testing.T) {
	prefix := runVcftoolsParity(t, "keep_info_flags.vcf", &Params{
		KeepINFO: "FLAG_A,FLAG_B",
		Recode:   true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "keep_info_or.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// -----------------------------------------------------------------------------
// Wave 18: --remove-INFO TAG (SITE FILTER)
// -----------------------------------------------------------------------------
//
// Uses the same keep_info_flags.vcf fixture as wave 17 (see header
// comment above the wave-17 block). Wave 18 polarity-inverts:
// --remove-INFO drops sites where the named Flag IS present, the
// complement of --keep-INFO.

// TestParity_RemoveINFO_SingleFlag — `--remove-INFO FLAG_A --recode`
// drops sites where FLAG_A is present (sites 100 and 300); sites 200
// and 400 survive. Site-set parity is byte-for-byte against upstream;
// the INFO column is "." on both sides because --recode-INFO-all is
// not set. Mirrors upstream's entry_filters.cpp:1068-1086.
func TestParity_RemoveINFO_SingleFlag(t *testing.T) {
	prefix := runVcftoolsParity(t, "keep_info_flags.vcf", &Params{
		RemoveINFO: "FLAG_A",
		Recode:     true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "remove_info_flagA.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_RemoveINFO_OR — multiple --remove-INFO tags veto via OR.
// Upstream entry_filters.cpp:1070-1083 iterates over flags_to_remove
// and sets `passed_filters = false` on the first present tag. With
// FLAG_A and FLAG_B as the remove set, only site 400 (neither flag)
// survives.
func TestParity_RemoveINFO_OR(t *testing.T) {
	prefix := runVcftoolsParity(t, "keep_info_flags.vcf", &Params{
		RemoveINFO: "FLAG_A,FLAG_B",
		Recode:     true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "remove_info_or.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_KeepAndRemoveINFO_Compose — upstream's keep-then-remove
// composition (entry_filters.cpp:1033-1086). With --keep-INFO FLAG_A
// and --remove-INFO FLAG_B, only site 100 (FLAG_A present AND FLAG_B
// absent) survives. Site-set parity is byte-for-byte; --recode-INFO-all
// is deliberately omitted so INFO collapses to "." on both sides (this
// side-steps the port-vs-upstream INFO-key-ordering quirk documented in
// tools/vcftools/cmd/vcftools/aliases_cli_test.go:138-144).
func TestParity_KeepAndRemoveINFO_Compose(t *testing.T) {
	prefix := runVcftoolsParity(t, "keep_info_flags.vcf", &Params{
		KeepINFO:   "FLAG_A",
		RemoveINFO: "FLAG_B",
		Recode:     true,
	})
	got := readFileBytes(t, prefix+".recode.vcf")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "keep_remove_info_compose.expected.recode.vcf"))
	if !bytes.Equal(got, want) {
		t.Errorf(".recode.vcf mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}
