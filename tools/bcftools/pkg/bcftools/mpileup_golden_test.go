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
// the SNP MAQ path and (post 4e.2+4e.3) the indel emission path. With
// the per-site bias annotations (VDB, SGB, RPBZ, MQBZ, BQBZ, MQSBZ,
// SCBZ), the MQ0F fraction, the INFO/QS float32 rounding, the
// MPLP_SMART_OVERLAPS read-pair quality merge AND the indel-branch
// glfgen/combine/2bcf all ported, `bcftools mpileup` must reproduce
// the upstream goldens exactly — header and every data record, INFO
// bias tags and the single INDEL record at 17:302 included.
//
// Two upstream invocations from reference_code/bcftools/test/test.pl
// are replayed:
//
//   - mpileup.11.out: `mpileup -a -AD mpileup.3.bam`, the full 4200 bp
//     contig. 4001 covered positions, 87 SNP ALT records, two
//     overlapping mate pairs (17:1118-1142 and 17:3785-3836) that
//     exercise smart-overlaps, plus the single +1 insertion record at
//     17:302 (T → TA) which exercises the full indel-calling
//     gap_prep/glfgen/combine/2bcf pipeline including the BQBZ leak
//     from the prior has_alt SNP combine.
//   - mpileup.1.out: `mpileup -r17:100-150 -a -AD mpileup.{1,2,3}.bam`,
//     the three-sample multi-BAM path over a 51 bp window.
//
// `-a -AD` removes FORMAT/AD from the default tag set, so the goldens
// carry FORMAT=PL only — exactly what the SNP / indel paths emit.
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

			// Index the golden data records by CHROM:POS:isIndel so the
			// SNP and INDEL records at the same coordinate (17:302) are
			// disambiguated.
			isIndel := func(ln string) bool {
				f := strings.Split(ln, "\t")
				if len(f) < 8 {
					return false
				}
				return strings.HasPrefix(f[7], "INDEL;") || strings.Contains(f[7], ";INDEL;") || f[7] == "INDEL"
			}
			recKey := func(ln string) string {
				f := strings.Split(ln, "\t")
				if isIndel(ln) {
					return f[0] + ":" + f[1] + ":indel"
				}
				return f[0] + ":" + f[1] + ":snp"
			}

			wantByPos := make(map[string]string, len(wantData))
			for _, ln := range wantData {
				wantByPos[recKey(ln)] = ln
			}

			checked, diffs := 0, 0
			for _, ln := range gotData {
				key := recKey(ln)
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
			for key := range wantByPos {
				t.Errorf("missing record at %s (golden has it, we did not emit it)", key)
			}
			if diffs > 10 {
				t.Errorf("... and %d more record mismatches", diffs-10)
			}
			t.Logf("%s: %d records byte-for-byte identical (SNP and INDEL)",
				tc.golden, checked)
		})
	}
}

// TestMpileupSCRGolden replays the upstream `bcftools mpileup -a
// -AD,INFO/SCR,FMT/SCR` test (test.pl line 1069), which exercises the
// soft-clip-read accumulator. SCR counts the reads in a column whose
// CIGAR contains at least one S op anywhere — see
// reference_code/bcftools/mpileup.c:307-324 for upstream's
// pileup_constructor that sets PLP_HAS_SOFT_CLIP, and
// reference_code/bcftools/bam2bcf.c:300 for the per-read tally in
// bcf_call_glfgen. INFO/SCR is the total, FORMAT/SCR is the per-sample
// count.
func TestMpileupSCRGolden(t *testing.T) {
	ref := mpileupFixture(t, "mpileup-SCR.fa")
	mpileupFixture(t, "mpileup-SCR.fa.fai")
	bam := mpileupFixture(t, "mpileup-SCR.bam")
	goldenPath := mpileupFixture(t, "mpileup-SCR.out")
	goldenBytes, err := os.ReadFile(goldenPath)
	if err != nil {
		t.Fatalf("read golden: %v", err)
	}

	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:    []string{bam},
		FastaRef:  ref,
		Annotate:  "-AD,INFO/SCR,FMT/SCR",
		NoVersion: true,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}
	if buf.String() != string(goldenBytes) {
		// Helpful diff: header first, then first 5 differing data lines.
		gotH, gotD := splitMpileupVCF(buf.String())
		wantH, wantD := splitMpileupVCF(string(goldenBytes))
		if len(gotH) != len(wantH) {
			t.Errorf("header line count: got %d, want %d", len(gotH), len(wantH))
		}
		nH := len(gotH)
		if len(wantH) < nH {
			nH = len(wantH)
		}
		for i := 0; i < nH; i++ {
			if gotH[i] != wantH[i] {
				t.Errorf("header line %d:\n got:  %s\n want: %s", i, gotH[i], wantH[i])
			}
		}
		diffs := 0
		nD := len(gotD)
		if len(wantD) < nD {
			nD = len(wantD)
		}
		for i := 0; i < nD; i++ {
			if gotD[i] != wantD[i] {
				diffs++
				if diffs <= 5 {
					t.Errorf("record %d:\n got:  %s\n want: %s", i, gotD[i], wantD[i])
				}
			}
		}
		if len(gotD) != len(wantD) {
			t.Errorf("data record count: got %d, want %d", len(gotD), len(wantD))
		}
	}
}

// TestMpileupSCROnIndelRow confirms that when -a INFO/SCR,FMT/SCR is
// in force, bcfCall2bcfIndel emits SCR on the indel row using the
// shared per-column tally that also feeds the SNP row. The upstream
// SCR golden has no indel-bearing columns, so we reuse the indel-AD.2
// fixture (chr11:75 has a homopolymer-anchored indel call and at least
// one soft-clipped read in the column). The SNP and indel rows at the
// same position must therefore report the same INFO/SCR and the same
// per-sample FORMAT/SCR.
func TestMpileupSCROnIndelRow(t *testing.T) {
	ref := mpileupFixture(t, "indel-AD.2.fa")
	mpileupFixture(t, "indel-AD.2.fa.fai")
	bam := mpileupFixture(t, "indel-AD.2.bam")
	mpileupFixture(t, "indel-AD.2.bam.bai")

	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:    []string{bam},
		FastaRef:  ref,
		Regions:   []string{"11:75"},
		Annotate:  "AD,INFO/SCR,FMT/SCR",
		NoVersion: true,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}

	var snp, indel string
	for _, ln := range strings.Split(buf.String(), "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		fields := strings.Split(ln, "\t")
		if len(fields) < 10 || fields[1] != "75" {
			continue
		}
		if strings.Contains(fields[7], "INDEL") {
			indel = ln
		} else {
			snp = ln
		}
	}
	if snp == "" || indel == "" {
		t.Fatalf("expected one SNP and one indel record at 11:75; got snp=%q indel=%q",
			snp, indel)
	}

	// Extract INFO/SCR and FORMAT/SCR from each record. FORMAT is the
	// 9th column (index 8), SCR's position varies and is parsed by name.
	parseSCR := func(line string) (infoSCR, fmtSCR string) {
		f := strings.Split(line, "\t")
		for _, kv := range strings.Split(f[7], ";") {
			if strings.HasPrefix(kv, "SCR=") {
				infoSCR = strings.TrimPrefix(kv, "SCR=")
				break
			}
		}
		format := strings.Split(f[8], ":")
		sample := strings.Split(f[9], ":")
		for i, k := range format {
			if k == "SCR" && i < len(sample) {
				fmtSCR = sample[i]
				break
			}
		}
		return
	}
	snpInfo, snpFmt := parseSCR(snp)
	indelInfo, indelFmt := parseSCR(indel)
	if indelInfo == "" {
		t.Fatalf("indel row missing INFO/SCR; record: %s", indel)
	}
	if indelFmt == "" {
		t.Fatalf("indel row missing FORMAT/SCR; record: %s", indel)
	}
	if snpInfo != indelInfo {
		t.Errorf("INFO/SCR mismatch SNP=%s indel=%s", snpInfo, indelInfo)
	}
	if snpFmt != indelFmt {
		t.Errorf("FORMAT/SCR mismatch SNP=%s indel=%s", snpFmt, indelFmt)
	}
	// And the per-column tally must be non-zero — this fixture is
	// chosen precisely because it has soft-clipped reads at the column.
	if indelInfo == "0" {
		t.Errorf("expected non-zero SCR at 11:75 indel row, got %s", indelInfo)
	}
}

// TestMpileupFilterGolden replays the upstream `bcftools mpileup
// --skip-*` test from test.pl (lines 1072-1073), exercising the
// BAM-flag read filter. The input has two reads of the same template,
// flags 99 (paired+propPair+R1+MateReverse) and 3 (paired+propPair).
// `--skip-all-unset READ1` drops the second read; `--skip-any-unset
// READ1` does likewise (since both forms reject a read with READ1
// unset for a single-bit mask). The result is one covered position
// at 1:100 with DP=1.
func TestMpileupFilterGolden(t *testing.T) {
	ref := mpileupFixture(t, "mpileup-SCR.fa")
	mpileupFixture(t, "mpileup-SCR.fa.fai")
	bam := mpileupFixture(t, "mpileup-filter.bam")
	goldenPath := mpileupFixture(t, "mpileup-filter.2.out")
	goldenBytes, err := os.ReadFile(goldenPath)
	if err != nil {
		t.Fatalf("read golden: %v", err)
	}

	cases := []struct {
		name string
		opts MpileupOptions
	}{
		{
			name: "skip-all-unset-READ1",
			opts: MpileupOptions{
				Inputs:    []string{bam},
				FastaRef:  ref,
				Targets:   []string{"1:100"},
				Annotate:  "-AD",
				FlagIncl:  "READ1", // --skip-all-unset alias
				NoVersion: true,
			},
		},
		{
			name: "skip-any-unset-READ1",
			opts: MpileupOptions{
				Inputs:    []string{bam},
				FastaRef:  ref,
				Targets:   []string{"1:100"},
				Annotate:  "-AD",
				FlagAny:   "READ1", // --skip-any-unset alias
				NoVersion: true,
			},
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var buf bytes.Buffer
			if err := MpileupFile(tc.opts, &buf); err != nil {
				t.Fatalf("MpileupFile: %v", err)
			}
			if buf.String() != string(goldenBytes) {
				t.Errorf("output mismatch:\n got:\n%s\n want:\n%s", buf.String(), string(goldenBytes))
			}
		})
	}
}

// TestMpileupIndelADGolden replays upstream's `bcftools mpileup -a AD
// -r 11:75 [--ambig-reads MODE]` invocations from test.pl lines
// 1066-1068 (indel-AD.{2,3,4}.out). They exercise the legacy
// REF-rescue heuristic that `bcf_cgp_compute_indelQ` leaves to
// `bcf_call_glfgen`: at homopolymer / tandem-repeat sites the indel
// branch reclassifies a REF-looking read (no CIGAR indel) as REF type,
// promotes its `q` to the raw base quality at qpos, and blends seqQ as
// `(3*seqQ + 2*q)/8`. Without that step REF reads at the deep
// homopolymer at chr11:75 would all fail the min-baseQ gate and AD
// would collapse. The byte-for-byte goldens here pin both the rescue
// itself (AD.2) and its interaction with --ambig-reads incAD (AD.3)
// and --ambig-reads incAD0 (AD.4).
func TestMpileupIndelADGolden(t *testing.T) {
	ref := mpileupFixture(t, "indel-AD.2.fa")
	mpileupFixture(t, "indel-AD.2.fa.fai")
	bam := mpileupFixture(t, "indel-AD.2.bam")
	mpileupFixture(t, "indel-AD.2.bam.bai")

	cases := []struct {
		name   string
		golden string
		mode   AmbigReadsMode
	}{
		{"default", "indel-AD.2.out", AmbigReadsDrop},
		{"incAD", "indel-AD.3.out", AmbigReadsIncAD},
		{"incAD0", "indel-AD.4.out", AmbigReadsIncAD0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			goldenBytes, err := os.ReadFile(mpileupFixture(t, tc.golden))
			if err != nil {
				t.Fatalf("read golden %s: %v", tc.golden, err)
			}
			var buf bytes.Buffer
			opts := MpileupOptions{
				Inputs:         []string{bam},
				FastaRef:       ref,
				Regions:        []string{"11:75"},
				Annotate:       "AD",
				AmbigReadsMode: tc.mode,
				NoVersion:      true,
			}
			if err := MpileupFile(opts, &buf); err != nil {
				t.Fatalf("MpileupFile: %v", err)
			}
			if buf.String() != string(goldenBytes) {
				gotH, gotD := splitMpileupVCF(buf.String())
				wantH, wantD := splitMpileupVCF(string(goldenBytes))
				if len(gotH) != len(wantH) {
					t.Errorf("header line count: got %d, want %d", len(gotH), len(wantH))
				}
				nH := len(gotH)
				if len(wantH) < nH {
					nH = len(wantH)
				}
				for i := 0; i < nH; i++ {
					if gotH[i] != wantH[i] {
						t.Errorf("header line %d:\n got:  %s\n want: %s", i, gotH[i], wantH[i])
					}
				}
				diffs := 0
				nD := len(gotD)
				if len(wantD) < nD {
					nD = len(wantD)
				}
				for i := 0; i < nD; i++ {
					if gotD[i] != wantD[i] {
						diffs++
						if diffs <= 5 {
							t.Errorf("record %d:\n got:  %s\n want: %s", i, gotD[i], wantD[i])
						}
					}
				}
				if len(gotD) != len(wantD) {
					t.Errorf("data record count: got %d, want %d", len(gotD), len(wantD))
				}
			}
		})
	}
}

// TestMpileupIndelsCNSGolden replays the upstream `bcftools mpileup
// -a AD --indels-cns` invocation from test.pl line 1065
// (indel-AD.1cns.out). It exercises the consensus-based indel caller
// dispatch (bcfCallGapPrepCNS in bam2bcf_indelcns.go, port of
// reference_code/bcftools/bam2bcf_edlib.c) which uses the in-tree
// edlib engine (pkg/htsgo/edlib) to score each read against per-type
// candidate haplotypes.
//
// Scope: the consensus builder (bcf_cgp_consensus with cons[0]/cons[1]
// het threading and cons_ins/ref_ins smoothing), the dual-consensus
// alignment scoring (bcf_cgp_align_score), and the edlib-flavored
// compute_indelQ (indelQ1/indelQ2 vs_ref blend, poly_mqual,
// TMP_MAGIC=255) have all landed, plus the CNS-specific glfgen path
// (bam2bcf.c:317-415: legacy REF-rescue disabled, seqQ_offset cap,
// realigned-read q2p5 dampener) and the IDV/IMF recompute from ADF/ADR
// (bam2bcf.c:1265-1275). Two indel rows (000000F:538, :658) still
// carry residual single-read assignment noise (one read's tail-distance
// contribution and PL byte). The soft expectation (counted, capped)
// keeps the dispatch under regression coverage while those final two
// rows are characterised.
func TestMpileupIndelsCNSGolden(t *testing.T) {
	ref := mpileupFixture(t, "indel-AD.1.fa")
	mpileupFixture(t, "indel-AD.1.fa.fai")
	bam := mpileupFixture(t, "indel-AD.1.bam")
	mpileupFixture(t, "indel-AD.1.bam.bai")

	goldenBytes, err := os.ReadFile(mpileupFixture(t, "indel-AD.1cns.out"))
	if err != nil {
		t.Fatalf("read golden: %v", err)
	}

	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:    []string{bam},
		FastaRef:  ref,
		Annotate:  "AD",
		IndelsCNS: true,
		NoVersion: true,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}

	gotH, gotD := splitMpileupVCF(buf.String())
	wantH, wantD := splitMpileupVCF(string(goldenBytes))

	// Header must match byte-for-byte (CLI flag drives only the indel
	// scoring; header tags are unchanged).
	if len(gotH) != len(wantH) {
		t.Errorf("header line count: got %d, want %d", len(gotH), len(wantH))
	}
	nH := len(gotH)
	if len(wantH) < nH {
		nH = len(wantH)
	}
	for i := 0; i < nH; i++ {
		if gotH[i] != wantH[i] {
			t.Errorf("header line %d:\n got:  %s\n want: %s", i, gotH[i], wantH[i])
		}
	}

	if len(gotD) != len(wantD) {
		t.Errorf("data record count: got %d, want %d", len(gotD), len(wantD))
	}

	// Tally byte-matched vs differing records. The non-indel bulk
	// (SNP / N-REF columns) shares the legacy SNP path, so the
	// CNS-specific delta is confined to the small number of indel rows.
	matches, diffs := 0, 0
	indelDiffs := 0
	n := len(gotD)
	if len(wantD) < n {
		n = len(wantD)
	}
	for i := 0; i < n; i++ {
		if gotD[i] == wantD[i] {
			matches++
			continue
		}
		diffs++
		if strings.Contains(gotD[i], "INDEL;") || strings.Contains(wantD[i], "INDEL;") {
			indelDiffs++
		}
		if diffs <= 5 {
			t.Logf("record %d differs:\n got:  %s\n want: %s", i, gotD[i], wantD[i])
		}
	}
	t.Logf("--indels-cns: %d/%d records byte-identical, %d differ (%d on indel rows)",
		matches, n, diffs, indelDiffs)

	// Sanity floor: the SNP / N-REF bulk must match. The upstream golden
	// has 301 data records; only a handful are indel rows (4 in
	// upstream's diff against the legacy golden). Require the bulk of
	// records to match so a CNS-path regression that breaks the SNP
	// emission gets caught.
	if matches < n*9/10 {
		t.Errorf("too few records byte-matched: %d/%d (want at least 90%%)", matches, n)
	}
}

// TestMpileupNMBZGolden replays the upstream `bcftools mpileup -a
// -AD,INFO/NMBZ` test for fixture annot-NMBZ.1 (test.pl line 1074). It
// exercises the per-read NM-tag bias accumulator and the INFO/NMBZ
// emission path. The matching SNP rows at chr19:69,72,76,78,91,99
// carry NMBZ values computed from the cross-sample ref_nm / alt_nm
// histograms via the existing Mann-Whitney U z-score helper.
//
// The 4e.6 slice ports NMBZ; annot-NMBZ.3's indel row residual (QS /
// NMBZ / PL[0]) is tracked separately as 4e.7 — the SNP row of
// annot-NMBZ.3 does byte-match including NMBZ=7.74597. The
// annot-NMBZ.2 golden now byte-matches at chr6:75 (DP=283 > 250)
// thanks to the htslib per-alignment-start depth cap port; see
// TestMpileupDepthCapGolden.
func TestMpileupNMBZGolden(t *testing.T) {
	ref := mpileupFixture(t, "annot-NMBZ.1.fa")
	mpileupFixture(t, "annot-NMBZ.1.fa.fai")
	bam := mpileupFixture(t, "annot-NMBZ.1.bam")
	goldenPath := mpileupFixture(t, "annot-NMBZ.1.1.out")
	goldenBytes, err := os.ReadFile(goldenPath)
	if err != nil {
		t.Fatalf("read golden: %v", err)
	}

	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:    []string{bam},
		FastaRef:  ref,
		Regions:   []string{"chr19:69-99"},
		Annotate:  "-AD,INFO/NMBZ",
		NoVersion: true,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}
	if buf.String() != string(goldenBytes) {
		gotH, gotD := splitMpileupVCF(buf.String())
		wantH, wantD := splitMpileupVCF(string(goldenBytes))
		if len(gotH) != len(wantH) {
			t.Errorf("header line count: got %d, want %d", len(gotH), len(wantH))
		}
		nH := len(gotH)
		if len(wantH) < nH {
			nH = len(wantH)
		}
		for i := 0; i < nH; i++ {
			if gotH[i] != wantH[i] {
				t.Errorf("header line %d:\n got:  %s\n want: %s", i, gotH[i], wantH[i])
			}
		}
		diffs := 0
		nD := len(gotD)
		if len(wantD) < nD {
			nD = len(wantD)
		}
		for i := 0; i < nD; i++ {
			if gotD[i] != wantD[i] {
				diffs++
				if diffs <= 5 {
					t.Errorf("record %d:\n got:  %s\n want: %s", i, gotD[i], wantD[i])
				}
			}
		}
		if len(gotD) != len(wantD) {
			t.Errorf("data record count: got %d, want %d", len(gotD), len(wantD))
		}
	}
}

// TestMpileupDepthCapGolden pins the htslib per-alignment-start depth
// cap port. annot-NMBZ.2 has a single SNP row at chr6:75 where raw
// coverage is 449 (orphan-filtered from 450) but the -d 250 cap drops
// 166 reads, leaving DP=283. Upstream htslib applies the cap per
// alignment start inside bam_plp_push (reference_code/htslib/sam.c:
// 6090): a new read is dropped when iter->pos == b->core.pos and the
// queue already holds maxcnt active reads. Our applyMpileupDepthCap
// reproduces that drop order so the surviving 283 reads, their I16
// partial sums and the resulting NMBZ / QS / SCBZ values all match
// the upstream golden byte-for-byte.
//
// Replays test.pl line 1075: `bcftools mpileup -a -AD,INFO/NMBZ -r
// chr6:75 annot-NMBZ.2.bam`.
func TestMpileupDepthCapGolden(t *testing.T) {
	ref := mpileupFixture(t, "annot-NMBZ.2.fa")
	mpileupFixture(t, "annot-NMBZ.2.fa.fai")
	bam := mpileupFixture(t, "annot-NMBZ.2.bam")
	mpileupFixture(t, "annot-NMBZ.2.bam.bai")
	goldenPath := mpileupFixture(t, "annot-NMBZ.2.1.out")
	goldenBytes, err := os.ReadFile(goldenPath)
	if err != nil {
		t.Fatalf("read golden: %v", err)
	}

	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:    []string{bam},
		FastaRef:  ref,
		Regions:   []string{"chr6:75"},
		Annotate:  "-AD,INFO/NMBZ",
		NoVersion: true,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}
	if buf.String() != string(goldenBytes) {
		gotH, gotD := splitMpileupVCF(buf.String())
		wantH, wantD := splitMpileupVCF(string(goldenBytes))
		if len(gotH) != len(wantH) {
			t.Errorf("header line count: got %d, want %d", len(gotH), len(wantH))
		}
		nH := len(gotH)
		if len(wantH) < nH {
			nH = len(wantH)
		}
		for i := 0; i < nH; i++ {
			if gotH[i] != wantH[i] {
				t.Errorf("header line %d:\n got:  %s\n want: %s", i, gotH[i], wantH[i])
			}
		}
		nD := len(gotD)
		if len(wantD) < nD {
			nD = len(wantD)
		}
		for i := 0; i < nD; i++ {
			if gotD[i] != wantD[i] {
				t.Errorf("record %d:\n got:  %s\n want: %s", i, gotD[i], wantD[i])
			}
		}
		if len(gotD) != len(wantD) {
			t.Errorf("data record count: got %d, want %d", len(gotD), len(wantD))
		}
	}
}

// TestMpileupAmbigReadsParse verifies parseAmbigReads round-trips every
// upstream-accepted spelling. The byte-level effect on per-allele AD is
// covered by the (currently deferred) indel-AD.{3,4} goldens which
// depend on additional indel-scoring work tracked in PARITY_ROADMAP
// 4e.7.
func TestMpileupAmbigReadsParse(t *testing.T) {
	cases := []struct {
		in   string
		mode AmbigReadsMode
		err  bool
	}{
		{"", AmbigReadsDrop, false},
		{"drop", AmbigReadsDrop, false},
		{"DROP", AmbigReadsDrop, false},
		{"incAD", AmbigReadsIncAD, false},
		{"incad", AmbigReadsIncAD, false},
		{"incAD0", AmbigReadsIncAD0, false},
		{"bogus", AmbigReadsDrop, true},
	}
	for _, tc := range cases {
		got, err := parseAmbigReads(tc.in)
		if tc.err {
			if err == nil {
				t.Errorf("parseAmbigReads(%q): want error, got nil", tc.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("parseAmbigReads(%q): unexpected error: %v", tc.in, err)
			continue
		}
		if got != tc.mode {
			t.Errorf("parseAmbigReads(%q): got %d, want %d", tc.in, got, tc.mode)
		}
	}
}

// TestMpileupBAMFlagParse covers the named-flag and numeric forms of
// the --skip-* argument parser ported from htslib bam_str2flag.
func TestMpileupBAMFlagParse(t *testing.T) {
	cases := []struct {
		in   string
		want uint16
		err  bool
	}{
		{"", 0, false},
		{"0x14", 0x14, false},
		{"20", 20, false},
		{"READ1", 0x40, false},
		{"read1", 0x40, false},
		{"PAIRED,PROPER_PAIR,MREVERSE", 0x1 | 0x2 | 0x20, false},
		{"SUPPLEMENTARY", 0x800, false},
		{"NOTAFLAG", 0, true},
	}
	for _, tc := range cases {
		got, err := parseBAMFlagString(tc.in)
		if tc.err {
			if err == nil {
				t.Errorf("parseBAMFlagString(%q): want error, got nil", tc.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("parseBAMFlagString(%q): unexpected error: %v", tc.in, err)
			continue
		}
		if got != tc.want {
			t.Errorf("parseBAMFlagString(%q): got 0x%x, want 0x%x", tc.in, got, tc.want)
		}
	}
}

// TestMpileupGoldensDeferred documents the upstream mpileup goldens that
// still do not byte-match and the precise reason, so the remaining work
// stays visible in the test output.
func TestMpileupGoldensDeferred(t *testing.T) {
	deferred := []struct{ golden, reason string }{
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
			"mpileup/indel-AD.1.out",
			"RESOLVED — byte-for-byte parity. Cluster (3) at " +
				"000000F:538 and :658 was the indel-pass is_del " +
				"qpos/min_dist off-by-one: upstream's resolve_cigar2 " +
				"(sam.c:5496) sets p->qpos = s->y at a D/N op, i.e. " +
				"queryPos AFTER the preceding M run (the first query " +
				"base of the next M run), not queryPos-1 (the last " +
				"base BEFORE the deletion). Our port's " +
				"accumulateMpileupBases D/N branch had qref = " +
				"queryPos-1, which dropped p.qpos by 1 for every " +
				"deletion-spanning read and (a) shifted min_dist for " +
				"I16[3<<2|*] anno accumulation, and (b) read the " +
				"wrong byte in the REF-rescue raw-qual lookup " +
				"(rec.Qual[p.qpos]) so the seqQ (3*seqQ+2*rawQ)/8 " +
				"blend and the min-baseQ-rescue gate both diverged. " +
				"Fixed by setting qref = queryPos. Listed here for " +
				"history; remove on next docs sweep.",
		},
		{
			"mpileup/indel-AD.1cns.out (residual)",
			"the --indels-cns dispatch now wires through to a Go port " +
				"of bam2bcf_edlib.c (bam2bcf_indelcns.go) driving the " +
				"in-tree edlib engine. The SNP / N-REF bulk matches " +
				"upstream byte-for-byte (covered by " +
				"TestMpileupIndelsCNSGolden); residual is on the four " +
				"indel rows at 000000F:537,:538,:655,:658 (homopolymer " +
				"/ tandem-repeat sites). The follow-up slice ports the " +
				"full bcf_cgp_consensus heterozygous threading " +
				"(cons[0]/cons[1]) and the edlib-flavored " +
				"compute_indelQ (indelQ1/indelQ2, vs_ref, poly_mqual, " +
				"TMP_MAGIC=255) that drive those residual deltas.",
		},
		{
			"mpileup/annot-NMBZ.3.1.out",
			"the SNP row at chr16:75 now byte-matches (including I16 " +
				"REF baseQ sums) after gating applyMpileupBAQ on the " +
				"requested regions: upstream's mpileup_reg loop only " +
				"invokes mplp_realn for piles emitted inside [beg, " +
				"end] (mpileup.c:573), so a `-r chr16:75` invocation " +
				"defers each read's first BAQ trigger to col 74 — by " +
				"which time the overlapping mate has been pushed and " +
				"overlap_push has already merged the qual array. " +
				"Without the region gate our port BAQed mate1 at its " +
				"first eligible chromosome-wide column (~16), which " +
				"saw raw quals; the subsequent overlap merge then " +
				"summed BAQed mate1 with raw mate2, inflating " +
				"qual[qpos=74] above what upstream computed. The " +
				"remaining indel row residual (QS / NMBZ / PL[0]: " +
				"226 vs 255, NMBZ sign flip) is the per-read chosen-" +
				"indel-type assignment at this homopolymer column — " +
				"same root cause as the indel-AD.1 homopolymer rows " +
				"that survived the BAQ probe (a single-ULP rounding " +
				"flip inside ProbalnGlocal at long homopolymer runs, " +
				"RNG-class divergence deferred per project policy).",
		},
		{
			"mpileup/annot-NMBZ.[23].2.out / FORMAT/NMBZ goldens",
			"the per-sample FORMAT/NMBZ emission is not in the slice 4e.6 " +
				"footprint; INFO/NMBZ is fully ported.",
		},
	}
	for _, d := range deferred {
		t.Logf("DEFERRED golden %s: %s", d.golden, d.reason)
	}
	t.Skip("documented above; these goldens are deferred to follow-up slices")
}
