package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// mpileupGoldenDir is where the upstream `bcftools mpileup` fixtures are
// vendored: the mpileup.{1,2,3}.bam read sets, the mpileup.ref.fa reference
// (with its .fai), and the various BAM/FASTA inputs. Expected VCF output is
// produced live by the upstream `bcftools mpileup` binary at test time, not
// from committed golden files.
const mpileupGoldenDir = "../../testdata/mpileup"

// mpileupFixture returns the absolute path to a vendored mpileup fixture.
// A missing fixture is fatal — these inputs are committed.
func mpileupFixture(t *testing.T, name string) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join(mpileupGoldenDir, name))
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	if _, err := os.Stat(abs); err != nil {
		t.Fatalf("vendored fixture %s missing: %v", name, err)
	}
	return abs
}

// indexedBAMCopy copies a BAM into a temp dir and runs `samtools index` so
// the upstream binary can honour region (`-r`) queries against fixtures
// that ship without a .bai/.csi sidecar. It returns the temp-dir BAM path.
// The samtools binary is the one built next to the vendored htslib (under
// reference_code/samtools); it is built on demand if absent.
func indexedBAMCopy(t *testing.T, bcftoolsBin, bam string) string {
	t.Helper()
	dir := t.TempDir()
	dst := filepath.Join(dir, filepath.Base(bam))
	if err := copyFile(bam, dst); err != nil {
		t.Fatalf("copy %s: %v", bam, err)
	}
	samtools := upstreamSamtoolsForBcftools(t, bcftoolsBin)
	cmd := exec.Command(samtools, "index", dst)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools index %s: %v\n%s", dst, err, out)
	}
	return dst
}

// upstreamSamtoolsForBcftools returns the path to a built `samtools` binary
// under reference_code/samtools, building it (plain make, after htslib) on
// demand. It is fatal — never skipped — if the build fails.
func upstreamSamtoolsForBcftools(t *testing.T, bcftoolsBin string) string {
	t.Helper()
	referenceCode := filepath.Dir(filepath.Dir(bcftoolsBin)) // .../reference_code
	root := filepath.Dir(referenceCode)
	samtoolsDir := filepath.Join(referenceCode, "samtools")
	bin := filepath.Join(samtoolsDir, "samtools")
	if fileExists(bin) {
		return bin
	}
	// Serialise with any concurrent samtools-package build into the shared
	// reference_code tree.
	unlock := lockBuild(root)
	defer unlock()
	if fileExists(bin) {
		return bin
	}
	if err := ensureHtslibBuilt(root, filepath.Join(referenceCode, "htslib")); err != nil {
		t.Fatalf("build htslib: %v", err)
	}
	if !fileExists(filepath.Join(samtoolsDir, "bamtk.c")) {
		if out, err := run(root, "git", "submodule", "update", "--init", "--recursive", "reference_code/samtools"); err != nil {
			t.Fatalf("git submodule samtools: %v\n%s", err, out)
		}
	}
	if out, err := run(samtoolsDir, "make", "-j4"); err != nil {
		t.Fatalf("make samtools: %v\n%s", err, out)
	}
	if !fileExists(bin) {
		t.Fatalf("samtools binary not produced at %s", bin)
	}
	return bin
}

// splitMpileupVCF splits VCF text into its header lines (the leading
// "##"/"#" block) and its data records, dropping the run-specific
// `##bcftools_*` and `##reference` header lines exactly as upstream's
// test.pl does (`grep -v ^##bcftools | grep -v ^##reference`).
func splitMpileupVCF(s string) (header, data []string) {
	for _, ln := range strings.Split(s, "\n") {
		if ln == "" {
			continue
		}
		if strings.HasPrefix(ln, "##bcftools") || strings.HasPrefix(ln, "##reference") {
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

// mpileupCase describes one live parity invocation: the same options drive
// both the Go MpileupFile port and the upstream `bcftools mpileup` binary.
type mpileupCase struct {
	name     string
	inputs   []string // BAM fixture basenames
	ref      string   // FASTA fixture basename
	regions  []string
	targets  []string
	annotate string
	flagIncl string // --skip-all-unset
	flagAny  string // --skip-any-unset
	ambig    AmbigReadsMode
}

// upstreamMpileupVCF runs the live `bcftools mpileup` binary for a case and
// returns its VCF output (with run-specific header lines stripped).
func upstreamMpileupVCF(t *testing.T, bin string, c mpileupCase) string {
	t.Helper()
	args := []string{"mpileup", "--no-version", "-f", mpileupFixture(t, c.ref)}
	if c.annotate != "" {
		args = append(args, "-a", c.annotate)
	}
	for _, r := range c.regions {
		args = append(args, "-r", r)
	}
	for _, tgt := range c.targets {
		args = append(args, "-t", tgt)
	}
	if c.flagIncl != "" {
		args = append(args, "--skip-all-unset", c.flagIncl)
	}
	if c.flagAny != "" {
		args = append(args, "--skip-any-unset", c.flagAny)
	}
	switch c.ambig {
	case AmbigReadsIncAD:
		args = append(args, "--ambig-reads", "incAD")
	case AmbigReadsIncAD0:
		args = append(args, "--ambig-reads", "incAD0")
	}
	args = append(args, "-Ov")
	for _, in := range c.inputs {
		bam := mpileupFixture(t, in)
		// Region (`-r`) queries require a BAM index. Most fixtures ship a
		// .bai sidecar; for those that do not, build an indexed copy in a
		// temp dir so the committed fixture tree stays index-free.
		if len(c.regions) > 0 && !fileExists(bam+".bai") && !fileExists(bam+".csi") {
			bam = indexedBAMCopy(t, bin, bam)
		}
		args = append(args, bam)
	}

	cmd := exec.Command(bin, args...)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools mpileup %v: %v\n%s", args, err, errBuf.String())
	}
	return out.String()
}

// goMpileupVCF runs the Go MpileupFile port for a case.
func goMpileupVCF(t *testing.T, c mpileupCase) string {
	t.Helper()
	inputs := make([]string, len(c.inputs))
	for i, in := range c.inputs {
		inputs[i] = mpileupFixture(t, in)
	}
	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:         inputs,
		FastaRef:       mpileupFixture(t, c.ref),
		Regions:        c.regions,
		Targets:        c.targets,
		Annotate:       c.annotate,
		FlagIncl:       c.flagIncl,
		FlagAny:        c.flagAny,
		AmbigReadsMode: c.ambig,
		NoVersion:      true,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}
	return buf.String()
}

// compareMpileupVCF compares the Go and upstream VCF outputs: the header
// must match line-for-line, and every data record must match byte-for-byte
// when indexed by CHROM:POS:{snp,indel} (so the SNP and INDEL rows at the
// same coordinate are disambiguated).
func compareMpileupVCF(t *testing.T, gotVCF, wantVCF string) {
	t.Helper()
	gotHeader, gotData := splitMpileupVCF(gotVCF)
	wantHeader, wantData := splitMpileupVCF(wantVCF)

	if len(gotHeader) != len(wantHeader) {
		t.Errorf("header line count: got %d, want %d\n got:  %v\n want: %v",
			len(gotHeader), len(wantHeader), gotHeader, wantHeader)
	}
	n := len(gotHeader)
	if len(wantHeader) < n {
		n = len(wantHeader)
	}
	for i := 0; i < n; i++ {
		if gotHeader[i] != wantHeader[i] {
			t.Errorf("header line %d:\n got:  %s\n want: %s", i, gotHeader[i], wantHeader[i])
		}
	}

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

	diffs := 0
	for _, ln := range gotData {
		key := recKey(ln)
		want, ok := wantByPos[key]
		if !ok {
			t.Errorf("emitted an extra record at %s:\n %s", key, ln)
			continue
		}
		delete(wantByPos, key)
		if ln != want {
			diffs++
			if diffs <= 10 {
				t.Errorf("record %s mismatch:\n got:  %s\n want: %s", key, ln, want)
			}
		}
	}
	for key := range wantByPos {
		t.Errorf("missing record at %s (upstream has it, we did not emit it)", key)
	}
	if diffs > 10 {
		t.Errorf("... and %d more record mismatches", diffs-10)
	}
}

// TestMpileup_UpstreamParity runs BOTH the live upstream `bcftools mpileup`
// binary and the Go MpileupFile port on the same fixtures and asserts the
// VCF output is identical (header line-for-line, data records byte-for-byte
// keyed by position+type). This replaces the former golden-file comparison;
// the upstream binary is built on demand and a build failure is fatal,
// never skipped.
//
// The invocations mirror reference_code/bcftools/test/test.pl:
//
//   - single-bam-full-contig: `mpileup -a -AD mpileup.3.bam`
//   - multi-bam-region:       `mpileup -r17:100-150 -a -AD mpileup.{1,2,3}.bam`
//   - SCR:                    `mpileup -a -AD,INFO/SCR,FMT/SCR mpileup-SCR.bam`
//   - filter:                 `mpileup --skip-{all,any}-unset READ1 ...`
//   - indel-AD:               `mpileup -a AD -r 11:75 [--ambig-reads MODE]`
//   - NMBZ:                    `mpileup -a -AD,INFO/NMBZ ...`
func TestMpileup_UpstreamParity(t *testing.T) {
	bin := upstreamBcftools(t)

	cases := []mpileupCase{
		{
			name:     "single-bam-full-contig",
			inputs:   []string{"mpileup.3.bam"},
			ref:      "mpileup.ref.fa",
			annotate: "-AD",
		},
		{
			name:     "multi-bam-region",
			inputs:   []string{"mpileup.1.bam", "mpileup.2.bam", "mpileup.3.bam"},
			ref:      "mpileup.ref.fa",
			regions:  []string{"17:100-150"},
			annotate: "-AD",
		},
		{
			name:     "scr",
			inputs:   []string{"mpileup-SCR.bam"},
			ref:      "mpileup-SCR.fa",
			annotate: "-AD,INFO/SCR,FMT/SCR",
		},
		{
			name:     "indel-AD-default",
			inputs:   []string{"indel-AD.2.bam"},
			ref:      "indel-AD.2.fa",
			regions:  []string{"11:75"},
			annotate: "AD",
			ambig:    AmbigReadsDrop,
		},
		{
			name:     "indel-AD-incAD",
			inputs:   []string{"indel-AD.2.bam"},
			ref:      "indel-AD.2.fa",
			regions:  []string{"11:75"},
			annotate: "AD",
			ambig:    AmbigReadsIncAD,
		},
		{
			name:     "indel-AD-incAD0",
			inputs:   []string{"indel-AD.2.bam"},
			ref:      "indel-AD.2.fa",
			regions:  []string{"11:75"},
			annotate: "AD",
			ambig:    AmbigReadsIncAD0,
		},
		{
			name:     "nmbz-1",
			inputs:   []string{"annot-NMBZ.1.bam"},
			ref:      "annot-NMBZ.1.fa",
			regions:  []string{"chr19:69-99"},
			annotate: "-AD,INFO/NMBZ",
		},
		{
			name:     "nmbz-2-depthcap",
			inputs:   []string{"annot-NMBZ.2.bam"},
			ref:      "annot-NMBZ.2.fa",
			regions:  []string{"chr6:75"},
			annotate: "-AD,INFO/NMBZ",
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// Ensure the .fai sidecar the FASTA reader needs is present.
			mpileupFixture(t, c.ref+".fai")
			want := upstreamMpileupVCF(t, bin, c)
			got := goMpileupVCF(t, c)
			compareMpileupVCF(t, got, want)
		})
	}
}

// TestMpileupSCROnIndelRow confirms that when -a INFO/SCR,FMT/SCR is in
// force, bcfCall2bcfIndel emits SCR on the indel row using the shared
// per-column tally that also feeds the SNP row. The SNP and indel rows at
// the same position must report the same INFO/SCR and FORMAT/SCR. This is a
// binary-free invariant check on the Go port's own output.
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
	if indelInfo == "0" {
		t.Errorf("expected non-zero SCR at 11:75 indel row, got %s", indelInfo)
	}
}

// TestMpileupFilterParity replays the upstream `bcftools mpileup --skip-*`
// flag filter live (test.pl lines 1072-1073) for both the all-unset and
// any-unset spellings, comparing Go output against the upstream binary.
func TestMpileupFilterParity(t *testing.T) {
	bin := upstreamBcftools(t)
	base := mpileupCase{
		inputs:   []string{"mpileup-filter.bam"},
		ref:      "mpileup-SCR.fa",
		targets:  []string{"1:100"},
		annotate: "-AD",
	}
	cases := []mpileupCase{
		func() mpileupCase { c := base; c.name = "skip-all-unset-READ1"; c.flagIncl = "READ1"; return c }(),
		func() mpileupCase { c := base; c.name = "skip-any-unset-READ1"; c.flagAny = "READ1"; return c }(),
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			mpileupFixture(t, c.ref+".fai")
			want := upstreamMpileupVCF(t, bin, c)
			got := goMpileupVCF(t, c)
			compareMpileupVCF(t, got, want)
		})
	}
}

// TestMpileupAmbigReadsParse verifies parseAmbigReads round-trips every
// upstream-accepted spelling.
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

// TestMpileupBAMFlagParse covers the named-flag and numeric forms of the
// --skip-* argument parser ported from htslib bam_str2flag.
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
