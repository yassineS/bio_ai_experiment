package tabix

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

// upstreamHtslibTabix builds (once) and returns absolute paths to the upstream
// htslib `tabix` and `bgzip` binaries from the vendored submodule under
// reference_code/htslib. It is the on-demand, memoised builder used by the
// live parity tests in this file. The build is performed at most once per test
// binary invocation via sync.Once; subsequent calls reuse the result.
//
// The helper deliberately t.Fatalf's (never t.Skip) when the submodule is
// present but the binaries cannot be produced, so a broken upstream build is
// surfaced as a test failure rather than silently skipped.
var (
	upstreamOnce   sync.Once
	upstreamTabix  string
	upstreamBgzip  string
	upstreamErr    error
	upstreamReason string
)

func upstreamHtslibTabix(t *testing.T) (tabixBin, bgzipBin string) {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamOnce.Do(func() {
		root, err := repoRoot()
		if err != nil {
			upstreamErr = err
			return
		}
		htslib := filepath.Join(root, "reference_code", "htslib")
		if _, statErr := os.Stat(filepath.Join(htslib, "tabix.c")); statErr != nil {
			upstreamReason = "reference_code/htslib submodule is not checked out " +
				"(run: git submodule update --init --recursive reference_code/htslib)"
			upstreamErr = statErr
			return
		}
		tabixBin = filepath.Join(htslib, "tabix")
		bgzipBin = filepath.Join(htslib, "bgzip")
		// Build the binaries if either is missing.
		if !isExecutable(tabixBin) || !isExecutable(bgzipBin) {
			if err := buildUpstream(htslib); err != nil {
				upstreamErr = err
				return
			}
		}
		if !isExecutable(tabixBin) || !isExecutable(bgzipBin) {
			upstreamReason = "upstream build did not produce tabix/bgzip binaries"
			upstreamErr = os.ErrNotExist
			return
		}
		upstreamTabix = tabixBin
		upstreamBgzip = bgzipBin
	})
	if upstreamErr != nil {
		t.Fatalf("upstream htslib tabix unavailable: %v (%s)", upstreamErr, upstreamReason)
	}
	return upstreamTabix, upstreamBgzip
}

// repoRoot walks up from this source file's directory to the module root
// (the directory containing go.mod).
func repoRoot() (string, error) {
	_, file, _, ok := runtime.Caller(0)
	if !ok {
		return "", os.ErrNotExist
	}
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", os.ErrNotExist
		}
		dir = parent
	}
}

func isExecutable(path string) bool {
	info, err := os.Stat(path)
	if err != nil {
		return false
	}
	return !info.IsDir() && info.Mode()&0o111 != 0
}

// buildUpstream runs autoreconf/configure/make to produce the tabix and bgzip
// binaries in htslib. It tolerates a pre-configured tree (skipping
// autoreconf/configure when configure has already run).
func buildUpstream(htslib string) error {
	runStep := func(name string, args ...string) error {
		cmd := exec.Command(name, args...)
		cmd.Dir = htslib
		out, err := cmd.CombinedOutput()
		if err != nil {
			return &buildError{step: name, output: string(out), err: err}
		}
		return nil
	}
	if _, err := os.Stat(filepath.Join(htslib, "configure")); err != nil {
		if err := runStep("autoreconf", "-i"); err != nil {
			return err
		}
	}
	if _, err := os.Stat(filepath.Join(htslib, "config.mk")); err != nil {
		if err := runStep("./configure"); err != nil {
			return err
		}
	}
	return runStep("make", "tabix", "bgzip")
}

type buildError struct {
	step   string
	output string
	err    error
}

func (e *buildError) Error() string {
	out := e.output
	if len(out) > 2000 {
		out = out[len(out)-2000:]
	}
	return "upstream build step " + e.step + " failed: " + e.err.Error() + "\n" + out
}

// bgzipCompress runs upstream bgzip to compress text into a fresh .gz file and
// returns its path.
func bgzipCompress(t *testing.T, bgzipBin, dir, name, text string) string {
	t.Helper()
	plain := filepath.Join(dir, strings.TrimSuffix(name, ".gz"))
	if err := os.WriteFile(plain, []byte(text), 0o644); err != nil {
		t.Fatalf("write %s: %v", plain, err)
	}
	gz := filepath.Join(dir, name)
	out, err := exec.Command(bgzipBin, "-c", plain).Output()
	if err != nil {
		t.Fatalf("bgzip compress: %v", err)
	}
	if err := os.WriteFile(gz, out, 0o644); err != nil {
		t.Fatalf("write %s: %v", gz, err)
	}
	return gz
}

// runUpstreamTabix runs the upstream tabix binary and returns its stdout.
func runUpstreamTabix(t *testing.T, tabixBin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(tabixBin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream tabix %v failed: %v\nstderr: %s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// TestTabix_TargetsStrictUpstreamParity compares the Go -T/--targets strict
// post-filter against upstream tabix -T byte-for-byte, across tab and BED
// coordinate conventions and boundary cases.
func TestTabix_TargetsStrictUpstreamParity(t *testing.T) {
	tabixBin, bgzipBin := upstreamHtslibTabix(t)
	dir := t.TempDir()

	const vcf = "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tT\t.\t.\t.\n" +
		"chr1\t150\t.\tC\tG\t.\t.\t.\n" +
		"chr1\t200\t.\tG\tA\t.\t.\t.\n" +
		"chr1\t250\t.\tT\tC\t.\t.\t.\n" +
		"chr2\t300\t.\tA\tG\t.\t.\t.\n"
	gz := bgzipCompress(t, bgzipBin, dir, "in.vcf.gz", vcf)
	// Build the .tbi with upstream so both implementations read the same index.
	runUpstreamTabix(t, tabixBin, "-p", "vcf", gz)

	cfg, err := PresetConfig(PresetVCF)
	if err != nil {
		t.Fatalf("preset: %v", err)
	}
	idx, err := ReadFile(gz + ".tbi")
	if err != nil {
		t.Fatalf("read index: %v", err)
	}

	cases := []struct {
		name    string
		tname   string // targets file name (suffix drives tab vs bed parsing)
		content string
		region  string // optional positional region; empty means whole-file
	}{
		{"tab_point", "t.txt", "chr1\t150\t150\n", ""},
		{"tab_range", "t.txt", "chr1\t100\t150\n", ""},
		{"tab_adjacent", "t.txt", "chr1\t100\t150\nchr1\t200\t250\n", ""},
		{"tab_gap_no_record", "t.txt", "chr1\t160\t190\n", ""},
		{"bed_range", "t.bed", "chr1\t149\t200\n", ""},
		{"chrom_only", "t.txt", "chr2\n", ""},
		{"with_region", "t.txt", "chr1\t140\t210\n", "chr1:100-160"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			tpath := filepath.Join(dir, c.name+"_"+c.tname)
			if err := os.WriteFile(tpath, []byte(c.content), 0o644); err != nil {
				t.Fatalf("write targets: %v", err)
			}
			upArgs := []string{"-T", tpath, gz}
			if c.region != "" {
				upArgs = append(upArgs, c.region)
			}
			want := runUpstreamTabix(t, tabixBin, upArgs...)

			got := goTargetsQuery(t, idx, cfg, gz, tpath, c.region)
			if !bytes.Equal(got, want) {
				t.Errorf("targets parity mismatch for %s\n go: %q\nwant: %q", c.name, got, want)
			}
		})
	}
}

// goTargetsQuery reproduces the CLI's -T behavior in-process: load the targets
// filter, query the region(s), and emit overlapping records as newline-
// terminated lines.
func goTargetsQuery(t *testing.T, idx *Index, cfg Config, gz, targetsPath, region string) []byte {
	t.Helper()
	tg, err := LoadTargets(targetsPath)
	if err != nil {
		t.Fatalf("LoadTargets: %v", err)
	}
	type qreg struct {
		chrom    string
		beg, end int
	}
	var regs []qreg
	if region == "" {
		for _, chrom := range idx.Chroms() {
			regs = append(regs, qreg{chrom, 0, 1 << 30})
		}
	} else {
		chrom, beg, end := parseRegionTest(t, region)
		regs = append(regs, qreg{chrom, beg, end})
	}
	var out bytes.Buffer
	for _, r := range regs {
		recs, err := idx.QueryRecords(gz, r.chrom, r.beg, r.end)
		if err != nil {
			t.Fatalf("QueryRecords: %v", err)
		}
		for _, rec := range recs {
			if !tg.Overlaps(r.chrom, rec.Beg, rec.End) {
				continue
			}
			out.Write(rec.Line)
			out.WriteByte('\n')
		}
	}
	return out.Bytes()
}

// parseRegionTest parses a CHROM:START-END region (1-based inclusive) into
// 0-based half-open coordinates for tests.
func parseRegionTest(t *testing.T, s string) (string, int, int) {
	t.Helper()
	colon := strings.IndexByte(s, ':')
	if colon < 0 {
		return s, 0, 1 << 30
	}
	chrom := s[:colon]
	rest := s[colon+1:]
	dash := strings.IndexByte(rest, '-')
	if dash < 0 {
		t.Fatalf("region without end: %q", s)
	}
	beg := atoiTest(t, rest[:dash])
	end := atoiTest(t, rest[dash+1:])
	return chrom, beg - 1, end
}

func atoiTest(t *testing.T, s string) int {
	t.Helper()
	n := 0
	for i := 0; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			t.Fatalf("bad number %q", s)
		}
		n = n*10 + int(s[i]-'0')
	}
	return n
}

// TestTabix_ReheaderUpstreamParity validates the Go --reheader output against
// upstream tabix as a live consumer. Upstream tabix's own reheader producer is
// broken at the vendored commit (it corrupts the body for typical inputs; see
// PARITY_ROADMAP.md), so it cannot serve as a byte-for-byte producer oracle.
// Instead this test confirms that the Go reheader output is a valid BGZF stream
// that upstream tabix can re-index and query, returning exactly the original
// records, and that upstream sees the replacement header verbatim.
func TestTabix_ReheaderUpstreamParity(t *testing.T) {
	tabixBin, bgzipBin := upstreamHtslibTabix(t)
	dir := t.TempDir()

	const body = "chr1\t100\t.\tA\tT\t.\t.\t.\n" +
		"chr1\t200\t.\tC\tG\t.\t.\t.\n" +
		"chr2\t150\t.\tT\tC\t.\t.\t.\n"
	const oldVCF = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + body
	gz := bgzipCompress(t, bgzipBin, dir, "in.vcf.gz", oldVCF)
	runUpstreamTabix(t, tabixBin, "-p", "vcf", gz)

	const newHdr = "##fileformat=VCFv4.2\n##NEWHDR=yes\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	hdrPath := filepath.Join(dir, "newhdr.txt")
	if err := os.WriteFile(hdrPath, []byte(newHdr), 0o644); err != nil {
		t.Fatalf("write header: %v", err)
	}

	// Produce the reheadered file with the Go implementation.
	var buf bytes.Buffer
	if err := Reheader(gz, hdrPath, '#', &buf); err != nil {
		t.Fatalf("Reheader: %v", err)
	}
	outGz := filepath.Join(dir, "out.vcf.gz")
	if err := os.WriteFile(outGz, buf.Bytes(), 0o644); err != nil {
		t.Fatalf("write reheadered: %v", err)
	}

	// Upstream must be able to index and query it.
	runUpstreamTabix(t, tabixBin, "-p", "vcf", outGz)

	// Upstream query of chr1 returns exactly the chr1 records.
	gotChr1 := runUpstreamTabix(t, tabixBin, outGz, "chr1")
	wantChr1 := "chr1\t100\t.\tA\tT\t.\t.\t.\nchr1\t200\t.\tC\tG\t.\t.\t.\n"
	if string(gotChr1) != wantChr1 {
		t.Errorf("upstream query of reheadered file:\n got %q\nwant %q", gotChr1, wantChr1)
	}

	// Upstream header view returns the replacement header verbatim.
	gotHdr := runUpstreamTabix(t, tabixBin, "-H", "-p", "vcf", outGz, "chr2")
	if !bytes.HasPrefix(gotHdr, []byte(newHdr)) {
		t.Errorf("upstream -H did not show the replacement header:\n%q", gotHdr)
	}

	// The full decompressed content must be exactly new header + original body.
	full := decodeBGZF(t, buf.Bytes())
	if full != newHdr+body {
		t.Errorf("decompressed reheader content:\n got %q\nwant %q", full, newHdr+body)
	}
}
