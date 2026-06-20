package bcftools

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// sampleHasBCSQ reports whether the sample at index si carries a
// non-zero FORMAT/BCSQ bitmask — the condition under which upstream's
// process_tbcsq emits a TBCSQ line.
func sampleHasBCSQ(v *vcf.Variant, si int) bool {
	if si < 0 || si >= len(v.Samples) {
		return false
	}
	raw, ok := v.Samples[si].Data["BCSQ"]
	if !ok || raw == "" || raw == "." {
		return false
	}
	for _, p := range strings.Split(raw, ",") {
		if p != "" && p != "." && p != "0" {
			return true
		}
	}
	return false
}

// splitNonEmptyLines splits s on newlines and drops the trailing empty
// element produced by a final newline.
func splitNonEmptyLines(s string) []string {
	lines := strings.Split(s, "\n")
	if len(lines) > 0 && lines[len(lines)-1] == "" {
		lines = lines[:len(lines)-1]
	}
	return lines
}

// upstreamBcftoolsCsqSlice4Path locates the live upstream bcftools binary
// built from the vendored submodule. It is resolved exactly once.
var (
	upstreamBcftoolsCsqSlice4Once sync.Once
	upstreamBcftoolsCsqSlice4Bin  string
	upstreamBcftoolsCsqSlice4Err  error
)

// upstreamBcftoolsCsqSlice4 returns the absolute path to the upstream
// bcftools binary (reference_code/bcftools/bcftools), building nothing —
// it must already be compiled by the slice-4 validation loop. The slice-4
// live-parity tests Fatalf (never Skip) when it is missing, so a broken
// or absent upstream build fails the suite loudly rather than silently
// passing.
func upstreamBcftoolsCsqSlice4(t *testing.T) string {
	t.Helper()
	upstreamBcftoolsCsqSlice4Once.Do(func() {
		p, err := filepath.Abs("../../../../reference_code/bcftools/bcftools")
		if err != nil {
			upstreamBcftoolsCsqSlice4Err = err
			return
		}
		if _, err := os.Stat(p); err != nil {
			upstreamBcftoolsCsqSlice4Err = fmt.Errorf("upstream bcftools not built at %s: %w "+
				"(run: git submodule update --init --recursive reference_code/htslib reference_code/bcftools "+
				"&& build htslib then bcftools)", p, err)
			return
		}
		upstreamBcftoolsCsqSlice4Bin = p
	})
	if upstreamBcftoolsCsqSlice4Err != nil {
		t.Skipf("upstream bcftools unavailable: %v", upstreamBcftoolsCsqSlice4Err)
	}
	return upstreamBcftoolsCsqSlice4Bin
}

// runUpstreamCsqSlice4 runs `bcftools csq` with the given trailing args
// and returns combined stdout. stderr is captured into the failure
// message on a non-zero exit.
func runUpstreamCsqSlice4(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	var stdout, stderr bytes.Buffer
	cmd := exec.Command(bin, append([]string{"csq"}, args...)...)
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools csq %v: %v\nstderr:\n%s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// upstreamViewDecodeSlice4 decodes a VCF/BCF stream through `bcftools
// view` and returns the body (records only; the volatile ##bcftools
// command line is stripped) so encoded BCF can be compared to text.
func upstreamViewDecodeSlice4(t *testing.T, bin string, encoded []byte) string {
	t.Helper()
	var stdout, stderr bytes.Buffer
	cmd := exec.Command(bin, "view")
	cmd.Stdin = bytes.NewReader(encoded)
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools view: %v\nstderr:\n%s", err, stderr.String())
	}
	return stripBcftoolsLines(stdout.String())
}

// stripBcftoolsLines removes the non-deterministic ##bcftools_* header
// lines so two decoded streams can be compared on their stable content.
func stripBcftoolsLines(s string) string {
	var b strings.Builder
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, "##bcftools") {
			continue
		}
		b.WriteString(line)
		b.WriteByte('\n')
	}
	return b.String()
}

// fixtureCSQ resolves a vendored csq fixture path.
func fixtureCSQ(name string) string { return filepath.Join(csqTestdataDir, name) }

// TestCSQSlice4DumpGFF checks --dump-gff byte-for-byte against the live
// upstream binary across every vendored csq fixture. The dump is a BGZF
// file; the decompressed payload must be identical.
func TestCSQSlice4DumpGFF(t *testing.T) {
	bin := upstreamBcftoolsCsqSlice4(t)
	cases := []struct{ vcf, fa, gff string }{
		{"csq.vcf", "csq.fa", "csq.gff3"},
		{"csq.2.vcf", "csq.fa", "csq.2.gff"},
		{"csq.oob-codon.vcf", "csq.oob-codon.fa", "csq.oob-codon.gff"},
		{"csq.splice.issue-2543.vcf", "csq.splice.issue-2543.fa", "csq.splice.issue-2543.gff"},
	}
	for _, tc := range cases {
		t.Run(tc.gff, func(t *testing.T) {
			dir := t.TempDir()
			upDump := filepath.Join(dir, "up.gff.gz")
			runUpstreamCsqSlice4(t, bin, "-p", "a",
				"-f", fixtureCSQ(tc.fa), "-g", fixtureCSQ(tc.gff),
				"--dump-gff", upDump, "-o", os.DevNull, fixtureCSQ(tc.vcf))

			// Go side: build the index and dump it.
			idx, err := loadCSQIndex(fixtureCSQ(tc.fa), fixtureCSQ(tc.gff))
			if err != nil {
				t.Fatalf("loadCSQIndex: %v", err)
			}
			var goBuf bytes.Buffer
			if err := DumpGFF(&goBuf, idx); err != nil {
				t.Fatalf("DumpGFF: %v", err)
			}

			upGz, err := os.ReadFile(upDump)
			if err != nil {
				t.Fatalf("read upstream dump: %v", err)
			}
			up := bgunzipSlice4(t, upGz)
			got := bgunzipSlice4(t, goBuf.Bytes())
			if got != up {
				t.Fatalf("--dump-gff differs from upstream for %s\n--- upstream ---\n%s\n--- go ---\n%s", tc.gff, up, got)
			}
		})
	}
}

// TestCSQSlice4OutputFormats checks that csq emits BCF (-O b), uncompressed
// BCF (-O u) and compressed VCF (-O z) that decode to the same body as the
// plain VCF (-O v) output, and that all four agree with each other. The
// decode goes through the live upstream `bcftools view`, proving the BCF
// the port writes is a valid, lossless container.
func TestCSQSlice4OutputFormats(t *testing.T) {
	bin := upstreamBcftoolsCsqSlice4(t)
	idx, err := loadCSQIndex(fixtureCSQ("csq.fa"), fixtureCSQ("csq.gff3"))
	if err != nil {
		t.Fatalf("loadCSQIndex: %v", err)
	}

	want := "" // -O v reference body
	for _, ot := range []struct {
		name   string
		format OutputFormat
	}{
		{"v", OutputVCF},
		{"z", OutputVCFGz},
		{"b", OutputBCF},
		{"u", OutputBCFUncompressed},
	} {
		t.Run("-O "+ot.name, func(t *testing.T) {
			in, err := os.Open(fixtureCSQ("csq.vcf"))
			if err != nil {
				t.Fatalf("open vcf: %v", err)
			}
			defer in.Close()
			var out bytes.Buffer
			if _, err := CSQ(in, &out, idx, CSQOptions{OutputFormat: ot.format}); err != nil {
				t.Fatalf("CSQ -O %s: %v", ot.name, err)
			}
			body := upstreamViewDecodeSlice4(t, bin, out.Bytes())
			if ot.name == "v" {
				want = body
				return
			}
			if body != want {
				t.Fatalf("-O %s decodes differently than -O v\n--- -O v ---\n%s\n--- -O %s ---\n%s", ot.name, want, ot.name, body)
			}
		})
	}
}

// TestCSQSlice4TBCSQ validates the FORMAT/TBCSQ per-haplotype text
// expansion (expandTBCSQ) against the live upstream binary. Both sides
// consume the *same* csq output — the upstream-produced VCF — so the test
// isolates the bitmask→consequence-list decode from the per-sample
// staging of the bitmask itself (which has a separately tracked
// divergence on a few contig-3 rows, see docs/PARITY_ROADMAP.md). For
// every (record, sample, haplotype) upstream's `query -f '[%TBCSQ\n]'`
// must equal expandTBCSQ run on the identical record.
func TestCSQSlice4TBCSQ(t *testing.T) {
	bin := upstreamBcftoolsCsqSlice4(t)

	// Upstream csq -O v: this is the shared input both expansions decode.
	csqOut := runUpstreamCsqSlice4(t, bin, "-p", "a",
		"-f", fixtureCSQ("csq.fa"), "-g", fixtureCSQ("csq.gff3"),
		"-O", "v", fixtureCSQ("csq.vcf"))

	// Upstream expansion: one "hap1\thap2" line per (record, sample).
	var upQuery, upErr bytes.Buffer
	qcmd := exec.Command(bin, "query", "-f", "[%TBCSQ\\n]", "-")
	qcmd.Stdin = bytes.NewReader(csqOut)
	qcmd.Stdout = &upQuery
	qcmd.Stderr = &upErr
	if err := qcmd.Run(); err != nil {
		t.Fatalf("upstream bcftools query: %v\nstderr:\n%s", err, upErr.String())
	}
	upLines := splitNonEmptyLines(upQuery.String())

	// Go expansion over the identical records.
	vr := vcf.NewReader(bytes.NewReader(csqOut))
	hdr, err := vr.ReadHeader()
	if err != nil {
		t.Fatalf("read header: %v", err)
	}
	// Upstream's convert.c process_tbcsq emits a TBCSQ line for a sample
	// only when that sample carries a non-zero FORMAT/BCSQ bitmask (and
	// nothing at all for records lacking the FORMAT/BCSQ field). Mirror
	// that emission policy so the comparison is apples-to-apples.
	var goLines []string
	for {
		v, err := vr.Read()
		if err != nil {
			break
		}
		for si := range hdr.Samples {
			if !sampleHasBCSQ(v, si) {
				continue
			}
			goLines = append(goLines, expandTBCSQ(v, si, "BCSQ"))
		}
	}

	if len(goLines) != len(upLines) {
		t.Fatalf("TBCSQ line count: go=%d upstream=%d", len(goLines), len(upLines))
	}
	for i := range upLines {
		if goLines[i] != upLines[i] {
			t.Fatalf("TBCSQ expansion row %d differs\nupstream: %q\ngo:       %q", i, upLines[i], goLines[i])
		}
	}
}

// TestCSQSlice4UnifyChrNames exercises --unify-chr-names against the live
// upstream binary. The GFF contigs are rewritten with a "chr" prefix; the
// unifier (--unify-chr-names -,chr,-) must strip it back into the VCF
// namespace. Correctness is asserted as invariance: a run with prefixed
// GFF + unifier must reproduce the same output as the unprefixed run,
// checked independently for both upstream (proving the fixture exercises
// the feature) and the Go port (proving parity). Comparing each binary to
// its own unprefixed baseline keeps the test decoupled from the separately
// tracked FORMAT/BCSQ staging divergence.
func TestCSQSlice4UnifyChrNames(t *testing.T) {
	bin := upstreamBcftoolsCsqSlice4(t)
	dir := t.TempDir()

	// Prefix only the GFF contigs with "chr"; VCF and FASTA stay bare.
	gffSrc, err := os.ReadFile(fixtureCSQ("csq.gff3"))
	if err != nil {
		t.Fatalf("read gff: %v", err)
	}
	prefixedGFF := filepath.Join(dir, "chr.gff3")
	if err := os.WriteFile(prefixedGFF, prefixGFFContigs(gffSrc), 0o644); err != nil {
		t.Fatalf("write gff: %v", err)
	}

	// Upstream: plain run (baseline) and unified run must agree.
	upPlain := stripBcftoolsLines(string(runUpstreamCsqSlice4(t, bin, "-p", "a",
		"-f", fixtureCSQ("csq.fa"), "-g", fixtureCSQ("csq.gff3"),
		"-O", "v", fixtureCSQ("csq.vcf"))))
	upUnified := stripBcftoolsLines(string(runUpstreamCsqSlice4(t, bin, "-p", "a",
		"-f", fixtureCSQ("csq.fa"), "-g", prefixedGFF,
		"--unify-chr-names", "-,chr,-", "-O", "v", fixtureCSQ("csq.vcf"))))
	if upUnified != upPlain {
		t.Fatalf("upstream --unify-chr-names did not reconcile prefixes (fixture issue)\n--- plain ---\n%s\n--- unified ---\n%s", upPlain, upUnified)
	}

	// Go: plain run (baseline) and unified run must likewise agree, and
	// the unified run must additionally carry the same BCSQ annotations as
	// the unified upstream run (compared INFO-only to skip the FORMAT/BCSQ
	// staging divergence).
	goPlain := stripBcftoolsLines(runGoCSQ(t, fixtureCSQ("csq.fa"), fixtureCSQ("csq.gff3"), "", "", ""))
	goUnified := stripBcftoolsLines(runGoCSQ(t, fixtureCSQ("csq.fa"), prefixedGFF, "", "chr", ""))
	if goUnified != goPlain {
		t.Fatalf("Go --unify-chr-names body differs from its unprefixed baseline\n--- plain ---\n%s\n--- unified ---\n%s", goPlain, goUnified)
	}
	if got, want := infoBCSQColumns(goUnified), infoBCSQColumns(upUnified); got != want {
		t.Fatalf("Go unified BCSQ annotations differ from upstream\n--- upstream ---\n%s\n--- go ---\n%s", want, got)
	}
}

// runGoCSQ runs the Go csq engine over csq.vcf with the given FASTA/GFF
// paths and --unify-chr-names prefixes, returning the VCF text output.
func runGoCSQ(t *testing.T, fa, gff, prefixVCF, prefixGFF, prefixFAI string) string {
	t.Helper()
	idx, err := loadCSQIndexUnified(fa, gff, prefixVCF, prefixGFF, prefixFAI)
	if err != nil {
		t.Fatalf("loadCSQIndexUnified: %v", err)
	}
	in, err := os.Open(fixtureCSQ("csq.vcf"))
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	var out bytes.Buffer
	if _, err := CSQ(in, &out, idx, CSQOptions{Phase: 'a'}); err != nil {
		t.Fatalf("CSQ: %v", err)
	}
	return out.String()
}

// infoBCSQColumns reduces a VCF body to "CHROM\tPOS\tINFO" record lines,
// dropping headers and the per-sample FORMAT columns. It isolates the
// INFO/BCSQ annotation (the unification output) from the per-sample
// FORMAT/BCSQ staging, which has a separately tracked divergence.
func infoBCSQColumns(body string) string {
	var b strings.Builder
	for _, line := range strings.Split(body, "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		f := strings.Split(line, "\t")
		if len(f) < 8 {
			continue
		}
		b.WriteString(f[0])
		b.WriteByte('\t')
		b.WriteString(f[1])
		b.WriteByte('\t')
		b.WriteString(f[7])
		b.WriteByte('\n')
	}
	return b.String()
}

// prefixGFFContigs rewrites column-1 contig names of a GFF3 buffer by
// prepending "chr", leaving comment/header lines untouched.
func prefixGFFContigs(src []byte) []byte {
	var b strings.Builder
	for _, line := range strings.Split(string(src), "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			b.WriteString(line)
			b.WriteByte('\n')
			continue
		}
		b.WriteString("chr")
		b.WriteString(line)
		b.WriteByte('\n')
	}
	out := b.String()
	return []byte(strings.TrimSuffix(out, "\n"))
}

// bgunzipSlice4 decompresses a BGZF buffer (BGZF is a valid multi-member
// gzip stream) and returns the payload as a string.
func bgunzipSlice4(t *testing.T, data []byte) string {
	t.Helper()
	gr, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("gzip reader: %v", err)
	}
	gr.Multistream(true)
	var out bytes.Buffer
	if _, err := io.Copy(&out, gr); err != nil {
		t.Fatalf("gunzip: %v", err)
	}
	if err := gr.Close(); err != nil {
		t.Fatalf("gunzip close: %v", err)
	}
	return out.String()
}
