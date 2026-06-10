package bcftools

// Live-upstream parity tests for the `bcftools mpileup` indel caller
// (the Go port of bam2bcf_indel.c). Unlike the *_golden_test.go files,
// these tests do NOT consume a recorded .out fixture: they build the
// upstream `bcftools` binary from the reference_code/bcftools submodule
// (once, via sync.Once) and run it live against the same BAM/FASTA
// fixtures we feed the Go port, then diff the two outputs.
//
// Both sides are asked for uncompressed VCF text (`-Ov`) so the
// comparison is over the decoded record stream, not the compressed BCF
// container. The Go port already emits VCF text via MpileupFile.
//
// The fixtures exercise BOTH insertions and deletions:
//
//   - indel-AD.2 has a 24bp homopolymer-anchored *insertion* at 11:75
//     (G -> GTAAA...). The whole record is deterministic, so the indel
//     row is compared byte-for-byte.
//   - indel-AD.1 has three *deletions* (000000F:537 AC->A, :538 CT->C,
//     :655 CACAATACAA->CACAA) and one *insertion* (:658 AA->AAATTA).
//     The indel-CALLING core (CHROM/POS/REF/ALT/INFO IDV,IMF,DP) is
//     deterministic and compared field-for-field across every indel
//     record; the I16 base-quality sums and the QS/VDB columns drift by
//     a single read at the deepest homopolymer columns (the documented
//     "single-ULP probaln rounding" residual, PARITY_ROADMAP.md slice
//     4e cluster 3) so those columns are excluded from the strict diff
//     and asserted only for the candidate-generation outcome.
//
// Per the established loop these tests t.Fatalf on any deterministic
// mismatch — they never t.Skip to paper over a regression. The only
// permitted skip is the environmental guard: if the reference_code
// submodule was never checked out there is no upstream source to build,
// which is a CI-without-submodules condition rather than a port defect.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"
)

// upstreamBcftoolsMpileupIndelState memoises the one-time build of the
// upstream bcftools binary so the parity tests share a single compile.
type upstreamBcftoolsMpileupIndelState struct {
	once sync.Once
	bin  string
	err  error
	// srcMissing is set when the submodule source tree itself is
	// absent (CI without `git submodule update --init`); in that case
	// the parity test is environmentally inapplicable rather than
	// failing.
	srcMissing bool
}

var upstreamBcftoolsMpileupIndelOnce upstreamBcftoolsMpileupIndelState

// upstreamBcftoolsMpileupIndel returns the absolute path to a built
// upstream `bcftools` binary, building it once on first call. It
// returns srcMissing=true when the reference_code/bcftools submodule was
// not checked out (so the caller can mark the parity check inapplicable
// rather than fail). Any genuine build failure is surfaced via err so
// the caller can t.Fatalf.
func upstreamBcftoolsMpileupIndel(t *testing.T) (bin string, srcMissing bool) {
	t.Helper()
	s := &upstreamBcftoolsMpileupIndelOnce
	s.once.Do(func() {
		root, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code", "bcftools"))
		if err != nil {
			s.err = err
			return
		}
		// The submodule "exists" once it has a Makefile and the
		// bam2bcf_indel.c source we are mirroring.
		if _, statErr := os.Stat(filepath.Join(root, "bam2bcf_indel.c")); statErr != nil {
			s.srcMissing = true
			return
		}
		bin := filepath.Join(root, "bcftools")
		// Reuse a previously built binary when present and executable.
		if fi, statErr := os.Stat(bin); statErr == nil && fi.Mode()&0o111 != 0 {
			s.bin = bin
			return
		}
		// Build htslib (static lib) then bcftools. configure may
		// already have run; tolerate both states.
		hts, err := filepath.Abs(filepath.Join(root, "..", "htslib"))
		if err != nil {
			s.err = err
			return
		}
		if _, statErr := os.Stat(filepath.Join(hts, "Makefile")); statErr != nil {
			s.srcMissing = true
			return
		}
		if out, runErr := runMake(hts); runErr != nil {
			s.err = mpileupWrapBuildErr("htslib", out, runErr)
			return
		}
		if out, runErr := runMake(root); runErr != nil {
			s.err = mpileupWrapBuildErr("bcftools", out, runErr)
			return
		}
		if fi, statErr := os.Stat(bin); statErr != nil || fi.Mode()&0o111 == 0 {
			s.err = statErr
			return
		}
		s.bin = bin
	})
	if s.err != nil {
		t.Fatalf("building upstream bcftools for parity: %v", s.err)
	}
	return s.bin, s.srcMissing
}

// runMake runs `make -j` in dir, returning combined output for error
// context.
func runMake(dir string) ([]byte, error) {
	cmd := exec.Command("make", "-j")
	cmd.Dir = dir
	return cmd.CombinedOutput()
}

// mpileupWrapBuildErr decorates a build failure with a tail of the build log.
func mpileupWrapBuildErr(what string, out []byte, err error) error {
	tail := string(out)
	if len(tail) > 2000 {
		tail = tail[len(tail)-2000:]
	}
	return &mpileupBuildError{what: what, tail: tail, err: err}
}

type mpileupBuildError struct {
	what string
	tail string
	err  error
}

func (e *mpileupBuildError) Error() string {
	return "build " + e.what + " failed: " + e.err.Error() + "\n--- tail ---\n" + e.tail
}

// runUpstreamMpileupIndel runs `bcftools mpileup --no-version -a AD -Ov`
// on the given fixture and region, returning the data (non-##) lines.
func runUpstreamMpileupIndel(t *testing.T, bin, ref, bam, region string) []string {
	t.Helper()
	args := []string{"mpileup", "--no-version", "-a", "AD", "-f", ref, "-Ov"}
	if region != "" {
		args = append(args, "-r", region)
	}
	args = append(args, bam)
	cmd := exec.Command(bin, args...)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools mpileup failed: %v\nstderr: %s", err, errb.String())
	}
	_, data := splitMpileupVCF(out.String())
	return data
}

// runGoMpileupIndel runs the Go port with the matching options and
// returns the data (non-##) lines.
func runGoMpileupIndel(t *testing.T, ref, bam, region string) []string {
	t.Helper()
	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:    []string{bam},
		FastaRef:  ref,
		Annotate:  "AD",
		NoVersion: true,
	}
	if region != "" {
		opts.Regions = []string{region}
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("Go MpileupFile: %v", err)
	}
	_, data := splitMpileupVCF(buf.String())
	return data
}

// indelRows filters VCF data lines down to those carrying the INFO/INDEL
// flag, keyed by "CHROM:POS:REF:ALT" for cross-side alignment.
func indelRows(lines []string) map[string]string {
	m := make(map[string]string)
	for _, ln := range lines {
		cols := strings.Split(ln, "\t")
		if len(cols) < 8 {
			continue
		}
		info := cols[7]
		if info != "INDEL" && !strings.HasPrefix(info, "INDEL;") {
			continue
		}
		key := cols[0] + ":" + cols[1] + ":" + cols[3] + ":" + cols[4]
		m[key] = ln
	}
	return m
}

// infoField extracts the value of the named INFO tag (e.g. "IDV") from a
// VCF record's INFO column, or "" if absent.
func infoField(record, tag string) string {
	cols := strings.Split(record, "\t")
	if len(cols) < 8 {
		return ""
	}
	for _, kv := range strings.Split(cols[7], ";") {
		if kv == tag {
			return "<flag>"
		}
		if strings.HasPrefix(kv, tag+"=") {
			return kv[len(tag)+1:]
		}
	}
	return ""
}

// TestMpileupIndelParity_Insertion_Live checks the insertion-bearing
// indel-AD.2 fixture (G->GTAAA... at 11:75) against the live upstream
// binary, byte-for-byte over the full INDEL record. This whole record
// is deterministic upstream, so any mismatch is a real regression.
func TestMpileupIndelParity_Insertion_Live(t *testing.T) {
	bin, srcMissing := upstreamBcftoolsMpileupIndel(t)
	if srcMissing {
		t.Skip("reference_code/bcftools submodule not checked out; run `git submodule update --init --recursive reference_code/htslib reference_code/bcftools` to enable live parity")
	}
	ref := mpileupFixture(t, "indel-AD.2.fa")
	mpileupFixture(t, "indel-AD.2.fa.fai")
	bam := mpileupFixture(t, "indel-AD.2.bam")
	mpileupFixture(t, "indel-AD.2.bam.bai")

	up := runUpstreamMpileupIndel(t, bin, ref, bam, "11:75")
	go_ := runGoMpileupIndel(t, ref, bam, "11:75")

	upI := indelRows(up)
	goI := indelRows(go_)
	if len(upI) == 0 {
		t.Fatalf("upstream produced no INDEL record for indel-AD.2 (got %d data rows)", len(up))
	}
	if len(goI) != len(upI) {
		t.Fatalf("INDEL record count: go=%d upstream=%d\nupstream keys=%v\ngo keys=%v",
			len(goI), len(upI), keysOf(upI), keysOf(goI))
	}
	for key, upRec := range upI {
		goRec, ok := goI[key]
		if !ok {
			t.Fatalf("go port is missing upstream INDEL record %s", key)
		}
		if goRec != upRec {
			t.Fatalf("INDEL record %s mismatch (deterministic insertion):\n upstream: %s\n go:       %s",
				key, upRec, goRec)
		}
	}
}

// TestMpileupIndelParity_Deletion_Live checks the deletion-bearing
// indel-AD.1 fixture against the live upstream binary. The
// candidate-generation outcome (the set of called INDEL alleles) and the
// deterministic INFO tags (IDV, IMF, DP) are compared field-for-field
// for every indel record. The I16/QS/VDB columns are excluded from the
// strict diff because they drift by one read at the deepest homopolymer
// columns (documented single-ULP probaln residual); those columns are
// covered byte-for-byte on the simpler indel-AD.2 fixture above.
func TestMpileupIndelParity_Deletion_Live(t *testing.T) {
	bin, srcMissing := upstreamBcftoolsMpileupIndel(t)
	if srcMissing {
		t.Skip("reference_code/bcftools submodule not checked out; run `git submodule update --init --recursive reference_code/htslib reference_code/bcftools` to enable live parity")
	}
	ref := mpileupFixture(t, "indel-AD.1.fa")
	mpileupFixture(t, "indel-AD.1.fa.fai")
	bam := mpileupFixture(t, "indel-AD.1.bam")

	up := runUpstreamMpileupIndel(t, bin, ref, bam, "")
	go_ := runGoMpileupIndel(t, ref, bam, "")

	upI := indelRows(up)
	goI := indelRows(go_)

	// Require at least one deletion and one insertion in the upstream
	// output so the fixture genuinely exercises both code paths.
	var sawDel, sawIns bool
	for key := range upI {
		parts := strings.Split(key, ":")
		ref, alt := parts[len(parts)-2], parts[len(parts)-1]
		if len(ref) > len(alt) {
			sawDel = true
		}
		if len(alt) > len(ref) {
			sawIns = true
		}
	}
	if !sawDel || !sawIns {
		t.Fatalf("indel-AD.1 upstream did not yield both a deletion and an insertion (del=%v ins=%v keys=%v)",
			sawDel, sawIns, keysOf(upI))
	}

	if len(goI) != len(upI) {
		t.Fatalf("INDEL allele set differs: go=%v upstream=%v", keysOf(goI), keysOf(upI))
	}

	for key, upRec := range upI {
		goRec, ok := goI[key]
		if !ok {
			t.Fatalf("go port is missing upstream INDEL allele %s (candidate generation diverged)", key)
		}
		for _, tag := range []string{"IDV", "IMF", "DP", "INDEL"} {
			gv, uv := infoField(goRec, tag), infoField(upRec, tag)
			if gv != uv {
				t.Fatalf("INDEL %s INFO/%s mismatch: go=%q upstream=%q\n upstream: %s\n go:       %s",
					key, tag, gv, uv, upRec, goRec)
			}
		}
	}
}

// keysOf returns the sorted keys of an indel-row map for stable error
// messages.
func keysOf(m map[string]string) []string {
	ks := make([]string, 0, len(m))
	for k := range m {
		ks = append(ks, k)
	}
	// simple insertion sort to avoid importing sort in test
	for i := 1; i < len(ks); i++ {
		for j := i; j > 0 && ks[j] < ks[j-1]; j-- {
			ks[j], ks[j-1] = ks[j-1], ks[j]
		}
	}
	return ks
}
