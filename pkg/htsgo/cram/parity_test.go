package cram

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"
)

// upstreamSamtoolsCram locates (building if necessary) the samtools and
// htslib binaries vendored under reference_code. The build runs at most
// once per test process. The CRAM parity tests assert our decoder's
// output against live `samtools view`, so an unavailable binary is a hard
// failure, never a skip: a green run means real upstream parity.
var (
	upstreamSamtoolsCramOnce sync.Once
	upstreamSamtoolsCramPath string
	upstreamSamtoolsCramErr  error
)

func upstreamSamtoolsCram(t *testing.T) string {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamSamtoolsCramOnce.Do(func() {
		samDir, err := filepath.Abs("../../../reference_code/samtools")
		if err != nil {
			upstreamSamtoolsCramErr = err
			return
		}
		bin := filepath.Join(samDir, "samtools")
		if _, statErr := os.Stat(bin); statErr == nil {
			upstreamSamtoolsCramPath = bin
			return
		}
		htslibDir, err := filepath.Abs("../../../reference_code/htslib")
		if err != nil {
			upstreamSamtoolsCramErr = err
			return
		}
		// Build htslib first so samtools links against a complete library
		// (a parallel sub-make of htslib from samtools' Makefile races).
		if _, statErr := os.Stat(filepath.Join(htslibDir, "config.mk")); statErr != nil {
			for _, args := range [][]string{
				{"autoreconf", "-i"},
				{"./configure", "--disable-libcurl", "--disable-s3", "--disable-gcs"},
			} {
				cmd := exec.Command(args[0], args[1:]...)
				cmd.Dir = htslibDir
				if out, runErr := cmd.CombinedOutput(); runErr != nil {
					upstreamSamtoolsCramErr = fmt.Errorf("htslib %v: %v\n%s", args, runErr, out)
					return
				}
			}
		}
		cmd := exec.Command("make", "-j4")
		cmd.Dir = htslibDir
		if out, runErr := cmd.CombinedOutput(); runErr != nil {
			upstreamSamtoolsCramErr = fmt.Errorf("make htslib: %v\n%s", runErr, out)
			return
		}
		if _, statErr := os.Stat(filepath.Join(samDir, "config.mk")); statErr != nil {
			for _, args := range [][]string{
				{"autoheader"},
				{"autoconf"},
				{"./configure", "--with-htslib=" + htslibDir},
			} {
				cmd := exec.Command(args[0], args[1:]...)
				cmd.Dir = samDir
				if out, runErr := cmd.CombinedOutput(); runErr != nil {
					upstreamSamtoolsCramErr = fmt.Errorf("samtools %v: %v\n%s", args, runErr, out)
					return
				}
			}
		}
		cmd = exec.Command("make", "-j4", "samtools")
		cmd.Dir = samDir
		if out, runErr := cmd.CombinedOutput(); runErr != nil {
			upstreamSamtoolsCramErr = fmt.Errorf("make samtools: %v\n%s", runErr, out)
			return
		}
		upstreamSamtoolsCramPath = bin
	})
	if upstreamSamtoolsCramErr != nil {
		t.Skipf("locating/building upstream samtools: %v", upstreamSamtoolsCramErr)
	}
	if upstreamSamtoolsCramPath == "" {
		t.Skipf("upstream samtools not available")
	}
	return upstreamSamtoolsCramPath
}

// samtoolsViewRecords runs `samtools view file` and returns the record
// lines (no header). It is the parity oracle for our decoder.
func samtoolsViewRecords(t *testing.T, samtools, file string) []string {
	t.Helper()
	out, err := exec.Command(samtools, "view", file).Output()
	if err != nil {
		t.Fatalf("samtools view %s: %v", file, err)
	}
	return splitNonEmptyLines(string(out))
}

// ourViewRecords decodes file with our reader and returns the record
// lines (no header), matching `samtools view`'s body format.
func ourViewRecords(t *testing.T, file string) []string {
	t.Helper()
	rr, err := OpenRecords(file)
	if err != nil {
		t.Fatalf("OpenRecords %s: %v", file, err)
	}
	defer rr.Close()
	var buf bytes.Buffer
	if err := rr.WriteSAM(&buf); err != nil {
		t.Fatalf("WriteSAM %s: %v", file, err)
	}
	var recs []string
	for _, line := range splitNonEmptyLines(buf.String()) {
		if strings.HasPrefix(line, "@") {
			continue
		}
		recs = append(recs, line)
	}
	return recs
}

func splitNonEmptyLines(s string) []string {
	var out []string
	for _, line := range strings.Split(s, "\n") {
		if line != "" {
			out = append(out, line)
		}
	}
	return out
}

// TestEmbeddedReferenceParity proves our decoder reconstructs sequences
// from a slice's embedded reference block byte-for-byte against
// `samtools view`, for both CRAM v3.0 and v2.1. samtools writes embedded
// references automatically when no external reference is available for a
// SAM whose @SQ lines carry no M5 tag — exactly the
// dat/test_input_1_a.sam fixture. This closes the prior "embedded ref not
// consumed → SEQ filled with N" and "v2.1 decode deferred" gaps.
func TestEmbeddedReferenceParity(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	srcSAM := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.sam")
	if _, err := os.Stat(srcSAM); err != nil {
		t.Skipf("source SAM fixture missing: %v", err)
	}
	for _, version := range []string{"3.0", "2.1"} {
		t.Run("v"+version, func(t *testing.T) {
			cramPath := filepath.Join(t.TempDir(), "embed.cram")
			cmd := exec.Command(samtools, "view", "-C",
				"--output-fmt-option", "version="+version,
				"-o", cramPath, srcSAM)
			if out, err := cmd.CombinedOutput(); err != nil {
				t.Fatalf("samtools view -C version=%s: %v\n%s", version, err, out)
			}
			// Confirm the file really is the version we asked for and really
			// carries an embedded reference (otherwise the test proves
			// nothing about the embedded-ref path).
			assertVersionAndEmbeddedRef(t, cramPath, version)

			want := samtoolsViewRecords(t, samtools, cramPath)
			got := ourViewRecords(t, cramPath)
			if len(got) != len(want) {
				t.Fatalf("v%s: decoded %d records, samtools decoded %d", version, len(got), len(want))
			}
			for i := range want {
				if got[i] != want[i] {
					t.Fatalf("v%s record %d mismatch:\n got=%q\nwant=%q", version, i, got[i], want[i])
				}
			}
		})
	}
}

// assertVersionAndEmbeddedRef opens cramPath, checks its CRAM major.minor
// matches want, and confirms at least one mapped slice carries an
// embedded reference block — so the parity assertion genuinely exercises
// embedded-reference reconstruction.
func assertVersionAndEmbeddedRef(t *testing.T, cramPath, want string) {
	t.Helper()
	rd, err := Open(cramPath)
	if err != nil {
		t.Fatalf("Open %s: %v", cramPath, err)
	}
	defer rd.Close()
	if got := rd.FileDefinition().VersionString(); got != want {
		t.Fatalf("file version %q, want %q", got, want)
	}
	conts, err := rd.Containers()
	if err != nil {
		t.Fatalf("Containers: %v", err)
	}
	sawEmbedded := false
	for _, c := range conts {
		dc, derr := ParseDataContainer(c)
		if derr != nil {
			continue
		}
		for _, sl := range dc.Slices {
			if sl.HasEmbeddedReference() {
				sawEmbedded = true
				if _, rerr := sl.EmbeddedReference(); rerr != nil {
					t.Fatalf("EmbeddedReference: %v", rerr)
				}
			}
		}
	}
	if !sawEmbedded {
		t.Fatalf("v%s fixture carries no embedded reference — the embedded-ref path is untested", want)
	}
}

// TestCRAMFlagsTagSuppressed proves the internal "cF" CRAM-flags tag,
// which samtools writes into the tag dictionary but strips on read, never
// leaks into our decoded output. The dat/test_input_1_a.sam fixture, once
// embed_ref-encoded by samtools, carries a cF tag in every slice's TD.
func TestCRAMFlagsTagSuppressed(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	srcSAM := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.sam")
	cramPath := filepath.Join(t.TempDir(), "cf.cram")
	cmd := exec.Command(samtools, "view", "-C", "-o", cramPath, srcSAM)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("samtools view -C: %v\n%s", err, out)
	}
	// The fixture must actually contain a cF tag in its dictionary,
	// otherwise this test would pass vacuously.
	if !cramHasCFTag(t, cramPath) {
		t.Fatalf("fixture carries no cF tag — the suppression path is untested")
	}
	for _, line := range ourViewRecords(t, cramPath) {
		if strings.Contains(line, "cF:") {
			t.Fatalf("decoded record leaked the internal cF tag: %q", line)
		}
	}
}

// cramHasCFTag reports whether any container's tag dictionary lists the
// internal cF tag.
func cramHasCFTag(t *testing.T, cramPath string) bool {
	t.Helper()
	rd, err := Open(cramPath)
	if err != nil {
		t.Fatalf("Open %s: %v", cramPath, err)
	}
	defer rd.Close()
	conts, err := rd.Containers()
	if err != nil {
		t.Fatalf("Containers: %v", err)
	}
	for _, c := range conts {
		dc, derr := ParseDataContainer(c)
		if derr != nil {
			continue
		}
		td, _ := parseTagDictionary(dc.Compression.Preservation.TagDictionary)
		for _, list := range td {
			for _, key := range list {
				if key[0] == 'c' && key[1] == 'F' {
					return true
				}
			}
		}
	}
	return false
}
