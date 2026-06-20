package conformance

import (
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// htslibSAMFixtures are SAM fixtures shipped by htslib's own test suite, each
// exercising a tricky parser path. "ref" names the FASTA in the same directory
// needed for CRAM (empty if CRAM is not attempted for this fixture).
var htslibSAMFixtures = []struct {
	name string // fixture basename without .sam
	ref  string // reference FASTA basename, "" if none
	why  string // what edge case it exercises
}{
	{"c1#noseq", "c1.fa", "records with no SEQ ('*') round-trip"},
	{"c1#pad1", "c1.fa", "padded CIGAR (P ops) variant 1"},
	{"c1#pad2", "c1.fa", "padded CIGAR variant 2"},
	{"c1#pad3", "c1.fa", "padded CIGAR variant 3"},
	{"c1#bounds", "c1.fa", "alignments at reference bounds"},
	{"c1#clip", "c1.fa", "hard/soft clipping at ends"},
	{"c1#unknown", "c1.fa", "reads vs unknown/extra reference"},
	{"c2#pad", "c2.fa", "padded CIGAR on second contig"},
	{"index_dos", "", "CRLF (DOS) line endings in SAM input"},
}

// TestHtslibSAM_RoundTripParity feeds each htslib SAM fixture through our
// samtools (SAM -> BAM -> SAM) and through upstream samtools, and asserts that
// the resulting record sets are identical (modulo @PG provenance and order).
// This proves we accept exactly the inputs upstream accepts and reproduce the
// same records, on htslib's own adversarial fixtures.
func TestHtslibSAM_RoundTripParity(t *testing.T) {
	dir := htslibTest(t)
	our := ourSamtools(t)
	up := upSamtools(t)

	for _, fx := range htslibSAMFixtures {
		fx := fx
		t.Run(fx.name, func(t *testing.T) {
			src := filepath.Join(dir, fx.name+".sam")
			if _, err := os.Stat(src); err != nil {
				t.Skipf("fixture %s missing: %v", src, err)
			}

			work := t.TempDir()
			ourBAM := filepath.Join(work, "our.bam")
			upBAM := filepath.Join(work, "up.bam")

			// Both producers: SAM -> BAM.
			if _, errOut, err := runCapture(t, our, "view", "-b", "-o", ourBAM, src); err != nil {
				t.Fatalf("our samtools view -b rejected %s (%s): %v\n%s", fx.name, fx.why, err, errOut)
			}
			if _, errOut, err := runCapture(t, up, "view", "-b", "-o", upBAM, src); err != nil {
				t.Skipf("upstream samtools rejected %s (%s): %v\n%s", fx.name, fx.why, err, errOut)
			}

			// BAM -> SAM (with header) for both.
			ourSAM, errOut, err := runCapture(t, our, "view", "-h", ourBAM)
			if err != nil {
				t.Fatalf("our samtools view -h failed on %s: %v\n%s", fx.name, err, errOut)
			}
			upSAM, errOut, err := runCapture(t, up, "view", "-h", upBAM)
			if err != nil {
				t.Fatalf("upstream samtools view -h failed on %s: %v\n%s", fx.name, err, errOut)
			}

			got := upstream.NormalizeSAM(ourSAM, true)
			want := upstream.NormalizeSAM(upSAM, true)
			if got != want {
				t.Errorf("%s (%s): our SAM round-trip diverges from upstream\n--- ours ---\n%s\n--- upstream ---\n%s",
					fx.name, fx.why, got, want)
			}
		})
	}
}

// seqOf returns the SEQ column (field 10) of every alignment record in sam,
// keyed nowhere — just the multiset of SEQ strings in record order, so a
// reference-relative codec that silently mis-decodes bases is caught.
func seqOf(sam string) []string {
	var out []string
	for _, ln := range strings.Split(sam, "\n") {
		if ln == "" || strings.HasPrefix(ln, "@") {
			continue
		}
		cols := strings.Split(ln, "\t")
		if len(cols) >= 10 {
			out = append(out, cols[0]+"\t"+cols[9])
		}
	}
	sort.Strings(out)
	return out
}

// TestHtslibSAM_CRAMRoundTrip encodes each reference-bearing SAM fixture to
// CRAM with our samtools using the fixture's external reference (-T), decodes
// it back, and asserts two things:
//
//  1. every read's SEQ survives the reference-relative encode/decode unchanged
//     (CRAM stores bases as deltas against the reference, so a reference-handling
//     bug shows up as silently corrupted bases — the highest-value failure to
//     guard against); and
//  2. our CRAM round-trip reproduces the SAME records as upstream samtools'
//     CRAM round-trip (CRAM auto-adds MD/NM tags, so the correct oracle is
//     upstream's own CRAM output, not the original SAM).
func TestHtslibSAM_CRAMRoundTrip(t *testing.T) {
	dir := htslibTest(t)
	our := ourSamtools(t)
	up := upSamtools(t)

	for _, fx := range htslibSAMFixtures {
		if fx.ref == "" {
			continue
		}
		fx := fx
		t.Run(fx.name, func(t *testing.T) {
			src := filepath.Join(dir, fx.name+".sam")
			ref := filepath.Join(dir, fx.ref)
			if _, err := os.Stat(src); err != nil {
				t.Skipf("fixture missing: %v", err)
			}
			if _, err := os.Stat(ref); err != nil {
				t.Skipf("reference %s missing: %v", ref, err)
			}

			work := t.TempDir()
			ourCRAM := filepath.Join(work, "our.cram")
			upCRAM := filepath.Join(work, "up.cram")

			// Establish the upstream oracle first: encode+decode CRAM with
			// upstream. If upstream itself can't handle this fixture on this
			// build, SKIP (don't assert against a moving target).
			if _, errOut, err := runCapture(t, up, "view", "-C", "-T", ref, "-o", upCRAM, src); err != nil {
				t.Skipf("upstream samtools view -C rejected %s: %v\n%s", fx.name, err, errOut)
			}
			upBack, errOut, err := runCapture(t, up, "view", "-h", "-T", ref, upCRAM)
			if err != nil {
				t.Skipf("upstream samtools CRAM decode failed on %s: %v\n%s", fx.name, err, errOut)
			}

			// Our SAM -> CRAM. Upstream accepted it, so any rejection here is a
			// genuine encoder parity gap (silent-corruption-adjacent: a refused
			// encode of valid input).
			if _, errOut, err := runCapture(t, our, "view", "-C", "-T", ref, "-o", ourCRAM, src); err != nil {
				t.Errorf("PARITY GAP: our samtools CRAM-encode rejected %s (%s) that upstream encodes fine.\n  our stderr: %s",
					fx.name, fx.why, strings.TrimSpace(errOut))
				return
			}
			ourBack, errOut, err := runCapture(t, our, "view", "-h", "-T", ref, ourCRAM)
			if err != nil {
				t.Errorf("PARITY GAP: our samtools CRAM-decode failed on %s (%s) that upstream round-trips fine.\n  our stderr: %s",
					fx.name, fx.why, strings.TrimSpace(errOut))
				return
			}

			// Source SAM read back through our samtools (baseline for SEQ).
			origSAM, errOut, err := runCapture(t, our, "view", "-h", src)
			if err != nil {
				t.Fatalf("our samtools view -h on source failed: %v\n%s", err, errOut)
			}

			// Assertion 1: SEQ bases survive the reference-relative round-trip
			// (the highest-value check: silent base corruption).
			if got, want := seqOf(ourBack), seqOf(origSAM); strings.Join(got, "\n") != strings.Join(want, "\n") {
				t.Errorf("%s: CRAM round-trip silently corrupted SEQ bases\n--- after CRAM ---\n%s\n--- original ---\n%s",
					fx.name, strings.Join(got, "\n"), strings.Join(want, "\n"))
			}

			// Assertion 2: our CRAM round-trip == upstream CRAM round-trip.
			got := upstream.NormalizeSAM(ourBack, true)
			want := upstream.NormalizeSAM(upBack, true)
			if got != want {
				t.Errorf("%s: our CRAM round-trip diverges from upstream's\n--- ours ---\n%s\n--- upstream ---\n%s",
					fx.name, got, want)
			}
		})
	}
}
