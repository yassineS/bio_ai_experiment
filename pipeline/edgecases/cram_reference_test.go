package edgecases

import (
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// cramRefFASTA is a tiny two-contig reference. The two contigs differ so that a
// reference mix-up between contigs (the RefSeqID=-2 multi-reference-slice bug
// class) would corrupt decoded bases detectably.
const cramRefFASTA = `>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA
`

// cramMultiRefSAM places reads on BOTH contigs. With the default slice size
// these reads share a single CRAM slice whose embedded reference id is the
// multi-reference sentinel (-2); decoding must select the per-record contig,
// not a single slice-wide reference. A read with a mismatch (X) against the
// reference exercises the substitution-feature path, where a wrong reference
// silently yields a wrong base.
const cramMultiRefSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:64
@SQ	SN:chr2	LN:64
r1	0	chr1	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF
r2	0	chr1	9	60	8M	*	0	0	ACGTACGT	FFFFFFFF
r3	0	chr1	17	60	8M	*	0	0	ACGTAAGT	FFFFFFFF
r4	0	chr2	1	60	8M	*	0	0	TTTTGGGG	FFFFFFFF
r5	0	chr2	9	60	8M	*	0	0	CCCCAAAA	FFFFFFFF
r6	0	chr2	17	60	8M	*	0	0	TTTTGCGG	FFFFFFFF
`

// seqByName extracts a sorted "QNAME\tSEQ" list from SAM text — the multiset of
// decoded bases, keyed by read, that a reference-handling bug would corrupt.
func seqByName(sam string) []string {
	var out []string
	for _, ln := range strings.Split(sam, "\n") {
		if ln == "" || strings.HasPrefix(ln, "@") {
			continue
		}
		c := strings.Split(ln, "\t")
		if len(c) >= 10 {
			out = append(out, c[0]+"\t"+c[9])
		}
	}
	sort.Strings(out)
	return out
}

// TestCRAMReferenceHandling is the highest-value silent-corruption guard:
// reference-relative CRAM encode/decode must reproduce every read's SEQ
// exactly. It covers (a) an external -T reference and (b) the
// multi-reference-slice case (reads on two contigs sharing one slice — the
// RefSeqID=-2 class fixed in commit 5f811f0). It asserts our decoded SEQ matches
// the source AND matches upstream samtools' own CRAM round-trip.
func TestCRAMReferenceHandling(t *testing.T) {
	our := ourBin(t, "samtools")
	up := upBin(t, "samtools")

	cases := []struct {
		name string
		sam  string
	}{
		{"external_reference_single_contig", `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:64
r1	0	chr1	1	60	8M	*	0	0	ACGTACGT	FFFFFFFF
r2	0	chr1	9	60	8M	*	0	0	ACGTAAGT	FFFFFFFF
r3	16	chr1	17	60	8M	*	0	0	ACGTACGT	FFFFFFFF
`},
		{"multi_reference_slice", cramMultiRefSAM},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			ref := writeFile(t, dir, "ref.fa", cramRefFASTA)
			// Generate the .fai with upstream samtools faidx (our samtools has
			// no faidx subcommand and our CRAM path does not require it, but
			// upstream's reader uses it). Skip if upstream faidx fails.
			if _, errOut, err := run(t, up, "faidx", ref); err != nil {
				t.Skipf("upstream samtools faidx failed: %v\n%s", err, errOut)
			}
			src := writeFile(t, dir, "in.sam", tc.sam)

			ourCRAM := filepath.Join(dir, "our.cram")
			upCRAM := filepath.Join(dir, "up.cram")

			// Our SAM -> CRAM -> SAM.
			mustRun(t, our, "view", "-C", "-T", ref, "-o", ourCRAM, src)
			ourBack := mustRun(t, our, "view", "-T", ref, ourCRAM)

			// Assertion 1: decoded SEQ == source SEQ (no silent base corruption).
			if got, want := seqByName(ourBack), seqByName(tc.sam); strings.Join(got, "\n") != strings.Join(want, "\n") {
				t.Fatalf("decoded SEQ != source SEQ (silent corruption!):\n got: %s\nwant: %s",
					strings.Join(got, "\n"), strings.Join(want, "\n"))
			}

			// Assertion 2: our decoded SEQ == upstream CRAM round-trip SEQ.
			mustRun(t, up, "view", "-C", "-T", ref, "-o", upCRAM, src)
			upBack := mustRun(t, up, "view", "-T", ref, upCRAM)
			if got, want := seqByName(ourBack), seqByName(upBack); strings.Join(got, "\n") != strings.Join(want, "\n") {
				t.Errorf("our CRAM SEQ diverges from upstream CRAM SEQ:\n ours: %s\nups: %s",
					strings.Join(got, "\n"), strings.Join(want, "\n"))
			}
		})
	}
}
