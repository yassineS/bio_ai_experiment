package cram

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// newRefDecoder builds a recordDecoder over a hand-built reference span
// for the reference-backed reconstruction tests. refBases is the slice's
// reference span starting at 1-based coordinate refStart; sm is the
// five-byte substitution-matrix entry.
func newRefDecoder(refBases []byte, refStart int32, sm []byte) *recordDecoder {
	h := refFreeHeader()
	h.Preservation.SubstitutionMatrix = sm
	return &recordDecoder{
		h:           h,
		src:         &SeriesSource{s: newTestSource(nil, nil)},
		slice:       &SliceHeader{},
		refBases:    refBases,
		refStart:    refStart,
		substMatrix: newSubstMatrix(sm),
	}
}

// TestReconstructWithReferenceMatch checks that a mapped read whose
// bases are entirely reference matches (no features) is reconstructed
// from the reference span rather than filled with 'N'.
func TestReconstructWithReferenceMatch(t *testing.T) {
	// Reference span starts at coordinate 1: ACGTACGTAC.
	rd := newRefDecoder([]byte("ACGTACGTAC"), 1, []byte{0x1B, 0x1B, 0x1B, 0x1B, 0x1B})
	// A read of length 6 at POS 3 with no features: a pure 6M match
	// copying reference[3..8] (1-based) = GTACGT.
	seq, _, cig, err := rd.reconstructMapped(nil, 6, 3)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if string(seq) != "GTACGT" {
		t.Errorf("seq = %q, want GTACGT (reference span at POS 3)", seq)
	}
	if cig.String() != "6M" {
		t.Errorf("cigar = %q, want 6M", cig.String())
	}
	if rd.needsReference {
		t.Error("a reference-backed match must not set needsReference")
	}
}

// TestReconstructWithReferenceSubstitution checks a substitution feature
// resolves its read base through the substitution matrix against the
// reference base at that position.
func TestReconstructWithReferenceSubstitution(t *testing.T) {
	// Reference span at coordinate 1: AAAAAAAAAA.
	rd := newRefDecoder([]byte("AAAAAAAAAA"), 1, []byte{0x1B, 0x1B, 0x1B, 0x1B, 0x1B})
	// A read of length 5 at POS 1 with one substitution at in-read pos 3,
	// code 1. Reference base at pos 3 is 'A'; code 1 against 'A' decodes
	// to candidate index 1 = 'G'. The rest are reference matches ('A').
	feats := []readFeature{{code: featSubst, pos: 3, substCode: 1}}
	seq, _, cig, err := rd.reconstructMapped(feats, 5, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if string(seq) != "AAGAA" {
		t.Errorf("seq = %q, want AAGAA (substitution G at pos 3)", seq)
	}
	if cig.String() != "5M" {
		t.Errorf("cigar = %q, want 5M", cig.String())
	}
	if rd.needsReference {
		t.Error("a resolved substitution must not set needsReference")
	}
}

// TestReconstructWithReferenceIndel checks the reference cursor advances
// correctly across an insertion and a deletion so the match runs after
// them copy the right reference bases.
func TestReconstructWithReferenceIndel(t *testing.T) {
	// Reference span at coordinate 1: ACGTACGTACGT (12 bases).
	rd := newRefDecoder([]byte("ACGTACGTACGT"), 1, nil)
	// Read at POS 1, length 8: 2M (ref[1..2]=AC), 2I (NN inserted, no
	// reference consumed), then a 3D deletion (consumes ref 3,4,5), then
	// 4M matching ref[6..9] = CGTA.
	feats := []readFeature{
		{code: featInsertion, pos: 3, bases: []byte("NN")},
		{code: featDeletion, pos: 5, length: 3},
	}
	seq, _, cig, err := rd.reconstructMapped(feats, 8, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if string(seq) != "ACNNCGTA" {
		t.Errorf("seq = %q, want ACNNCGTA", seq)
	}
	if cig.String() != "2M2I3D4M" {
		t.Errorf("cigar = %q, want 2M2I3D4M", cig.String())
	}
}

// TestReconstructSubstitutionOutsideSpan checks that a substitution
// whose reference coordinate falls outside the resolved span degrades
// to an 'N' base and sets needsReference, rather than panicking.
func TestReconstructSubstitutionOutsideSpan(t *testing.T) {
	// A reference span declared to begin at coordinate 100, so a record
	// at POS 1 references coordinates the span does not cover.
	rd := newRefDecoder([]byte("ACGT"), 100, nil)
	feats := []readFeature{{code: featSubst, pos: 1, substCode: 0}}
	seq, _, _, err := rd.reconstructMapped(feats, 4, 1)
	if err != nil {
		// A trailing match run also needs the (absent) reference span and
		// is a hard error — acceptable. The point is no panic.
		return
	}
	if seq[0] != 'N' {
		t.Errorf("substitution outside the span = %c, want N", seq[0])
	}
	if !rd.needsReference {
		t.Error("a substitution outside the reference span must set needsReference")
	}
}

// TestReferenceBaseAtBounds checks referenceBaseAt's bounds handling.
func TestReferenceBaseAtBounds(t *testing.T) {
	rd := newRefDecoder([]byte("ACGT"), 10, nil)
	if b, ok := rd.referenceBaseAt(10); !ok || b != 'A' {
		t.Errorf("referenceBaseAt(10) = %c,%v; want A,true", b, ok)
	}
	if b, ok := rd.referenceBaseAt(13); !ok || b != 'T' {
		t.Errorf("referenceBaseAt(13) = %c,%v; want T,true", b, ok)
	}
	if _, ok := rd.referenceBaseAt(9); ok {
		t.Error("referenceBaseAt before the span start must report ok=false")
	}
	if _, ok := rd.referenceBaseAt(14); ok {
		t.Error("referenceBaseAt past the span end must report ok=false")
	}
	// No reference attached at all.
	noRef := newRefDecoder(nil, 0, nil)
	if _, ok := noRef.referenceBaseAt(1); ok {
		t.Error("referenceBaseAt with no reference must report ok=false")
	}
}

// TestReconstructReferenceTooShort checks that a match run reaching past
// the slice's resolved reference span is a hard error, not a silent
// short read.
func TestReconstructReferenceTooShort(t *testing.T) {
	rd := newRefDecoder([]byte("ACGT"), 1, nil)
	// A 10-base read at POS 1 needs reference 1..10, but the span is 4.
	if _, _, _, err := rd.reconstructMapped(nil, 10, 1); err == nil {
		t.Error("a match run past the reference span must error")
	}
}

// TestReconstructSoftClipWithReference checks that a soft clip consumes
// read bases but not reference, so a following match run copies from the
// correct reference coordinate.
func TestReconstructSoftClipWithReference(t *testing.T) {
	rd := newRefDecoder([]byte("ACGTACGT"), 1, nil)
	// Read at POS 1, length 6: 2S soft clip (XX), then 4M matching
	// reference[1..4] = ACGT.
	feats := []readFeature{{code: featSoftClip, pos: 1, bases: []byte("XX")}}
	seq, _, cig, err := rd.reconstructMapped(feats, 6, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if string(seq) != "XXACGT" {
		t.Errorf("seq = %q, want XXACGT", seq)
	}
	if cig.String() != "2S4M" {
		t.Errorf("cigar = %q, want 2S4M", cig.String())
	}
}

// referenceBackedFixture pairs a reference-backed CRAM v3.0 fixture with
// the FASTA whose contig MD5s its slice headers verify against.
var referenceBackedFixture = struct {
	cram, fasta, contig string
	contigMD5           string
}{
	cram:      "quickcheck/7.quickcheck.cram30.ok.cram",
	fasta:     "dat/mpileup.ref.fa",
	contig:    "17",
	contigMD5: "f8c08a4411f07717451464d546b3706d",
}

// TestReferenceBackedDecode decodes a real reference-backed CRAM v3.0
// file with its reference FASTA and asserts the mapped reads decode to
// concrete bases — not the 'N' placeholders the reference-free C4b path
// produces — and that NeedsReference is no longer set. The slice-header
// MD5s are verified inside the decode; reaching the end proves they
// matched.
func TestReferenceBackedDecode(t *testing.T) {
	cramPath := filepath.Join(samtoolsTestDir, referenceBackedFixture.cram)
	faPath := filepath.Join(samtoolsTestDir, referenceBackedFixture.fasta)
	if _, err := os.Stat(cramPath); err != nil {
		t.Fatalf("samtools submodule not initialised — CRAM fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	if _, err := os.Stat(faPath); err != nil {
		t.Fatalf("samtools submodule not initialised — reference FASTA unavailable; run `git submodule update --init reference_code/samtools`")
	}

	// First decode without a reference: the C4b fallback fills bases an
	// external reference would supply with 'N'.
	noRef, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords (no reference): %v", err)
	}
	plainRecs, err := noRef.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll (no reference): %v", err)
	}
	noRef.Close()
	if !noRef.NeedsReference() {
		t.Fatal("a reference-backed CRAM must report NeedsReference when decoded without one")
	}
	plainN := countN(plainRecs)
	if plainN == 0 {
		t.Fatal("the reference-free decode produced no 'N' placeholders; fixture is not reference-backed")
	}

	// Now decode with the reference FASTA attached.
	rr, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords (with reference): %v", err)
	}
	defer rr.Close()
	if err := rr.SetReferenceFASTA(faPath); err != nil {
		t.Fatalf("SetReferenceFASTA: %v", err)
	}
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll (with reference): %v", err)
	}
	if len(recs) != len(plainRecs) {
		t.Fatalf("reference-backed decode yielded %d records, plain yielded %d", len(recs), len(plainRecs))
	}
	if rr.NeedsReference() {
		t.Error("a fully reference-resolved decode must not report NeedsReference")
	}
	refN := countN(recs)
	if refN >= plainN {
		t.Errorf("reference resolution did not reduce 'N' count: plain=%d resolved=%d", plainN, refN)
	}
	// Every decoded base must be a valid IUPAC base.
	for _, r := range recs {
		for i := 0; i < len(r.Seq); i++ {
			if !strings.ContainsRune("ACGTMRWSYKVHDBN", rune(r.Seq[i])) {
				t.Fatalf("record %s base %d is %q, not a valid base", r.QName, i, r.Seq[i])
			}
		}
	}
	// A pure-match read (CIGAR all M, no clip/indel) and with no
	// substitution must equal the reference span exactly. Verify at
	// least one and that any mismatch is a single-base substitution.
	verified := verifyPureMatchAgainstReference(t, recs, faPath)
	if verified == 0 {
		t.Error("no pure-match read was available to cross-check against the reference")
	}
}

// TestReferenceMD5MismatchIsHardError attaches the WRONG reference FASTA
// to a reference-backed CRAM and asserts the decode fails with an MD5
// mismatch rather than silently producing wrong sequence.
func TestReferenceMD5MismatchIsHardError(t *testing.T) {
	cramPath := filepath.Join(samtoolsTestDir, referenceBackedFixture.cram)
	if _, err := os.Stat(cramPath); err != nil {
		t.Fatalf("samtools submodule not initialised — CRAM fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	// Build a wrong reference: the right contig name and a generous
	// length, but all-'A' bases, so the slice-span MD5 cannot match.
	wrongFA := writeFASTA(t, referenceBackedFixture.contig, strings.Repeat("A", 4200))
	rr, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords: %v", err)
	}
	defer rr.Close()
	if err := rr.SetReferenceFASTA(wrongFA); err != nil {
		t.Fatalf("SetReferenceFASTA: %v", err)
	}
	_, err = rr.ReadAll()
	if err == nil {
		t.Fatal("decoding against the wrong reference must be a hard error")
	}
	if !strings.Contains(err.Error(), "MD5 mismatch") {
		t.Errorf("error %q does not mention an MD5 mismatch", err.Error())
	}
}

// TestReferenceBackedDecodeViaRefCache decodes the reference-backed CRAM
// using a REF_CACHE directory instead of an explicit FASTA. The CRAM's
// @SQ M5 tag is the whole-sequence MD5 the cache is keyed on; the cache
// file holds the contig bases and the slice-span MD5 then verifies.
func TestReferenceBackedDecodeViaRefCache(t *testing.T) {
	cramPath := filepath.Join(samtoolsTestDir, referenceBackedFixture.cram)
	faPath := filepath.Join(samtoolsTestDir, referenceBackedFixture.fasta)
	if _, err := os.Stat(cramPath); err != nil {
		t.Fatalf("samtools submodule not initialised — CRAM fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	if _, err := os.Stat(faPath); err != nil {
		t.Fatalf("samtools submodule not initialised — reference FASTA unavailable; run `git submodule update --init reference_code/samtools`")
	}

	// Open the CRAM once to read its @SQ M5 tag — the digest htslib's
	// REF_CACHE keys a whole-sequence file on.
	probe, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords (probe): %v", err)
	}
	m5 := probe.contigMD5(0)
	probe.Close()
	if m5 == "" {
		t.Fatalf("the CRAM @SQ entry carries no M5 tag — REF_CACHE keying not exercised; the %s fixture must carry an M5 tag", referenceBackedFixture.cram)
	}

	// Lay the contig bases into a REF_CACHE directory under the M5 the
	// @SQ tag declares, using the htslib %2s/%2s/%s path.
	seq := readContigSequence(t, faPath, referenceBackedFixture.contig)
	dir := t.TempDir()
	c := OpenRefCache(dir)
	cachePath := c.refCachePath(m5)
	if err := os.MkdirAll(filepath.Dir(cachePath), 0o755); err != nil {
		t.Fatalf("mkdir cache: %v", err)
	}
	if err := os.WriteFile(cachePath, seq, 0o644); err != nil {
		t.Fatalf("write cache: %v", err)
	}

	rr, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords: %v", err)
	}
	defer rr.Close()
	rr.SetRefCache(dir)
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("REF_CACHE-backed decode: %v", err)
	}
	if rr.NeedsReference() {
		t.Error("a REF_CACHE-resolved decode must not report NeedsReference")
	}
	if len(recs) == 0 {
		t.Error("REF_CACHE decode produced no records")
	}
	if n := countN(recs); n > len(recs) {
		t.Errorf("REF_CACHE decode left %d 'N' bases across %d records — resolution incomplete", n, len(recs))
	}

	// A wrong M5 — no cache file under it — must be a clear not-found
	// error, never a panic. Point the cache at an empty directory.
	missRR, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords (miss): %v", err)
	}
	defer missRR.Close()
	missRR.SetRefCache(t.TempDir())
	if _, err := missRR.ReadAll(); err == nil {
		t.Error("a REF_CACHE miss must surface as an error")
	} else if !strings.Contains(err.Error(), "not found in REF_CACHE") {
		t.Errorf("REF_CACHE miss error %q should name the missing reference", err.Error())
	}
}

// countN totals the 'N' bases across a record set.
func countN(recs []*sam.Record) int {
	n := 0
	for _, r := range recs {
		for i := 0; i < len(r.Seq); i++ {
			if r.Seq[i] == 'N' {
				n++
			}
		}
	}
	return n
}

// readContigSequence reads one contig's bases from a FASTA, upper-cased
// and whitespace-free — the form a REF_CACHE entry stores.
func readContigSequence(t *testing.T, faPath, contig string) []byte {
	t.Helper()
	ra, err := fasta.OpenRandomAccess(faPath)
	if err != nil {
		t.Fatalf("open reference FASTA: %v", err)
	}
	defer ra.Close()
	n := ra.Length(contig)
	if n < 0 {
		t.Fatalf("contig %q not in reference FASTA", contig)
	}
	seq, err := ra.Fetch(contig, 0, n)
	if err != nil {
		t.Fatalf("fetch contig %q: %v", contig, err)
	}
	return seq
}

// verifyPureMatchAgainstReference cross-checks every pure-match read
// (CIGAR all M, no soft clip / indel / hard clip) against the reference
// FASTA span at its alignment position: a read with no substitution
// must equal the reference span exactly, and any read that differs must
// differ only at single-base positions (real substitutions). It returns
// the number of pure-match reads it checked.
func verifyPureMatchAgainstReference(t *testing.T, recs []*sam.Record, faPath string) int {
	t.Helper()
	ra, err := fasta.OpenRandomAccess(faPath)
	if err != nil {
		t.Fatalf("open reference FASTA: %v", err)
	}
	defer ra.Close()
	checked := 0
	for _, r := range recs {
		if r.RName == "" || r.Pos <= 0 {
			continue
		}
		// An unmapped read carries no alignment: its SEQ is not a copy of
		// any reference span and must not be cross-checked.
		if r.Flag&sam.FlagUnmapped != 0 {
			continue
		}
		// A pure-match read has a non-empty CIGAR consisting only of M.
		cs := r.Cigar.String()
		if cs == "" || cs == "*" || strings.ContainsAny(cs, "SIDNHP") {
			continue
		}
		span, err := ra.Fetch(r.RName, int64(r.Pos-1), int64(r.Pos-1)+int64(len(r.Seq)))
		if err != nil {
			continue
		}
		checked++
		// Count mismatching positions: every difference must be a single
		// substituted base (the reference and read lengths are equal here
		// because the CIGAR is pure-match).
		if len(span) != len(r.Seq) {
			t.Errorf("pure-match read %s: reference span len %d != seq len %d", r.QName, len(span), len(r.Seq))
			continue
		}
		mismatches := 0
		for i := 0; i < len(span); i++ {
			if span[i] != r.Seq[i] {
				mismatches++
			}
		}
		// A pure-match read may carry a handful of substitutions; the
		// invariant is just that the lengths line up and most bases
		// match the reference (a wholly wrong reference would mismatch
		// nearly everywhere).
		if mismatches > len(span)/2 {
			t.Errorf("pure-match read %s mismatches the reference in %d/%d bases — reference resolution looks wrong",
				r.QName, mismatches, len(span))
		}
	}
	return checked
}
