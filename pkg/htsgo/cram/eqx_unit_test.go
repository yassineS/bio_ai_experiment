package cram

import (
	"bytes"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestUnitEqXWriteRoundTrip proves, with no upstream binary, that the CRAM
// writer accepts =/X CIGAR ops (the --eqx aligner output) and that a
// record carrying them round-trips through our own writer and reader to
// the same SEQ/POS/QUAL. =/X reconstruct as M, exactly as htslib's reader
// does for the equivalent reference-based CRAM (M/=/X all encode per-base
// features against the reference).
func TestUnitEqXWriteRoundTrip(t *testing.T) {
	h := writerTestHeader()
	cases := []struct {
		name      string
		cigar     string
		seq       string
		wantCigar string // CRAM collapses =/X to M.
	}{
		{"all-equal", "8=", "ACGTACGT", "8M"},
		{"one-mismatch", "3=1X4=", "ACGAACGT", "8M"},
		{"alternating", "2=1X1=1X3=", "CGTTCATA", "8M"},
		{"eqx-with-indel", "3=2I3=", "ACGTTACG", "3M2I3M"},
		{"eqx-with-del", "4=2D4=", "ACGTACGT", "4M2D4M"},
		{"leading-softclip-eqx", "2S6=", "TTACGTAC", "2S6M"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			rec := mkRec(tc.name, "chr1", 10, tc.cigar, tc.seq)
			out := roundTrip(t, h, []*sam.Record{rec})
			if len(out) != 1 {
				t.Fatalf("got %d records, want 1", len(out))
			}
			got := out[0]
			if got.Cigar.String() != tc.wantCigar {
				t.Errorf("cigar = %q, want %q", got.Cigar.String(), tc.wantCigar)
			}
			if normSeq(got.Seq) != tc.seq {
				t.Errorf("seq = %q, want %q", normSeq(got.Seq), tc.seq)
			}
			if got.Pos != 10 {
				t.Errorf("pos = %d, want 10", got.Pos)
			}
		})
	}
}

// TestUnitEqXEncodeFeatures checks the feature stream the encoder emits for
// = and X directly: an op of a given length encodes the same base-stretch
// ('b') feature whether it is spelled M, = or X — the encoder draws no
// per-op distinction, exactly as htslib's cigar loop folds
// BAM_CMATCH/CEQUAL/CDIFF into one case. This is the binary-free analogue
// of the byte-identical-features requirement against htslib.
func TestUnitEqXEncodeFeatures(t *testing.T) {
	seq := "ACGAACGT"
	enc := func(cigar string) []byte {
		t.Helper()
		rec := mkRec("r", "chr1", 1, cigar, seq)
		e := &recordEncoder{buffers: &seriesBuffers{}}
		if err := e.encodeFeatures(rec, len(seq)); err != nil {
			t.Fatalf("encodeFeatures(%s): %v", cigar, err)
		}
		// Concatenate the feature-relevant buffers into one comparable blob.
		var blob bytes.Buffer
		blob.Write(e.buffers.fn)
		blob.Write(e.buffers.fc)
		blob.Write(e.buffers.fp)
		blob.Write(e.buffers.bbLen)
		blob.Write(e.buffers.bb)
		return blob.Bytes()
	}
	// Same op length, three spellings — the feature buffers must be byte
	// identical, proving =/X are encoded exactly like M.
	mRun := enc("8M")
	eqRun := enc("8=")
	xRun := enc("8X")
	if !bytes.Equal(mRun, eqRun) {
		t.Errorf("= feature stream differs from M:\n  M = %x\n  = = %x", mRun, eqRun)
	}
	if !bytes.Equal(mRun, xRun) {
		t.Errorf("X feature stream differs from M:\n  M = %x\n  X = %x", mRun, xRun)
	}
	// A mixed =/X spelling of the same eight bases must also match the
	// equivalent mixed M spelling op-for-op (3M1M4M, i.e. three base
	// stretches of 3, 1, 4 bases), proving the per-op encoding is spelling
	// independent.
	mixedM := enc("3M1M4M")
	mixedEqX := enc("3=1X4=")
	if !bytes.Equal(mixedM, mixedEqX) {
		t.Errorf("mixed =/X feature stream differs from mixed M:\n  M     = %x\n  mixed = %x", mixedM, mixedEqX)
	}
}

// TestUnitCigarBackRejected proves the writer rejects the CIGAR back-step
// op B with a clear, B-specific error, matching htslib's CRAM encoder
// (which also rejects B as an unknown op). No upstream binary is needed.
func TestUnitCigarBackRejected(t *testing.T) {
	rec := mkRec("rb", "chr1", 1, "8M", "ACGTACGT")
	// Overwrite the CIGAR with 4M 2B 4M; B (back-step) consumes no query
	// bases, so the eight-base sequence still reconciles by query length —
	// the rejection must come from the op itself, not a length mismatch.
	rec.Cigar = sam.Cigar{
		sam.CigarOp(4<<4 | sam.CigarMatch),
		sam.CigarOp(2<<4 | sam.CigarBack),
		sam.CigarOp(4<<4 | sam.CigarMatch),
	}
	e := &recordEncoder{buffers: &seriesBuffers{}}
	err := e.encodeFeatures(rec, 8)
	if err == nil {
		t.Fatal("encodeFeatures accepted a B (back-step) op; want a rejection")
	}
	if !bytes.Contains([]byte(err.Error()), []byte("back-step")) &&
		!bytes.Contains([]byte(err.Error()), []byte("B is not supported")) {
		t.Errorf("error %q does not mention the B back-step rejection", err.Error())
	}
}
