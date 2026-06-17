package fastp

// FASTQ-similarity helper used to validate the genuinely heuristic /
// sampling-dependent fastp paths (SE adapter auto-detection) against the
// upstream binary.
//
// Deterministic transforms (poly-G/poly-X trimming, sliding-window quality
// trimming) are validated by strict byte-parity elsewhere in this package.
// Adapter auto-detection is a kmer/nucleotide-tree heuristic whose exact
// per-read trim length can differ from upstream by a base or two when the
// detected adapter string differs by a trailing base; for that path the
// contract is a documented SIMILARITY bound, computed here:
//
//   - LengthAgreement: fraction of reads whose trimmed length is IDENTICAL
//     between the two outputs (matched by read ID).
//   - MaxLenDelta:     the largest |len_go - len_up| over all matched reads.
//   - BaseIdentity:    fraction of overlapping bases (min of the two lengths,
//     per read) that are identical between the two outputs.
//
// These three numbers together prove the Go output is "sensible and
// comparable to upstream" without demanding byte-equality on a heuristic.

import (
	"strings"
	"testing"
)

// fastqSimilarity holds the aggregate similarity metrics between two FASTQ
// outputs, keyed by read ID.
type fastqSimilarity struct {
	Matched         int     // reads present (by ID) in BOTH outputs
	OnlyA           int     // reads present only in the first output
	OnlyB           int     // reads present only in the second output
	LengthAgreement float64 // fraction of matched reads with identical length
	MaxLenDelta     int     // largest |lenA-lenB| across matched reads
	BaseIdentity    float64 // identical overlapping bases / total overlapping bases
}

// parsedFastqRead is one record extracted from a FASTQ blob.
type parsedFastqRead struct {
	id  string
	seq string
}

// parseFastqBlob splits raw FASTQ bytes into records, returning a map keyed
// by the read ID (the first whitespace-delimited token of the header, sans
// the leading '@'). Quality and '+' lines are ignored — similarity is
// computed on the sequence. Malformed trailing groups are skipped.
func parseFastqBlob(b []byte) map[string]parsedFastqRead {
	out := map[string]parsedFastqRead{}
	lines := strings.Split(string(b), "\n")
	for i := 0; i+3 < len(lines) || (i+1 < len(lines) && strings.HasPrefix(lines[i], "@")); i += 4 {
		if i+1 >= len(lines) {
			break
		}
		header := lines[i]
		if !strings.HasPrefix(header, "@") {
			// Resync: not a record boundary (shouldn't happen for fastp out).
			break
		}
		id := strings.Fields(strings.TrimPrefix(header, "@"))
		var key string
		if len(id) > 0 {
			key = id[0]
		}
		out[key] = parsedFastqRead{id: key, seq: lines[i+1]}
	}
	return out
}

// compareFastqSimilarity computes the similarity metrics between two FASTQ
// blobs (typically the Go port's output and the upstream binary's output).
func compareFastqSimilarity(a, b []byte) fastqSimilarity {
	ma := parseFastqBlob(a)
	mb := parseFastqBlob(b)

	var sim fastqSimilarity
	identicalLen := 0
	var overlapTotal, overlapMatch int
	for id, ra := range ma {
		rb, ok := mb[id]
		if !ok {
			sim.OnlyA++
			continue
		}
		sim.Matched++
		la, lb := len(ra.seq), len(rb.seq)
		if la == lb {
			identicalLen++
		}
		delta := la - lb
		if delta < 0 {
			delta = -delta
		}
		if delta > sim.MaxLenDelta {
			sim.MaxLenDelta = delta
		}
		n := la
		if lb < n {
			n = lb
		}
		for i := 0; i < n; i++ {
			overlapTotal++
			if ra.seq[i] == rb.seq[i] {
				overlapMatch++
			}
		}
	}
	for id := range mb {
		if _, ok := ma[id]; !ok {
			sim.OnlyB++
		}
	}
	if sim.Matched > 0 {
		sim.LengthAgreement = float64(identicalLen) / float64(sim.Matched)
	}
	if overlapTotal > 0 {
		sim.BaseIdentity = float64(overlapMatch) / float64(overlapTotal)
	}
	return sim
}

// TestUnitFastqSimilarityHelper exercises the FASTQ-similarity helper itself
// with hand-built blobs so its behaviour is pinned WITHOUT needing the
// upstream binary. It checks the identical case, a single-base-shorter read,
// and an id-mismatch case.
func TestUnitFastqSimilarityHelper(t *testing.T) {
	identical := []byte("@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTTTTGGGG\n+\nIIIIIIII\n")
	sim := compareFastqSimilarity(identical, identical)
	if sim.Matched != 2 {
		t.Fatalf("matched=%d, want 2", sim.Matched)
	}
	if sim.LengthAgreement != 1.0 {
		t.Fatalf("LengthAgreement=%v, want 1.0", sim.LengthAgreement)
	}
	if sim.BaseIdentity != 1.0 {
		t.Fatalf("BaseIdentity=%v, want 1.0", sim.BaseIdentity)
	}
	if sim.MaxLenDelta != 0 {
		t.Fatalf("MaxLenDelta=%d, want 0", sim.MaxLenDelta)
	}

	// r2 is one base shorter in B and r1 has one differing base.
	a := []byte("@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nTTTTGGGG\n+\nIIIIIIII\n")
	b := []byte("@r1\nACGTACGA\n+\nIIIIIIII\n@r2\nTTTTGGG\n+\nIIIIIII\n")
	sim = compareFastqSimilarity(a, b)
	if sim.Matched != 2 {
		t.Fatalf("matched=%d, want 2", sim.Matched)
	}
	if sim.MaxLenDelta != 1 {
		t.Fatalf("MaxLenDelta=%d, want 1 (r2 short by 1)", sim.MaxLenDelta)
	}
	// 1 of 2 reads has identical length.
	if sim.LengthAgreement != 0.5 {
		t.Fatalf("LengthAgreement=%v, want 0.5", sim.LengthAgreement)
	}
	// Overlap bases: r1=8 (7 match), r2=min(8,7)=7 (all match) => 14/15.
	wantBI := 14.0 / 15.0
	if sim.BaseIdentity < wantBI-1e-9 || sim.BaseIdentity > wantBI+1e-9 {
		t.Fatalf("BaseIdentity=%v, want %v", sim.BaseIdentity, wantBI)
	}

	// id-mismatch: A has r3, B has r4.
	a3 := []byte("@r3\nACGT\n+\nIIII\n")
	b4 := []byte("@r4\nACGT\n+\nIIII\n")
	sim = compareFastqSimilarity(a3, b4)
	if sim.Matched != 0 || sim.OnlyA != 1 || sim.OnlyB != 1 {
		t.Fatalf("got matched=%d onlyA=%d onlyB=%d, want 0/1/1", sim.Matched, sim.OnlyA, sim.OnlyB)
	}
}
