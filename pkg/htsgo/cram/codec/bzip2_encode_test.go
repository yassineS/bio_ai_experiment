package codec

import (
	"bytes"
	"compress/bzip2"
	"io"
	"math/rand"
	"os/exec"
	"testing"
)

// bzPayloads is the shared set of round-trip inputs spanning the edge
// cases that exercise each pipeline stage: empty, tiny, highly
// repetitive (RLE1 + long MTF zero runs), random (worst case for BWT and
// Huffman), multi-block (> one 900k block) and binary.
func bzPayloads() map[string][]byte {
	rng := rand.New(rand.NewSource(1))
	random := make([]byte, 5000)
	rng.Read(random)
	big := make([]byte, 2_100_000) // forces >2 blocks at level 9 (900k)
	for i := range big {
		big[i] = byte(rng.Intn(7)) // small alphabet, some runs
	}
	binary := make([]byte, 4096)
	for i := range binary {
		binary[i] = byte(i * 7)
	}
	return map[string][]byte{
		"empty":      {},
		"single":     {0x42},
		"repeated":   bytes.Repeat([]byte{'A'}, 100000),
		"acgt":       bytes.Repeat([]byte("ACGT"), 50000),
		"text":       []byte("the quick brown fox jumps over the lazy dog. "),
		"random":     random,
		"multiblock": big,
		"binary":     binary,
		"allbytes": func() []byte {
			b := make([]byte, 256)
			for i := range b {
				b[i] = byte(i)
			}
			return b
		}(),
		"runs_only":     bytes.Repeat([]byte{0x00}, 1234567),
		"alternating":   bytes.Repeat([]byte{0xAA, 0x55}, 30000),
		"near_boundary": bytes.Repeat([]byte("X"), 900000-25),
	}
}

// TestBzip2EncodeStdlibRoundTrip is the primary correctness gate: our
// encoder output must decode under Go's compress/bzip2 back to the exact
// input, for every payload shape.
func TestBzip2EncodeStdlibRoundTrip(t *testing.T) {
	for name, in := range bzPayloads() {
		t.Run(name, func(t *testing.T) {
			enc, err := bzip2Encode(in, 9)
			if err != nil {
				t.Fatalf("bzip2Encode: %v", err)
			}
			if len(in) > 0 {
				// Stream must start with the BZh9 magic.
				if !bytes.HasPrefix(enc, []byte("BZh9")) {
					t.Fatalf("output missing BZh9 magic, got %x", enc[:4])
				}
			}
			got, err := io.ReadAll(bzip2.NewReader(bytes.NewReader(enc)))
			if err != nil {
				t.Fatalf("stdlib bzip2 decode: %v", err)
			}
			if !bytes.Equal(got, in) {
				t.Fatalf("round-trip mismatch: decoded %d bytes, want %d (first diff at %d)",
					len(got), len(in), firstDiff(got, in))
			}
		})
	}
}

// TestBzip2EncodeLevels round-trips a moderately sized payload at every
// block-size level, exercising single- and multi-block carving.
func TestBzip2EncodeLevels(t *testing.T) {
	in := bytes.Repeat([]byte("bioinformatics pipelines compress well-ish "), 8000)
	for level := 1; level <= 9; level++ {
		enc, err := bzip2Encode(in, level)
		if err != nil {
			t.Fatalf("level %d: encode: %v", level, err)
		}
		got, err := io.ReadAll(bzip2.NewReader(bytes.NewReader(enc)))
		if err != nil {
			t.Fatalf("level %d: decode: %v", level, err)
		}
		if !bytes.Equal(got, in) {
			t.Fatalf("level %d: round-trip mismatch", level)
		}
	}
}

// TestBzip2EncodeBadLevel pins the level validation.
func TestBzip2EncodeBadLevel(t *testing.T) {
	for _, lvl := range []int{-1, 0, 10, 100} {
		if _, err := bzip2Encode([]byte("x"), lvl); err == nil {
			t.Errorf("level %d: expected error", lvl)
		}
	}
}

// TestBzip2EncodeSystemBzip2 decompresses our output with the system
// bzip2 -d when available. It is a cross-tool gate: never skipped when
// the binary exists.
func TestBzip2EncodeSystemBzip2(t *testing.T) {
	path, err := exec.LookPath("bzip2")
	if err != nil {
		t.Fatalf("system bzip2 not available; install the bzip2 package to run the bzip2 cross-tool gate")
	}
	for name, in := range bzPayloads() {
		if len(in) == 0 {
			continue // empty stream has no block; bzip2 -d emits nothing — fine but trivial
		}
		t.Run(name, func(t *testing.T) {
			enc, err := bzip2Encode(in, 9)
			if err != nil {
				t.Fatalf("encode: %v", err)
			}
			cmd := exec.Command(path, "-d", "-c")
			cmd.Stdin = bytes.NewReader(enc)
			var out, errb bytes.Buffer
			cmd.Stdout = &out
			cmd.Stderr = &errb
			if err := cmd.Run(); err != nil {
				t.Fatalf("bzip2 -d failed: %v (stderr: %s)", err, errb.String())
			}
			if !bytes.Equal(out.Bytes(), in) {
				t.Fatalf("system bzip2 -d mismatch: got %d bytes, want %d (first diff %d)",
					out.Len(), len(in), firstDiff(out.Bytes(), in))
			}
		})
	}
}

// --- pipeline-stage unit tests ----------------------------------------

// TestRLE1RoundTrip checks rle1Encode reverses with a reference RLE1
// decoder for run-heavy and run-free inputs.
func TestRLE1RoundTrip(t *testing.T) {
	cases := [][]byte{
		{},
		{1},
		bytes.Repeat([]byte{7}, 3),   // below the 4-byte run threshold
		bytes.Repeat([]byte{7}, 4),   // exactly the threshold
		bytes.Repeat([]byte{7}, 255), // max single run
		bytes.Repeat([]byte{7}, 300), // run longer than 255 restarts
		[]byte("aaaabbbbcccc"),
		[]byte("no runs here at all!"),
	}
	for i, in := range cases {
		enc := rle1Encode(in)
		dec := rle1Decode(enc)
		if !bytes.Equal(dec, in) {
			t.Errorf("case %d: RLE1 round-trip failed: got %v want %v", i, dec, in)
		}
	}
}

// rle1Decode reverses rle1Encode: four equal bytes signal a run whose
// length is the following count byte plus four.
func rle1Decode(in []byte) []byte {
	var out []byte
	i := 0
	for i < len(in) {
		out = append(out, in[i])
		if i >= 3 && in[i] == in[i-1] && in[i-1] == in[i-2] && in[i-2] == in[i-3] {
			// We just appended the 4th of a run; next byte is the count.
			i++
			if i < len(in) {
				count := int(in[i])
				for k := 0; k < count; k++ {
					out = append(out, in[i-1])
				}
			}
		}
		i++
	}
	return out
}

// TestBWTForwardInverse checks the forward BWT and a reference inverse
// agree (the inverse reconstructs the original), and the origin pointer
// is correct.
func TestBWTForwardInverse(t *testing.T) {
	cases := [][]byte{
		[]byte("banana"),
		[]byte("a"),
		[]byte("abracadabra"),
		[]byte("mississippi"),
		bytes.Repeat([]byte{0}, 50),
		[]byte("the quick brown fox"),
	}
	for _, in := range cases {
		bwt, ptr := bwtEncode(in)
		got := bwtInverse(bwt, ptr)
		if !bytes.Equal(got, in) {
			t.Errorf("BWT inverse mismatch for %q: got %q", in, got)
		}
	}
}

// bwtInverse reconstructs the original string from a BWT last-column and
// origin pointer using the standard LF-mapping.
func bwtInverse(bwt []byte, origPtr int) []byte {
	n := len(bwt)
	if n == 0 {
		return nil
	}
	// Counting sort to build the first column and the next-array.
	var count [256]int
	for _, b := range bwt {
		count[b]++
	}
	var total [256]int
	sum := 0
	for i := 0; i < 256; i++ {
		total[i] = sum
		sum += count[i]
	}
	next := make([]int, n)
	var seen [256]int
	for i, b := range bwt {
		next[total[b]+seen[b]] = i
		seen[b]++
	}
	out := make([]byte, n)
	p := next[origPtr]
	for i := 0; i < n; i++ {
		out[i] = bwt[p]
		p = next[p]
	}
	return out
}

// TestMTFAndRLE2 checks that the MTF/RLE2 transform produces a stream the
// inverse can reverse to the original BWT bytes.
func TestMTFAndRLE2(t *testing.T) {
	cases := [][]byte{
		[]byte("banana"),
		bytes.Repeat([]byte{0x10}, 100), // long zero-run after MTF
		[]byte("ACGTACGTACGT"),
		[]byte{0, 1, 2, 3, 4, 5, 0, 0, 0, 5, 5, 1},
	}
	for _, bwt := range cases {
		var inUse [256]bool
		for _, b := range bwt {
			inUse[b] = true
		}
		var seqToUnseq []byte
		for i := 0; i < 256; i++ {
			if inUse[i] {
				seqToUnseq = append(seqToUnseq, byte(i))
			}
		}
		nInUse := len(seqToUnseq)
		mtfv := mtfAndRLE2(bwt, seqToUnseq, nInUse)
		got := mtfRLE2Inverse(mtfv, seqToUnseq, nInUse)
		if !bytes.Equal(got, bwt) {
			t.Errorf("MTF/RLE2 round-trip mismatch for %v: got %v", bwt, got)
		}
	}
}

// mtfRLE2Inverse reverses mtfAndRLE2.
func mtfRLE2Inverse(mtfv []uint16, seqToUnseq []byte, nInUse int) []byte {
	mtf := make([]byte, nInUse)
	copy(mtf, seqToUnseq)
	eob := uint16(nInUse + 1)
	var out []byte
	runBit := 0
	runVal := 0
	flush := func() {
		for k := 0; k < runVal; k++ {
			out = append(out, mtf[0])
		}
		runVal = 0
		runBit = 0
	}
	for _, s := range mtfv {
		switch {
		case s == bzRunA:
			runVal += 1 << runBit
			runBit++
		case s == bzRunB:
			runVal += 2 << runBit
			runBit++
		case s == eob:
			flush()
		default:
			flush()
			j := int(s) - 1 // literal rank
			v := mtf[j]
			copy(mtf[1:j+1], mtf[0:j])
			mtf[0] = v
			out = append(out, v)
		}
	}
	return out
}

// TestHuffmanLengthsValid checks that bzHuffmanLengths yields a valid
// prefix code (Kraft sum <= 1) within the length cap for assorted
// frequency distributions.
func TestHuffmanLengthsValid(t *testing.T) {
	cases := [][]int{
		{1, 1, 1, 1},
		{100, 1, 1, 1, 1, 1},
		{1, 2, 4, 8, 16, 32, 64, 128, 256},
		make([]int, 50), // all zero -> all weight 1
		{1000000, 1},
	}
	for ci, freq := range cases {
		lengths := make([]uint8, len(freq))
		bzHuffmanLengths(freq, lengths, bzMaxCodeLen)
		var kraft float64
		for i, l := range lengths {
			if l < 1 || l > bzMaxCodeLen {
				t.Errorf("case %d sym %d: length %d out of range", ci, i, l)
			}
			kraft += 1.0 / float64(uint64(1)<<l)
		}
		if kraft > 1.0000001 {
			t.Errorf("case %d: Kraft sum %g exceeds 1 (not a valid prefix code)", ci, kraft)
		}
	}
}

// TestCanonicalCodesPrefixFree checks the canonical code assignment
// produces distinct, prefix-free codes consistent with the lengths.
func TestCanonicalCodesPrefixFree(t *testing.T) {
	lengths := []uint8{2, 2, 3, 3, 3, 3}
	codes := make([]uint32, len(lengths))
	assignCanonicalCodes(lengths, codes)
	// Verify no code is a prefix of another.
	for i := range lengths {
		for j := range lengths {
			if i == j {
				continue
			}
			li, lj := lengths[i], lengths[j]
			if li > lj {
				continue
			}
			// Is code[i] (length li) a prefix of code[j]?
			if codes[j]>>(lj-li) == codes[i] {
				t.Errorf("code %d is a prefix of code %d", i, j)
			}
		}
	}
}

// TestBzCRC32 checks the bzip2 CRC against a known value: the CRC of the
// empty string is 0, and a round-trip via our encoder is validated by the
// stdlib decoder elsewhere (which checks the CRC).
func TestBzCRC32(t *testing.T) {
	if got := bzCRC32(nil); got != 0 {
		t.Errorf("CRC of empty = %08x, want 0", got)
	}
	// Determinism / sensitivity.
	a := bzCRC32([]byte("hello"))
	b := bzCRC32([]byte("hellp"))
	if a == b {
		t.Errorf("CRC collision on single-bit change")
	}
}
