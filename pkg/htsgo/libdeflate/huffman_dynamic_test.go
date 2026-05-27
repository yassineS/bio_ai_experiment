package libdeflate

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"testing"
)

// parseHuffmanLens extracts the BTRACE HUFFMAN_LEN events for the
// litlen and offset tables in source order. The trace harness emits
// these immediately after libdeflate's make_huffman_codes call for
// the dynamic block, so the BTRACE order matches the canonical
// (symbol-index) order.
func parseHuffmanLens(t *testing.T, path string) (litlen []uint8, offset []uint8) {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open trace: %v", err)
	}
	defer f.Close()

	litlen = make([]uint8, numLitlenSyms)
	offset = make([]uint8, numOffsetSyms)

	s := bufio.NewScanner(f)
	s.Buffer(make([]byte, 1<<20), 1<<20)
	// Skip the initial HUFFMAN_LEN events that record the static
	// Huffman codes (emitted by deflate_init_static_codes). Only the
	// HUFFMAN_LEN events emitted after BLOCK_BEGIN describe the
	// dynamic code we want to verify, and the apply_btrace patch
	// only emits a length entry when codes->lens[sym] != 0.
	seenBlockBegin := false
	for s.Scan() {
		line := s.Text()
		if strings.HasPrefix(line, "BTRACE BLOCK_BEGIN ") {
			seenBlockBegin = true
			continue
		}
		if !seenBlockBegin {
			continue
		}
		if !strings.HasPrefix(line, "BTRACE HUFFMAN_LEN ") {
			continue
		}
		fields := strings.Fields(line)
		kv := map[string]string{}
		for _, f := range fields[2:] {
			eq := strings.IndexByte(f, '=')
			if eq < 0 {
				continue
			}
			kv[f[:eq]] = f[eq+1:]
		}
		sym, err := parseUintAny(kv["sym"])
		if err != nil {
			t.Fatalf("parse sym in %q: %v", line, err)
		}
		clen, err := parseUintAny(kv["codelen"])
		if err != nil {
			t.Fatalf("parse codelen in %q: %v", line, err)
		}
		switch kv["table"] {
		case "litlen":
			litlen[sym] = uint8(clen)
		case "offset":
			offset[sym] = uint8(clen)
		}
	}
	if err := s.Err(); err != nil {
		t.Fatalf("scan trace: %v", err)
	}
	return litlen, offset
}

// runHuffmanOracleTest builds the dynamic Huffman code from the
// freqs implied by the matchfinder output for the named fixture
// (with the +1 EOB increment) and verifies the per-symbol code
// lengths byte-match the BTRACE HUFFMAN_LEN events.
func runHuffmanOracleTest(t *testing.T, name string) {
	t.Helper()
	binPath := fmt.Sprintf("testdata/oracle/%s.bin", name)
	tracePath := fmt.Sprintf("testdata/oracle/%s.trace", name)

	data, err := os.ReadFile(binPath)
	if err != nil {
		t.Fatalf("read input: %v", err)
	}
	blocks := lazyEmitBlocks(data, 6)
	if len(blocks) != 1 {
		t.Fatalf("expected 1 block, got %d", len(blocks))
	}
	blk := blocks[0]
	blk.freqs.litlen[endOfBlock]++
	d := buildDynamicCode(&blk.freqs)

	wantLitlen, wantOffset := parseHuffmanLens(t, tracePath)
	for sym := 0; sym < numLitlenSyms; sym++ {
		if d.litlenLens[sym] != wantLitlen[sym] {
			t.Fatalf("%s litlen sym=%d: got len=%d want=%d",
				name, sym, d.litlenLens[sym], wantLitlen[sym])
		}
	}
	for sym := 0; sym < numOffsetSyms; sym++ {
		if d.offsetLens[sym] != wantOffset[sym] {
			t.Fatalf("%s offset sym=%d: got len=%d want=%d",
				name, sym, d.offsetLens[sym], wantOffset[sym])
		}
	}
}

func TestDynamicHuffmanOracleRandom64K(t *testing.T) {
	runHuffmanOracleTest(t, "random_64k")
}

func TestDynamicHuffmanOracleBGZFPayload(t *testing.T) {
	runHuffmanOracleTest(t, "bgzf_payload")
}
