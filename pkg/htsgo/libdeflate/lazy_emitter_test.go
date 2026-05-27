package libdeflate

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
	"testing"
)

// traceItem is the parsed form of a BTRACE LIT or MATCH event.
type traceItem struct {
	pos  uint32
	kind string // "LIT" or "MATCH"
	a    uint32 // byte value for LIT, len for MATCH
	b    uint32 // unused for LIT, offset for MATCH
}

// parseTrace extracts the LIT / MATCH event sequence from a BTRACE
// trace file in source order. The BTRACE harness emits these exactly
// where libdeflate's lazy parser calls deflate_choose_literal /
// deflate_choose_match, so the event order is the matchfinder's
// authoritative output order.
func parseTrace(t *testing.T, path string) []traceItem {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open trace: %v", err)
	}
	defer f.Close()

	var out []traceItem
	s := bufio.NewScanner(f)
	s.Buffer(make([]byte, 1<<20), 1<<20)
	for s.Scan() {
		line := s.Text()
		if !strings.HasPrefix(line, "BTRACE ") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}
		kind := fields[1]
		if kind != "LIT" && kind != "MATCH" {
			continue
		}
		kv := map[string]string{}
		for _, f := range fields[2:] {
			eq := strings.IndexByte(f, '=')
			if eq < 0 {
				continue
			}
			kv[f[:eq]] = f[eq+1:]
		}
		pos, err := parseUintAny(kv["pos"])
		if err != nil {
			t.Fatalf("parse pos in %q: %v", line, err)
		}
		ti := traceItem{pos: uint32(pos), kind: kind}
		switch kind {
		case "LIT":
			byteVal, err := parseUintAny(kv["byte"])
			if err != nil {
				t.Fatalf("parse byte in %q: %v", line, err)
			}
			ti.a = uint32(byteVal)
		case "MATCH":
			ln, err := parseUintAny(kv["len"])
			if err != nil {
				t.Fatalf("parse len in %q: %v", line, err)
			}
			off, err := parseUintAny(kv["offset"])
			if err != nil {
				t.Fatalf("parse offset in %q: %v", line, err)
			}
			ti.a = uint32(ln)
			ti.b = uint32(off)
		}
		out = append(out, ti)
	}
	if err := s.Err(); err != nil {
		t.Fatalf("scan trace: %v", err)
	}
	return out
}

// parseUintAny accepts decimal or 0x-prefixed hex.
func parseUintAny(s string) (uint64, error) {
	if strings.HasPrefix(s, "0x") || strings.HasPrefix(s, "0X") {
		return strconv.ParseUint(s[2:], 16, 64)
	}
	return strconv.ParseUint(s, 10, 64)
}

// runOracleStreamTest verifies that lazyEmit's item stream for the
// given fixture matches the BTRACE oracle event stream.
func runOracleStreamTest(t *testing.T, name string) {
	t.Helper()
	binPath := fmt.Sprintf("testdata/oracle/%s.bin", name)
	tracePath := fmt.Sprintf("testdata/oracle/%s.trace", name)

	data, err := os.ReadFile(binPath)
	if err != nil {
		t.Fatalf("read input: %v", err)
	}
	want := parseTrace(t, tracePath)

	got := lazyEmit(data, 6)

	// Replay each item, tracking its absolute position in `data`,
	// so we can compare against the trace's `pos=` field directly.
	var pos uint32
	if len(got) != len(want) {
		// Find the first divergence to make the failure useful.
		n := len(got)
		if len(want) < n {
			n = len(want)
		}
		for i := 0; i < n; i++ {
			gi := got[i]
			wi := want[i]
			if itemMatches(gi, wi, pos) {
				if gi.isLiteral() {
					pos++
				} else {
					pos += uint32(gi.length)
				}
				continue
			}
			t.Fatalf("%s: item %d diverges at pos=%d\n  got:  %s\n  want: %s",
				name, i, pos, formatItem(gi), formatTrace(wi))
		}
		t.Fatalf("%s: item count mismatch got=%d want=%d (matched first %d)",
			name, len(got), len(want), n)
	}

	for i, gi := range got {
		wi := want[i]
		if !itemMatches(gi, wi, pos) {
			t.Fatalf("%s: item %d diverges at pos=%d\n  got:  %s\n  want: %s",
				name, i, pos, formatItem(gi), formatTrace(wi))
		}
		if gi.isLiteral() {
			pos++
		} else {
			pos += uint32(gi.length)
		}
	}
	if pos != uint32(len(data)) {
		t.Fatalf("%s: total replayed bytes %d != input length %d", name, pos, len(data))
	}
}

func itemMatches(gi item, wi traceItem, pos uint32) bool {
	if wi.pos != pos {
		return false
	}
	if gi.isLiteral() {
		return wi.kind == "LIT" && uint32(gi.literal) == wi.a
	}
	return wi.kind == "MATCH" && uint32(gi.length) == wi.a && uint32(gi.offset) == wi.b
}

func formatItem(it item) string {
	if it.isLiteral() {
		return fmt.Sprintf("LIT byte=0x%02x", it.literal)
	}
	return fmt.Sprintf("MATCH len=%d offset=%d", it.length, it.offset)
}

func formatTrace(ti traceItem) string {
	if ti.kind == "LIT" {
		return fmt.Sprintf("LIT byte=0x%02x", ti.a)
	}
	return fmt.Sprintf("MATCH len=%d offset=%d", ti.a, ti.b)
}

func TestLazyEmitterOracleRandom64K(t *testing.T) {
	runOracleStreamTest(t, "random_64k")
}

func TestLazyEmitterOracleBGZFPayload(t *testing.T) {
	runOracleStreamTest(t, "bgzf_payload")
}
