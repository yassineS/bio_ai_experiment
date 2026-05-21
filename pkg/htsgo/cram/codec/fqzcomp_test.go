package codec

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// loadFQZCorpus loads an htscodecs q-file and applies the same transform
// the fqzcomp_qual test tool's parse_lines does: it keeps only the first
// tab-separated column of each line, strips the newlines, and subtracts
// 33 (ASCII phred -> raw quality). The result is the binary quality
// buffer that fqz_decompress reconstructs — the decode oracle.
func loadFQZCorpus(t *testing.T, name string) ([]byte, bool) {
	t.Helper()
	raw, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", name))
	if err != nil {
		return nil, false
	}
	var out []byte
	for _, line := range bytes.Split(raw, []byte{'\n'}) {
		if tab := bytes.IndexByte(line, '\t'); tab >= 0 {
			line = line[:tab]
		}
		for _, c := range line {
			out = append(out, c-33)
		}
	}
	return out, true
}

// fqzQFiles is the set of q-files the htscodecs fqzcomp corpus ships
// vectors for; each has strategy levels 0..3.
var fqzQFiles = []string{"q4", "q8", "qvar", "q40+dir"}

// loadFQZSlice loads an htscodecs q-file and reproduces the fqz_slice
// the fqzcomp_qual test tool builds via count_lines + parse_lines: each
// line is one read; column 1 is the phred-encoded quality string, the
// optional column 2 is the read2 marker, the optional column 3 is the
// selector. It returns the binary quality buffer (column 1, minus 33)
// and the matching slice. The precompressed dat/fqzcomp vectors were
// produced from the full files (q40+dir carries a read2 column), so the
// columns must be honoured for byte-exact encode comparison.
func loadFQZSlice(t *testing.T, name string) ([]byte, *fqzSlice, bool) {
	t.Helper()
	raw, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", name))
	if err != nil {
		return nil, nil, false
	}
	atoiPrefix := func(b []byte) int {
		n, sign, started := 0, 1, false
		for i := 0; i < len(b); i++ {
			c := b[i]
			if !started && (c == '+' || c == '-') {
				if c == '-' {
					sign = -1
				}
				started = true
				continue
			}
			if c < '0' || c > '9' {
				break
			}
			n = n*10 + int(c-'0')
			started = true
		}
		return n * sign
	}
	isSpace := func(c byte) bool {
		return c == ' ' || c == '\t' || c == '\n' || c == '\r' ||
			c == '\v' || c == '\f'
	}

	var qual []byte
	var lens, flags []uint32
	start := 0
	for i := 0; i < len(raw); i++ {
		c := raw[i]
		if c == '\n' || c == ' ' || c == '\t' {
			lens = append(lens, uint32(i-start))
			r2, sel := 0, 0
			j := i
			for j < len(raw) && raw[j] != '\n' && isSpace(raw[j]) {
				j++
			}
			if j < len(raw) && raw[j] != '\n' {
				r2 = atoiPrefix(raw[j:])
			}
			for j < len(raw) && !isSpace(raw[j]) {
				j++
			}
			for j < len(raw) && raw[j] != '\n' && isSpace(raw[j]) {
				j++
			}
			if j < len(raw) && raw[j] != '\n' {
				sel = atoiPrefix(raw[j:])
			}
			for j < len(raw) && raw[j] != '\n' {
				j++
			}
			flags = append(flags, uint32(r2)*fqzFRead2|uint32(sel)<<16)
			i = j
			start = i + 1
		} else {
			qual = append(qual, c-33)
		}
	}
	s := &fqzSlice{
		numRecords: len(lens),
		length:     lens,
		flags:      flags,
	}
	return qual, s, true
}

// TestFQZComp_ComplianceVectors decodes every pre-compressed htscodecs
// fqzcomp vector (q{4,8,var,40+dir}.{0,1,2,3}) and asserts byte-for-byte
// against the expected raw quality buffer. This is the compliance
// oracle: it proves the decoder matches the reference C on-wire format.
func TestFQZComp_ComplianceVectors(t *testing.T) {
	ran := 0
	for _, qfile := range fqzQFiles {
		want, ok := loadFQZCorpus(t, qfile)
		if !ok {
			continue
		}
		for s := 0; s < 4; s++ {
			comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "fqzcomp",
				qfile+"."+itoa(s)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".s"+itoa(s), func(t *testing.T) {
				got, err := FQZCompDecode(comp)
				if err != nil {
					t.Fatalf("decode: %v", err)
				}
				if !bytes.Equal(got, want) {
					t.Fatalf("decoded %d bytes, want %d; first mismatch at %d",
						len(got), len(want), firstDiff(got, want))
				}
			})
		}
	}
	if ran == 0 {
		t.Skip("htscodecs submodule not initialised — fqzcomp vectors unavailable")
	}
}

// TestFQZComp_EncodeMatchesHTScodecs checks the encoder reproduces the
// htscodecs fqzcomp vectors byte-for-byte. fqzcomp parameter selection
// is fully deterministic given (strategy, input, per-read lengths), so
// a byte-exact match is the strongest possible encoder test. Any vector
// that does not match is reported (not silently skipped) so the
// encoder's per-strategy status stays visible.
func TestFQZComp_EncodeMatchesHTScodecs(t *testing.T) {
	ran := 0
	for _, qfile := range fqzQFiles {
		qual, s, ok := loadFQZSlice(t, qfile)
		if !ok {
			continue
		}
		for st := 0; st < 4; st++ {
			comp, err := os.ReadFile(filepath.Join(htscodecsDir, "dat", "fqzcomp",
				qfile+"."+itoa(st)))
			if err != nil {
				continue
			}
			ran++
			t.Run(qfile+".s"+itoa(st), func(t *testing.T) {
				got, err := FQZCompEncode(qual, st, s)
				if err != nil {
					t.Fatalf("encode: %v", err)
				}
				if !bytes.Equal(got, comp) {
					t.Errorf("encoded %d bytes, htscodecs vector is %d; first mismatch at %d",
						len(got), len(comp), firstDiff(got, comp))
				}
			})
		}
	}
	if ran == 0 {
		t.Skip("htscodecs submodule not initialised — fqzcomp vectors unavailable")
	}
}

// TestFQZComp_RoundTrip exercises encode->decode as the identity across
// every q-file and strategy, with no dependency on the precompressed
// vectors. It proves the encoder produces a stream the decoder can
// reconstruct exactly even where it is not byte-identical to htscodecs.
func TestFQZComp_RoundTrip(t *testing.T) {
	ran := 0
	for _, qfile := range fqzQFiles {
		qual, s, ok := loadFQZSlice(t, qfile)
		if !ok {
			continue
		}
		for st := 0; st < 4; st++ {
			ran++
			t.Run(qfile+".s"+itoa(st), func(t *testing.T) {
				comp, err := FQZCompEncode(qual, st, s)
				if err != nil {
					t.Fatalf("encode: %v", err)
				}
				got, err := FQZCompDecode(comp)
				if err != nil {
					t.Fatalf("decode: %v", err)
				}
				if !bytes.Equal(got, qual) {
					t.Fatalf("round trip differs: %d bytes, want %d; first mismatch at %d",
						len(got), len(qual), firstDiff(got, qual))
				}
			})
		}
	}
	if ran == 0 {
		t.Skip("htscodecs submodule not initialised — fqzcomp corpus unavailable")
	}
}

// FuzzFQZCompDecode feeds arbitrary bytes to the decoder; it must never
// panic — only return data or an error. Seeds are kept small so short
// fuzz runs explore the header/parameter parser quickly: a small
// self-built stream plus a few degenerate inputs. (The full corpus
// vectors are large and the CTX_SIZE model set is expensive to build,
// which would otherwise dominate a short run.)
func FuzzFQZCompDecode(f *testing.F) {
	if seed, err := FQZCompEncode([]byte{30, 31, 32, 2, 30, 31}, 0, nil); err == nil {
		f.Add(seed)
	}
	if seed, err := FQZCompEncode([]byte{5, 5, 5, 5, 5, 5, 5, 5}, 3, nil); err == nil {
		f.Add(seed)
	}
	f.Add([]byte{})
	f.Add([]byte{0})
	f.Add([]byte{1, fqzVers, 0})
	f.Add(make([]byte, 32))
	f.Fuzz(func(t *testing.T, data []byte) {
		_, _ = FQZCompDecode(data)
	})
}
