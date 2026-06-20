// Package difffuzz is a differential-fuzzing harness for the ported
// bioinformatics tools. For each "target" (a tool/subcommand plus a CLI
// template, e.g. `bcftools view`, `samtools flagstat`, `bgzip -d`) it
//
//  1. generates fuzzed inputs three ways — by MUTATING a valid seed fixture,
//     by STRUCTURED random generation of the relevant text format, and by
//     RAW random bytes (to exercise the never-crash / error-parity path);
//  2. runs BOTH our binary and the vendored upstream binary on each input with
//     identical arguments;
//  3. diffs stdout, stderr, AND exit code, classifying every divergence
//     (stdout-differs / stderr-differs / exitcode-differs / one-crashed /
//     both-crashed) after applying the SAME provenance normalization the
//     parity harness uses, so a reported divergence is a real behavioral
//     difference rather than a benign version stamp;
//  4. minimizes a divergence-triggering input by delta-debugging so the report
//     carries small reproducers;
//  5. optionally captures Go branch coverage of our tool over the run.
//
// All randomness flows from a caller-supplied seed so a run is reproducible.
//
// The harness reuses pipeline/internal/upstream for binary resolution and
// pipeline/runner for the normalization; it adds NO third-party dependencies.
package difffuzz

import (
	"fmt"
	"math/rand"
)

// Format names the input format a target consumes. It selects which structured
// generator the fuzzer uses and which seed fixtures it mutates.
type Format string

// The input formats the fuzzer understands. RawBytes targets receive only
// raw-random and mutation inputs (no structured generator).
const (
	FormatVCF      Format = "vcf"
	FormatSAM      Format = "sam"
	FormatBED      Format = "bed"
	FormatFASTA    Format = "fasta"
	FormatFASTQ    Format = "fastq"
	FormatGzip     Format = "gzip"     // bgzip/gzip member stream (binary)
	FormatRawBytes Format = "rawbytes" // no structured generator
)

// Origin records how a particular fuzzed input was produced; it is reported per
// divergence so a reader knows which generation strategy found the bug.
type Origin string

// The three generation strategies (the manuscript's a/b/c).
const (
	OriginMutation   Origin = "mutation"   // (a) mutated valid fixture
	OriginStructured Origin = "structured" // (b) structured random of the format
	OriginRaw        Origin = "raw"        // (c) raw random bytes
)

// Input is one fuzzed input: its bytes plus how it was produced.
type Input struct {
	Data   []byte
	Origin Origin
}

// Mutator generates fuzzed inputs for a single format from a seedable RNG.
// It is deterministic: the same seed and seed-fixture bytes always yield the
// same sequence of inputs.
type Mutator struct {
	rng    *rand.Rand
	format Format
	// seed is the valid fixture bytes that mutation operators perturb. May be
	// nil/empty (e.g. for a RawBytes target), in which case mutation falls back
	// to perturbing raw random bytes.
	seed []byte
}

// NewMutator returns a Mutator for format seeded by seedBytes (the bytes of a
// valid fixture) and the given RNG seed. seedBytes may be nil.
func NewMutator(format Format, seedBytes []byte, rngSeed int64) *Mutator {
	cp := make([]byte, len(seedBytes))
	copy(cp, seedBytes)
	return &Mutator{
		rng:    rand.New(rand.NewSource(rngSeed)),
		format: format,
		seed:   cp,
	}
}

// Next returns the next fuzzed input. The origin is chosen by rotating through
// the three strategies weighted toward mutation (the most fruitful) but always
// exercising structured and raw generation too. idx is the global input index
// (used only to vary the rotation deterministically).
func (m *Mutator) Next(idx int) Input {
	switch idx % 4 {
	case 0, 1:
		return Input{Data: m.mutate(), Origin: OriginMutation}
	case 2:
		if m.format == FormatRawBytes {
			return Input{Data: m.rawRandom(), Origin: OriginRaw}
		}
		return Input{Data: m.structured(), Origin: OriginStructured}
	default:
		return Input{Data: m.rawRandom(), Origin: OriginRaw}
	}
}

// --- (a) mutation operators -------------------------------------------------

// mutate applies one randomly chosen mutation operator to a copy of the seed
// bytes. When the seed is empty it perturbs a fresh raw-random buffer so the
// operator always has something to chew on.
func (m *Mutator) mutate() []byte {
	base := m.seed
	if len(base) == 0 {
		base = m.rawRandom()
	}
	b := make([]byte, len(base))
	copy(b, base)
	switch m.rng.Intn(7) {
	case 0:
		return m.bitFlip(b)
	case 1:
		return m.byteFlip(b)
	case 2:
		return m.truncate(b)
	case 3:
		return m.duplicateRecords(b)
	case 4:
		return m.reorderRecords(b)
	case 5:
		return m.boundaryValue(b)
	default:
		return m.insertBytes(b)
	}
}

// bitFlip flips a single random bit.
func (m *Mutator) bitFlip(b []byte) []byte {
	if len(b) == 0 {
		return b
	}
	i := m.rng.Intn(len(b))
	b[i] ^= 1 << uint(m.rng.Intn(8))
	return b
}

// byteFlip overwrites a small random run of bytes with random values.
func (m *Mutator) byteFlip(b []byte) []byte {
	if len(b) == 0 {
		return b
	}
	n := 1 + m.rng.Intn(4)
	start := m.rng.Intn(len(b))
	for k := 0; k < n && start+k < len(b); k++ {
		b[start+k] = byte(m.rng.Intn(256))
	}
	return b
}

// truncate cuts the buffer at a random offset (exercising premature-EOF paths).
func (m *Mutator) truncate(b []byte) []byte {
	if len(b) == 0 {
		return b
	}
	return b[:m.rng.Intn(len(b))]
}

// insertBytes splices a short random run into the buffer.
func (m *Mutator) insertBytes(b []byte) []byte {
	n := 1 + m.rng.Intn(8)
	ins := make([]byte, n)
	for i := range ins {
		ins[i] = byte(m.rng.Intn(256))
	}
	at := 0
	if len(b) > 0 {
		at = m.rng.Intn(len(b) + 1)
	}
	out := make([]byte, 0, len(b)+n)
	out = append(out, b[:at]...)
	out = append(out, ins...)
	out = append(out, b[at:]...)
	return out
}

// duplicateRecords duplicates a random line ("record" for text formats),
// stressing dedup/sort/merge logic.
func (m *Mutator) duplicateRecords(b []byte) []byte {
	lines := splitKeepEmpty(b)
	if len(lines) == 0 {
		return b
	}
	i := m.rng.Intn(len(lines))
	dup := append([][]byte{}, lines[:i+1]...)
	dup = append(dup, lines[i])
	dup = append(dup, lines[i+1:]...)
	return joinLines(dup)
}

// reorderRecords swaps two random lines, stressing order-dependence and
// sortedness assumptions.
func (m *Mutator) reorderRecords(b []byte) []byte {
	lines := splitKeepEmpty(b)
	if len(lines) < 2 {
		return b
	}
	i, j := m.rng.Intn(len(lines)), m.rng.Intn(len(lines))
	lines[i], lines[j] = lines[j], lines[i]
	return joinLines(lines)
}

// boundaryValue replaces a numeric-looking token on a random non-header data
// line with a boundary value (0, -1, huge, etc.), exercising coordinate and
// range edge cases.
func (m *Mutator) boundaryValue(b []byte) []byte {
	boundaries := []string{"0", "-1", "1", "2147483647", "2147483648",
		"9223372036854775807", "4294967296", "-2147483648", "1000000000000"}
	lines := splitKeepEmpty(b)
	if len(lines) == 0 {
		return b
	}
	// Prefer a data line (skip header markers common to the bio formats).
	order := m.rng.Perm(len(lines))
	for _, li := range order {
		ln := lines[li]
		if len(ln) == 0 || isHeaderLine(ln) {
			continue
		}
		fields := splitFields(ln)
		if len(fields) == 0 {
			continue
		}
		// Find a numeric field to overwrite.
		fOrder := m.rng.Perm(len(fields))
		for _, fi := range fOrder {
			if looksNumeric(fields[fi]) {
				fields[fi] = []byte(boundaries[m.rng.Intn(len(boundaries))])
				lines[li] = joinFields(fields)
				return joinLines(lines)
			}
		}
	}
	return joinLines(lines)
}

// --- (c) raw random ---------------------------------------------------------

// rawRandom returns a buffer of uniformly random bytes of random short length.
// This exercises the never-crash / clean-error-parity path: neither binary
// should panic; both should reject (or both accept) the same garbage.
func (m *Mutator) rawRandom() []byte {
	n := m.rng.Intn(256)
	b := make([]byte, n)
	for i := range b {
		b[i] = byte(m.rng.Intn(256))
	}
	return b
}

// --- (b) structured random --------------------------------------------------

// structured generates a syntactically plausible (but randomly populated)
// document in the target's format. These inputs parse far enough to reach the
// tool's logic, so divergences here tend to be semantic rather than parse-level.
func (m *Mutator) structured() []byte {
	switch m.format {
	case FormatVCF:
		return m.structuredVCF()
	case FormatSAM:
		return m.structuredSAM()
	case FormatBED:
		return m.structuredBED()
	case FormatFASTA:
		return m.structuredFASTA()
	case FormatFASTQ:
		return m.structuredFASTQ()
	default:
		return m.rawRandom()
	}
}

func (m *Mutator) structuredVCF() []byte {
	var b []byte
	b = append(b, "##fileformat=VCFv4.2\n"...)
	b = append(b, "##contig=<ID=chr1,length=100000>\n"...)
	b = append(b, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n"...)
	b = append(b, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"...)
	n := 1 + m.rng.Intn(8)
	bases := []byte("ACGT")
	for i := 0; i < n; i++ {
		pos := 1 + m.rng.Intn(100000)
		ref := bases[m.rng.Intn(4)]
		alt := bases[m.rng.Intn(4)]
		qual := m.rng.Intn(100)
		b = append(b, []byte(fmt.Sprintf("chr1\t%d\t.\t%c\t%c\t%d\tPASS\tDP=%d\n",
			pos, ref, alt, qual, m.rng.Intn(200)))...)
	}
	return b
}

func (m *Mutator) structuredSAM() []byte {
	var b []byte
	b = append(b, "@HD\tVN:1.6\tSO:coordinate\n"...)
	b = append(b, "@SQ\tSN:chr1\tLN:100000\n"...)
	n := 1 + m.rng.Intn(8)
	bases := []byte("ACGT")
	for i := 0; i < n; i++ {
		l := 5 + m.rng.Intn(20)
		seq := make([]byte, l)
		qual := make([]byte, l)
		for j := range seq {
			seq[j] = bases[m.rng.Intn(4)]
			qual[j] = byte('!' + m.rng.Intn(40))
		}
		pos := 1 + m.rng.Intn(100000)
		flag := []int{0, 4, 16, 99, 147}[m.rng.Intn(5)]
		b = append(b, []byte(fmt.Sprintf("r%d\t%d\tchr1\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\n",
			i, flag, pos, l, seq, qual))...)
	}
	return b
}

func (m *Mutator) structuredBED() []byte {
	var b []byte
	n := 1 + m.rng.Intn(10)
	for i := 0; i < n; i++ {
		start := m.rng.Intn(99000)
		end := start + 1 + m.rng.Intn(1000)
		b = append(b, []byte(fmt.Sprintf("chr1\t%d\t%d\tf%d\t%d\t%c\n",
			start, end, i, m.rng.Intn(1000), "+-"[m.rng.Intn(2)]))...)
	}
	return b
}

func (m *Mutator) structuredFASTA() []byte {
	var b []byte
	n := 1 + m.rng.Intn(4)
	bases := []byte("ACGTN")
	for i := 0; i < n; i++ {
		b = append(b, []byte(fmt.Sprintf(">seq%d desc\n", i))...)
		l := 10 + m.rng.Intn(60)
		seq := make([]byte, l)
		for j := range seq {
			seq[j] = bases[m.rng.Intn(len(bases))]
		}
		b = append(b, seq...)
		b = append(b, '\n')
	}
	return b
}

func (m *Mutator) structuredFASTQ() []byte {
	var b []byte
	n := 1 + m.rng.Intn(6)
	bases := []byte("ACGTN")
	for i := 0; i < n; i++ {
		l := 10 + m.rng.Intn(40)
		seq := make([]byte, l)
		qual := make([]byte, l)
		for j := range seq {
			seq[j] = bases[m.rng.Intn(len(bases))]
			qual[j] = byte('!' + m.rng.Intn(40))
		}
		b = append(b, []byte(fmt.Sprintf("@r%d\n%s\n+\n%s\n", i, seq, qual))...)
	}
	return b
}

// --- small line/field helpers (no bytes-package gymnastics needed) ----------

// splitKeepEmpty splits b on '\n', keeping a trailing empty element so that
// joinLines round-trips the trailing newline.
func splitKeepEmpty(b []byte) [][]byte {
	var out [][]byte
	start := 0
	for i := 0; i < len(b); i++ {
		if b[i] == '\n' {
			out = append(out, b[start:i])
			start = i + 1
		}
	}
	out = append(out, b[start:])
	return out
}

// joinLines rejoins lines with '\n' (inverse of splitKeepEmpty).
func joinLines(lines [][]byte) []byte {
	var out []byte
	for i, ln := range lines {
		if i > 0 {
			out = append(out, '\n')
		}
		out = append(out, ln...)
	}
	return out
}

func splitFields(ln []byte) [][]byte {
	var out [][]byte
	start := 0
	for i := 0; i < len(ln); i++ {
		if ln[i] == '\t' {
			out = append(out, ln[start:i])
			start = i + 1
		}
	}
	out = append(out, ln[start:])
	return out
}

func joinFields(fields [][]byte) []byte {
	var out []byte
	for i, f := range fields {
		if i > 0 {
			out = append(out, '\t')
		}
		out = append(out, f...)
	}
	return out
}

func isHeaderLine(ln []byte) bool {
	return len(ln) > 0 && (ln[0] == '#' || ln[0] == '@' || ln[0] == '>' || ln[0] == '+')
}

func looksNumeric(f []byte) bool {
	if len(f) == 0 {
		return false
	}
	for i, c := range f {
		if c >= '0' && c <= '9' {
			continue
		}
		if i == 0 && (c == '-' || c == '+') {
			continue
		}
		return false
	}
	return true
}
