// Package edlib is a pure-Go port of Martin Sosic's edlib bit-parallel
// edit-distance library (https://github.com/Martinsos/edlib), tracking the
// cut-down C variant vendored inside upstream bcftools at
// reference_code/bcftools/edlib.c. The port is byte-faithful to that
// reference where bit-level behaviour matters; line cross-references in the
// implementation files point at the original C.
//
// Scope of this slice: the API and the bit-parallel Myers semi-global
// algorithm (HW and SHW modes) that bcftools mpileup --indels-cns drives.
// The classical Needleman-Wunsch (NW) global mode is supported via a
// straight dynamic-programming fallback because the upstream cut-down C
// removed it (see edlib.c:152-172 — myersCalcEditDistanceNW is commented
// out). NW here is used mostly by tests; the production --indels-cns path
// uses HW only.
//
// Task.Path returns the alignment as the same opcode stream the original
// edlib emits: 0 = match, 1 = insertion to target (= deletion from query),
// 2 = insertion to query (= deletion from target), 3 = mismatch. Path is
// reconstructed by running a small banded DP between the discovered
// start/end locations rather than by the (commented-out) bit-parallel
// hint-recovery from the C original; that keeps the port simple and
// matches what the upstream consensus caller consumes downstream.
package edlib

import "errors"

// Mode selects how gaps at the query boundaries are scored.
//
// Maps to edlib.h:36-62 EdlibAlignMode.
type Mode int

const (
	// ModeNW is global alignment (Needleman-Wunsch): gaps at either end
	// of the query are penalised. Equivalent of EDLIB_MODE_NW.
	ModeNW Mode = iota
	// ModeSHW is prefix alignment: gaps at the end of the query (i.e.
	// trailing target suffix) are free. Equivalent of EDLIB_MODE_SHW.
	ModeSHW
	// ModeHW is infix alignment: gaps at both ends of the query (i.e.
	// surrounding target context) are free. Equivalent of EDLIB_MODE_HW.
	// This is the mode bcftools mpileup --indels-cns drives.
	ModeHW
)

// Task selects how much information the aligner returns. Each higher task
// implies more work; cheaper tasks should be preferred when possible.
//
// Maps to edlib.h:66-71 EdlibAlignTask.
type Task int

const (
	// TaskDistance asks only for the edit distance plus the end
	// locations of optimal alignments in target. Equivalent of
	// EDLIB_TASK_DISTANCE.
	TaskDistance Task = iota
	// TaskLoc additionally computes the start locations in target.
	// Equivalent of EDLIB_TASK_LOC.
	TaskLoc
	// TaskPath additionally reconstructs the alignment opcode sequence
	// for the first (start, end) pair. Equivalent of EDLIB_TASK_PATH.
	TaskPath
)

// Edit operations, matching edlib.h:84-87.
const (
	OpMatch    = 0 // Match
	OpInsert   = 1 // Insertion to target (= deletion from query)
	OpDelete   = 2 // Insertion to query (= deletion from target)
	OpMismatch = 3 // Mismatch
)

// Config controls a single alignment call. Maps to the
// EdlibAlignConfig struct at edlib.h:100-140. Equality extensions
// (additionalEqualities) are not yet ported — bcftools --indels-cns does
// not use them.
type Config struct {
	// K is the search-distance cap. If non-negative and the edit
	// distance exceeds K, Align returns EditDistance == -1. If negative
	// (the typical setting), the aligner doubles K internally until a
	// solution is found.
	K int
	// Mode selects the alignment boundary scoring.
	Mode Mode
	// Task selects how much output is requested.
	Task Task
}

// DefaultConfig returns the same default as edlibDefaultAlignConfig()
// (edlib.c:654-656): K=-1, ModeNW, TaskDistance.
func DefaultConfig() Config {
	return Config{K: -1, Mode: ModeNW, Task: TaskDistance}
}

// Result mirrors EdlibAlignResult at edlib.h:162-218.
//
// EndLocations and StartLocations are 0-based positions in target. For
// ModeHW the alignment spans target[StartLocations[i] : EndLocations[i]+1].
// EditDistance is -1 if no alignment with cost <= Config.K exists.
//
// Alignment is non-nil only when Config.Task == TaskPath; opcodes are the
// OpMatch / OpInsert / OpDelete / OpMismatch constants above and apply to
// the first (start, end) pair.
type Result struct {
	EditDistance   int
	EndLocations   []int
	StartLocations []int
	Alignment      []byte
	AlphabetLength int
}

// ErrInvalidMode is returned for Align with both sequences empty and an
// unknown Mode value. Real edlib reports this via EDLIB_STATUS_ERROR.
var ErrInvalidMode = errors.New("edlib: invalid alignment mode")

// Align computes the edit distance between query and target subject to
// cfg.Mode / cfg.Task / cfg.K. The returned Result mirrors the C
// EdlibAlignResult fields most callers need. Inputs are treated as byte
// sequences and compared by exact byte equality (no IUPAC wildcards yet).
//
// See edlib.c:100-243 edlibAlign for the reference flow.
func Align(query, target []byte, cfg Config) (Result, error) {
	res := Result{EditDistance: -1}

	qLen := len(query)
	tLen := len(target)

	// Transform sequences (edlib.c:598-639). We need a compact alphabet
	// so the Peq table stays small. Last symbol of the alphabet acts as
	// a wildcard with bits all set, matching edlib.c:273-279; we don't
	// expose IUPAC equalities here so the wildcard column is unused but
	// kept for API parity.
	q, t, alpha := transformSequences(query, target)
	res.AlphabetLength = alpha

	// Special-case: at least one input empty (edlib.c:120-140).
	if qLen == 0 || tLen == 0 {
		switch cfg.Mode {
		case ModeNW:
			res.EditDistance = maxInt(qLen, tLen)
			res.EndLocations = []int{tLen - 1}
		case ModeSHW, ModeHW:
			res.EditDistance = qLen
			res.EndLocations = []int{-1}
		default:
			return res, ErrInvalidMode
		}
		if cfg.Task == TaskLoc || cfg.Task == TaskPath {
			res.StartLocations = []int{0}
		}
		if cfg.Task == TaskPath {
			res.Alignment = emptyAlignment(qLen, tLen, cfg.Mode)
		}
		return res, nil
	}

	switch cfg.Mode {
	case ModeNW:
		// The bcftools cut-down edlib does not implement NW under the
		// bit-parallel path (edlib.c:167-172). For parity-relevant
		// callers we run a straight banded DP — slow but correct and
		// only used by tests for typical sequence lengths.
		return nwAlign(q, t, alpha, cfg)
	case ModeHW, ModeSHW:
		return semiGlobalAlign(q, t, alpha, cfg)
	default:
		return res, ErrInvalidMode
	}
}

// transformSequences mirrors edlib.c:598-639. It returns the transformed
// query and target plus the discovered alphabet size. The "alphabet" itself
// is implicit — we keep a 256-entry lookup table because Go has no MAX_UCHAR
// optimisation worth replicating.
func transformSequences(query, target []byte) (q, t []byte, alphabetSize int) {
	var letterIdx [256]byte
	var inAlphabet [256]bool

	q = make([]byte, len(query))
	t = make([]byte, len(target))
	next := byte(0)
	for i, c := range query {
		if !inAlphabet[c] {
			inAlphabet[c] = true
			letterIdx[c] = next
			next++
		}
		q[i] = letterIdx[c]
	}
	for i, c := range target {
		if !inAlphabet[c] {
			inAlphabet[c] = true
			letterIdx[c] = next
			next++
		}
		t[i] = letterIdx[c]
	}
	return q, t, int(next)
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// emptyAlignment fabricates an opcode list for the degenerate case where
// one input is empty. For NW we must walk all of target as deletions or
// all of query as insertions; for SHW/HW the free-gap modes have no
// opcodes to emit (the gap is not part of the alignment, mirroring
// edlibAlign's note at edlib.c:121-134 where alignment stays NULL).
func emptyAlignment(qLen, tLen int, mode Mode) []byte {
	if mode != ModeNW {
		// Free trailing/leading gap — no opcodes recorded.
		return nil
	}
	if qLen == 0 {
		out := make([]byte, tLen)
		for i := range out {
			out[i] = OpInsert // insertion-to-target = consume target only
		}
		return out
	}
	out := make([]byte, qLen)
	for i := range out {
		out[i] = OpDelete // insertion-to-query = consume query only
	}
	return out
}
