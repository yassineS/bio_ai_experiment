package mosdepth

import (
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// maxQuantizeBound mirrors Nim's high(int) sentinel that upstream mosdepth
// appends when the quantize spec ends in ':'. It marks the open-ended top
// bin whose label renders as "N:inf". On a 64-bit platform Nim's int is
// 64-bit so this matches math.MaxInt64; any observed depth is far below it,
// so it behaves as "+inf" for the purpose of bin lookup.
const maxQuantizeBound = math.MaxInt64

// ParseQuantize parses a ':'-separated quantize segment specification into
// the sorted slice of integer bin boundaries upstream mosdepth uses.
//
// It mirrors upstream's get_quantize_args exactly:
//   - a leading ':' (e.g. ":1:4") is treated as "0:1:4" — a 0 lower bound is
//     prepended.
//   - a trailing ':' (e.g. "1:4:") appends the high(int) sentinel, giving an
//     open-ended top bin "N:inf".
//   - the resulting integers are sorted ascending.
//
// An empty spec yields a nil slice (quantize disabled).
func ParseQuantize(spec string) ([]int, error) {
	if spec == "" {
		return nil, nil
	}
	a := spec
	if a[0] == ':' {
		a = "0" + a
	}
	if a[len(a)-1] == ':' {
		a = a + strconv.FormatInt(maxQuantizeBound, 10)
	}
	parts := strings.Split(a, ":")
	out := make([]int, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "" {
			return nil, fmt.Errorf("mosdepth: bad --quantize segment in %q (empty value)", spec)
		}
		v, err := strconv.Atoi(p)
		if err != nil {
			return nil, fmt.Errorf("mosdepth: bad --quantize value %q: %w", p, err)
		}
		out = append(out, v)
	}
	if len(out) < 2 {
		return nil, fmt.Errorf("mosdepth: --quantize needs at least two boundaries, got %q", spec)
	}
	sort.Ints(out)
	return out, nil
}

// quantizeLabels builds the per-bin label slice for the parsed quantize
// boundaries, honouring the MOSDEPTH_Q<i> environment overrides exactly like
// upstream's make_lookup. There is one label per bin, i.e. len(quants)-1
// labels. For bin i the default label is "quants[i]:quants[i+1]", except the
// open-ended top bin (quants[i+1] == high(int)) which renders as
// "quants[i]:inf". Setting MOSDEPTH_Q<i> overrides bin i's label verbatim.
func quantizeLabels(quants []int) []string {
	n := len(quants) - 1
	if n < 1 {
		return nil
	}
	labels := make([]string, n)
	for i := 0; i < n; i++ {
		if env := os.Getenv("MOSDEPTH_Q" + strconv.Itoa(i)); env != "" {
			labels[i] = env
			continue
		}
		if quants[i+1] == maxQuantizeBound {
			labels[i] = strconv.Itoa(quants[i]) + ":inf"
		} else {
			labels[i] = strconv.Itoa(quants[i]) + ":" + strconv.Itoa(quants[i+1])
		}
	}
	return labels
}

// quantizeBin returns the bin index that depth d falls into, mirroring
// upstream's linear_search. Bins are the half-open intervals [quants[i],
// quants[i+1]) for i in 0..len(quants)-2, so a valid index is always in
// range for the len(quants)-1 labels. A depth that lands in no bin — below
// quants[0], above quants[high], or exactly equal to the finite top boundary
// quants[high] (which no half-open bin covers) — returns -1. Such positions
// produce no quantized output line, leaving a gap in the BED, exactly as
// upstream does.
func quantizeBin(d int, quants []int) int {
	if d < quants[0] || d >= quants[len(quants)-1] {
		return -1
	}
	for i := 0; i+1 < len(quants); i++ {
		if d < quants[i+1] {
			return i
		}
	}
	// Unreachable: the guard above already excludes d >= quants[high].
	return -1
}

// emitQuantized walks the accumulator's per-base depth for one chromosome and
// writes collapsed quantized runs to w. Each contiguous run of positions that
// map to the same bin label is emitted as a single BED record
// "chrom\tstart\tend\tlabel". Positions whose depth maps to bin -1 (outside
// the quantize range) are skipped, matching upstream's gen_quantized.
func emitQuantized(w *bedGzWriter, chrom string, a *covAccum, quants []int, labels []string) error {
	var runStart = -1
	var runBin = -1
	var emitErr error
	flush := func(end int) {
		if emitErr != nil {
			return
		}
		if runStart >= 0 && runBin >= 0 && runBin < len(labels) {
			emitErr = w.writeBED(chrom, runStart, end, labels[runBin])
		}
	}
	a.emit(func(pos int, depth int32) {
		if emitErr != nil {
			return
		}
		b := quantizeBin(int(depth), quants)
		if b == runBin {
			return
		}
		flush(pos)
		runStart = pos
		runBin = b
	})
	if emitErr != nil {
		return emitErr
	}
	// Flush the trailing run out to the reference length.
	end := a.refLen
	if end <= runStart {
		end = runStart + 1
	}
	flush(end)
	return emitErr
}
