// Shared float-binning helper for native bcftools plugins, porting bin.c
// (bin_init/bin_get_idx/bin_get_value/bin_get_size). It bins float values into
// predefined half-open intervals exactly as upstream does, including the
// min/max boundary fix-up performed when min!=max. The af-dist plugin uses it
// to mirror upstream's distribution tables byte-for-byte.
package bcftools

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// bins holds the parsed boundary values, mirroring struct _bin_t in bin.c. All
// values are stored as 32-bit floats so the binary search and boundary fix-ups
// match upstream's C `float` arithmetic.
type bins struct {
	vals []float32
}

// binInit parses a comma-separated boundary list and applies the min/max
// boundary fix-up used by upstream when min!=max. The list must contain at
// least one comma (the file-name form of bin_init is not used by af-dist's
// default bins and is reported as unsupported by the caller). It ports
// bin_init() from bin.c.
func binInit(listDef string, min, max float32) (*bins, error) {
	if !strings.Contains(listDef, ",") {
		return nil, fmt.Errorf("bin: reading bins from a file is not supported: %s", listDef)
	}
	parts := strings.Split(listDef, ",")
	b := &bins{}
	for _, p := range parts {
		f, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
		if err != nil {
			return nil, fmt.Errorf("bin: could not parse %s: %s", listDef, p)
		}
		v := float32(f)
		if min != max && (v < min || v > max) {
			return nil, fmt.Errorf("bin: expected values from the interval [%f,%f], found %s", min, max, p)
		}
		b.vals = append(b.vals, v)
	}

	if min != max {
		if len(b.vals) <= 1 {
			return nil, fmt.Errorf("bin: at least two boundaries required")
		}
		maxErr := (b.vals[1] - b.vals[0]) * 1e-6
		if abs32(b.vals[0]-min) > maxErr {
			b.vals = append([]float32{min}, b.vals...)
		}
		if abs32(b.vals[len(b.vals)-1]-max) > maxErr {
			b.vals = append(b.vals, max)
		}
	}
	return b, nil
}

// size returns the number of boundaries; subtract 1 for the number of bins.
// It ports bin_get_size().
func (b *bins) size() int { return len(b.vals) }

// value returns the i-th boundary value. It ports bin_get_value().
func (b *bins) value(i int) float32 { return b.vals[i] }

// idx finds the bin index for value via the same binary search as
// bin_get_idx(): 0 <= idx <= size-2 for in-range values, with size-1 returned
// for values above the last boundary.
func (b *bins) idx(value float32) int {
	n := len(b.vals)
	if b.vals[n-1] < value {
		return n - 1
	}
	imin, imax := 0, n-2
	for imin < imax {
		i := (imin + imax) / 2
		if value < b.vals[i] {
			imax = i - 1
		} else if value > b.vals[i] {
			imin = i + 1
		} else {
			return i
		}
	}
	if b.vals[imax] <= value {
		return imax
	}
	if imin > 0 {
		return imin - 1
	}
	return 0
}

// abs32 is the 32-bit absolute value used by binInit's boundary fix-up.
func abs32(v float32) float32 { return float32(math.Abs(float64(v))) }
