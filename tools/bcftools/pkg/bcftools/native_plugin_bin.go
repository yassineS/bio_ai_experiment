// Shared float-binning helper for native bcftools plugins, porting bin.c
// (bin_init/bin_get_idx/bin_get_value/bin_get_size). It bins float values into
// predefined half-open intervals exactly as upstream does, including the
// min/max boundary fix-up performed when min!=max. The af-dist plugin uses it
// to mirror upstream's distribution tables byte-for-byte.
//
// bin_init accepts the boundary list either inline (a comma-separated string
// such as "0,0.25,0.5,0.75,1") or as a file path (no comma), one boundary value
// per line — exactly as upstream's hts_readlist decides (a comma indicates a
// list, otherwise a file). The file form is what `bcftools +af-dist -p
// bins.txt` / `-d bins.txt` uses.
package bcftools

import (
	"bufio"
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// bins holds the parsed boundary values, mirroring struct _bin_t in bin.c. All
// values are stored as 32-bit floats so the binary search and boundary fix-ups
// match upstream's C `float` arithmetic.
type bins struct {
	vals []float32
}

// binInit parses a boundary list and applies the min/max boundary fix-up used
// by upstream when min!=max. The list is taken inline when it contains a comma
// ("0,0.25,...") or read from a file (one boundary per line) otherwise, exactly
// as upstream's bin_init decides via hts_readlist. It ports bin_init() from
// bin.c.
func binInit(listDef string, min, max float32) (*bins, error) {
	// A comma indicates an inline list, otherwise a file (bin.c: is_file =
	// strchr(list_def,',') ? 0 : 1).
	var tokens []string
	if strings.Contains(listDef, ",") {
		tokens = strings.Split(listDef, ",")
	} else {
		var err error
		tokens, err = readBinFile(listDef)
		if err != nil {
			return nil, fmt.Errorf("bin: failed to read %s", listDef)
		}
	}

	b := &bins{}
	for _, p := range tokens {
		// Upstream calls strtod(list[i],&tmp) and errors if *tmp is non-empty,
		// i.e. the token must be a clean float with no trailing characters. It
		// does NOT trim whitespace for inline tokens; readBinFile reproduces the
		// verbatim-line behaviour of hts_readlist for the file form.
		f, err := strconv.ParseFloat(p, 64)
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

// readBinFile reads a bin-boundary file, returning one token per non-empty
// line, mirroring htslib's hts_readlist(is_file=1): each non-empty line becomes
// one verbatim list item (no trimming, no comment handling), and the file may
// be gzip/bgzf-compressed (handled transparently by iohelper). Empty lines are
// skipped, matching hts_readlist's `if (str.l == 0) continue`.
func readBinFile(path string) ([]string, error) {
	f, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	var out []string
	for sc.Scan() {
		// hts_readlist splits on '\n' only; a trailing '\r' (CRLF input) would
		// remain part of the line and make strtod fail, exactly as upstream. We
		// strip only the '\n' that bufio already removed, leaving any '\r'.
		line := sc.Text()
		if line == "" {
			continue
		}
		out = append(out, line)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}
