package prinseq

// pngconvert.go ports the histogram->matrix data-shaping helpers of
// upstream prinseq-graphs.pl (convertOdToBinMatrix, convertToBoxValues,
// convertToBarValues, convertOdToStackBinMatrix). These produce the
// exact numeric series the plots draw; they are unit-tested for
// value-level parity independent of the rendered pixels.

import "sort"

// convertOdToBinMatrix ports prinseq-graphs.pl convertOdToBinMatrix
// (lines 917-956). It bins a value->count histogram into up to 100
// bins and returns (matrix, xmax, ymax). `min` is the inclusive
// lower bound (0 or 1); `nonice` suppresses the ymax rounding.
func convertOdToBinMatrix(data map[int]int, min int, nonice bool) ([]int, int, int) {
	xmax := 0
	for k := range data {
		if k > xmax {
			xmax = k
		}
	}
	bin := getBinVal(xmax)
	xmax = bin * 100
	xmin := min

	ymax := 0
	tmp := 0
	tmpbin := bin
	var matrix []int
	for i := xmin; i <= xmax; i++ {
		if v, ok := data[i]; ok {
			tmp += v
		}
		tmpbin--
		if tmpbin <= 0 {
			tmpbin = bin
			if tmp > ymax {
				ymax = tmp
			}
			matrix = append(matrix, tmp)
			tmp = 0
		}
	}

	if !nonice {
		ymax = niceCeil4(ymax)
	}
	return matrix, xmax, ymax
}

// boxRow is one boxplot column: position plus the five-number
// summary (min, p25, median, p75, max).
type boxRow struct {
	pos    int
	min    int
	p25    int
	median int
	p75    int
	max    int
}

// convertToBoxValues ports prinseq-graphs.pl convertToBoxValues
// (lines 1582-1597). It flattens a position->GDStat map into an
// ordered slice of boxRows and returns (matrix, ymax). `niceval`
// rounds ymax up to the next multiple (4 for the quality plots).
func convertToBoxValues(data map[int]GDStat, niceval int) ([]boxRow, int) {
	keys := make([]int, 0, len(data))
	for k := range data {
		keys = append(keys, k)
	}
	sort.Ints(keys)

	ymax := 0
	matrix := make([]boxRow, 0, len(keys))
	for _, k := range keys {
		st := data[k]
		matrix = append(matrix, boxRow{
			pos:    k,
			min:    st.Min,
			p25:    st.P25,
			median: st.Median,
			p75:    st.P75,
			max:    st.Max,
		})
		if st.Max > ymax {
			ymax = st.Max
		}
	}
	if niceval > 0 && ymax%niceval != 0 {
		ymax = (ymax/niceval + 1) * niceval
	}
	return matrix, ymax
}

// convertToBarValues ports prinseq-graphs.pl convertToBarValues
// (lines 1802-1829). It expands a sparse value->count map into a
// dense slice [start..xmax] and returns (matrix, xmax, ymax).
func convertToBarValues(data map[int]int, niceval, start int) ([]int, int, int) {
	xmax := 0
	for k := range data {
		if k > xmax {
			xmax = k
		}
	}
	if niceval > 0 && xmax%niceval != 0 {
		xmax = (xmax/niceval + 1) * niceval
	}

	ymax := 0
	matrix := make([]int, 0, xmax-start+1)
	for q := start; q <= xmax; q++ {
		v := data[q]
		if v > ymax {
			ymax = v
		}
		matrix = append(matrix, v)
	}
	ymax = niceCeil4(ymax)
	return matrix, xmax, ymax
}

// convertOdToStackBinMatrix ports prinseq-graphs.pl
// convertOdToStackBinMatrix (lines 1968-2015). It bins a
// value->(stack->count) table into up to 100 bins for `stacks`
// layers, returning (matrix[stack][bin], xmax, ymax). `max` (>0)
// fixes the xmax; otherwise it is derived from the data.
func convertOdToStackBinMatrix(data map[int]map[int]int, stacks, min, max int, nonice bool) ([][]int, int, int) {
	xmax := max
	if xmax == 0 {
		for k := range data {
			if k > xmax {
				xmax = k
			}
		}
	}
	bin := getBinVal(xmax)
	xmax = bin * 100
	xmin := min

	matrix := make([][]int, stacks)
	sums := make([]int, stacks)
	ymax := 0
	sum := 0
	tmpbin := bin
	for i := xmin; i <= xmax; i++ {
		if inner, ok := data[i]; ok {
			for s := 0; s < stacks; s++ {
				if v, ok := inner[s]; ok {
					sums[s] += v
					sum += v
				}
			}
		}
		tmpbin--
		if tmpbin <= 0 {
			tmpbin = bin
			if sum > ymax {
				ymax = sum
			}
			sum = 0
			for s := 0; s < stacks; s++ {
				matrix[s] = append(matrix[s], sums[s])
				sums[s] = 0
			}
		}
	}

	if !nonice {
		ymax = niceCeil4(ymax)
	}
	return matrix, xmax, ymax
}
