package prinseq

// pngplots.go contains the per-graph drawing routines ported from
// prinseq-graphs.pl (createAnnotBarPlot, createBoxPlot, createBarPlot,
// createStackBarPlot, createOddsRatioPlot). Geometry (margins, plot
// height, bar width) follows upstream so the plots are visually
// equivalent; the rasterisation differs (stdlib vs Cairo) so pixels
// are not byte-identical, by design.

import (
	"image/color"
	"sort"
)

// drawAnnotBarPlot ports createAnnotBarPlot (lines 984-1211): a
// single-series bar plot with an optional mean ± std annotation band.
// `annot` is nil to omit the annotation (used by tail/Ns plots).
func drawAnnotBarPlot(matrix []int, xmax, ymax int, annot *GDStat, xlab, ylab string, zero bool, add string) *canvas {
	bin := 1
	if xmax > 100 {
		bin = xmax / 100
		xmax = 100
	}
	if ymax <= 0 {
		ymax = 1
	}

	const (
		size   = 6
		offset = 20
		left   = 40
		bottom = 15
		top    = 20
		height = 200
	)
	zeroI := boolToInt(zero)
	width := left + offset*2 + (xmax+zeroI)*size
	totalH := bottom + top + offset*2 + height
	c := newCanvas(width, totalH)

	plotX := left + offset
	plotY := top + offset
	plotW := (xmax+zeroI)*size - 1

	// plot background
	c.fillRect(plotX, plotY, plotW, height, colBackplot)

	// x-axis ticks
	startI := 1
	if zero {
		startI = 0
	}
	for i := startI; i <= xmax; i++ {
		tx := left + offset + size/2 + size*i - boolToInt(!zero)*size - 1
		if i%5 == 0 && i > 1 && i < xmax {
			c.vLine(tx, plotY+height, plotY+height+3, colTick)
		} else {
			c.vLine(tx, plotY+height, plotY+height+1, colTick)
		}
	}
	// y-axis ticks
	c.hLine(plotX-3, plotX, plotY, colTick)
	c.hLine(plotX-3, plotX, plotY+height-1, colTick)

	// helplines
	for j := 1; j <= 3; j++ {
		hy := plotY + height*j/4 - 1
		c.hLine(plotX, plotX+(xmax+zeroI)*size, hy, colHelpline)
	}

	// y tick labels
	c.textRight(plotX-5, plotY-3, addCommas(ymax), colTick)
	c.textRight(plotX-5, plotY+height-3, "0", colTick)
	// x tick labels (every 10)
	for i := startI; i <= xmax; i++ {
		if i%10 == 0 && i > 1 && i < xmax {
			tx := left + offset + size/2 + size*i - boolToInt(!zero)*size - 1
			c.textCentered(tx, plotY+height+4, pngItoa(i*bin), colTick)
		}
	}

	// axis labels
	xl := xlab
	if bin > 1 {
		xl = xlab + " (Bin size: " + pngItoa(bin) + add + ")"
	}
	c.textCentered(plotX+(xmax+zeroI)*size/2, plotY+height+14, xl, colLabel)
	yl := ylab
	if bin > 1 {
		yl = ylab + " (per bin)"
	}
	c.textVertical(offset-4, plotY+height/2+textWidth(yl)/2, yl, colLabel)

	// annotation band (mean ± std)
	if annot != nil {
		drawAnnotBand(c, annot, bin, size, plotX, plotY, height)
	}

	// bars
	for pos := 0; pos < len(matrix); pos++ {
		v := matrix[pos]
		if v == 0 {
			continue
		}
		frac := float64(v) / float64(ymax)
		barH := int(frac * float64(height))
		if barH > 0 {
			c.fillRect(plotX+pos*size, plotY+height-barH, size-1, barH, colBar)
		}
	}
	return c
}

// drawAnnotBand draws the mean line and 1/2-SD shaded boxes/lines.
func drawAnnotBand(c *canvas, annot *GDStat, bin, size, plotX, plotY, height int) {
	mean := int(annot.Mean)
	std := int(annot.Std)
	std1l := mean - std
	std2l := mean - 2*std
	std1r := mean + std
	std2r := mean + 2*std
	if std1l == std1r {
		return
	}
	clampL := func(v int) int {
		if v < 0 {
			return 0
		}
		return v / bin
	}
	clampR := func(v int) int {
		if v/bin > 100 {
			return 100
		}
		return v / bin
	}
	s1l := clampL(std1l)
	s2l := clampL(std2l)
	s1r := clampR(std1r)
	s2r := clampR(std2r)

	c.fillRect(plotX+s2l*size+2, plotY, (s2r-s2l)*size, height, colStd2)
	c.fillRect(plotX+s1l*size+2, plotY, (s1r-s1l)*size, height, colStd1)

	meanX := plotX + (mean/bin)*size + 2
	c.vLine(meanX, plotY-5, plotY+height, colMean)
	if s1l > 0 {
		c.vLine(plotX+s1l*size+2, plotY-5, plotY+height, colStdLine)
	}
	if s2l > 0 {
		c.vLine(plotX+s2l*size+2, plotY-5, plotY+height, colStdLine)
	}
	if s1r < 100 {
		c.vLine(plotX+s1r*size+2, plotY-5, plotY+height, colStdLine)
	}
	if s2r < 100 {
		c.vLine(plotX+s2r*size+2, plotY-5, plotY+height, colStdLine)
	}
	c.textCentered(meanX, plotY-12, "M", colTick)
	if s1l > 0 {
		c.textCentered(plotX+s1l*size+2, plotY-12, "1SD", colTick)
	}
	if s2l > 0 {
		c.textCentered(plotX+s2l*size+3, plotY-12, "2SD", colTick)
	}
	if s1r < 100 {
		c.textCentered(plotX+s1r*size+2, plotY-12, "1SD", colTick)
	}
	if s2r < 100 {
		c.textCentered(plotX+s2r*size+3, plotY-12, "2SD", colTick)
	}
}

// drawBoxPlot ports createBoxPlot (lines 1599-1800). Each column is a
// box-and-whisker; the x-axis is always 100 wide (relative %) or
// binned bp positions.
func drawBoxPlot(matrix []boxRow, ymax int, xlab, ylab string, bin int, add string) *canvas {
	if bin <= 0 {
		bin = 1
	}
	if ymax <= 0 {
		ymax = 1
	}
	const (
		size   = 6
		offset = 20
		left   = 25
		bottom = 25
		top    = 5
		height = 300
		xmax   = 100
	)
	width := left + offset*2 + xmax*size
	totalH := bottom + offset*2 + height
	c := newCanvas(width, totalH)

	plotX := left + offset
	plotY := top + offset
	plotW := xmax*size - 1

	c.fillRect(plotX, plotY, plotW, height, colBackplot)

	// legend
	lx := plotX + size*50
	legend := []struct {
		col   color.RGBA
		label string
	}{
		{colWhisker, "Min/Max value"},
		{colBox, "25th to 75th percentile"},
		{colMedian, "Median"},
	}
	for _, l := range legend {
		c.fillRect(lx, top+5, 10, 10, l.col)
		lx += 15
		c.text(lx, top+5, l.label, colTick)
		lx += textWidth(l.label) + 15
	}

	// x ticks (every 10)
	for i := 1; i <= 9; i++ {
		tx := plotX + size/2 + size*10*i - size - 1
		c.vLine(tx, plotY+height, plotY+height+3, colTick)
	}
	// y ticks (4 divisions)
	for j := 0; j <= 4; j++ {
		ty := plotY + height*j/4
		c.hLine(plotX-3, plotX, ty, colTick)
	}
	// helplines
	for j := 1; j <= 3; j++ {
		hy := plotY + height*j/4 - 1
		c.hLine(plotX, plotX+xmax*size, hy, colHelpline)
	}

	// x tick labels
	for i := 1; i <= 9; i++ {
		tx := plotX + size/2 + size*10*i - size - 1
		c.textCentered(tx, plotY+height+4, pngItoa(i*10*bin), colTick)
	}
	// y tick labels
	for j := 0; j <= 4; j++ {
		ty := plotY + height*(4-j)/4
		c.textRight(plotX-5, ty-3, addCommas(ymax*j/4), colTick)
	}

	// axis labels
	xl := xlab
	if bin > 1 {
		xl = xlab + " (Bin size: " + pngItoa(bin) + add + ")"
	}
	c.textCentered(plotX+xmax*size/2, plotY+height+14, xl, colLabel)
	c.textVertical(offset-4, plotY+height/2+textWidth(ylab)/2, ylab, colLabel)

	// boxes
	factor := float64(height) / float64(ymax)
	for _, v := range matrix {
		bx := plotX + size*v.pos
		// whiskers
		if v.min != v.p25 {
			y := plotY + height - int(float64(v.min)*factor) - 1
			c.hLine(bx+1, bx+size-2, y, colWhisker)
		}
		if v.p75 != v.max {
			y := plotY + height - int(float64(v.max)*factor)
			c.hLine(bx+1, bx+size-2, y, colWhisker)
		}
		// whisker stems
		if v.min != v.p25 {
			c.vLine(bx+size/2-1, plotY+height-int(float64(v.p25)*factor), plotY+height-int(float64(v.min)*factor), colWhisker)
		}
		if v.p75 != v.max {
			c.vLine(bx+size/2-1, plotY+height-int(float64(v.max)*factor), plotY+height-int(float64(v.p75)*factor)-1, colWhisker)
		}
		// box
		if v.p25 != v.median || v.p75 != v.median {
			boxTop := plotY + height - int(float64(v.p75)*factor)
			boxH := int(float64(v.p75-v.p25) * factor)
			c.fillRect(bx, boxTop, size-1, boxH, colWhisker)
			drawRectOutline(c, bx, boxTop, size-2, boxH, colBox)
		} else {
			y := plotY + height - int(float64(v.median)*factor)
			c.hLine(bx, bx+size-1, y, colBox)
		}
		// median
		my := plotY + height - int(float64(v.median)*factor)
		c.hLine(bx+1, bx+size-2, my, colMedian)
	}
	return c
}

// drawBarPlot ports createBarPlot (lines 1831-1966): a simple
// single-series bar plot (no annotation), used for the per-read mean
// quality histogram.
func drawBarPlot(matrix []int, xmax, ymax int, xlab, ylab string) *canvas {
	if ymax <= 0 {
		ymax = 1
	}
	size := 10
	if xmax > 100 {
		size = 3
	} else if xmax > 50 {
		size = 6
	}
	const (
		offset = 20
		left   = 25
		bottom = 15
		top    = 0
		height = 200
	)
	width := left + offset*2 + xmax*size
	totalH := bottom + offset*2 + height
	c := newCanvas(width, totalH)

	plotX := left + offset
	plotY := top + offset
	plotW := xmax*size - 1

	c.fillRect(plotX, plotY, plotW, height, colBackplot)

	for i := 1; i <= xmax; i++ {
		tx := plotX + size/2 + size*i - size - 1
		if i%5 == 0 && i > 1 && i < xmax {
			c.vLine(tx, plotY+height, plotY+height+3, colTick)
		} else {
			c.vLine(tx, plotY+height, plotY+height+1, colTick)
		}
	}
	c.hLine(plotX-3, plotX, plotY, colTick)
	c.hLine(plotX-3, plotX, plotY+height-1, colTick)
	for j := 1; j <= 3; j++ {
		hy := plotY + height*j/4 - 1
		c.hLine(plotX, plotX+xmax*size, hy, colHelpline)
	}

	for i := 1; i <= xmax; i++ {
		if i%5 == 0 && i > 1 && i < xmax {
			tx := plotX + size/2 + size*i - size - 1
			c.textCentered(tx, plotY+height+4, pngItoa(i), colTick)
		}
	}
	c.textRight(plotX-5, plotY-3, addCommas(ymax), colTick)
	c.textRight(plotX-5, plotY+height-3, "0", colTick)

	c.textCentered(plotX+xmax*size/2, plotY+height+14, xlab, colLabel)
	c.textVertical(offset-4, plotY+height/2+textWidth(ylab)/2, ylab, colLabel)

	// bars: matrix index i corresponds to value (start+i)=1+i; pos maps
	// to (pos+1)*size, matching upstream's "$pos+1" offset for zero=0.
	for pos := 0; pos < xmax; pos++ {
		idx := pos + 1
		if idx >= len(matrix) {
			break
		}
		v := matrix[idx]
		if v == 0 {
			continue
		}
		frac := float64(v) / float64(ymax)
		barH := int(frac * float64(height))
		if barH > 0 {
			c.fillRect(plotX+idx*size, plotY+height-barH, size-1, barH, colBar)
		}
	}
	return c
}

// drawStackBarPlot ports createStackBarPlot (lines 2017-2183): a
// stacked bar plot over up to five duplicate-level layers with a
// colour legend.
func drawStackBarPlot(matrix [][]int, xmax, ymax int, xlab, ylab, add string) *canvas {
	bin := 1
	if xmax > 100 {
		bin = xmax / 100
		xmax = 100
	}
	if ymax <= 0 {
		ymax = 1
	}
	const (
		size   = 6
		offset = 20
		left   = 40
		bottom = 15
		top    = 20
		height = 200
	)
	width := left + offset*2 + xmax*size
	totalH := bottom + top + offset*2 + height
	c := newCanvas(width, totalH)

	plotX := left + offset
	plotY := top + offset
	plotW := xmax*size - 1

	c.fillRect(plotX, plotY, plotW, height, colBackplot)

	// legend (right-aligned, drawn in reverse)
	x := plotX + size*100 - 5
	for i := len(stackLegend) - 1; i >= 0; i-- {
		x -= textWidth(stackLegend[i])
		c.text(x, top+5, stackLegend[i], colTick)
		x -= 15
		c.fillRect(x, top+5, 10, 10, stackColors[i])
		x -= 15
	}

	for i := 1; i <= xmax; i++ {
		tx := plotX + size/2 + size*i - size - 1
		if i%5 == 0 && i > 1 && i < xmax {
			c.vLine(tx, plotY+height, plotY+height+3, colTick)
		} else {
			c.vLine(tx, plotY+height, plotY+height+1, colTick)
		}
	}
	c.hLine(plotX-3, plotX, plotY, colTick)
	c.hLine(plotX-3, plotX, plotY+height-1, colTick)
	for j := 1; j <= 3; j++ {
		hy := plotY + height*j/4 - 1
		c.hLine(plotX, plotX+xmax*size, hy, colHelpline)
	}

	for i := 1; i <= xmax; i++ {
		if i%10 == 0 && i > 1 && i < xmax {
			tx := plotX + size/2 + size*i - size - 1
			c.textCentered(tx, plotY+height+4, pngItoa(i*bin), colTick)
		}
	}
	c.textRight(plotX-5, plotY-3, addCommas(ymax), colTick)
	c.textRight(plotX-5, plotY+height-3, "0", colTick)

	xl := xlab
	if bin > 1 {
		xl = xlab + " (Bin size: " + pngItoa(bin) + add + ")"
	}
	c.textCentered(plotX+xmax*size/2, plotY+height+14, xl, colLabel)
	yl := ylab
	if bin > 1 {
		yl = ylab + " (per bin)"
	}
	c.textVertical(offset-4, plotY+height/2+textWidth(yl)/2, yl, colLabel)

	// stacked bars
	stacks := len(matrix)
	for pos := 0; pos < xmax; pos++ {
		acc := 0.0
		for s := 0; s < stacks; s++ {
			if pos >= len(matrix[s]) {
				continue
			}
			v := matrix[s][pos]
			if v == 0 {
				continue
			}
			cur := float64(v) / float64(ymax)
			y := plotY + height - int(acc*float64(height))
			barH := int(cur * float64(height))
			if barH > 0 {
				c.fillRect(plotX+pos*size, y-barH, size-1, barH, stackColors[s])
			}
			acc += cur
		}
	}
	return c
}

// drawOddsRatioPlot ports createOddsRatioPlot (lines 1427-1580): a
// scatter of dinucleotide odds ratios on a fixed 0.5..1.5 y-scale.
func drawOddsRatioPlot(odds map[string]float64) *canvas {
	yvalues := []float64{0.5, 0.78, 1.00, 1.23, 1.5}
	keys := make([]string, 0, len(odds))
	for k := range odds {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	num := len(keys)

	const (
		size   = 40
		offset = 20
		left   = 35
		right  = 90
		bottom = 20
		top    = 0
		height = 100
	)
	plotW := size * 10
	width := left + offset*2 + plotW + right
	totalH := bottom + offset*2 + height
	c := newCanvas(width, totalH)

	plotX := left + offset
	plotY := top + offset

	c.fillRect(plotX, plotY, plotW, height, colBackplot)

	// right-side over/under marks
	h77 := scaleH(0.77/2, height)
	h78 := scaleH(0.78/2, height)
	c.fillRect(plotX+plotW+8, plotY, 3, h77, color.RGBA{255, 127, 127, 153})
	c.fillRect(plotX+plotW+8, plotY+height-h78, 3, h78, color.RGBA{255, 127, 127, 153})

	// x ticks
	for i := 0; i < num; i++ {
		tx := plotX + size/2 + i*size
		c.vLine(tx, plotY+height, plotY+height+3, colTick)
	}
	// y ticks + helplines
	for _, yv := range yvalues {
		ty := plotY + height - int(yv/2*float64(height))
		c.hLine(plotX-3, plotX, ty, colTick)
		alpha := uint8(77)
		if yv == 0.5 || yv == 1.00 || yv == 1.50 {
			alpha = 26
		}
		c.hLine(plotX, plotX+plotW, ty, color.RGBA{0, 0, 0, alpha})
	}

	// x tick labels (split dinucleotide pairs with "/")
	for i, k := range keys {
		label := splitPairs(k)
		tx := plotX + size/2 + size*i
		c.textCentered(tx-1, plotY+height+4, label, colTick)
	}
	// y tick labels
	for _, yv := range yvalues {
		ty := plotY + height - int(yv/2*float64(height))
		c.textRight(plotX-5, ty-3, format2(yv), colTick)
	}
	// right-side labels
	c.text(plotX+plotW+15, plotY+height-scaleH(1.6/2, height)-3, "Over-represented", colLabel)
	c.text(plotX+plotW+15, plotY+height-scaleH(0.4/2, height)-3, "Under-represented", colLabel)

	// axis labels
	c.textCentered(plotX+plotW/2, plotY+height+14, "Dinucleotide", colLabel)
	c.textVertical(offset-4, plotY+height/2+textWidth("Odds ratio")/2, "Odds ratio", colLabel)

	// dots
	for i, k := range keys {
		v := odds[k]
		col := colOddsUnder
		if v > 1.23 || v < 0.78 {
			col = colOddsOver
		}
		cx := plotX + size/2 + size*i
		cy := plotY + height - int(v/2*float64(height))
		c.fillCircle(cx, cy, 5, col)
	}
	return c
}

// drawRectOutline strokes a rectangle border (1px).
func drawRectOutline(c *canvas, x, y, w, h int, col color.RGBA) {
	c.hLine(x, x+w, y, col)
	c.hLine(x, x+w, y+h, col)
	c.vLine(x, y, y+h, col)
	c.vLine(x+w, y, y+h, col)
}

// scaleH returns int(frac*height), evaluated at runtime to avoid Go's
// untyped-constant float->int conversion restriction.
func scaleH(frac float64, height int) int {
	return int(frac * float64(height))
}

func boolToInt(b bool) int {
	if b {
		return 1
	}
	return 0
}

// splitPairs turns a 4-char dinucleotide pair key (e.g. "AATT") into
// "AA/TT", matching upstream `join("/",(m/../g))`.
func splitPairs(s string) string {
	var parts []string
	for i := 0; i+1 < len(s); i += 2 {
		parts = append(parts, s[i:i+2])
	}
	out := ""
	for i, p := range parts {
		if i > 0 {
			out += "/"
		}
		out += p
	}
	if out == "" {
		return s
	}
	return out
}
