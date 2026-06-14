package prinseq

// pngcanvas.go is a minimal pure-stdlib drawing surface used by the
// PNG graph renderer. It wraps an *image.RGBA and offers the small
// set of primitives the prinseq graphs need: filled rectangles,
// axis-aligned and arbitrary lines, circles, and bitmap text. No
// third-party graphics library is used (none is sanctioned in this
// repo); anti-aliasing is intentionally omitted — upstream
// prinseq-graphs.pl itself draws most plot elements with
// antialiasing disabled.

import (
	"image"
	"image/color"
	"math"
)

// canvas is a drawable RGBA surface.
type canvas struct {
	img *image.RGBA
	w   int
	h   int
}

// newCanvas returns a w×h canvas filled white (the prinseq graphs
// use a white page background).
func newCanvas(w, h int) *canvas {
	if w < 1 {
		w = 1
	}
	if h < 1 {
		h = 1
	}
	img := image.NewRGBA(image.Rect(0, 0, w, h))
	c := &canvas{img: img, w: w, h: h}
	c.fillRect(0, 0, w, h, color.RGBA{255, 255, 255, 255})
	return c
}

// blend draws src over the existing pixel at (x,y) using src's alpha
// (simple source-over compositing).
func (c *canvas) blend(x, y int, src color.RGBA) {
	if x < 0 || y < 0 || x >= c.w || y >= c.h {
		return
	}
	if src.A == 255 {
		c.img.SetRGBA(x, y, src)
		return
	}
	if src.A == 0 {
		return
	}
	dst := c.img.RGBAAt(x, y)
	a := float64(src.A) / 255.0
	out := color.RGBA{
		R: uint8(float64(src.R)*a + float64(dst.R)*(1-a)),
		G: uint8(float64(src.G)*a + float64(dst.G)*(1-a)),
		B: uint8(float64(src.B)*a + float64(dst.B)*(1-a)),
		A: 255,
	}
	c.img.SetRGBA(x, y, out)
}

// fillRect fills the rectangle [x,x+w) × [y,y+h) with col.
func (c *canvas) fillRect(x, y, w, h int, col color.RGBA) {
	for j := y; j < y+h; j++ {
		for i := x; i < x+w; i++ {
			c.blend(i, j, col)
		}
	}
}

// hLine draws a horizontal line from x1 to x2 (inclusive) at y.
func (c *canvas) hLine(x1, x2, y int, col color.RGBA) {
	if x2 < x1 {
		x1, x2 = x2, x1
	}
	for x := x1; x <= x2; x++ {
		c.blend(x, y, col)
	}
}

// vLine draws a vertical line from y1 to y2 (inclusive) at x.
func (c *canvas) vLine(x, y1, y2 int, col color.RGBA) {
	if y2 < y1 {
		y1, y2 = y2, y1
	}
	for y := y1; y <= y2; y++ {
		c.blend(x, y, col)
	}
}

// fillCircle fills a disc of radius r centred at (cx,cy).
func (c *canvas) fillCircle(cx, cy, r int, col color.RGBA) {
	r2 := r * r
	for dy := -r; dy <= r; dy++ {
		for dx := -r; dx <= r; dx++ {
			if dx*dx+dy*dy <= r2 {
				c.blend(cx+dx, cy+dy, col)
			}
		}
	}
}

// text draws s with the 5x7 bitmap font, top-left at (x,y).
func (c *canvas) text(x, y int, s string, col color.RGBA) {
	cursor := x
	for _, r := range s {
		glyph, ok := font5x7[r]
		if !ok {
			// Unknown rune: leave a blank cell and advance.
			cursor += glyphW + glyphGapX
			continue
		}
		for row := 0; row < glyphH; row++ {
			bits := glyph[row]
			for col2 := 0; col2 < glyphW; col2++ {
				if bits&(1<<uint(glyphW-1-col2)) != 0 {
					c.blend(cursor+col2, y+row, col)
				}
			}
		}
		cursor += glyphW + glyphGapX
	}
}

// textCentered draws s centred horizontally at cx, with the top at y.
func (c *canvas) textCentered(cx, y int, s string, col color.RGBA) {
	c.text(cx-textWidth(s)/2, y, s, col)
}

// textRight draws s right-aligned so its right edge is at x.
func (c *canvas) textRight(x, y int, s string, col color.RGBA) {
	c.text(x-textWidth(s), y, s, col)
}

// textVertical draws s rotated 90° counter-clockwise (reading bottom
// to top), used for y-axis labels. (x,y) is the bottom-left corner
// of the resulting vertical strip.
func (c *canvas) textVertical(x, y int, s string, col color.RGBA) {
	cursor := y
	for _, r := range s {
		glyph, ok := font5x7[r]
		if !ok {
			cursor -= glyphW + glyphGapX
			continue
		}
		for row := 0; row < glyphH; row++ {
			bits := glyph[row]
			for col2 := 0; col2 < glyphW; col2++ {
				if bits&(1<<uint(glyphW-1-col2)) != 0 {
					// Rotate: glyph column -> upward, glyph row -> rightward.
					c.blend(x+row, cursor-col2, col)
				}
			}
		}
		cursor -= glyphW + glyphGapX
	}
}

// addCommas formats an integer with thousands separators, matching
// upstream prinseq-graphs.pl addCommas (line 602).
func addCommas(n int) string {
	neg := n < 0
	if neg {
		n = -n
	}
	s := pngItoa(n)
	if len(s) <= 3 {
		if neg {
			return "-" + s
		}
		return s
	}
	var out []byte
	pre := len(s) % 3
	if pre > 0 {
		out = append(out, s[:pre]...)
	}
	for i := pre; i < len(s); i += 3 {
		if len(out) > 0 {
			out = append(out, ',')
		}
		out = append(out, s[i:i+3]...)
	}
	if neg {
		return "-" + string(out)
	}
	return string(out)
}

func pngItoa(n int) string {
	if n == 0 {
		return "0"
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	return string(buf[i:])
}

// niceCeil4 rounds v up to the next multiple of 4, matching the
// upstream `sprintf("%d",($ymax/4)+1)*4 if($ymax % 4)` idiom.
func niceCeil4(v int) int {
	if v%4 != 0 {
		return (v/4 + 1) * 4
	}
	return v
}

// format2 formats a float with two decimals (Perl sprintf %.2f).
func format2(f float64) string {
	return strconvFormat(f, 2)
}

// strconvFormat is a tiny %.<prec>f formatter avoiding -0.00.
func strconvFormat(f float64, prec int) string {
	if math.Abs(f) < 0.5/pow10(prec) {
		f = 0
	}
	scale := pow10(prec)
	scaled := math.Round(f * scale)
	intPart := int(math.Abs(scaled) / scale)
	frac := int(math.Mod(math.Abs(scaled), scale))
	sign := ""
	if scaled < 0 {
		sign = "-"
	}
	fracStr := pngItoa(frac)
	for len(fracStr) < prec {
		fracStr = "0" + fracStr
	}
	return sign + pngItoa(intPart) + "." + fracStr
}

func pow10(n int) float64 {
	r := 1.0
	for i := 0; i < n; i++ {
		r *= 10
	}
	return r
}
