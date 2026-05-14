// HTML report generator. The output is a single self-contained file:
// embedded CSS, inline SVG charts, no JavaScript, no CDN fetches. The
// report is composed via html/template so user-controlled strings (tool
// version, etc.) are escaped automatically.

package fastp

import (
	"bytes"
	"fmt"
	"html/template"
	"os"
	"sort"
	"strings"
	"time"
)

// WriteHTMLReport renders an HTML report from stats and writes it to
// path. The output is a single self-contained file (embedded CSS, inline
// SVG) and contains no JavaScript or CDN references.
func WriteHTMLReport(path string, stats *ProcessStats) error {
	if stats == nil {
		return fmt.Errorf("stats is nil")
	}
	data := buildHTMLData(stats)
	tpl, err := template.New("fastp-report").Funcs(template.FuncMap{
		"safeSVG": func(s string) template.HTML { return template.HTML(s) },
		"mul100":  mul100,
	}).Parse(htmlTemplate)
	if err != nil {
		return fmt.Errorf("parse template: %w", err)
	}
	var buf bytes.Buffer
	if err := tpl.Execute(&buf, data); err != nil {
		return fmt.Errorf("render template: %w", err)
	}
	f, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create HTML report %q: %w", path, err)
	}
	defer f.Close()
	if _, err := f.Write(buf.Bytes()); err != nil {
		return fmt.Errorf("write HTML report: %w", err)
	}
	return nil
}

// htmlData is the value passed to the template.
type htmlData struct {
	Generated      string
	Version        string
	PairedEnd      bool
	BeforeReads    int64
	AfterReads     int64
	BeforeBases    int64
	AfterBases     int64
	BeforeQ20Rate  float64
	BeforeQ30Rate  float64
	AfterQ20Rate   float64
	AfterQ30Rate   float64
	BeforeGC       float64
	AfterGC        float64
	FilteredPassed int64
	FilteredLowQ   int64
	FilteredManyN  int64
	FilteredShort  int64
	FilteredLong   int64
	AdapterReads   int64
	AdapterBases   int64
	DetectedR1     string
	DetectedR2     string
	QualitySVGR1   string
	QualitySVGR2   string
	CompositionR1  string
	CompositionR2  string
	LengthSVG      string
	FilterReasons  []reasonRow
	AdapterRows    []adapterRow

	// Duplication section. Visible only when ShowDuplication is true; the
	// CLI sets that whenever --dup_calc_accuracy >= 1 (or --dedup) was on
	// for the run, so empty Duplication panels don't bloat reports of
	// runs that didn't ask for it.
	ShowDuplication bool
	DupRate         float64
	DupTotal        int64
	DedupDropped    int64
	DupHistRows     []dupHistRow
	DupHistSVG      string
}

type dupHistRow struct {
	Count int
	Reads int64
}

type reasonRow struct {
	Reason string
	Count  int64
}

type adapterRow struct {
	Label string
	Value string
}

// buildHTMLData converts a ProcessStats into the value plugged into the
// template.
func buildHTMLData(s *ProcessStats) htmlData {
	totalBefore := s.TotalBasesR1 + s.TotalBasesR2
	if totalBefore == 0 {
		totalBefore = s.TotalBases
	}
	totalAfter := s.CleanBasesR1 + s.CleanBasesR2
	if totalAfter == 0 {
		totalAfter = s.CleanBases
	}
	d := htmlData{
		Generated:      time.Now().UTC().Format("2006-01-02 15:04:05 UTC"),
		Version:        ToolVersion,
		PairedEnd:      s.TotalReadsR2 > 0,
		BeforeReads:    int64(s.TotalReads),
		AfterReads:     int64(s.CleanReads),
		BeforeBases:    totalBefore,
		AfterBases:     totalAfter,
		BeforeQ20Rate:  safeDiv(s.Q20BasesBefore, totalBefore),
		BeforeQ30Rate:  safeDiv(s.Q30BasesBefore, totalBefore),
		AfterQ20Rate:   safeDiv(s.Q20BasesAfter, totalAfter),
		AfterQ30Rate:   safeDiv(s.Q30BasesAfter, totalAfter),
		BeforeGC:       safeDiv(s.GCBasesBefore, totalBefore),
		AfterGC:        safeDiv(s.GCBasesAfter, totalAfter),
		FilteredPassed: int64(s.CleanReads),
		FilteredLowQ:   int64(s.LowQualityReads),
		FilteredManyN:  int64(s.TooManyNReads),
		FilteredShort:  int64(s.TooShortReads),
		FilteredLong:   int64(s.TooLongReads),
		AdapterReads:   int64(s.AdapterTrimmedReads),
		AdapterBases:   s.AdapterTrimmedBases,
		DetectedR1:     s.DetectedAdapterR1,
		DetectedR2:     s.DetectedAdapterR2,
	}

	d.QualitySVGR1 = renderQualitySVG(s, 0, "Read 1 per-base mean quality")
	if d.PairedEnd {
		d.QualitySVGR2 = renderQualitySVG(s, 1, "Read 2 per-base mean quality")
	}
	d.CompositionR1 = renderCompositionSVG(s, 0, "Read 1 per-base composition")
	if d.PairedEnd {
		d.CompositionR2 = renderCompositionSVG(s, 1, "Read 2 per-base composition")
	}
	d.LengthSVG = renderLengthSVG(s, "Length distribution (after filtering)")

	d.FilterReasons = []reasonRow{
		{"Passed filter", int64(s.CleanReads)},
		{"Low quality", int64(s.LowQualityReads)},
		{"Too many N", int64(s.TooManyNReads)},
		{"Too short", int64(s.TooShortReads)},
		{"Too long", int64(s.TooLongReads)},
	}
	d.AdapterRows = []adapterRow{
		{"Adapter-trimmed reads", fmt.Sprintf("%d", s.AdapterTrimmedReads)},
		{"Adapter-trimmed bases", fmt.Sprintf("%d", s.AdapterTrimmedBases)},
	}
	if s.DetectedAdapterR1 != "" {
		d.AdapterRows = append(d.AdapterRows, adapterRow{"Detected R1 adapter", s.DetectedAdapterR1})
	}
	if s.DetectedAdapterR2 != "" {
		d.AdapterRows = append(d.AdapterRows, adapterRow{"Detected R2 adapter", s.DetectedAdapterR2})
	}

	if s.DupTotal > 0 || len(s.DupHist) > 0 || s.DedupDropped > 0 {
		d.ShowDuplication = true
		d.DupRate = s.DupRate
		d.DupTotal = s.DupTotal
		d.DedupDropped = int64(s.DedupDropped)
		// Build a table sorted by count (ascending). Cap the level
		// histogram at 20 so very high duplication tails don't blow up
		// the table; we collapse the tail into a single "20+" bucket.
		keys := make([]int, 0, len(s.DupHist))
		for k := range s.DupHist {
			keys = append(keys, k)
		}
		sort.Ints(keys)
		var tailTotal int64
		for _, k := range keys {
			if k >= 20 {
				tailTotal += s.DupHist[k]
				continue
			}
			d.DupHistRows = append(d.DupHistRows, dupHistRow{Count: k, Reads: s.DupHist[k]})
		}
		if tailTotal > 0 {
			d.DupHistRows = append(d.DupHistRows, dupHistRow{Count: 20, Reads: tailTotal})
		}
		d.DupHistSVG = renderDupHistSVG(d.DupHistRows, "Duplication levels")
	}
	return d
}

// renderDupHistSVG renders a bar plot of the duplication-level
// histogram. rows is expected to be sorted ascending by Count.
func renderDupHistSVG(rows []dupHistRow, title string) string {
	if len(rows) == 0 {
		return emptyChart(title, "no data")
	}
	xs := make([]float64, len(rows))
	ys := make([]float64, len(rows))
	maxY := 1.0
	for i, r := range rows {
		xs[i] = float64(r.Count)
		ys[i] = float64(r.Reads)
		if ys[i] > maxY {
			maxY = ys[i]
		}
	}
	return barPlot(title, "occurrence count", "reads", xs, ys, maxY, "#9c27b0")
}

// renderQualitySVG returns an inline SVG plot of per-cycle mean quality
// for the given read index.
func renderQualitySVG(s *ProcessStats, readIdx int, title string) string {
	n := len(s.QualSumByCycle[readIdx])
	if n == 0 {
		return emptyChart(title, "no data")
	}
	mean := make([]float64, n)
	maxQ := 1.0
	for i := 0; i < n; i++ {
		c := s.QualCountByCycle[readIdx][i]
		if c > 0 {
			mean[i] = float64(s.QualSumByCycle[readIdx][i]) / float64(c)
			if mean[i] > maxQ {
				maxQ = mean[i]
			}
		}
	}
	// Cap the y-axis at maxQ rounded up to nearest 5, but at least 40.
	yMax := 40.0
	if maxQ > yMax {
		yMax = (float64(int(maxQ/5)) + 1) * 5
	}
	return linePlot(title, "cycle", "mean Q", mean, yMax, "#4CAF50")
}

// renderCompositionSVG returns an inline SVG plot of per-cycle base
// composition (A/C/G/T/N fractions) for the given read index. Each base
// is a separate colored line.
func renderCompositionSVG(s *ProcessStats, readIdx int, title string) string {
	n := len(s.QualSumByCycle[readIdx])
	if n == 0 {
		return emptyChart(title, "no data")
	}
	series := make(map[string][]float64, 5)
	colors := map[string]string{
		"A": "#1f77b4", "C": "#ff7f0e", "G": "#2ca02c", "T": "#d62728", "N": "#7f7f7f",
	}
	order := []string{"A", "C", "G", "T", "N"}
	for b, name := range order {
		row := make([]float64, n)
		for i := 0; i < n; i++ {
			c := s.QualCountByCycle[readIdx][i]
			if c > 0 {
				row[i] = float64(s.BaseCountByCycle[readIdx][b][i]) / float64(c)
			}
		}
		series[name] = row
	}
	return multiLinePlot(title, "cycle", "fraction", series, order, colors, 1.0)
}

// renderLengthSVG renders the AFTER-filtering length distribution as a
// histogram. R1 + R2 are combined into a single distribution.
func renderLengthSVG(s *ProcessStats, title string) string {
	combined := map[int]int64{}
	for k, v := range s.LengthHistAfter[0] {
		combined[k] += v
	}
	for k, v := range s.LengthHistAfter[1] {
		combined[k] += v
	}
	if len(combined) == 0 {
		return emptyChart(title, "no data")
	}
	type lc struct {
		l int
		c int64
	}
	all := make([]lc, 0, len(combined))
	for k, v := range combined {
		all = append(all, lc{k, v})
	}
	sort.Slice(all, func(i, j int) bool { return all[i].l < all[j].l })
	xs := make([]float64, len(all))
	ys := make([]float64, len(all))
	maxY := 1.0
	for i, p := range all {
		xs[i] = float64(p.l)
		ys[i] = float64(p.c)
		if ys[i] > maxY {
			maxY = ys[i]
		}
	}
	return barPlot(title, "length", "count", xs, ys, maxY, "#3f51b5")
}

// emptyChart returns an SVG placeholder for an empty dataset.
func emptyChart(title, msg string) string {
	return fmt.Sprintf(`<svg class="chart" viewBox="0 0 600 60" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="%s"><text x="10" y="35" font-family="sans-serif" font-size="14" fill="#777">%s: %s</text></svg>`,
		htmlEscape(title), htmlEscape(title), htmlEscape(msg))
}

// linePlot renders a single-series line plot as an inline SVG.
func linePlot(title, xLabel, yLabel string, ys []float64, yMax float64, color string) string {
	const w, h = 600, 220
	const padL, padR, padT, padB = 50, 20, 30, 40
	plotW := w - padL - padR
	plotH := h - padT - padB
	n := len(ys)
	if n == 0 || yMax <= 0 {
		return emptyChart(title, "no data")
	}
	var b strings.Builder
	fmt.Fprintf(&b, `<svg class="chart" viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="%s">`, w, h, htmlEscape(title))
	fmt.Fprintf(&b, `<rect x="0" y="0" width="%d" height="%d" fill="white"/>`, w, h)
	fmt.Fprintf(&b, `<text x="%d" y="20" font-family="sans-serif" font-size="14" font-weight="bold" fill="#333">%s</text>`, padL, htmlEscape(title))
	// axes
	fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#888" stroke-width="1"/>`, padL, padT, padL, padT+plotH)
	fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#888" stroke-width="1"/>`, padL, padT+plotH, padL+plotW, padT+plotH)
	// y-axis ticks
	for i := 0; i <= 4; i++ {
		val := yMax * float64(i) / 4
		y := padT + plotH - int(float64(plotH)*float64(i)/4)
		fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#eee"/>`, padL, y, padL+plotW, y)
		fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="10" fill="#666" text-anchor="end">%.0f</text>`, padL-4, y+3, val)
	}
	// path
	var path strings.Builder
	denom := float64(n - 1)
	if denom <= 0 {
		denom = 1
	}
	for i, v := range ys {
		x := padL + int(float64(plotW)*float64(i)/denom)
		y := padT + plotH - int(float64(plotH)*v/yMax)
		if i == 0 {
			fmt.Fprintf(&path, "M %d %d", x, y)
		} else {
			fmt.Fprintf(&path, " L %d %d", x, y)
		}
	}
	fmt.Fprintf(&b, `<path d="%s" fill="none" stroke="%s" stroke-width="1.5"/>`, path.String(), color)
	// axis labels
	fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="11" fill="#555" text-anchor="middle">%s</text>`, padL+plotW/2, h-10, htmlEscape(xLabel))
	fmt.Fprintf(&b, `<text x="15" y="%d" font-family="sans-serif" font-size="11" fill="#555" transform="rotate(-90 15,%d)">%s</text>`, padT+plotH/2, padT+plotH/2, htmlEscape(yLabel))
	b.WriteString(`</svg>`)
	return b.String()
}

// multiLinePlot renders multiple series on the same axes, with a legend.
func multiLinePlot(title, xLabel, yLabel string, series map[string][]float64, order []string, colors map[string]string, yMax float64) string {
	const w, h = 600, 240
	const padL, padR, padT, padB = 50, 100, 30, 40
	plotW := w - padL - padR
	plotH := h - padT - padB
	var n int
	for _, name := range order {
		if len(series[name]) > n {
			n = len(series[name])
		}
	}
	if n == 0 || yMax <= 0 {
		return emptyChart(title, "no data")
	}
	var b strings.Builder
	fmt.Fprintf(&b, `<svg class="chart" viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="%s">`, w, h, htmlEscape(title))
	fmt.Fprintf(&b, `<rect x="0" y="0" width="%d" height="%d" fill="white"/>`, w, h)
	fmt.Fprintf(&b, `<text x="%d" y="20" font-family="sans-serif" font-size="14" font-weight="bold" fill="#333">%s</text>`, padL, htmlEscape(title))
	fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#888" stroke-width="1"/>`, padL, padT, padL, padT+plotH)
	fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#888" stroke-width="1"/>`, padL, padT+plotH, padL+plotW, padT+plotH)
	for i := 0; i <= 4; i++ {
		val := yMax * float64(i) / 4
		y := padT + plotH - int(float64(plotH)*float64(i)/4)
		fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#eee"/>`, padL, y, padL+plotW, y)
		fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="10" fill="#666" text-anchor="end">%.2f</text>`, padL-4, y+3, val)
	}
	denom := float64(n - 1)
	if denom <= 0 {
		denom = 1
	}
	for _, name := range order {
		col := colors[name]
		if col == "" {
			col = "#000"
		}
		var path strings.Builder
		for i, v := range series[name] {
			x := padL + int(float64(plotW)*float64(i)/denom)
			y := padT + plotH - int(float64(plotH)*v/yMax)
			if i == 0 {
				fmt.Fprintf(&path, "M %d %d", x, y)
			} else {
				fmt.Fprintf(&path, " L %d %d", x, y)
			}
		}
		fmt.Fprintf(&b, `<path d="%s" fill="none" stroke="%s" stroke-width="1.2"/>`, path.String(), col)
	}
	// legend
	for i, name := range order {
		col := colors[name]
		if col == "" {
			col = "#000"
		}
		y := padT + 10 + i*16
		fmt.Fprintf(&b, `<rect x="%d" y="%d" width="10" height="10" fill="%s"/>`, padL+plotW+10, y, col)
		fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="11" fill="#333">%s</text>`, padL+plotW+25, y+9, htmlEscape(name))
	}
	fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="11" fill="#555" text-anchor="middle">%s</text>`, padL+plotW/2, h-10, htmlEscape(xLabel))
	fmt.Fprintf(&b, `<text x="15" y="%d" font-family="sans-serif" font-size="11" fill="#555" transform="rotate(-90 15,%d)">%s</text>`, padT+plotH/2, padT+plotH/2, htmlEscape(yLabel))
	b.WriteString(`</svg>`)
	return b.String()
}

// barPlot renders a simple bar plot.
func barPlot(title, xLabel, yLabel string, xs, ys []float64, yMax float64, color string) string {
	const w, h = 600, 220
	const padL, padR, padT, padB = 50, 20, 30, 40
	plotW := w - padL - padR
	plotH := h - padT - padB
	n := len(xs)
	if n == 0 || yMax <= 0 {
		return emptyChart(title, "no data")
	}
	var b strings.Builder
	fmt.Fprintf(&b, `<svg class="chart" viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="%s">`, w, h, htmlEscape(title))
	fmt.Fprintf(&b, `<rect x="0" y="0" width="%d" height="%d" fill="white"/>`, w, h)
	fmt.Fprintf(&b, `<text x="%d" y="20" font-family="sans-serif" font-size="14" font-weight="bold" fill="#333">%s</text>`, padL, htmlEscape(title))
	fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#888" stroke-width="1"/>`, padL, padT, padL, padT+plotH)
	fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#888" stroke-width="1"/>`, padL, padT+plotH, padL+plotW, padT+plotH)
	for i := 0; i <= 4; i++ {
		val := yMax * float64(i) / 4
		y := padT + plotH - int(float64(plotH)*float64(i)/4)
		fmt.Fprintf(&b, `<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#eee"/>`, padL, y, padL+plotW, y)
		fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="10" fill="#666" text-anchor="end">%.0f</text>`, padL-4, y+3, val)
	}
	barW := plotW / n
	if barW < 1 {
		barW = 1
	}
	for i := 0; i < n; i++ {
		x := padL + i*plotW/n
		y := padT + plotH - int(float64(plotH)*ys[i]/yMax)
		bh := padT + plotH - y
		fmt.Fprintf(&b, `<rect x="%d" y="%d" width="%d" height="%d" fill="%s"/>`, x, y, barW-1, bh, color)
	}
	// x-axis ticks: at most 5 evenly spaced labels.
	ticks := 5
	if n < ticks {
		ticks = n
	}
	if ticks <= 1 {
		fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="10" fill="#666" text-anchor="middle">%.0f</text>`, padL+plotW/2, padT+plotH+12, xs[0])
	} else {
		for i := 0; i < ticks; i++ {
			idx := i * (n - 1) / (ticks - 1)
			x := padL + idx*plotW/n + barW/2
			fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="10" fill="#666" text-anchor="middle">%.0f</text>`, x, padT+plotH+12, xs[idx])
		}
	}
	fmt.Fprintf(&b, `<text x="%d" y="%d" font-family="sans-serif" font-size="11" fill="#555" text-anchor="middle">%s</text>`, padL+plotW/2, h-10, htmlEscape(xLabel))
	fmt.Fprintf(&b, `<text x="15" y="%d" font-family="sans-serif" font-size="11" fill="#555" transform="rotate(-90 15,%d)">%s</text>`, padT+plotH/2, padT+plotH/2, htmlEscape(yLabel))
	b.WriteString(`</svg>`)
	return b.String()
}

// htmlEscape escapes characters that would break SVG/HTML attribute
// values. It is intentionally limited to the small set used by our
// generators; user-controlled strings reach the template (not these
// helpers) and are escaped there.
func htmlEscape(s string) string {
	r := strings.NewReplacer("&", "&amp;", "<", "&lt;", ">", "&gt;", `"`, "&quot;")
	return r.Replace(s)
}

const htmlTemplate = `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>fastp report</title>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif; margin: 20px; background: #fafafa; color: #222; }
.container { max-width: 1100px; margin: 0 auto; background: #fff; padding: 24px 30px; border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.08); }
h1 { border-bottom: 3px solid #4CAF50; padding-bottom: 8px; margin-top: 0; }
h2 { color: #444; border-bottom: 1px solid #ddd; padding-bottom: 6px; margin-top: 32px; }
table { width: 100%; border-collapse: collapse; margin: 12px 0 20px; }
th, td { border: 1px solid #e0e0e0; padding: 8px 12px; text-align: left; }
th { background: #f1f8e9; color: #2e7d32; }
tr:nth-child(even) td { background: #fafafa; }
.chart { display: block; margin: 12px 0; max-width: 100%; height: auto; border: 1px solid #eee; border-radius: 4px; background: #fff; }
.summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; margin: 16px 0; }
.kpi { background: #e8f5e9; border-left: 4px solid #4CAF50; padding: 10px 14px; border-radius: 4px; }
.kpi .label { font-size: 12px; color: #555; text-transform: uppercase; letter-spacing: 0.05em; }
.kpi .value { font-size: 22px; font-weight: 600; color: #1b5e20; }
.footer { margin-top: 32px; padding-top: 12px; border-top: 1px solid #e0e0e0; color: #777; font-size: 12px; }
code { background: #f5f5f5; padding: 1px 4px; border-radius: 3px; }
</style>
</head>
<body>
<div class="container">
<h1>fastp report</h1>
<p>Generated {{.Generated}} by fastp v{{.Version}} (Go implementation).</p>

<h2>Summary</h2>
<div class="summary-grid">
  <div class="kpi"><div class="label">Reads (before)</div><div class="value">{{.BeforeReads}}</div></div>
  <div class="kpi"><div class="label">Reads (after)</div><div class="value">{{.AfterReads}}</div></div>
  <div class="kpi"><div class="label">Bases (before)</div><div class="value">{{.BeforeBases}}</div></div>
  <div class="kpi"><div class="label">Bases (after)</div><div class="value">{{.AfterBases}}</div></div>
</div>
<table>
  <tr><th>Metric</th><th>Before filtering</th><th>After filtering</th></tr>
  <tr><td>Q20 rate</td><td>{{printf "%.2f%%" (mul100 .BeforeQ20Rate)}}</td><td>{{printf "%.2f%%" (mul100 .AfterQ20Rate)}}</td></tr>
  <tr><td>Q30 rate</td><td>{{printf "%.2f%%" (mul100 .BeforeQ30Rate)}}</td><td>{{printf "%.2f%%" (mul100 .AfterQ30Rate)}}</td></tr>
  <tr><td>GC content</td><td>{{printf "%.2f%%" (mul100 .BeforeGC)}}</td><td>{{printf "%.2f%%" (mul100 .AfterGC)}}</td></tr>
</table>

<h2>Per-base quality</h2>
{{safeSVG .QualitySVGR1}}
{{if .PairedEnd}}{{safeSVG .QualitySVGR2}}{{end}}

<h2>Per-base composition</h2>
{{safeSVG .CompositionR1}}
{{if .PairedEnd}}{{safeSVG .CompositionR2}}{{end}}

<h2>Length distribution</h2>
{{safeSVG .LengthSVG}}

<h2>Filtering reasons</h2>
<table>
  <tr><th>Reason</th><th>Reads</th></tr>
  {{range .FilterReasons}}<tr><td>{{.Reason}}</td><td>{{.Count}}</td></tr>{{end}}
</table>

<h2>Adapter trimming</h2>
<table>
  <tr><th>Field</th><th>Value</th></tr>
  {{range .AdapterRows}}<tr><td>{{.Label}}</td><td><code>{{.Value}}</code></td></tr>{{end}}
</table>

{{if .ShowDuplication}}
<h2>Duplication</h2>
<table>
  <tr><th>Field</th><th>Value</th></tr>
  <tr><td>Duplication rate</td><td>{{printf "%.2f%%" (mul100 .DupRate)}}</td></tr>
  <tr><td>Reads scanned</td><td>{{.DupTotal}}</td></tr>
  {{if .DedupDropped}}<tr><td>Reads dropped by --dedup</td><td>{{.DedupDropped}}</td></tr>{{end}}
</table>
{{if .DupHistRows}}
{{safeSVG .DupHistSVG}}
<table>
  <tr><th>Occurrence count</th><th>Reads</th></tr>
  {{range .DupHistRows}}<tr><td>{{.Count}}</td><td>{{.Reads}}</td></tr>{{end}}
</table>
{{end}}
{{end}}

<div class="footer">fastp (Go) v{{.Version}}. Self-contained report &mdash; no scripts, no external resources.</div>
</div>
</body>
</html>
`

// mul100 is registered as a template function via the FuncMap in
// WriteHTMLReport. It scales a fraction (0..1) to a percentage (0..100)
// so the template can write {{printf "%.2f%%" (mul100 .Rate)}}.
func mul100(f float64) float64 { return f * 100 }
