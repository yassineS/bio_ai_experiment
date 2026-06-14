package prinseq

// pnggraphs.go ports the PNG graphical-report flow of upstream
// `prinseq-graphs.pl` (vendored at reference_code/prinseq/) to pure
// Go. Upstream reads a `-graph_data` (.gd) file and renders a set of
// PNG plots (and an optional HTML index) using the Perl Cairo
// bindings. This port reads the same .gd (see gdparse.go) and draws
// the same set of plots onto an *image.RGBA, encoded with the stdlib
// image/png. A tiny hand-rolled bitmap font (pngfont.go) supplies
// axis/legend text since `golang.org/x/image/font` is not a
// sanctioned dependency.
//
// PARITY NOTE: byte-for-byte PNG identity with upstream's Cairo/GD
// output is NOT achievable (different rasteriser, font, and
// antialiasing model) and is explicitly NOT a goal — this is the
// standard "logical parity where byte-identity needs a foreign
// library" pattern used elsewhere in the repo (cf. mosdepth plots,
// vcftools --pca). What IS asserted: the exact set of graph PNGs
// upstream emits (same filenames), and the plotted data series
// extracted from the .gd (unit-tested for exact values).

import (
	"fmt"
	"image/color"
	"image/png"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// colour palette (mirrors prinseq-graphs.pl RGBA tuples).
var (
	colBar       = color.RGBA{127, 127, 255, 255} // b2b2ff bars
	colMean      = color.RGBA{255, 127, 127, 255} // ffb2b2 mean line
	colBackplot  = color.RGBA{242, 242, 242, 255} // 0.95 grey plot bg
	colTick      = color.RGBA{0, 0, 0, 204}       // 0.8 alpha ticks
	colLabel     = color.RGBA{0, 0, 0, 255}
	colHelpline  = color.RGBA{255, 255, 255, 230} // 0.9 white helplines
	colWhisker   = color.RGBA{178, 178, 255, 230} // b2b2ff 0.9
	colBox       = color.RGBA{127, 127, 255, 255} // 7f7fff
	colMedian    = color.RGBA{0, 0, 0, 128}       // 0.5 black
	colStd1      = color.RGBA{0, 0, 0, 10}        // ~0.04 alpha
	colStd2      = color.RGBA{0, 0, 0, 8}         // ~0.03 alpha
	colStdLine   = color.RGBA{178, 178, 255, 204} // 0.8
	colOddsUnder = color.RGBA{127, 127, 255, 255}
	colOddsOver  = color.RGBA{255, 127, 127, 255}
)

// stackColors mirrors the five duplicate-level colours
// (prinseq-graphs.pl createStackBarPlot @cols).
var stackColors = []color.RGBA{
	{69, 114, 167, 255},
	{137, 165, 78, 255},
	{170, 70, 67, 255},
	{147, 169, 207, 255},
	{51, 102, 102, 255},
}

// stackLegend mirrors the duplicate-level legend labels.
var stackLegend = []string{
	"Exact dupl.", "5' dupl.", "3' dupl.",
	"Rev. compl. exact dupl.", "Rev. compl. 5'/3' dupl.",
}

// GraphPNG is a single rendered graph: a filename suffix (e.g.
// "_ld.png"), a human title, and the encoded image.
type GraphPNG struct {
	Suffix string
	Title  string
	canvas *canvas
}

// RenderGDGraphs is the Go equivalent of prinseq-graphs.pl
// generateGraphs: it inspects which sub-tables the .gd contains and
// renders the corresponding set of plots. The returned slice is in a
// deterministic order matching upstream's emission order.
func RenderGDGraphs(d *GDData) []GraphPNG {
	var graphs []GraphPNG
	add := func(suffix, title string, cv *canvas) {
		if cv != nil {
			graphs = append(graphs, GraphPNG{Suffix: suffix, Title: title, canvas: cv})
		}
	}

	// Length distribution.
	if c, ok := d.Counts["length"]; ok && len(c) > 0 {
		st := d.Stats["length"]
		matrix, xmax, ymax := convertOdToBinMatrix(c, 1, false)
		add("_ld.png", "Length Distribution",
			drawAnnotBarPlot(matrix, xmax, ymax, &st, "Read Length in bp", "# Sequences", false, " bp"))
	}

	// Poly-A/T tail distributions (5' and 3').
	if d.Tail != 0 {
		if c, ok := d.Counts["tail5"]; ok && len(c) > 0 {
			matrix, xmax, ymax := convertOdToBinMatrix(c, 1, false)
			add("_td5.png", "Poly-A/T Tail Distribution (5')",
				drawAnnotBarPlot(matrix, xmax, ymax, nil, "5' Tail Length in bp", "# Sequences", false, " bp"))
		}
		if c, ok := d.Counts["tail3"]; ok && len(c) > 0 {
			matrix, xmax, ymax := convertOdToBinMatrix(c, 1, false)
			add("_td3.png", "Poly-A/T Tail Distribution (3')",
				drawAnnotBarPlot(matrix, xmax, ymax, nil, "3' Tail Length in bp", "# Sequences", false, " bp"))
		}
	}

	// Percentage of N's.
	if c, ok := d.Counts["ns"]; ok && len(c) > 0 {
		matrix, xmax, ymax := convertOdToBinMatrix(c, 1, false)
		add("_ns.png", "Percentage of N's",
			drawAnnotBarPlot(matrix, xmax, ymax, nil, "Percentage of N's per Read (1-100%)", "# Sequences", false, ""))
	}

	// GC content distribution.
	if c, ok := d.Counts["gc"]; ok && len(c) > 0 {
		st := d.Stats["gc"]
		matrix, xmax, ymax := convertOdToBinMatrix(c, 0, false)
		add("_gc.png", "GC Content Distribution",
			drawAnnotBarPlot(matrix, xmax, ymax, &st, "GC Content (0-100%)", "Number of Sequences", true, ""))
	}

	// Sequence complexity - DUST.
	if len(d.ComplDust) > 0 {
		matrix, xmax, ymax := convertOdToBinMatrix(d.ComplDust, 0, false)
		add("_cd.png", "Sequence complexity distribution (DUST)",
			drawAnnotBarPlot(matrix, xmax, ymax, nil, "Mean sequence complexity (DUST scores)", "Number of sequences", true, ""))
	}

	// Sequence complexity - entropy.
	if len(d.ComplEntropy) > 0 {
		matrix, xmax, ymax := convertOdToBinMatrix(d.ComplEntropy, 0, false)
		add("_ce.png", "Sequence complexity distribution (Entropy)",
			drawAnnotBarPlot(matrix, xmax, ymax, nil, "Mean sequence complexity (Entropy values)", "Number of sequences", true, ""))
	}

	// Dinucleotide odds-ratio PCA (microbial / viral) + odds ratio.
	if len(d.DinucOdds) > 0 {
		row := dinucOddsRow(d.DinucOdds)
		add("_pm.png", "PCA (microbial)", drawPCAPlot(row, "m"))
		add("_pv.png", "PCA (viral)", drawPCAPlot(row, "v"))
		add("_or.png", "Odds ratios", drawOddsRatioPlot(d.DinucOdds))
	}

	// Base quality distribution (relative).
	if len(d.Quals) > 0 {
		matrix, ymax := convertToBoxValues(d.Quals, 4)
		add("_qd.png", "Base Quality Distribution",
			drawBoxPlot(matrix, ymax, "Read position in %", "Quality score", 1, ""))
	}
	// Base quality distribution (binned by bp).
	if len(d.QualsBin) > 0 {
		matrix, ymax := convertToBoxValues(d.QualsBin, 4)
		add("_qd2.png", "Base Quality Distribution",
			drawBoxPlot(matrix, ymax, "Read position in bp", "Quality score", d.BinVal, "bp"))
	}
	// Sequence quality distribution (mean per read).
	if len(d.QualsMean) > 0 {
		matrix, xmax, ymax := convertToBarValues(d.QualsMean, 5, 1)
		add("_qd3.png", "Sequence Quality Distribution",
			drawBarPlot(matrix, xmax, ymax, "Mean of quality scores per sequence", "Number of sequences"))
	}

	// Sequence duplication levels.
	if len(d.DubsCounts) > 0 {
		matrix, xmax, ymax := convertOdToStackBinMatrix(d.DubsCounts, 5, 1, 100, false)
		add("_df.png", "Sequence duplication level",
			drawStackBarPlot(matrix, xmax, ymax, "Number of duplicates", "Number of sequences", ""))
	}
	if len(d.DubsLength) > 0 {
		matrix, xmax, ymax := convertOdToStackBinMatrix(d.DubsLength, 5, 1, 0, false)
		add("_dl.png", "Sequence duplication level",
			drawStackBarPlot(matrix, xmax, ymax, "Read Length in bp", "Number of duplicates", " bp"))
	}
	if len(d.DubsCounts) > 0 {
		dubsmax := dubsMaxTable(d.DubsCounts)
		matrix, xmax, ymax := convertOdToStackBinMatrix(dubsmax, 5, 1, 100, false)
		add("_dm.png", "Sequence duplication level",
			drawStackBarPlot(matrix, xmax, ymax, "Sequence", "Number of duplicates", ""))
	}

	return graphs
}

// WriteGDGraphs renders the .gd graph set and writes one PNG per
// graph to disk, named "<prefix><suffix>". When writeHTML is true an
// HTML index ("<prefix>.html") embedding the PNGs by filename is
// also written. The returned slice lists the written file paths in
// emission order. This is the Go equivalent of the prinseq-graphs.pl
// `-png_all` (and `-html_all`) driver.
func WriteGDGraphs(d *GDData, prefix string, writeHTML bool) ([]string, error) {
	graphs := RenderGDGraphs(d)
	var written []string
	for _, g := range graphs {
		path := prefix + g.Suffix
		f, err := os.Create(path)
		if err != nil {
			return written, err
		}
		if err := png.Encode(f, g.canvas.img); err != nil {
			f.Close()
			return written, err
		}
		if err := f.Close(); err != nil {
			return written, err
		}
		written = append(written, path)
	}
	if writeHTML {
		htmlPath := prefix + ".html"
		if err := writeGDHTML(d, prefix, graphs, htmlPath); err != nil {
			return written, err
		}
		written = append(written, htmlPath)
	}
	return written, nil
}

// writeGDHTML emits a minimal HTML report linking the PNG files,
// mirroring the structure of prinseq-graphs.pl generateHtml (input
// info header + one section per graph). Images are referenced by
// their basename so the HTML and PNGs are portable as a set.
func writeGDHTML(d *GDData, prefix string, graphs []GraphPNG, htmlPath string) error {
	var b strings.Builder
	b.WriteString("<!DOCTYPE html>\n<html>\n<head>\n")
	b.WriteString(`<meta http-equiv="Content-Type" content="text/html; charset=UTF-8">` + "\n")
	b.WriteString("<title>PRINSEQ-graphs Report</title>\n</head>\n<body>\n<center>\n")
	b.WriteString("<h2>PRINSEQ-graphs HTML Report</h2>\n")
	b.WriteString("<div style=\"text-align:left;width:740px;margin:auto;\">\n")
	b.WriteString("<h3>Input Information</h3>\n<table>\n")
	fmt.Fprintf(&b, "<tr><td># Sequences:</td><td>%s</td></tr>\n", addCommas(d.NumSeqs))
	fmt.Fprintf(&b, "<tr><td>Total bases:</td><td>%s</td></tr>\n", addCommas(d.NumBases))
	if d.Format1 != "" {
		fmt.Fprintf(&b, "<tr><td>Input format:</td><td>%s</td></tr>\n", strings.ToUpper(d.Format1))
	}
	b.WriteString("</table>\n<hr>\n")
	base := filepath.Base(prefix)
	for _, g := range graphs {
		fmt.Fprintf(&b, "<h3>%s</h3>\n", g.Title)
		fmt.Fprintf(&b, "<img src=\"%s%s\" alt=\"%s\"><hr>\n", base, g.Suffix, g.Title)
	}
	b.WriteString("</div>\n</center>\n</body>\n</html>\n")
	return os.WriteFile(htmlPath, []byte(b.String()), 0o644)
}

// dinucOddsRow returns the 10 dinucodds values in sorted-key order,
// matching upstream `map {...} sort keys %{$data->{dinucodds}}`.
func dinucOddsRow(odds map[string]float64) []float64 {
	keys := make([]string, 0, len(odds))
	for k := range odds {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	row := make([]float64, len(keys))
	for i, k := range keys {
		row[i] = odds[k]
	}
	return row
}

// dubsMaxTable reproduces upstream's `%dubsmax` construction
// (prinseq-graphs.pl lines 896-912): walk dup precounts descending,
// expanding into per-rank (1..100) entries of {type -> precount}.
func dubsMaxTable(dubs map[int]map[int]int) map[int]map[int]int {
	out := map[int]map[int]int{}
	counts := make([]int, 0, len(dubs))
	for n := range dubs {
		counts = append(counts, n)
	}
	sort.Sort(sort.Reverse(sort.IntSlice(counts)))
	rank := 1
	for _, n := range counts {
		inner := dubs[n]
		// Iterate dup-types in ascending order for determinism.
		types := make([]int, 0, len(inner))
		for s := range inner {
			types = append(types, s)
		}
		sort.Ints(types)
		stop := false
		for _, s := range types {
			for i := 0; i < inner[s]; i++ {
				m, ok := out[rank]
				if !ok {
					m = map[int]int{}
					out[rank] = m
				}
				m[s] = n
				rank++
				if rank > 100 {
					stop = true
					break
				}
			}
			if stop {
				break
			}
		}
		if stop {
			break
		}
	}
	return out
}
