package prinseq

import (
	"fmt"
	"io"
	"math"
	"sort"
)

// GraphType represents different graph types
type GraphType string

const (
	GraphTypeLength     GraphType = "length"
	GraphTypeGC         GraphType = "gc"
	GraphTypeQuality    GraphType = "quality"
	GraphTypeDinuc      GraphType = "dinucleotides"
	GraphTypePositional GraphType = "positional_quality"
)

// GenerateGraph generates a simple ASCII/SVG graph from statistics
func GenerateGraph(stats *Stats, graphType GraphType, writer io.Writer) error {
	switch graphType {
	case GraphTypeLength:
		return generateLengthGraph(stats, writer)
	case GraphTypeGC:
		return generateGCGraph(stats, writer)
	case GraphTypeQuality:
		return generateQualityGraph(stats, writer)
	case GraphTypeDinuc:
		return generateDinucGraph(stats, writer)
	case GraphTypePositional:
		return generatePositionalQualityGraph(stats, writer)
	default:
		return fmt.Errorf("unknown graph type: %s", graphType)
	}
}

// GenerateSVG generates SVG graphs for all available statistics
func GenerateSVG(stats *Stats, writer io.Writer) error {
	_, err := fmt.Fprintf(writer, `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="800" viewBox="0 0 1200 800">
  <style>
    .title { font: bold 16px sans-serif; }
    .axis { font: 12px sans-serif; }
    .bar { fill: steelblue; }
    .line { fill: none; stroke: steelblue; stroke-width: 2; }
    .grid { stroke: #ccc; stroke-width: 0.5; }
  </style>
`)
	if err != nil {
		return err
	}

	// Generate length distribution graph
	if len(stats.LengthDistribution) > 0 {
		if err := generateLengthSVG(stats, writer, 50, 50); err != nil {
			return err
		}
	}

	// Generate quality distribution graph (if available)
	if len(stats.QualityDistribution) > 0 {
		if err := generateQualitySVG(stats, writer, 650, 50); err != nil {
			return err
		}
	}

	// Generate positional quality graph (if available)
	if len(stats.PositionalQuality) > 0 {
		if err := generatePositionalQualitySVG(stats, writer, 50, 450); err != nil {
			return err
		}
	}

	_, err = fmt.Fprintf(writer, "</svg>\n")
	return err
}

func generateLengthGraph(stats *Stats, writer io.Writer) error {
	if len(stats.LengthDistribution) == 0 {
		return fmt.Errorf("no length distribution data available")
	}

	// Get sorted lengths
	lengths := make([]int, 0, len(stats.LengthDistribution))
	for length := range stats.LengthDistribution {
		lengths = append(lengths, length)
	}
	sort.Ints(lengths)

	// Find max count for scaling
	maxCount := 0
	for _, count := range stats.LengthDistribution {
		if count > maxCount {
			maxCount = count
		}
	}

	fmt.Fprintf(writer, "Length Distribution:\n")
	fmt.Fprintf(writer, "%-10s | %-10s | %s\n", "Length", "Count", "Histogram")
	fmt.Fprintf(writer, "%s\n", "------------------------------------------------------------")

	for _, length := range lengths {
		count := stats.LengthDistribution[length]
		barWidth := int(float64(count) / float64(maxCount) * 40)
		bar := ""
		for i := 0; i < barWidth; i++ {
			bar += "█"
		}
		fmt.Fprintf(writer, "%-10d | %-10d | %s\n", length, count, bar)
	}

	return nil
}

func generateGCGraph(stats *Stats, writer io.Writer) error {
	fmt.Fprintf(writer, "GC Content: %.2f%%\n", stats.GCContent)
	
	// Create a simple bar representation
	gcPercent := int(stats.GCContent)
	atPercent := 100 - gcPercent
	
	fmt.Fprintf(writer, "\n")
	fmt.Fprintf(writer, "GC: [")
	for i := 0; i < gcPercent/2; i++ {
		fmt.Fprintf(writer, "█")
	}
	fmt.Fprintf(writer, "] %.1f%%\n", stats.GCContent)
	
	fmt.Fprintf(writer, "AT: [")
	for i := 0; i < atPercent/2; i++ {
		fmt.Fprintf(writer, "█")
	}
	fmt.Fprintf(writer, "] %.1f%%\n", 100.0-stats.GCContent)
	
	return nil
}

func generateQualityGraph(stats *Stats, writer io.Writer) error {
	if len(stats.QualityDistribution) == 0 {
		return fmt.Errorf("no quality distribution data available")
	}

	// Get sorted quality scores
	qualities := make([]int, 0, len(stats.QualityDistribution))
	for qual := range stats.QualityDistribution {
		qualities = append(qualities, qual)
	}
	sort.Ints(qualities)

	// Find max count for scaling
	maxCount := 0
	for _, count := range stats.QualityDistribution {
		if count > maxCount {
			maxCount = count
		}
	}

	fmt.Fprintf(writer, "Quality Score Distribution:\n")
	fmt.Fprintf(writer, "%-10s | %-10s | %s\n", "Quality", "Count", "Histogram")
	fmt.Fprintf(writer, "%s\n", "------------------------------------------------------------")

	for _, qual := range qualities {
		count := stats.QualityDistribution[qual]
		barWidth := int(float64(count) / float64(maxCount) * 40)
		bar := ""
		for i := 0; i < barWidth; i++ {
			bar += "█"
		}
		fmt.Fprintf(writer, "%-10d | %-10d | %s\n", qual, count, bar)
	}

	return nil
}

func generateDinucGraph(stats *Stats, writer io.Writer) error {
	if len(stats.Dinucleotides) == 0 {
		return fmt.Errorf("no dinucleotide data available")
	}

	// Get sorted dinucleotides
	dinucs := make([]string, 0, len(stats.Dinucleotides))
	for dinuc := range stats.Dinucleotides {
		dinucs = append(dinucs, dinuc)
	}
	sort.Strings(dinucs)

	// Find max count for scaling
	maxCount := 0
	totalCount := 0
	for _, count := range stats.Dinucleotides {
		if count > maxCount {
			maxCount = count
		}
		totalCount += count
	}

	fmt.Fprintf(writer, "Dinucleotide Frequencies:\n")
	fmt.Fprintf(writer, "%-6s | %-10s | %-8s | %s\n", "Dinuc", "Count", "Percent", "Histogram")
	fmt.Fprintf(writer, "%s\n", "------------------------------------------------------------")

	for _, dinuc := range dinucs {
		count := stats.Dinucleotides[dinuc]
		percent := float64(count) / float64(totalCount) * 100.0
		barWidth := int(float64(count) / float64(maxCount) * 30)
		bar := ""
		for i := 0; i < barWidth; i++ {
			bar += "█"
		}
		fmt.Fprintf(writer, "%-6s | %-10d | %6.2f%% | %s\n", dinuc, count, percent, bar)
	}

	return nil
}

func generatePositionalQualityGraph(stats *Stats, writer io.Writer) error {
	if len(stats.PositionalQuality) == 0 {
		return fmt.Errorf("no positional quality data available")
	}

	// Find min/max for scaling
	minQual, maxQual := math.MaxFloat64, -math.MaxFloat64
	for _, qual := range stats.PositionalQuality {
		if qual < minQual {
			minQual = qual
		}
		if qual > maxQual {
			maxQual = qual
		}
	}

	fmt.Fprintf(writer, "Positional Quality Scores:\n")
	fmt.Fprintf(writer, "%-10s | %-10s | %s\n", "Position", "Avg Quality", "Graph")
	fmt.Fprintf(writer, "%s\n", "------------------------------------------------------------")

	for i, qual := range stats.PositionalQuality {
		barWidth := int((qual - minQual) / (maxQual - minQual) * 40)
		bar := ""
		for j := 0; j < barWidth; j++ {
			bar += "█"
		}
		fmt.Fprintf(writer, "%-10d | %10.2f | %s\n", i+1, qual, bar)
	}

	return nil
}

// SVG generation functions
func generateLengthSVG(stats *Stats, writer io.Writer, x, y int) error {
	lengths := make([]int, 0, len(stats.LengthDistribution))
	for length := range stats.LengthDistribution {
		lengths = append(lengths, length)
	}
	sort.Ints(lengths)

	maxCount := 0
	for _, count := range stats.LengthDistribution {
		if count > maxCount {
			maxCount = count
		}
	}

	width := 500
	height := 300
	
	fmt.Fprintf(writer, `  <g transform="translate(%d,%d)">
    <text class="title" x="%d" y="0">Length Distribution</text>
`, x, y, width/2)

	// Draw axes
	fmt.Fprintf(writer, `    <line x1="0" y1="%d" x2="%d" y2="%d" stroke="black" stroke-width="1"/>
    <line x1="0" y1="0" x2="0" y2="%d" stroke="black" stroke-width="1"/>
`, height, width, height, height)

	// Draw bars
	if len(lengths) > 0 {
		barWidth := float64(width) / float64(len(lengths))
		for i, length := range lengths {
			count := stats.LengthDistribution[length]
			barHeight := float64(count) / float64(maxCount) * float64(height-20)
			fmt.Fprintf(writer, `    <rect class="bar" x="%.2f" y="%.2f" width="%.2f" height="%.2f"/>
`, float64(i)*barWidth, float64(height)-barHeight, barWidth*0.8, barHeight)
		}
	}

	fmt.Fprintf(writer, "  </g>\n")
	return nil
}

func generateQualitySVG(stats *Stats, writer io.Writer, x, y int) error {
	qualities := make([]int, 0, len(stats.QualityDistribution))
	for qual := range stats.QualityDistribution {
		qualities = append(qualities, qual)
	}
	sort.Ints(qualities)

	maxCount := 0
	for _, count := range stats.QualityDistribution {
		if count > maxCount {
			maxCount = count
		}
	}

	width := 500
	height := 300
	
	fmt.Fprintf(writer, `  <g transform="translate(%d,%d)">
    <text class="title" x="%d" y="0">Quality Score Distribution</text>
`, x, y, width/2)

	// Draw axes
	fmt.Fprintf(writer, `    <line x1="0" y1="%d" x2="%d" y2="%d" stroke="black" stroke-width="1"/>
    <line x1="0" y1="0" x2="0" y2="%d" stroke="black" stroke-width="1"/>
`, height, width, height, height)

	// Draw bars
	if len(qualities) > 0 {
		barWidth := float64(width) / float64(len(qualities))
		for i, qual := range qualities {
			count := stats.QualityDistribution[qual]
			barHeight := float64(count) / float64(maxCount) * float64(height-20)
			fmt.Fprintf(writer, `    <rect class="bar" x="%.2f" y="%.2f" width="%.2f" height="%.2f"/>
`, float64(i)*barWidth, float64(height)-barHeight, barWidth*0.8, barHeight)
		}
	}

	fmt.Fprintf(writer, "  </g>\n")
	return nil
}

func generatePositionalQualitySVG(stats *Stats, writer io.Writer, x, y int) error {
	width := 1100
	height := 300
	
	fmt.Fprintf(writer, `  <g transform="translate(%d,%d)">
    <text class="title" x="%d" y="0">Positional Quality Scores</text>
`, x, y, width/2)

	// Draw axes
	fmt.Fprintf(writer, `    <line x1="0" y1="%d" x2="%d" y2="%d" stroke="black" stroke-width="1"/>
    <line x1="0" y1="0" x2="0" y2="%d" stroke="black" stroke-width="1"/>
`, height, width, height, height)

	// Find min/max for scaling
	minQual, maxQual := math.MaxFloat64, -math.MaxFloat64
	for _, qual := range stats.PositionalQuality {
		if qual < minQual {
			minQual = qual
		}
		if qual > maxQual {
			maxQual = qual
		}
	}

	// Draw line graph
	if len(stats.PositionalQuality) > 0 {
		points := "M"
		stepX := float64(width) / float64(len(stats.PositionalQuality))
		for i, qual := range stats.PositionalQuality {
			y := float64(height) - ((qual - minQual) / (maxQual - minQual) * float64(height-20))
			if i == 0 {
				fmt.Fprintf(writer, "%s%.2f,%.2f", points, float64(i)*stepX, y)
			} else {
				fmt.Fprintf(writer, " L%.2f,%.2f", float64(i)*stepX, y)
			}
		}
		fmt.Fprintf(writer, `" class="line"/>
`)
	}

	fmt.Fprintf(writer, "  </g>\n")
	return nil
}
