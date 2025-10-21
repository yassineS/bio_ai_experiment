package fastp

import (
	"fmt"
	"os"
	"time"
)

// GenerateHTMLReport generates an HTML report from processing statistics
func GenerateHTMLReport(stats *ProcessStats, opts ProcessOptions, outputPath string) error {
	f, err := os.Create(outputPath)
	if err != nil {
		return fmt.Errorf("failed to create HTML report: %w", err)
	}
	defer f.Close()

	html := generateHTMLContent(stats, opts)
	_, err = f.WriteString(html)
	return err
}

func generateHTMLContent(stats *ProcessStats, opts ProcessOptions) string {
	timestamp := time.Now().Format("2006-01-02 15:04:05")

	// Calculate percentages
	cleanPercent := 0.0
	if stats.TotalReads > 0 {
		cleanPercent = 100.0 * float64(stats.CleanReads) / float64(stats.TotalReads)
	}

	adapterPercent := 0.0
	if stats.TotalReads > 0 && stats.AdapterTrimmedReads > 0 {
		adapterPercent = 100.0 * float64(stats.AdapterTrimmedReads) / float64(stats.TotalReads)
	}

	polyGPercent := 0.0
	if stats.TotalReads > 0 && stats.PolyGTrimmedReads > 0 {
		polyGPercent = 100.0 * float64(stats.PolyGTrimmedReads) / float64(stats.TotalReads)
	}

	filteredPercent := 0.0
	if stats.TotalReads > 0 {
		filtered := stats.TotalReads - stats.CleanReads
		filteredPercent = 100.0 * float64(filtered) / float64(stats.TotalReads)
	}

	html := fmt.Sprintf(`<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Fastp Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            color: #333;
            border-bottom: 3px solid #4CAF50;
            padding-bottom: 10px;
        }
        h2 {
            color: #555;
            margin-top: 30px;
            border-bottom: 2px solid #ddd;
            padding-bottom: 8px;
        }
        .summary {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        .stat-box {
            background: linear-gradient(135deg, #667eea 0%%, #764ba2 100%%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        .stat-box.green {
            background: linear-gradient(135deg, #11998e 0%%, #38ef7d 100%%);
        }
        .stat-box.blue {
            background: linear-gradient(135deg, #4facfe 0%%, #00f2fe 100%%);
        }
        .stat-box.orange {
            background: linear-gradient(135deg, #fa709a 0%%, #fee140 100%%);
        }
        .stat-label {
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
        }
        .stat-value {
            font-size: 28px;
            font-weight: bold;
        }
        .stat-percent {
            font-size: 16px;
            opacity: 0.9;
            margin-top: 5px;
        }
        table {
            width: 100%%;
            border-collapse: collapse;
            margin: 20px 0;
        }
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #4CAF50;
            color: white;
            font-weight: bold;
        }
        tr:hover {
            background-color: #f5f5f5;
        }
        .footer {
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #777;
            font-size: 14px;
        }
        .progress-bar {
            width: 100%%;
            height: 25px;
            background-color: #e0e0e0;
            border-radius: 12px;
            overflow: hidden;
            margin: 10px 0;
        }
        .progress-fill {
            height: 100%%;
            background: linear-gradient(90deg, #4CAF50 0%%, #45a049 100%%);
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            font-size: 14px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Fastp Processing Report</h1>
        <p><strong>Generated:</strong> %s</p>
        
        <h2>Summary Statistics</h2>
        <div class="summary">
            <div class="stat-box green">
                <div class="stat-label">Total Reads</div>
                <div class="stat-value">%d</div>
            </div>
            <div class="stat-box blue">
                <div class="stat-label">Clean Reads</div>
                <div class="stat-value">%d</div>
                <div class="stat-percent">%.2f%%</div>
            </div>
            <div class="stat-box orange">
                <div class="stat-label">Filtered Reads</div>
                <div class="stat-value">%d</div>
                <div class="stat-percent">%.2f%%</div>
            </div>
            <div class="stat-box">
                <div class="stat-label">Total Bases</div>
                <div class="stat-value">%d</div>
            </div>
        </div>
        
        <h2>Read Quality</h2>
        <div class="progress-bar">
            <div class="progress-fill" style="width: %.2f%%">%.2f%% Clean</div>
        </div>
        
        <h2>Processing Details</h2>
        <table>
            <tr>
                <th>Metric</th>
                <th>Count</th>
                <th>Percentage</th>
            </tr>
            <tr>
                <td>Total reads processed</td>
                <td>%d</td>
                <td>100.00%%</td>
            </tr>
            <tr>
                <td>Clean reads (passed filters)</td>
                <td>%d</td>
                <td>%.2f%%</td>
            </tr>
            <tr>
                <td>Low quality reads filtered</td>
                <td>%d</td>
                <td>%.2f%%</td>
            </tr>
            <tr>
                <td>Too short reads filtered</td>
                <td>%d</td>
                <td>%.2f%%</td>
            </tr>
            <tr>
                <td>Too long reads filtered</td>
                <td>%d</td>
                <td>%.2f%%</td>
            </tr>
            <tr>
                <td>Too many N reads filtered</td>
                <td>%d</td>
                <td>%.2f%%</td>
            </tr>
        </table>
        
        <h2>Adapter and Trimming</h2>
        <table>
            <tr>
                <th>Operation</th>
                <th>Reads Affected</th>
                <th>Bases Removed</th>
                <th>Percentage</th>
            </tr>`,
		timestamp,
		stats.TotalReads,
		stats.CleanReads, cleanPercent,
		stats.TotalReads-stats.CleanReads, filteredPercent,
		stats.TotalBases,
		cleanPercent, cleanPercent,
		stats.TotalReads,
		stats.CleanReads, cleanPercent,
		stats.LowQualityReads, 100.0*float64(stats.LowQualityReads)/float64(stats.TotalReads),
		stats.TooShortReads, 100.0*float64(stats.TooShortReads)/float64(stats.TotalReads),
		stats.TooLongReads, 100.0*float64(stats.TooLongReads)/float64(stats.TotalReads),
		stats.TooManyNReads, 100.0*float64(stats.TooManyNReads)/float64(stats.TotalReads),
	)

	// Adapter trimming row
	if stats.AdapterTrimmedReads > 0 {
		html += fmt.Sprintf(`
            <tr>
                <td>Adapter trimming</td>
                <td>%d</td>
                <td>%d</td>
                <td>%.2f%%</td>
            </tr>`,
			stats.AdapterTrimmedReads,
			stats.AdapterTrimmedBases,
			adapterPercent,
		)
	}

	// Poly-G trimming row
	if stats.PolyGTrimmedReads > 0 {
		html += fmt.Sprintf(`
            <tr>
                <td>Poly-G tail trimming</td>
                <td>%d</td>
                <td>%d</td>
                <td>%.2f%%</td>
            </tr>`,
			stats.PolyGTrimmedReads,
			stats.PolyGTrimmedBases,
			polyGPercent,
		)
	}

	html += `
        </table>`

	// Add new features section
	if stats.DetectedAdapter != "" || stats.UMIExtracted > 0 || stats.BasesCorrected > 0 || stats.MergedReads > 0 {
		html += `
        <h2>Advanced Features</h2>
        <table>
            <tr>
                <th>Feature</th>
                <th>Value</th>
            </tr>`

		if stats.DetectedAdapter != "" {
			html += fmt.Sprintf(`
            <tr>
                <td>Auto-detected adapter</td>
                <td>%s</td>
            </tr>`, stats.DetectedAdapter)
		}

		if stats.UMIExtracted > 0 {
			html += fmt.Sprintf(`
            <tr>
                <td>UMIs extracted</td>
                <td>%d reads</td>
            </tr>`, stats.UMIExtracted)
		}

		if stats.BasesCorrected > 0 {
			html += fmt.Sprintf(`
            <tr>
                <td>Bases corrected</td>
                <td>%d bases</td>
            </tr>`, stats.BasesCorrected)
		}

		if stats.MergedReads > 0 {
			html += fmt.Sprintf(`
            <tr>
                <td>Overlapping reads merged</td>
                <td>%d pairs (%.2f%%)</td>
            </tr>`, stats.MergedReads, 100.0*float64(stats.MergedReads)/float64(stats.TotalReads/2))
		}

		html += `
        </table>`
	}

	// Configuration section
	html += `
        <h2>Configuration</h2>
        <table>
            <tr>
                <th>Parameter</th>
                <th>Value</th>
            </tr>`

	if opts.Adapter3 != "" {
		html += fmt.Sprintf(`
            <tr>
                <td>3' Adapter</td>
                <td>%s</td>
            </tr>`, opts.Adapter3)
	}

	if opts.Adapter5 != "" {
		html += fmt.Sprintf(`
            <tr>
                <td>5' Adapter</td>
                <td>%s</td>
            </tr>`, opts.Adapter5)
	}

	html += fmt.Sprintf(`
            <tr>
                <td>Quality threshold</td>
                <td>%d</td>
            </tr>
            <tr>
                <td>Minimum length</td>
                <td>%d</td>
            </tr>
            <tr>
                <td>Maximum N count</td>
                <td>%d</td>
            </tr>
            <tr>
                <td>Threads</td>
                <td>%d</td>
            </tr>`,
		opts.QualThreshold,
		opts.MinLength,
		opts.MaxNCount,
		opts.Threads,
	)

	if opts.TrimPolyG {
		html += `
            <tr>
                <td>Poly-G trimming</td>
                <td>Enabled</td>
            </tr>`
	}

	if opts.LowComplexity {
		html += fmt.Sprintf(`
            <tr>
                <td>Complexity filtering</td>
                <td>Enabled (threshold: %.2f)</td>
            </tr>`, opts.ComplexityThreshold)
	}

	if opts.BaseCorrection {
		html += `
            <tr>
                <td>Base correction</td>
                <td>Enabled</td>
            </tr>`
	}

	if opts.MergeOverlap {
		html += fmt.Sprintf(`
            <tr>
                <td>Overlap merging</td>
                <td>Enabled (min overlap: %d, max mismatch: %d)</td>
            </tr>`, opts.MinOverlap, opts.MaxMismatch)
	}

	html += `
        </table>
        
        <div class="footer">
            <p>Generated by <strong>fastp</strong> (Go implementation) - An all-in-one FASTQ preprocessor</p>
            <p>Visit <a href="https://github.com/yassineS/bio_ai_experiment">github.com/yassineS/bio_ai_experiment</a> for more information</p>
        </div>
    </div>
</body>
</html>`

	return html
}
