package prinseq

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"time"
)

// GenerateHTMLReport generates an HTML report with embedded graphs
func GenerateHTMLReport(stats *Stats, writer io.Writer) error {
	// Generate SVG graphs
	var svgBuf bytes.Buffer
	if err := GenerateSVG(stats, &svgBuf); err != nil {
		return fmt.Errorf("error generating SVG: %w", err)
	}

	// Convert stats to JSON for embedding
	statsJSON, err := json.MarshalIndent(stats, "    ", "  ")
	if err != nil {
		return fmt.Errorf("error marshaling stats: %w", err)
	}

	// Write HTML
	_, err = fmt.Fprintf(writer, `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PRINSEQ Quality Control Report</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
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
            padding-bottom: 5px;
        }
        .summary {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        .stat-card {
            background: linear-gradient(135deg, #667eea 0%%, #764ba2 100%%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .stat-card.green {
            background: linear-gradient(135deg, #4CAF50 0%%, #45a049 100%%);
        }
        .stat-card.blue {
            background: linear-gradient(135deg, #2196F3 0%%, #1976D2 100%%);
        }
        .stat-card.orange {
            background: linear-gradient(135deg, #FF9800 0%%, #F57C00 100%%);
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
        .graphs {
            margin: 30px 0;
        }
        .graph-container {
            margin: 20px 0;
            text-align: center;
        }
        .json-view {
            background: #f8f8f8;
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 15px;
            overflow-x: auto;
        }
        pre {
            margin: 0;
            font-family: 'Courier New', monospace;
            font-size: 12px;
        }
        .footer {
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            text-align: center;
            color: #777;
            font-size: 14px;
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
        }
        tr:hover {
            background-color: #f5f5f5;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 PRINSEQ Quality Control Report</h1>
        
        <div class="summary">
            <div class="stat-card green">
                <div class="stat-label">Total Reads</div>
                <div class="stat-value">%d</div>
            </div>
            <div class="stat-card blue">
                <div class="stat-label">Total Bases</div>
                <div class="stat-value">%d</div>
            </div>
            <div class="stat-card orange">
                <div class="stat-label">Average Length</div>
                <div class="stat-value">%.1f</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">GC Content</div>
                <div class="stat-value">%.1f%%</div>
            </div>
        </div>

        <h2>📊 Summary Statistics</h2>
        <table>
            <tr>
                <th>Metric</th>
                <th>Value</th>
            </tr>
            <tr>
                <td>Number of Reads</td>
                <td>%d</td>
            </tr>
            <tr>
                <td>Total Bases</td>
                <td>%d</td>
            </tr>
            <tr>
                <td>Minimum Length</td>
                <td>%d</td>
            </tr>
            <tr>
                <td>Maximum Length</td>
                <td>%d</td>
            </tr>
            <tr>
                <td>Average Length</td>
                <td>%.2f</td>
            </tr>
            <tr>
                <td>GC Content</td>
                <td>%.2f%%</td>
            </tr>
            <tr>
                <td>Number of Ns</td>
                <td>%d</td>
            </tr>
`,
		stats.NumReads,
		stats.TotalBases,
		stats.AvgLength,
		stats.GCContent,
		stats.NumReads,
		stats.TotalBases,
		stats.MinLength,
		stats.MaxLength,
		stats.AvgLength,
		stats.GCContent,
		stats.NumNs,
	)

	if err != nil {
		return err
	}

	if stats.AvgQuality > 0 {
		_, err = fmt.Fprintf(writer, `            <tr>
                <td>Average Quality Score</td>
                <td>%.2f</td>
            </tr>
`, stats.AvgQuality)
		if err != nil {
			return err
		}
	}

	_, err = fmt.Fprintf(writer, `        </table>

        <h2>📈 Graphs and Distributions</h2>
        <div class="graphs">
            <div class="graph-container">
                %s
            </div>
        </div>

        <h2>📄 Detailed Statistics (JSON)</h2>
        <div class="json-view">
            <pre>%s</pre>
        </div>

        <div class="footer">
            Generated by PRINSEQ v1.1.0 on %s
        </div>
    </div>
</body>
</html>
`, svgBuf.String(), string(statsJSON), time.Now().Format("2006-01-02 15:04:05"))

	return err
}
