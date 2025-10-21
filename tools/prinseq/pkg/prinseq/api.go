package prinseq

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
)

// APIServer provides HTTP endpoints for PRINSEQ operations
type APIServer struct {
	addr string
}

// NewAPIServer creates a new API server
func NewAPIServer(addr string) *APIServer {
	return &APIServer{addr: addr}
}

// Start starts the API server
func (s *APIServer) Start() error {
	http.HandleFunc("/api/stats", s.handleStats)
	http.HandleFunc("/api/filter", s.handleFilter)
	http.HandleFunc("/api/benchmark", s.handleBenchmark)
	http.HandleFunc("/api/report", s.handleReport)
	http.HandleFunc("/api/graph", s.handleGraph)
	http.HandleFunc("/health", s.handleHealth)
	http.HandleFunc("/", s.handleIndex)

	fmt.Printf("Starting PRINSEQ API server on %s\n", s.addr)
	return http.ListenAndServe(s.addr, nil)
}

// handleStats processes statistics requests
func (s *APIServer) handleStats(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	// Read request body
	data, err := io.ReadAll(r.Body)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error reading request: %v", err), http.StatusBadRequest)
		return
	}
	defer r.Body.Close()

	// Determine format from query parameter
	isFastq := r.URL.Query().Get("format") == "fastq"
	enhanced := r.URL.Query().Get("enhanced") == "true"

	var stats *Stats
	if enhanced {
		stats, err = CalculateEnhancedStats(bytes.NewReader(data), isFastq)
	} else {
		stats, err = CalculateStats(bytes.NewReader(data), isFastq)
	}

	if err != nil {
		http.Error(w, fmt.Sprintf("Error calculating stats: %v", err), http.StatusInternalServerError)
		return
	}

	// Return JSON response
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(stats)
}

// handleFilter processes filtering requests
func (s *APIServer) handleFilter(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	// Parse filter options from query parameters
	opts := FilterOptions{
		MinLen:      parseIntQuery(r, "min_len", 0),
		MaxLen:      parseIntQuery(r, "max_len", 0),
		MinGC:       parseFloat64Query(r, "min_gc", 0),
		MaxGC:       parseFloat64Query(r, "max_gc", 0),
		MinQualMean: parseFloat64Query(r, "min_qual", 0),
		MaxNsP:      parseFloat64Query(r, "max_ns_p", 0),
		MaxNsN:      parseIntQuery(r, "max_ns_n", 0),
	}

	// Read request body
	data, err := io.ReadAll(r.Body)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error reading request: %v", err), http.StatusBadRequest)
		return
	}
	defer r.Body.Close()

	// Determine format
	isFastq := r.URL.Query().Get("format") == "fastq"

	// Filter sequences
	var output bytes.Buffer
	err = Filter(bytes.NewReader(data), &output, isFastq, opts)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error filtering: %v", err), http.StatusInternalServerError)
		return
	}

	// Return filtered sequences
	if isFastq {
		w.Header().Set("Content-Type", "text/plain")
	} else {
		w.Header().Set("Content-Type", "text/plain")
	}
	w.Write(output.Bytes())
}

// handleBenchmark runs benchmark suite
func (s *APIServer) handleBenchmark(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	// Read request body
	data, err := io.ReadAll(r.Body)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error reading request: %v", err), http.StatusBadRequest)
		return
	}
	defer r.Body.Close()

	isFastq := r.URL.Query().Get("format") == "fastq"

	// Run benchmark suite
	results, err := RunBenchmarkSuite(bytes.NewReader(data), isFastq)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error running benchmark: %v", err), http.StatusInternalServerError)
		return
	}

	// Return JSON response
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(results)
}

// handleReport generates HTML report
func (s *APIServer) handleReport(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	// Read request body
	data, err := io.ReadAll(r.Body)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error reading request: %v", err), http.StatusBadRequest)
		return
	}
	defer r.Body.Close()

	isFastq := r.URL.Query().Get("format") == "fastq"

	// Calculate enhanced stats
	stats, err := CalculateEnhancedStats(bytes.NewReader(data), isFastq)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error calculating stats: %v", err), http.StatusInternalServerError)
		return
	}

	// Generate HTML report
	var report bytes.Buffer
	err = GenerateHTMLReport(stats, &report)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error generating report: %v", err), http.StatusInternalServerError)
		return
	}

	// Return HTML
	w.Header().Set("Content-Type", "text/html")
	w.Write(report.Bytes())
}

// handleGraph generates graph
func (s *APIServer) handleGraph(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "Method not allowed", http.StatusMethodNotAllowed)
		return
	}

	// Read request body
	data, err := io.ReadAll(r.Body)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error reading request: %v", err), http.StatusBadRequest)
		return
	}
	defer r.Body.Close()

	isFastq := r.URL.Query().Get("format") == "fastq"
	graphType := GraphType(r.URL.Query().Get("type"))
	if graphType == "" {
		graphType = GraphTypeLength
	}

	// Calculate enhanced stats
	stats, err := CalculateEnhancedStats(bytes.NewReader(data), isFastq)
	if err != nil {
		http.Error(w, fmt.Sprintf("Error calculating stats: %v", err), http.StatusInternalServerError)
		return
	}

	// Generate graph
	var graph bytes.Buffer
	if r.URL.Query().Get("svg") == "true" {
		err = GenerateSVG(stats, &graph)
		w.Header().Set("Content-Type", "image/svg+xml")
	} else {
		err = GenerateGraph(stats, graphType, &graph)
		w.Header().Set("Content-Type", "text/plain")
	}

	if err != nil {
		http.Error(w, fmt.Sprintf("Error generating graph: %v", err), http.StatusInternalServerError)
		return
	}

	w.Write(graph.Bytes())
}

// handleHealth returns health status
func (s *APIServer) handleHealth(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(map[string]string{
		"status": "healthy",
		"version": "1.1.0",
	})
}

// handleIndex returns API documentation
func (s *APIServer) handleIndex(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "text/html")
	fmt.Fprintf(w, `<!DOCTYPE html>
<html>
<head>
    <title>PRINSEQ API</title>
    <style>
        body { font-family: sans-serif; max-width: 800px; margin: 40px auto; padding: 20px; }
        h1 { color: #333; }
        code { background: #f4f4f4; padding: 2px 6px; border-radius: 3px; }
        pre { background: #f4f4f4; padding: 15px; border-radius: 5px; overflow-x: auto; }
        .endpoint { margin: 20px 0; padding: 15px; border-left: 3px solid #4CAF50; background: #f9f9f9; }
    </style>
</head>
<body>
    <h1>🧬 PRINSEQ API Server</h1>
    <p>Version 1.1.0</p>
    
    <h2>Endpoints</h2>
    
    <div class="endpoint">
        <h3>POST /api/stats</h3>
        <p>Calculate sequence statistics</p>
        <p>Query parameters:</p>
        <ul>
            <li><code>format</code>: "fastq" or "fasta"</li>
            <li><code>enhanced</code>: "true" for detailed statistics</li>
        </ul>
        <pre>curl -X POST -d @sequences.fastq "http://localhost:8080/api/stats?format=fastq&enhanced=true"</pre>
    </div>
    
    <div class="endpoint">
        <h3>POST /api/filter</h3>
        <p>Filter sequences based on criteria</p>
        <p>Query parameters:</p>
        <ul>
            <li><code>format</code>: "fastq" or "fasta"</li>
            <li><code>min_len</code>: minimum length</li>
            <li><code>max_len</code>: maximum length</li>
            <li><code>min_gc</code>: minimum GC content</li>
            <li><code>max_gc</code>: maximum GC content</li>
            <li><code>min_qual</code>: minimum quality score</li>
        </ul>
        <pre>curl -X POST -d @sequences.fastq "http://localhost:8080/api/filter?format=fastq&min_len=100&min_qual=20"</pre>
    </div>
    
    <div class="endpoint">
        <h3>POST /api/benchmark</h3>
        <p>Run performance benchmark</p>
        <p>Query parameters:</p>
        <ul>
            <li><code>format</code>: "fastq" or "fasta"</li>
        </ul>
        <pre>curl -X POST -d @sequences.fastq "http://localhost:8080/api/benchmark?format=fastq"</pre>
    </div>
    
    <div class="endpoint">
        <h3>POST /api/report</h3>
        <p>Generate HTML quality report</p>
        <p>Query parameters:</p>
        <ul>
            <li><code>format</code>: "fastq" or "fasta"</li>
        </ul>
        <pre>curl -X POST -d @sequences.fastq "http://localhost:8080/api/report?format=fastq" > report.html</pre>
    </div>
    
    <div class="endpoint">
        <h3>POST /api/graph</h3>
        <p>Generate quality graphs</p>
        <p>Query parameters:</p>
        <ul>
            <li><code>format</code>: "fastq" or "fasta"</li>
            <li><code>type</code>: "length", "gc", "quality", "dinucleotides", "positional_quality"</li>
            <li><code>svg</code>: "true" for SVG output</li>
        </ul>
        <pre>curl -X POST -d @sequences.fastq "http://localhost:8080/api/graph?format=fastq&svg=true" > graph.svg</pre>
    </div>
    
    <div class="endpoint">
        <h3>GET /health</h3>
        <p>Health check endpoint</p>
        <pre>curl http://localhost:8080/health</pre>
    </div>
</body>
</html>`)
}

// Helper functions for parsing query parameters
func parseIntQuery(r *http.Request, key string, defaultValue int) int {
	val := r.URL.Query().Get(key)
	if val == "" {
		return defaultValue
	}
	var result int
	fmt.Sscanf(val, "%d", &result)
	return result
}

func parseFloat64Query(r *http.Request, key string, defaultValue float64) float64 {
	val := r.URL.Query().Get(key)
	if val == "" {
		return defaultValue
	}
	var result float64
	fmt.Sscanf(val, "%f", &result)
	return result
}
