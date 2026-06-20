package giab

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// Engine identifies the GA4GH benchmarking engine used for biological
// concordance against the GIAB truth set.
type Engine string

const (
	// EngineHappy is Illumina's hap.py (typically run via the pkrusche/hap.py
	// container). It writes a <prefix>.summary.csv.
	EngineHappy Engine = "hap.py"
	// EngineVcfeval is RTG Tools' vcfeval. It writes a summary.txt and, with
	// --output-mode=ga4gh / rtg vcfeval ... a roc/summary; the harness parses
	// its fixed-width summary.txt.
	EngineVcfeval Engine = "vcfeval"
)

// BenchMetrics holds precision/recall/F1 for one variant type in one stratum
// for one call set. Counts are carried through when the engine reports them so
// downstream tooling can recompute or weight.
type BenchMetrics struct {
	VarType   string  `json:"var_type"` // "SNP"/"SNV" or "INDEL"
	Stratum   string  `json:"stratum"`  // "*" (whole high-conf) or a stratification name
	Recall    float64 `json:"recall"`
	Precision float64 `json:"precision"`
	F1        float64 `json:"f1"`
	TruthTP   int     `json:"truth_tp,omitempty"`
	TruthFN   int     `json:"truth_fn,omitempty"`
	QueryFP   int     `json:"query_fp,omitempty"`
}

// ParseHappySummary parses a hap.py <prefix>.summary.csv stream into per-type
// metrics for the whole region (stratum "*"). hap.py emits one row per
// (Type, Filter) pair; the harness keeps the PASS rows (Filter == "PASS"),
// falling back to ALL if no PASS row exists.
//
// The relevant columns are: Type, Filter, METRIC.Recall, METRIC.Precision,
// METRIC.F1_Score, TRUTH.TP, TRUTH.FN, QUERY.FP.
func ParseHappySummary(r io.Reader) ([]BenchMetrics, error) {
	cr := csv.NewReader(r)
	cr.FieldsPerRecord = -1
	rows, err := cr.ReadAll()
	if err != nil {
		return nil, err
	}
	if len(rows) < 2 {
		return nil, fmt.Errorf("hap.py summary has no data rows")
	}
	idx := map[string]int{}
	for i, h := range rows[0] {
		idx[strings.TrimSpace(h)] = i
	}
	need := []string{"Type", "Filter", "METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score"}
	for _, c := range need {
		if _, ok := idx[c]; !ok {
			return nil, fmt.Errorf("hap.py summary missing column %q", c)
		}
	}
	get := func(row []string, col string) string {
		i, ok := idx[col]
		if !ok || i >= len(row) {
			return ""
		}
		return strings.TrimSpace(row[i])
	}

	// Prefer PASS rows; remember ALL rows as a fallback per type.
	pass := map[string]BenchMetrics{}
	all := map[string]BenchMetrics{}
	for _, row := range rows[1:] {
		typ := normType(get(row, "Type"))
		if typ == "" {
			continue
		}
		m := BenchMetrics{
			VarType:   typ,
			Stratum:   "*",
			Recall:    parseFloat(get(row, "METRIC.Recall")),
			Precision: parseFloat(get(row, "METRIC.Precision")),
			F1:        parseFloat(get(row, "METRIC.F1_Score")),
			TruthTP:   parseInt(get(row, "TRUTH.TP")),
			TruthFN:   parseInt(get(row, "TRUTH.FN")),
			QueryFP:   parseInt(get(row, "QUERY.FP")),
		}
		switch get(row, "Filter") {
		case "PASS":
			pass[typ] = m
		case "ALL":
			all[typ] = m
		default:
			if _, ok := all[typ]; !ok {
				all[typ] = m
			}
		}
	}
	var out []BenchMetrics
	for _, typ := range []string{"SNP", "INDEL"} {
		if m, ok := pass[typ]; ok {
			out = append(out, m)
		} else if m, ok := all[typ]; ok {
			out = append(out, m)
		}
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("hap.py summary contained no SNP/INDEL rows")
	}
	return out, nil
}

// ParseVcfevalSummary parses an RTG vcfeval summary.txt stream. vcfeval writes
// a fixed-width table with a header, a separator line of dashes, the per-
// threshold rows, and a final "None"/"Total" row. The harness reads the final
// (whole-callset) row and returns a single combined metric (vcfeval does not
// split SNP/indel in summary.txt — that requires --output-mode=ga4gh, which
// hap.py-style parsing covers). The returned slice has one entry with VarType
// "ALL".
//
// Columns (RTG Tools): Threshold, True-pos-baseline, True-pos-call,
// False-pos, False-neg, Precision, Sensitivity, F-measure. "Sensitivity" is
// recall.
func ParseVcfevalSummary(r io.Reader) ([]BenchMetrics, error) {
	data, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(string(data), "\n")
	var header []string
	var dataRows [][]string
	for _, ln := range lines {
		ln = strings.TrimRight(ln, "\r")
		if strings.TrimSpace(ln) == "" {
			continue
		}
		if strings.HasPrefix(strings.TrimSpace(ln), "---") {
			continue
		}
		fields := strings.Fields(ln)
		if header == nil {
			header = fields
			continue
		}
		dataRows = append(dataRows, fields)
	}
	if header == nil || len(dataRows) == 0 {
		return nil, fmt.Errorf("vcfeval summary has no data rows")
	}
	idx := map[string]int{}
	for i, h := range header {
		idx[h] = i
	}
	// The summary's last row is the aggregate (Threshold "None").
	row := dataRows[len(dataRows)-1]
	get := func(col string) string {
		i, ok := idx[col]
		if !ok || i >= len(row) {
			return ""
		}
		return row[i]
	}
	m := BenchMetrics{
		VarType:   "ALL",
		Stratum:   "*",
		Precision: parseFloat(get("Precision")),
		Recall:    parseFloat(get("Sensitivity")),
		F1:        parseFloat(get("F-measure")),
		TruthTP:   parseInt(get("True-pos-baseline")),
		QueryFP:   parseInt(get("False-pos")),
		TruthFN:   parseInt(get("False-neg")),
	}
	if m.Recall == 0 && get("Sensitivity") == "" {
		// Try alternative RTG header spelling.
		m.Recall = parseFloat(get("Recall"))
	}
	return []BenchMetrics{m}, nil
}

// ParseHappySummaryFile parses a hap.py summary.csv file by path.
func ParseHappySummaryFile(path string) ([]BenchMetrics, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ParseHappySummary(f)
}

// ParseVcfevalSummaryFile parses an RTG vcfeval summary.txt file by path.
func ParseVcfevalSummaryFile(path string) ([]BenchMetrics, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ParseVcfevalSummary(f)
}

// normType maps hap.py's "SNP"/"INDEL" type tokens to our canonical labels,
// returning "" for anything else.
func normType(s string) string {
	switch strings.ToUpper(strings.TrimSpace(s)) {
	case "SNP", "SNV", "SNPS":
		return "SNP"
	case "INDEL", "INDELS":
		return "INDEL"
	default:
		return ""
	}
}

func parseFloat(s string) float64 {
	s = strings.TrimSpace(s)
	if s == "" || s == "." || strings.EqualFold(s, "NA") {
		return 0
	}
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0
	}
	return v
}

func parseInt(s string) int {
	s = strings.TrimSpace(s)
	if s == "" || s == "." {
		return 0
	}
	// Some engines print counts as floats (e.g. "123.0").
	if f, err := strconv.ParseFloat(s, 64); err == nil {
		return int(f)
	}
	return 0
}
