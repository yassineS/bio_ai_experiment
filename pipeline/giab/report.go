package giab

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

// WriteReports writes giab_concordance.json and giab_concordance.md into dir,
// creating dir if needed. It returns the two written paths.
func WriteReports(rep *Report, dir string) (jsonPath, mdPath string, err error) {
	if dir == "" {
		dir = "."
	}
	if err = os.MkdirAll(dir, 0o755); err != nil {
		return "", "", err
	}
	jsonPath = filepath.Join(dir, "giab_concordance.json")
	mdPath = filepath.Join(dir, "giab_concordance.md")

	b, err := json.MarshalIndent(rep, "", "  ")
	if err != nil {
		return "", "", err
	}
	b = append(b, '\n')
	if err := os.WriteFile(jsonPath, b, 0o644); err != nil {
		return "", "", err
	}
	if err := os.WriteFile(mdPath, []byte(rep.Markdown()), 0o644); err != nil {
		return "", "", err
	}
	return jsonPath, mdPath, nil
}

// Markdown renders the human-readable report.
func (rep *Report) Markdown() string {
	var b strings.Builder
	w := func(format string, args ...any) { fmt.Fprintf(&b, format+"\n", args...) }

	w("# GIAB variant-calling concordance: %s", orNone(rep.Sample))
	w("")
	if rep.Build != "" {
		w("- Reference build: `%s`", rep.Build)
	}
	w("- Generated: %s", rep.Generated.Format("2006-01-02T15:04:05Z"))
	w("- Call sets: ours=`%s` upstream=`%s`", orNone(rep.OursVCF), orNone(rep.UpVCF))
	w("")
	w("| Stage | Status | Detail |")
	w("|-------|--------|--------|")
	w("| Produce call sets | %s | %s |", rep.CallStatus, mdCell(rep.CallDetail))
	w("| Ours-vs-upstream concordance | %s | %s |", rep.ConcordanceStatus, mdCell(rep.ConcordanceDetail))
	w("| Biological concordance (truth) | %s | %s |", rep.BiologicalStatus, mdCell(rep.BiologicalDetail))
	w("")

	// Ours-vs-upstream concordance table + ULP-flip result.
	w("## Ours vs upstream (record-exact within the high-confidence BED)")
	w("")
	if rep.Concordance == nil {
		w("_Not run (%s)._", mdCell(rep.ConcordanceDetail))
		w("")
	} else {
		c := rep.Concordance
		w("| Metric | Count |")
		w("|--------|-------|")
		w("| Common sites | %d |", c.Common)
		w("| Identical | %d |", c.Identical)
		w("| Differ | %d |", c.Differ)
		w("| ... differ only at QUAL/PL ULP floor | %d |", c.QualULPOnly)
		w("| ... flip a genotype or PASS/FAIL | **%d** |", c.GenotypeOrFilterFlips)
		w("| Only in ours | %d |", c.OnlyOurs)
		w("| Only in upstream | %d |", c.OnlyUp)
		w("")
		w("**ULP-flip result:** %d site(s) differ in QUAL/PL only at the ULP/Phred "+
			"floor, %d of which flip a genotype or PASS/FAIL verdict.",
			c.QualULPOnly, c.GenotypeOrFilterFlips)
		w("")
		if len(c.Diffs) > 0 {
			w("### Differing sites (first %d)", len(c.Diffs))
			w("")
			w("| CHROM | POS | ours QUAL | up QUAL | ours GT | up GT | ours PASS | up PASS | ULP-only | flips | note |")
			w("|-------|-----|-----------|---------|---------|-------|-----------|---------|----------|-------|------|")
			for _, d := range c.Diffs {
				w("| %s | %d | %s | %s | %s | %s | %t | %t | %t | %t | %s |",
					d.Chrom, d.Pos, mdCell(d.OursQual), mdCell(d.UpQual),
					mdCell(d.OursGT), mdCell(d.UpGT), d.OursPass, d.UpPass,
					d.QualULP, d.Flips(), mdCell(d.Note))
			}
			w("")
		}
	}

	// Biological concordance per stratum, ours vs upstream side by side.
	w("## Biological concordance vs GIAB truth (ours vs upstream)")
	w("")
	if len(rep.Biological) == 0 {
		w("_Not run (%s)._", mdCell(rep.BiologicalDetail))
		w("")
	} else {
		for _, er := range rep.Biological {
			w("### Stratum: `%s` (engine: %s)", er.Stratum, er.Engine)
			w("")
			if er.Status != StatusPass {
				w("_%s: %s_", er.Status, mdCell(er.Detail))
				w("")
				continue
			}
			w("| Var type | Recall (ours / up) | Precision (ours / up) | F1 (ours / up) |")
			w("|----------|--------------------|------------------------|----------------|")
			byType := pairByType(er.Ours, er.Up)
			for _, vt := range orderedTypes(byType) {
				p := byType[vt]
				w("| %s | %s / %s | %s / %s | %s / %s |", vt,
					f3(p.ours.Recall), f3(p.up.Recall),
					f3(p.ours.Precision), f3(p.up.Precision),
					f3(p.ours.F1), f3(p.up.F1))
			}
			w("")
		}
	}

	if rep.Failed() {
		w("> **FAILED:** a substantive discrepancy was found (see the flip count " +
			"and/or biological stage detail above).")
	} else {
		w("> All stages PASS or SKIP; no genotype/PASS-FAIL flip and no benchmarking error.")
	}
	w("")
	return b.String()
}

type metricPair struct{ ours, up BenchMetrics }

// pairByType joins ours/up metric slices by var type.
func pairByType(ours, up []BenchMetrics) map[string]metricPair {
	m := map[string]metricPair{}
	for _, o := range ours {
		p := m[o.VarType]
		p.ours = o
		m[o.VarType] = p
	}
	for _, u := range up {
		p := m[u.VarType]
		p.up = u
		m[u.VarType] = p
	}
	return m
}

// orderedTypes returns var types in a stable, human-friendly order.
func orderedTypes(m map[string]metricPair) []string {
	var out []string
	for _, t := range []string{"SNP", "INDEL", "ALL"} {
		if _, ok := m[t]; ok {
			out = append(out, t)
		}
	}
	for t := range m {
		if t != "SNP" && t != "INDEL" && t != "ALL" {
			out = append(out, t)
		}
	}
	return out
}

// f3 formats a metric to three decimals.
func f3(v float64) string { return fmt.Sprintf("%.4f", v) }

// mdCell escapes a string for safe inclusion in a markdown table cell and
// renders empty as a dash.
func mdCell(s string) string {
	s = strings.TrimSpace(s)
	if s == "" {
		return "-"
	}
	s = strings.ReplaceAll(s, "|", "\\|")
	s = strings.ReplaceAll(s, "\n", " ")
	return s
}
