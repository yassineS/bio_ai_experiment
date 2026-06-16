// split-vep feature machinery ported from split-vep.c: the EXPRESSION transcript
// selector (init_select_tr_expr / get_matching_transcript), the -g/--gene-list
// restriction and prioritisation (init_gene_list / restrict_csqs_to_genes), the
// file-based severity-scale override (the args->severity branch of init_data),
// and the --columns-types regex table (init_column2type / get_column_type). See
// native_plugin_splitvep.go for the plugin lifecycle and option parsing.
package bcftools

import (
	"fmt"
	"os"
	"regexp"
	"strings"
)

// svDefaultSeverityText is the default consequence severity scale, byte-for-byte
// from default_severity() in split-vep.c. It is the text printed by -S - / -S ?
// and the source parsed by initSeverityScale when no -S FILE is given.
const svDefaultSeverityText = "# Default consequence substrings ordered in ascending order by severity.\n" +
	"# Consequences with the same severity can be put on the same line in arbitrary order.\n" +
	"# See also https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html\n" +
	"intergenic\n" +
	"feature_truncation feature_elongation\n" +
	"regulatory\n" +
	"TF_binding_site TFBS\n" +
	"downstream upstream\n" +
	"non_coding_transcript non_coding\n" +
	"intron NMD_transcript\n" +
	"non_coding_transcript_exon\n" +
	"5_prime_utr 3_prime_utr\n" +
	"coding_sequence mature_miRNA\n" +
	"stop_retained start_retained synonymous\n" +
	"incomplete_terminal_codon\n" +
	"splice_region\n" +
	"missense inframe protein_altering\n" +
	"transcript_amplification\n" +
	"exon_loss\n" +
	"disruptive\n" +
	"start_lost stop_lost stop_gained frameshift\n" +
	"splice_acceptor splice_donor\n" +
	"transcript_ablation\n"

// svDefaultColumnTypesText is the default CSQ subfield type table, byte-for-byte
// from default_column_types() in split-vep.c. It is the text printed by
// --columns-types - and the source parsed by initColumn2Type when no override
// FILE is given. The leading column is a regex (anchored with ^...$) and the
// second column the type.
const svDefaultColumnTypesText = "# Default CSQ subfield types, unlisted fields are type String.\n" +
	"# Note that the name search is done using regular expressions, with\n" +
	"# \"^\" and \"$\" appended automatically\n" +
	"DISTANCE                   Integer\n" +
	"STRAND                     Integer\n" +
	"TSL                        Integer\n" +
	"GENE_PHENO                 Integer\n" +
	"HGVS_OFFSET                Integer\n" +
	".*_POPS                    String\n" +
	"AF                         Float\n" +
	".*_AF                      Float\n" +
	"MAX_AF_.*                  Float\n" +
	"MOTIF_POS                  Integer\n" +
	"MOTIF_SCORE_CHANGE         Float\n" +
	"existing_InFrame_oORFs     Integer\n" +
	"existing_OutOfFrame_oORFs  Integer\n" +
	"existing_uORFs             Integer\n" +
	"SpliceAI_pred_DP_.*        Integer\n" +
	"SpliceAI_pred_DS_.*        Float\n"

// initSelectTrExpr parses a --select EXPRESSION (primary/pick/mane expand to
// CANONICAL=YES / PICK=1 / MANE_SELECT!="" before this is called), porting
// init_select_tr_expr. The field is resolved through the annot-prefix and must
// exist in INFO/<tag>; the operator is one of =,!=,~,!~, with the right-hand
// value optionally double-quoted. Regex operators compile the value as a POSIX
// extended regular expression, matching upstream's regcomp.
func (p *splitVepPlugin) initSelectTrExpr(expr string) error {
	var field, value string
	var op int
	switch {
	case strings.Contains(expr, "!="):
		i := strings.Index(expr, "!=")
		field, value, op = expr[:i], expr[i+2:], svTrOpNe
	case strings.Contains(expr, "!~"):
		i := strings.Index(expr, "!~")
		field, value, op = expr[:i], expr[i+2:], svTrOpNr
	case strings.IndexByte(expr, '=') >= 0:
		i := strings.IndexByte(expr, '=')
		field, value, op = expr[:i], expr[i+1:], svTrOpEq
	case strings.IndexByte(expr, '~') >= 0:
		i := strings.IndexByte(expr, '~')
		field, value, op = expr[:i], expr[i+1:], svTrOpRe
	default:
		return fmt.Errorf("Could not parse the expression: -s %s", expr)
	}
	value = svUnquote(value)
	field = p.sanitizeField(strings.TrimPrefix(field, p.annotPrefix))
	idx, ok := p.field2idx[field]
	if !ok {
		return fmt.Errorf("The field \"%s\" was requested via \"%s\" but it is not present in INFO/%s", field, expr, p.vepTag)
	}
	spec := svTrExprSpec{field: field, idx: idx, op: op, value: value}
	if op == svTrOpRe || op == svTrOpNr {
		re, err := regexp.Compile(value)
		if err != nil {
			return fmt.Errorf("Error: fail to compile the regular expression \"%s\"", value)
		}
		spec.regex = re
	}
	p.trExpr = spec
	p.selTr = svTrExpr
	return nil
}

// svUnquote strips one matching pair of surrounding double quotes, mirroring the
// is_quoted handling in init_select_tr_expr.
func svUnquote(s string) string {
	if len(s) >= 2 && s[0] == '"' && s[len(s)-1] == '"' {
		return s[1 : len(s)-1]
	}
	return s
}

// matchingTranscripts returns the indices of the transcripts that satisfy the
// EXPRESSION selector, porting get_matching_transcript. A transcript with too
// few columns is an error, matching upstream.
func (p *splitVepPlugin) matchingTranscripts(transcripts [][]string) ([]int, error) {
	var out []int
	for i, tr := range transcripts {
		if p.trExpr.idx >= len(tr) {
			return nil, fmt.Errorf("Too few columns: %d (field %s) >= %d", p.trExpr.idx, p.trExpr.field, len(tr))
		}
		val := tr[p.trExpr.idx]
		var match bool
		switch p.trExpr.op {
		case svTrOpEq:
			match = p.trExpr.value == val
		case svTrOpNe:
			match = p.trExpr.value != val
		case svTrOpRe:
			match = p.trExpr.regex.MatchString(val)
		case svTrOpNr:
			match = !p.trExpr.regex.MatchString(val)
		}
		if match {
			out = append(out, i)
		}
	}
	return out, nil
}

// initGeneList loads the -g/--gene-list file (one gene per line), records the
// restrict-vs-prioritise mode (a leading '+' on the filename selects prioritise),
// and resolves --gene-list-fields (default SYMBOL,Gene,gene) to CSQ subfield
// indices. Ports init_gene_list.
func (p *splitVepPlugin) initGeneList() error {
	spec := p.genesFname
	p.genesMode = svGenesRestrict
	if strings.HasPrefix(spec, "+") {
		p.genesMode = svGenesPrioritize
		spec = spec[1:]
	}
	data, err := os.ReadFile(spec)
	if err != nil {
		return fmt.Errorf("Could not read the file provided with --gene-list %s", spec)
	}
	p.genes = map[string]bool{}
	for _, line := range strings.Split(string(data), "\n") {
		line = strings.TrimRight(line, "\r")
		if line == "" {
			continue
		}
		p.genes[line] = true
	}
	if len(p.genes) == 0 {
		return fmt.Errorf("Could not read the file provided with --gene-list %s", spec)
	}

	fieldsStr := p.geneFields
	if fieldsStr == "" {
		fieldsStr = "SYMBOL,Gene,gene"
	}
	p.geneFieldIdx = nil
	for _, raw := range strings.Split(fieldsStr, ",") {
		raw = strings.TrimSpace(raw)
		if raw == "" {
			continue
		}
		name := p.sanitizeField(raw)
		if idx, ok := p.field2idx[name]; ok {
			p.geneFieldIdx = append(p.geneFieldIdx, idx)
		}
	}
	if len(p.geneFieldIdx) == 0 {
		return fmt.Errorf("None of the \"%s\" fields is present in INFO/%s", fieldsStr, p.vepTag)
	}
	return nil
}

// restrictCsqsToGenes filters/reorders the transcript list by the -g gene set,
// porting restrict_csqs_to_genes. In restrict mode only transcripts whose gene
// matches survive (none => empty); in prioritise mode all survive but matched
// transcripts are partitioned to the front using the same two-pointer swap as
// upstream. It returns the (possibly reordered) transcript slice and its
// retained length.
func (p *splitVepPlugin) restrictCsqsToGenes(transcripts [][]string) ([][]string, int) {
	n := len(transcripts)
	hit := make([]bool, n)
	nhit := 0
	for i, tr := range transcripts {
		matched := false
		for _, idx := range p.geneFieldIdx {
			if idx >= len(tr) {
				continue
			}
			if p.genes[tr[idx]] {
				matched = true
				break
			}
		}
		hit[i] = matched
		if matched {
			nhit++
		}
	}
	if nhit == 0 {
		if p.genesMode == svGenesRestrict {
			return transcripts, 0 // no gene of interest
		}
		return transcripts, n
	}
	// Partition the hits to the front with upstream's two-pointer swap.
	i, j := 0, n-1
	for i < j {
		if hit[i] {
			i++
			continue
		}
		if !hit[j] {
			j--
			continue
		}
		transcripts[i], transcripts[j] = transcripts[j], transcripts[i]
		hit[i], hit[j] = hit[j], hit[i]
		i++
		j--
	}
	return transcripts, nhit
}

// initColumn2Type builds the regex column-type table, porting init_column2type.
// The source is the --columns-types FILE when given, otherwise the built-in
// default_column_types text. Each non-comment line is "<regex><whitespace><Type>";
// the regex is anchored with ^...$ and the type must be one of
// Float/Integer/Flag/String.
func (p *splitVepPlugin) initColumn2Type() error {
	var text string
	if p.columnTypes != "" && p.columnTypes != "-" {
		data, err := os.ReadFile(p.columnTypes)
		if err != nil {
			return fmt.Errorf("Cannot read %s", p.columnTypes)
		}
		text = string(data)
	} else {
		text = svDefaultColumnTypesText
	}
	var table []svCol2Type
	for _, line := range strings.Split(text, "\n") {
		if line == "" || line[0] == '#' {
			continue
		}
		// Split into the leading non-space token (the regex) and the trailing
		// non-space token (the type), exactly as upstream scans.
		var pat, typ string
		k := 0
		for k < len(line) && !isSpaceByte(line[k]) {
			k++
		}
		pat = line[:k]
		for k < len(line) && isSpaceByte(line[k]) {
			k++
		}
		typ = strings.TrimRight(line[k:], " \t\r")
		if pat == "" || typ == "" {
			return fmt.Errorf("Error: failed to parse the column type \"%s\"", line)
		}
		t, ok := svColumnTypeCode(typ)
		if !ok {
			return fmt.Errorf("Error: the column type \"%s\" is not supported: %s", typ, line)
		}
		re, err := regexp.Compile("^" + pat + "$")
		if err != nil {
			return fmt.Errorf("Error: fail to compile the column type regular expression \"^%s$\": %s", pat, line)
		}
		table = append(table, svCol2Type{regex: re, typ: t})
	}
	if len(table) == 0 {
		return fmt.Errorf("Failed to parse the column types")
	}
	p.column2type = table
	return nil
}

// svColumnTypeCode maps a column-type spelling to a value type. Flag maps to the
// String renderer here (the plugin emits a Flag header but treats the value as a
// string for output, matching upstream's handling for split-vep output where
// only Integer/Float drive numeric re-parsing).
func svColumnTypeCode(s string) (int, bool) {
	switch s {
	case "Float":
		return svTypeReal, true
	case "Integer":
		return svTypeInt, true
	case "String":
		return svTypeStr, true
	case "Flag":
		return svTypeStr, true
	}
	return 0, false
}

// isSpaceByte reports whether c is an ASCII whitespace byte, matching isspace_c.
func isSpaceByte(c byte) bool {
	return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f'
}
