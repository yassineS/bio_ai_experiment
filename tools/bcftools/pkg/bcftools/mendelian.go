package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// Trio identifies one (child, father, mother) sample tuple for the
// mendelian-error scan.
type Trio struct {
	Child  string
	Father string
	Mother string
}

// MendelianMode captures the -m/--mode behaviour. The upstream plugin
// exposes a comma-separated set; we accept a single character per call
// in v1 — the most common choice. Modes that affect output content
// (annotate / delete / x-chrom) are implemented; others are accepted
// for CLI compatibility.
type MendelianMode int

const (
	// MendelianAnnotate adds INFO/MERR to every record. The default.
	MendelianAnnotate MendelianMode = iota
	// MendelianCount emits a tab-separated trio-level summary instead of VCF.
	MendelianCount
	// MendelianDelete drops records where at least one trio has a
	// Mendel error.
	MendelianDelete
	// MendelianXChrom flags X-chromosome handling (paternal allele
	// ignored on chrX outside PAR for male children). The v1 port
	// detects the contig name heuristically (`X`, `chrX`) and treats
	// the father as haploid there; the rest of the pipeline is
	// unchanged.
	MendelianXChrom
	// MendelianPlusPG is upstream's `+` mode: keep all records, also
	// append a ##bcftools_PG header line. v1 treats it as a synonym
	// for MendelianAnnotate (the PG line is provenance, not parity-
	// critical).
	MendelianPlusPG
)

// ParseMendelianMode parses the one-letter mode string from -m/--mode.
func ParseMendelianMode(s string) (MendelianMode, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "", "a":
		return MendelianAnnotate, nil
	case "c":
		return MendelianCount, nil
	case "x":
		return MendelianXChrom, nil
	case "d":
		return MendelianDelete, nil
	case "+":
		return MendelianPlusPG, nil
	}
	return 0, fmt.Errorf("bcftools mendelian: unknown -m mode %q (accept a|c|x|d|+)", s)
}

// MendelianOptions controls the behaviour of MendelianFile.
type MendelianOptions struct {
	// Trios are the (child, father, mother) tuples to check. At least
	// one must be supplied.
	Trios []Trio
	// TrioFile is read in addition to Trios; lines are comma- or
	// tab-separated CHILD,FATHER,MOTHER triplets. Blank lines and
	// '#' comments are ignored.
	TrioFile string
	// Mode controls how errors surface in the output stream.
	Mode MendelianMode
	// Count is the legacy short-flag form of MendelianCount (-c).
	// When both Count and Mode are set, Count wins.
	Count bool
	// Delete is the legacy short-flag form of MendelianDelete (-d).
	// When both Delete and Mode are set, Delete wins.
	Delete bool
	// RulesFile, when set, is read for per-contig ploidy overrides
	// (CHROM <tab> from <tab> to <tab> sex <tab> ploidy). v1
	// accepts the path but only honours the heuristic X-chrom
	// behaviour for the male child case (matching MendelianXChrom).
	RulesFile string
	// OutputFormat selects between VCF/VCF.gz/BCF for the streaming
	// output paths (a/d/x/+). Ignored when Mode is MendelianCount.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
}

// MendelianFile is the file-aware entry point used by the CLI. It opens
// path through iohelper (transparent gzip + BCF auto-detect), walks
// every record, evaluates the per-trio consistency rule, and emits the
// requested output.
func MendelianFile(path string, out io.Writer, opts MendelianOptions) (MendelianSummary, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return MendelianSummary{}, fmt.Errorf("bcftools mendelian: open %s: %w", path, err)
	}
	defer in.Close()
	return Mendelian(in, out, opts)
}

// MendelianSummary is the per-trio + per-input rollup returned by every
// Mendelian* entry point. It is the data backing -c summary output.
type MendelianSummary struct {
	// Trios contains the per-trio error counts in the order supplied
	// in opts.Trios + opts.TrioFile.
	Trios []TrioStats
	// TotalRecords is the number of variant records seen.
	TotalRecords int
	// RecordsWithError is the number of records where at least one
	// trio had a Mendel error.
	RecordsWithError int
}

// TrioStats holds the rollup numbers for a single (child, father,
// mother) triple.
type TrioStats struct {
	Trio
	NTested  int // records where all three samples had complete diploid GTs
	NError   int // records that fail Mendelian inheritance for this trio
	NMissing int // records skipped because at least one sample was missing/uncalled
}

// Mendelian streams VCF/BCF input through the trio rule and writes the
// requested output to w. The streaming-vs-buffered split mirrors
// Convert / ConcatFiles for consistency.
func Mendelian(in io.Reader, out io.Writer, opts MendelianOptions) (MendelianSummary, error) {
	// Legacy short-flag overrides win when set.
	if opts.Count {
		opts.Mode = MendelianCount
	} else if opts.Delete {
		opts.Mode = MendelianDelete
	}

	trios, err := loadTrios(opts.Trios, opts.TrioFile)
	if err != nil {
		return MendelianSummary{}, fmt.Errorf("bcftools mendelian: %w", err)
	}
	if len(trios) == 0 {
		return MendelianSummary{}, fmt.Errorf("bcftools mendelian: at least one trio (-t or -T) is required")
	}

	hdr, variants, err := readAllVariants(in)
	if err != nil {
		return MendelianSummary{}, err
	}

	// Resolve trio sample-name → header index once per file.
	indices, err := resolveTrioIndices(hdr, trios)
	if err != nil {
		return MendelianSummary{}, fmt.Errorf("bcftools mendelian: %w", err)
	}

	summary := MendelianSummary{
		Trios: make([]TrioStats, len(trios)),
	}
	for i := range trios {
		summary.Trios[i].Trio = trios[i]
	}

	// Inject the INFO/MERR meta line when we emit VCF.
	annotatedHdr := hdr
	if opts.Mode != MendelianCount {
		annotatedHdr = withMERRHeader(hdr)
	}

	var writer variantWriter
	var finish func()
	if opts.Mode != MendelianCount {
		writer, finish, err = openOutput(out, ViewOptions{
			OutputFormat:  opts.OutputFormat,
			CompressLevel: opts.CompressLevel,
		}, annotatedHdr)
		if err != nil {
			return summary, fmt.Errorf("bcftools mendelian: %w", err)
		}
		defer finish()
		if err := writer.WriteHeader(); err != nil {
			return summary, err
		}
	}

	for _, v := range variants {
		summary.TotalRecords++
		totalErrors := 0
		for ti, idx := range indices {
			child, father, mother, complete := readTrioGenotypes(v, idx)
			if !complete {
				summary.Trios[ti].NMissing++
				continue
			}
			summary.Trios[ti].NTested++

			// X-chromosome handling: if mode is MendelianXChrom and
			// chrom looks like chrX/X, treat the father as haploid.
			haploidFather := opts.Mode == MendelianXChrom && isXChrom(v.Chrom)
			if !mendelianConsistent(child, father, mother, haploidFather) {
				summary.Trios[ti].NError++
				totalErrors++
			}
		}
		if totalErrors > 0 {
			summary.RecordsWithError++
		}

		switch opts.Mode {
		case MendelianCount:
			// no output per record; we'll emit a summary at the end.
		case MendelianDelete:
			if totalErrors == 0 {
				if err := writer.Write(v); err != nil {
					return summary, err
				}
			}
		default:
			// Annotate / X-chrom / Plus-PG all write every record,
			// optionally with INFO/MERR set.
			if v.Info == nil {
				v.Info = make(map[string]string)
			}
			if _, exists := v.Info["MERR"]; !exists {
				v.InfoOrder = append(v.InfoOrder, "MERR")
			}
			v.Info["MERR"] = strconv.Itoa(totalErrors)
			if err := writer.Write(v); err != nil {
				return summary, err
			}
		}
	}

	if writer != nil {
		if err := writer.Flush(); err != nil {
			return summary, err
		}
	}

	if opts.Mode == MendelianCount {
		if err := writeMendelianSummary(out, summary); err != nil {
			return summary, err
		}
	}
	return summary, nil
}

// loadTrios merges the in-options trio slice with the -T trio-file
// contents and validates them.
func loadTrios(inline []Trio, path string) ([]Trio, error) {
	out := append([]Trio{}, inline...)
	if path == "" {
		return out, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		// Accept comma-separated or whitespace-separated triplets.
		var fields []string
		if strings.Contains(line, ",") {
			fields = strings.Split(line, ",")
		} else {
			fields = strings.Fields(line)
		}
		if len(fields) < 3 {
			return nil, fmt.Errorf("trio-file %s: malformed line %q", path, line)
		}
		out = append(out, Trio{
			Child:  strings.TrimSpace(fields[0]),
			Father: strings.TrimSpace(fields[1]),
			Mother: strings.TrimSpace(fields[2]),
		})
	}
	return out, sc.Err()
}

// ParseTrioFlag splits a single -t/--trio CHILD,FATHER,MOTHER argument.
func ParseTrioFlag(s string) (Trio, error) {
	parts := strings.Split(s, ",")
	if len(parts) != 3 {
		return Trio{}, fmt.Errorf("bcftools mendelian: bad -t %q (need CHILD,FATHER,MOTHER)", s)
	}
	t := Trio{
		Child:  strings.TrimSpace(parts[0]),
		Father: strings.TrimSpace(parts[1]),
		Mother: strings.TrimSpace(parts[2]),
	}
	if t.Child == "" || t.Father == "" || t.Mother == "" {
		return Trio{}, fmt.Errorf("bcftools mendelian: bad -t %q (empty member)", s)
	}
	return t, nil
}

// trioIndex holds the resolved sample positions for one trio.
type trioIndex struct {
	child, father, mother int
}

// resolveTrioIndices maps each trio's sample names to their position in
// the input header. A missing name produces a descriptive error.
func resolveTrioIndices(hdr *vcf.Header, trios []Trio) ([]trioIndex, error) {
	byName := make(map[string]int, len(hdr.Samples))
	for i, s := range hdr.Samples {
		byName[s] = i
	}
	out := make([]trioIndex, len(trios))
	for i, t := range trios {
		cidx, ok := byName[t.Child]
		if !ok {
			return nil, fmt.Errorf("trio %d: child %q not in input", i+1, t.Child)
		}
		fidx, ok := byName[t.Father]
		if !ok {
			return nil, fmt.Errorf("trio %d: father %q not in input", i+1, t.Father)
		}
		midx, ok := byName[t.Mother]
		if !ok {
			return nil, fmt.Errorf("trio %d: mother %q not in input", i+1, t.Mother)
		}
		out[i] = trioIndex{child: cidx, father: fidx, mother: midx}
	}
	return out, nil
}

// readTrioGenotypes pulls the three GTs for the trio at variant v. The
// returned slices are 2-element diploid genotypes (paternal/maternal
// allele indexes); missing or no-call samples return complete == false.
func readTrioGenotypes(v *vcf.Variant, idx trioIndex) (child, father, mother []int, complete bool) {
	c, ok := parseTrioGT(sampleData(v, idx.child))
	if !ok {
		return nil, nil, nil, false
	}
	f, ok := parseTrioGT(sampleData(v, idx.father))
	if !ok {
		return nil, nil, nil, false
	}
	m, ok := parseTrioGT(sampleData(v, idx.mother))
	if !ok {
		return nil, nil, nil, false
	}
	return c, f, m, true
}

// sampleData returns the GT field of the i-th sample, or "" when the
// index is out of range.
func sampleData(v *vcf.Variant, i int) string {
	if i < 0 || i >= len(v.Samples) {
		return ""
	}
	gt, ok := v.Samples[i].Data["GT"]
	if !ok {
		return ""
	}
	return gt
}

// parseTrioGT splits "0/1" / "0|1" / "1/." into two integer allele
// indexes. Returns ok == false if the genotype contains a missing
// allele (".") or is not diploid.
func parseTrioGT(gt string) ([]int, bool) {
	if gt == "" || gt == "." {
		return nil, false
	}
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	if len(parts) != 2 {
		return nil, false
	}
	out := make([]int, 2)
	for i, p := range parts {
		if p == "" || p == "." {
			return nil, false
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, false
		}
		out[i] = n
	}
	return out, true
}

// mendelianConsistent returns true if the child's two alleles can each
// be sourced from exactly one parent. When haploidFather is true (X-
// chromosome mode for male children) the father contributes only one
// allele, so we only require the maternal side to match.
func mendelianConsistent(child, father, mother []int, haploidFather bool) bool {
	if len(child) != 2 || len(father) != 2 || len(mother) != 2 {
		return false
	}
	// Mendelian rule: for some ordering of the child's alleles, one
	// allele equals one of the father's and the other equals one of
	// the mother's.
	fset := map[int]bool{father[0]: true, father[1]: true}
	mset := map[int]bool{mother[0]: true, mother[1]: true}
	if haploidFather {
		// On chrX a male child inherits an X from the mother only;
		// the paternal contribution is the Y (i.e. effectively
		// nothing testable). Treat the child as if it were a single
		// allele drawn from the mother.
		return mset[child[0]] && mset[child[1]]
	}
	// Either (c0 from father AND c1 from mother) OR (c0 from mother AND c1 from father).
	if fset[child[0]] && mset[child[1]] {
		return true
	}
	if mset[child[0]] && fset[child[1]] {
		return true
	}
	return false
}

// isXChrom returns true for contigs that look like the human X
// chromosome.
func isXChrom(name string) bool {
	switch strings.ToUpper(name) {
	case "X", "CHRX":
		return true
	}
	return false
}

// withMERRHeader returns a copy of hdr that contains an
// ##INFO=<ID=MERR,...> line if it doesn't already.
func withMERRHeader(hdr *vcf.Header) *vcf.Header {
	if hdr == nil {
		return hdr
	}
	for _, m := range hdr.MetaInfo {
		if strings.HasPrefix(m, "##INFO=<ID=MERR,") {
			return hdr
		}
	}
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	out.MetaInfo = append(out.MetaInfo,
		`##INFO=<ID=MERR,Number=1,Type=Integer,Description="Number of trios with a Mendelian inheritance error at this site (bcftools mendelian)">`)
	return out
}

// writeMendelianSummary emits the human-readable rollup used by the
// -c/--count mode. Format mirrors the upstream plugin's "summary" mode
// (TSV with a leading header row).
func writeMendelianSummary(out io.Writer, s MendelianSummary) error {
	w := bufio.NewWriter(out)
	defer w.Flush()
	if _, err := fmt.Fprintf(w, "# trio\tchild\tfather\tmother\tn_tested\tn_error\tn_missing\n"); err != nil {
		return err
	}
	for i, t := range s.Trios {
		if _, err := fmt.Fprintf(w, "%d\t%s\t%s\t%s\t%d\t%d\t%d\n",
			i+1, t.Child, t.Father, t.Mother, t.NTested, t.NError, t.NMissing); err != nil {
			return err
		}
	}
	if _, err := fmt.Fprintf(w, "# totals\trecords=%d\twith_error=%d\n",
		s.TotalRecords, s.RecordsWithError); err != nil {
		return err
	}
	return nil
}
