// Package giab implements the GIAB variant-calling concordance harness
// (manuscript experiment C2/P1). It drives two call sets — one produced with
// OUR bcftools, one with the UPSTREAM bcftools — from the same reference and
// reads, then reports:
//
//  1. Ours-vs-upstream concordance, record-by-record, restricted to the GIAB
//     high-confidence BED. This is the byte-/record-exact claim. Crucially it
//     includes a ULP-flip detector: where QUAL or PL differ by only the last
//     unit in the last place (Phred ULP), it verifies the genotype (GT) and the
//     FILTER (PASS/FAIL) are unchanged, so the libm floor is shown to be a
//     non-issue.
//  2. Biological concordance vs the GIAB truth set, via hap.py or RTG vcfeval
//     when one is available, stratified by the GA4GH/GIAB stratification BEDs.
//
// Every external prerequisite (GIAB VCF/BED, reads BAM, hap.py/vcfeval, our and
// upstream bcftools) is checked up front; a missing prerequisite produces a
// clear SKIP with the reason and a pointer to docs/GIAB_CONCORDANCE.md rather
// than a hard error, so the harness runs (doing nothing) on a machine with no
// GIAB data. This mirrors how pipeline/internal/upstream and pipeline/roundtrip
// degrade.
package giab

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

// VCFRecord is a minimal parsed VCF data row — only the fields the concordance
// comparator needs. INFO and the raw line are retained verbatim so callers can
// surface exact context.
type VCFRecord struct {
	Chrom  string
	Pos    int
	ID     string
	Ref    string
	Alt    string
	Qual   string // kept as the raw token; "." means missing
	Filter string
	Info   string
	Format string   // the FORMAT field (e.g. "GT:PL:DP")
	Sample []string // per-sample columns (FORMAT-keyed values)
	Line   string   // the original tab-joined data line
}

// Key returns the locus identity (CHROM:POS:REF:ALT) used to pair records
// across the two call sets.
func (r VCFRecord) Key() string {
	return r.Chrom + ":" + strconv.Itoa(r.Pos) + ":" + r.Ref + ":" + r.Alt
}

// IsSNV reports whether the record is a single-nucleotide variant (REF and ALT
// both single bases, no symbolic ALT). Multi-allelic ALTs (comma-separated) are
// treated as indels for stratification convenience only by this helper; callers
// needing exact typing should split first.
func (r VCFRecord) IsSNV() bool {
	if len(r.Ref) != 1 {
		return false
	}
	for _, a := range strings.Split(r.Alt, ",") {
		if len(a) != 1 || a == "*" || strings.ContainsAny(a, "<>[]") {
			return false
		}
	}
	return true
}

// gtField returns the value for a FORMAT subfield key (e.g. "GT" or "PL") in
// the first sample column, or "" if absent. GIAB single-sample call sets have
// exactly one sample; the first column is the one that matters.
func (r VCFRecord) gtField(key string) string {
	if len(r.Sample) == 0 || r.Format == "" {
		return ""
	}
	keys := strings.Split(r.Format, ":")
	vals := strings.Split(r.Sample[0], ":")
	for i, k := range keys {
		if k == key {
			if i < len(vals) {
				return vals[i]
			}
			return ""
		}
	}
	return ""
}

// GT returns the genotype call (the GT subfield), normalised so phasing does
// not spuriously register as a difference: "0|1" and "0/1" compare equal, and
// the allele indices are sorted. An unphased het and a phased het with the same
// alleles are the same genotype for concordance purposes.
func (r VCFRecord) GT() string {
	raw := r.gtField("GT")
	if raw == "" {
		return ""
	}
	alleles := strings.FieldsFunc(raw, func(c rune) bool { return c == '|' || c == '/' })
	sort.Strings(alleles)
	return strings.Join(alleles, "/")
}

// PL returns the PL subfield (phred-scaled genotype likelihoods), or "".
func (r VCFRecord) PL() string { return r.gtField("PL") }

// PassFail collapses the FILTER column to a boolean PASS verdict. "PASS" and
// "." (unfiltered, conventionally treated as pass by GIAB tooling) are PASS;
// anything else is a fail. The raw filter is preserved on the record.
func (r VCFRecord) PassFail() bool {
	f := strings.TrimSpace(r.Filter)
	return f == "PASS" || f == "." || f == ""
}

// ParseVCF reads a (plain-text) VCF stream and returns its data records. Header
// lines are skipped. It tolerates trailing whitespace and blank lines. The
// caller is responsible for decompression (feed a decompressed stream).
func ParseVCF(r io.Reader) ([]VCFRecord, error) {
	var recs []VCFRecord
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 1<<20), 1<<26)
	for sc.Scan() {
		line := sc.Text()
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) < 8 {
			// Not a well-formed data row; skip rather than fail the whole parse.
			continue
		}
		pos, err := strconv.Atoi(cols[1])
		if err != nil {
			return nil, fmt.Errorf("bad POS %q: %w", cols[1], err)
		}
		rec := VCFRecord{
			Chrom:  cols[0],
			Pos:    pos,
			ID:     cols[2],
			Ref:    cols[3],
			Alt:    cols[4],
			Qual:   cols[5],
			Filter: cols[6],
			Info:   cols[7],
			Line:   line,
		}
		if len(cols) >= 9 {
			rec.Format = cols[8]
			rec.Sample = cols[9:]
		}
		recs = append(recs, rec)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return recs, nil
}

// ParseVCFFile parses a plain-text VCF file by path.
func ParseVCFFile(path string) ([]VCFRecord, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ParseVCF(f)
}

// qualToFloat parses a QUAL token, returning (value, ok). "." (missing) is ok
// with value NaN handled by the caller.
func qualToFloat(s string) (float64, bool) {
	s = strings.TrimSpace(s)
	if s == "" || s == "." {
		return 0, false
	}
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		return 0, false
	}
	return v, true
}

// floatsClose reports whether two floats are within a tolerance that captures a
// last-place (ULP/Phred) difference. QUAL is conventionally emitted to two
// decimals by bcftools; a one-ULP Phred difference shows up as a sub-0.5
// wobble. We treat anything within absTol (default below) OR within a relative
// epsilon as "ULP-level".
func floatsClose(a, b, absTol float64) bool {
	if a == b {
		return true
	}
	diff := math.Abs(a - b)
	if diff <= absTol {
		return true
	}
	// Relative tolerance for larger magnitudes.
	scale := math.Max(math.Abs(a), math.Abs(b))
	return diff <= scale*1e-6
}
