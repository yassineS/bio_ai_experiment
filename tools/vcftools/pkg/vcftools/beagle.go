// BEAGLE genotype-likelihood output for vcftools.
//
// vcftools supports two flavours of BEAGLE output:
//
//	--BEAGLE-GL → <prefix>.BEAGLE.GL — log10-scale GL triplets per sample,
//	               derived from the VCF FORMAT/PL field
//	--BEAGLE-PL → <prefix>.BEAGLE.PL — raw Phred PL triplets per sample
//
// In both cases the output uses BEAGLE's "marker" layout: a single header
// line, then one row per biallelic SNP with the format
//
//	marker allele1 allele2 <sample1_v1 sample1_v2 sample1_v3> ...
//
// where marker is "<CHR>:<POS>" (matching upstream vcftools), allele1 is
// REF, allele2 is the first ALT, and the three per-sample values are the
// genotype likelihoods for the AA / AB / BB genotypes.
//
// Records without a PL FORMAT field are skipped; we warn to stderr once at
// the start of the run (when at least one site lacked PL) rather than
// repeating the warning per site.
//
// PL → GL conversion follows the VCF spec: GL = -PL / 10 (so GL is a
// log10-probability with maximum 0). Missing per-sample PL values are
// emitted as the conventional placeholder triplet
//
//	GL: -0.481 -0.481 -0.481   (≈ log10(1/3); equally-likely genotypes)
//	PL:  0      0      0       (raw missing-encoded zeros)
//
// Sites with more than one ALT allele are skipped — BEAGLE format is
// biallelic. Indels and symbolic ALTs are also skipped. These choices match
// upstream vcftools' "snps with PL only" filter.
package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// log10OneThird is log10(1/3), used as the equally-likely-genotype GL
// placeholder for samples missing PL.
const log10OneThird = -0.47712125471966244

// beagleMode picks between GL (log10) and PL (raw Phred) output.
type beagleMode int

const (
	beagleGL beagleMode = iota
	beaglePL
)

// beagleWriter incrementally writes a BEAGLE genotype-likelihood file.
type beagleWriter struct {
	mode beagleMode
	f    io.WriteCloser
	w    *bufio.Writer

	// warnedNoPL is set once we've emitted the "no PL" warning to stderr.
	warnedNoPL bool

	// headerWritten is true after the first call to write() finishes.
	headerWritten bool
}

// newBEAGLEWriter creates the output file and prepares it for writes. It does
// not yet emit the header — that happens lazily on the first record so we can
// include the actual sample names from the (possibly filtered) header.
func newBEAGLEWriter(prefix string, mode beagleMode) (*beagleWriter, error) {
	suffix := ".BEAGLE.GL"
	if mode == beaglePL {
		suffix = ".BEAGLE.PL"
	}
	path := prefix + suffix
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return nil, fmt.Errorf("opening %s: %w", path, err)
	}
	return &beagleWriter{
		mode: mode,
		f:    f,
		w:    bufio.NewWriter(f),
	}, nil
}

// writeHeader emits the BEAGLE header row. samples are written in their
// header order; each sample contributes three columns (one per genotype).
func (b *beagleWriter) writeHeader(samples []string) error {
	if b.headerWritten {
		return nil
	}
	b.headerWritten = true
	fields := make([]string, 0, 3+3*len(samples))
	fields = append(fields, "marker", "allele1", "allele2")
	for _, s := range samples {
		fields = append(fields, s, s, s)
	}
	if _, err := fmt.Fprintln(b.w, strings.Join(fields, " ")); err != nil {
		return err
	}
	return nil
}

// write consumes a single variant. Sites without PL or with more than one
// ALT / non-SNP REF/ALT are skipped (the first PL-less site triggers a
// stderr warning once per run).
func (b *beagleWriter) write(v *vcf.Variant, samples []string) error {
	if err := b.writeHeader(samples); err != nil {
		return err
	}
	if len(v.Alt) != 1 {
		return nil
	}
	if !isSimpleSNP(v.Ref, v.Alt[0]) {
		return nil
	}
	if !formatContains(v.Format, "PL") {
		if !b.warnedNoPL {
			fmt.Fprintln(os.Stderr, "warning: VCF lacks a PL FORMAT field; --BEAGLE-GL/--BEAGLE-PL will skip sites without PL")
			b.warnedNoPL = true
		}
		return nil
	}

	marker := fmt.Sprintf("%s:%d", v.Chrom, v.Pos)
	fields := make([]string, 0, 3+3*len(samples))
	fields = append(fields, marker, v.Ref, v.Alt[0])

	// Index samples by name so we can emit them in header order even if the
	// variant's sample slice has been reordered.
	byName := make(map[string]*vcf.Sample, len(v.Samples))
	for i := range v.Samples {
		byName[v.Samples[i].Name] = &v.Samples[i]
	}

	for _, name := range samples {
		sample := byName[name]
		var pl [3]float64
		var have bool
		if sample != nil {
			pl, have = parseBiallelicPL(sample.Data["PL"])
		}
		if !have {
			// Missing PL → equally-likely placeholder.
			if b.mode == beagleGL {
				fields = append(fields, formatGL(log10OneThird), formatGL(log10OneThird), formatGL(log10OneThird))
			} else {
				fields = append(fields, "0", "0", "0")
			}
			continue
		}
		if b.mode == beagleGL {
			fields = append(fields, formatGL(-pl[0]/10), formatGL(-pl[1]/10), formatGL(-pl[2]/10))
		} else {
			fields = append(fields, formatPL(pl[0]), formatPL(pl[1]), formatPL(pl[2]))
		}
	}
	if _, err := fmt.Fprintln(b.w, strings.Join(fields, " ")); err != nil {
		return err
	}
	return nil
}

// close flushes the bufio.Writer and closes the underlying file. It is safe
// to call on a nil receiver to simplify the caller's deferred cleanup.
func (b *beagleWriter) close() error {
	if b == nil {
		return nil
	}
	var firstErr error
	if b.w != nil {
		// Lazy header: emit it even for empty outputs so downstream parsers
		// don't choke. We don't know samples here, so this is a no-op if
		// writeHeader was never reached — that's fine; an empty file is
		// also acceptable.
		if err := b.w.Flush(); err != nil {
			firstErr = err
		}
	}
	if b.f != nil {
		if err := b.f.Close(); err != nil && firstErr == nil {
			firstErr = err
		}
	}
	return firstErr
}

// parseBiallelicPL parses a comma-separated PL string. For a biallelic site
// we expect three values (AA, AB, BB); anything else (including "." or an
// empty string) is treated as missing.
func parseBiallelicPL(s string) ([3]float64, bool) {
	var zero [3]float64
	if s == "" || s == "." {
		return zero, false
	}
	parts := strings.Split(s, ",")
	if len(parts) != 3 {
		return zero, false
	}
	var out [3]float64
	for i, p := range parts {
		if p == "" || p == "." {
			return zero, false
		}
		v, err := strconv.ParseFloat(p, 64)
		if err != nil {
			return zero, false
		}
		out[i] = v
	}
	return out, true
}

// formatContains reports whether tag is present in the FORMAT list.
func formatContains(format []string, tag string) bool {
	for _, f := range format {
		if f == tag {
			return true
		}
	}
	return false
}

// isSimpleSNP returns true when REF and ALT are both single canonical bases
// (A/C/G/T/N). Indels, symbolic alleles, and breakends are excluded.
func isSimpleSNP(ref, alt string) bool {
	if len(ref) != 1 || len(alt) != 1 {
		return false
	}
	return isACGTN(ref[0]) && isACGTN(alt[0])
}

func isACGTN(b byte) bool {
	switch b {
	case 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n':
		return true
	}
	return false
}

// formatGL formats a log10-probability with vcftools-compatible precision
// (six decimals; bare zero rendered as "0").
func formatGL(v float64) string {
	if v == 0 {
		return "0"
	}
	return strconv.FormatFloat(v, 'f', 6, 64)
}

// formatPL formats a raw PL value. Integer values render without a trailing
// ".0" to mirror upstream vcftools.
func formatPL(v float64) string {
	if v == float64(int64(v)) {
		return strconv.FormatInt(int64(v), 10)
	}
	return strconv.FormatFloat(v, 'f', -1, 64)
}
