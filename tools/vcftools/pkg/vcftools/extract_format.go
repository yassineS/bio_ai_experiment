package vcftools

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// extractFormatRunner streams variants and writes
// `<prefix>.<NAME>.FORMAT` with one row per kept site: CHROM, POS, then
// one column per kept sample carrying that sample's value for the named
// FORMAT field. Ported from upstream
// variant_file_format_convert.cpp:1204-1263 (output_FORMAT_information):
//
//   - Sites whose FORMAT column does not list NAME are skipped entirely
//     (upstream line 1247-1248).
//   - For sites that do carry NAME, each sample's value is taken from
//     its colon-separated value vector at NAME's position. A sample
//     whose vector is too short to reach that index emits a literal "."
//     (upstream vcf_entry.cpp:618 + the early `break` at line 637).
//   - The header row is `CHROM\tPOS` followed by one column per sample
//     (the post-filter sample list, matching `include_indv` at
//     upstream line 1227).
type extractFormatRunner struct {
	w        *bufio.Writer
	f        io.WriteCloser
	formatID string
	samples  []string
}

// newExtractFormatRunner opens `<prefix>.<formatID>.FORMAT` and writes
// the header row. The caller must ensure formatID is non-empty.
func newExtractFormatRunner(prefix, formatID string, samples []string) (*extractFormatRunner, error) {
	if formatID == "" {
		return nil, fmt.Errorf("--extract-FORMAT-info: empty FORMAT name")
	}
	f, err := iohelper.OpenWriter(prefix + "." + formatID + ".FORMAT")
	if err != nil {
		return nil, fmt.Errorf("opening %s.%s.FORMAT: %w", prefix, formatID, err)
	}
	w := bufio.NewWriter(f)
	if _, err := w.WriteString("CHROM\tPOS"); err != nil {
		_ = f.Close()
		return nil, err
	}
	for _, s := range samples {
		if _, err := w.WriteString("\t" + s); err != nil {
			_ = f.Close()
			return nil, err
		}
	}
	if _, err := w.WriteString("\n"); err != nil {
		_ = f.Close()
		return nil, err
	}
	return &extractFormatRunner{
		w:        w,
		f:        f,
		formatID: formatID,
		samples:  append([]string(nil), samples...),
	}, nil
}

// addVariant emits one row for v when its FORMAT column lists the
// runner's tag. Otherwise the variant is skipped (matches upstream's
// `if (e->FORMAT_id_exists(FORMAT_id) == false) continue;` at
// variant_file_format_convert.cpp:1247).
func (r *extractFormatRunner) addVariant(v *vcf.Variant) error {
	if r == nil {
		return nil
	}
	if !formatContains(v.Format, r.formatID) {
		return nil
	}
	if _, err := fmt.Fprintf(r.w, "%s\t%d", v.Chrom, v.Pos); err != nil {
		return err
	}
	for _, sample := range v.Samples {
		// Sample.Data is populated only for FORMAT keys whose colon-
		// separated value index existed in the sample string. A
		// missing key therefore means "value vector too short" — match
		// upstream's "." sentinel (vcf_entry.cpp:618).
		val, ok := sample.Data[r.formatID]
		if !ok || val == "" {
			val = "."
		}
		if _, err := r.w.WriteString("\t" + val); err != nil {
			return err
		}
	}
	if _, err := r.w.WriteString("\n"); err != nil {
		return err
	}
	return nil
}

// close flushes the buffered writer and closes the underlying file.
// Safe to call on a nil runner.
func (r *extractFormatRunner) close() error {
	if r == nil {
		return nil
	}
	if err := r.w.Flush(); err != nil {
		_ = r.f.Close()
		return err
	}
	return r.f.Close()
}

// formatContains is defined in beagle.go (shared helper).
