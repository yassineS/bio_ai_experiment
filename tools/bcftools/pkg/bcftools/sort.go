// Package bcftools — `bcftools sort` subcommand.
//
// `bcftools sort` re-orders the records of a VCF/BCF by (CHROM, POS), where
// the contig order is the one declared in the merged `##contig` lines. The
// upstream binary supports an external-merge sort with a configurable memory
// budget; for the v1 port we use an in-memory sort because the typical
// per-chromosome shard already fits comfortably in RAM. The CLI still accepts
// `-m/--max-mem` and `-T/--tmpdir` so existing scripts work; behaviour
// reduces to "load everything, sort, emit". The deviation is documented in
// docs/PARITY_ROADMAP.md.
package bcftools

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// SortOptions controls the `bcftools sort` subcommand. All fields are
// optional; defaults match the upstream binary.
type SortOptions struct {
	// OutputFormat selects the output encoding. Defaults to OutputVCF.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output (negative means
	// gzip's default).
	CompressLevel int
	// MaxMem is the upstream `-m/--max-mem` value (e.g. "768M"). Currently
	// accepted but not enforced; v1 always sorts in-memory.
	MaxMem string
	// TmpDir is the upstream `-T/--tmpdir` value. Currently accepted but
	// not used; v1 has no on-disk merge step.
	TmpDir string
}

// SortFile is the file-aware entry point for `bcftools sort`. It opens path
// through iohelper, fully reads the records, sorts them by (contig-order,
// POS, REF, ALT), and writes the result to out using the requested output
// format.
func SortFile(path string, out io.Writer, opts SortOptions) (int, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("bcftools sort: open %s: %w", path, err)
	}
	defer in.Close()
	return Sort(in, out, opts)
}

// Sort reads every record from in, sorts by (contig-order-from-header, POS,
// REF, ALT), and writes the sorted stream to out.
func Sort(in io.Reader, out io.Writer, opts SortOptions) (int, error) {
	hdr, recs, err := readAllVariants(in)
	if err != nil {
		return 0, fmt.Errorf("bcftools sort: %w", err)
	}
	order := contigOrder(hdr)
	sortVariantsForSort(recs, order)

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
	}, hdr)
	if err != nil {
		return 0, err
	}
	defer finish()
	if err := w.WriteHeader(); err != nil {
		return 0, err
	}
	for _, v := range recs {
		if err := w.Write(v); err != nil {
			return 0, err
		}
	}
	return len(recs), w.Flush()
}

// sortVariantsForSort performs an in-place stable sort on recs using the
// same ordering rule that `concat -a` uses (lessVariant). A stable sort is
// used so that records with identical keys preserve their original input
// order — this matches upstream's external-merge behaviour for ties.
func sortVariantsForSort(recs []*vcf.Variant, order map[string]int) {
	sort.SliceStable(recs, func(i, j int) bool {
		return lessVariant(recs[i], recs[j], order)
	})
}
