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
	"bufio"
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// SortOptions controls the `bcftools sort` subcommand. All fields are
// optional; defaults match the upstream binary.
type SortOptions struct {
	// OutputFormat selects the output encoding. Defaults to OutputVCF.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output (negative means
	// gzip's default).
	CompressLevel int
	// Threads is the -@/--threads value; >1 enables parallel BGZF compression
	// of -O z and -O b output via bgzf.MultiWriter (see ViewOptions.Threads).
	Threads int
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
	var hdr *vcf.Header
	var recs []*vcf.Variant
	var err error

	// RAW-LINE FAST PATH. For uncompressed VCF text output (-O v) the sort never
	// mutates a record, so each VCF data line is captured verbatim, sorted on the
	// light parsed key (CHROM/POS/REF/ALT) and re-emitted unchanged — skipping the
	// per-record INFO/sample-map parse→re-encode round-trip that dominated peak
	// RSS. The captured line re-emits byte-identically to a full parse+re-encode
	// for a well-formed record. The raw path is confined to -O v output; BCF and
	// compressed output keep the full parse path so their binary/BGZF encoders run
	// as before. For BCF *input* (which has no raw text line) the reader
	// transparently falls back to the full parse, so the records still sort and
	// emit correctly under -O v.
	if opts.OutputFormat == OutputVCF {
		hdr, recs, err = readAllVariantsRawLine(in)
	} else {
		hdr, recs, err = readAllVariants(in)
	}
	if err != nil {
		return 0, fmt.Errorf("bcftools sort: %w", err)
	}
	order := contigOrder(hdr)
	sortVariantsForSort(recs, order)

	w, finish, err := openOutput(out, ViewOptions{
		OutputFormat:  opts.OutputFormat,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
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

// readAllVariantsRawLine reads every record for the sort -O v fast path. For a
// VCF text stream it captures each data line verbatim (vcf.KeepRawLine),
// parsing only the CHROM/POS/REF/ALT sort key and re-emitting the line unchanged
// on write. For a BCF stream — which has no raw text line — it transparently
// falls back to the full parse (readAllBCF), returning ordinary parsed records
// whose empty RawLine makes the writer re-encode them normally. Either way the
// records sort by the same key and produce byte-identical -O v output.
func readAllVariantsRawLine(in io.Reader) (*vcf.Header, []*vcf.Variant, error) {
	br := bufio.NewReader(in)
	head, perr := br.Peek(5)
	if perr != nil && perr != io.EOF {
		return nil, nil, perr
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		return readAllBCF(br)
	}
	r := vcf.NewReader(br)
	r.KeepRawLine(true)
	hdr, err := r.ReadHeader()
	if err != nil {
		return nil, nil, err
	}
	var out []*vcf.Variant
	for {
		// Read (owned) so each Variant — including its RawLine — owns its bytes and
		// is safe to buffer across the full input.
		v, rerr := r.Read()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return hdr, out, rerr
		}
		out = append(out, v)
	}
	return hdr, out, nil
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
