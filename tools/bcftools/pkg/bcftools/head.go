// Package bcftools — `bcftools head` subcommand.
//
// `bcftools head` emits the VCF/BCF header (or a slice of it). It is the
// counterpart to `bcftools view -h`, but with finer control: callers can ask
// for the first N header lines, the sample list pulled from `#CHROM`, or
// both. Compared with `view -h` it never opens a writer for records, so it
// is a fast no-allocation way to introspect a file.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// HeadOptions controls the `bcftools head` subcommand.
//
// All fields are optional; when every field is zero/false the full header is
// emitted (matching upstream's default behaviour).
type HeadOptions struct {
	// NumLines, when > 0, caps the number of *header* lines (counting both
	// `##` meta lines and the `#CHROM` line). When NumLines == 0 the entire
	// header is emitted.
	NumLines int
	// SamplesOnly causes the command to print one sample-name per line
	// (parsed from the `#CHROM` line) and skip every other meta line.
	SamplesOnly bool
	// NumRecords, when > 0, appends that many variant records after the
	// full header. This mirrors upstream's `-n/--records INT`.
	NumRecords int
}

// HeadFile is the file-aware entry point for `bcftools head`. It opens path
// through iohelper, parses the header (handling both VCF and BCF), and writes
// the requested slice of it to out.
func HeadFile(path string, out io.Writer, opts HeadOptions) error {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return err
	}
	defer in.Close()
	return Head(in, out, opts)
}

// Head reads a VCF or BCF stream from in, extracts the header, and emits the
// requested slice. It is the streaming entry point used by HeadFile and by
// callers wanting to avoid opening a file twice.
func Head(in io.Reader, out io.Writer, opts HeadOptions) error {
	bw := bufio.NewWriter(out)
	defer bw.Flush()

	// When `-n/--records` is requested we need a full streaming reader so
	// we can decode records after the header. Otherwise we can take the
	// fast header-only path.
	if opts.NumRecords > 0 {
		return headWithRecords(in, bw, opts)
	}

	hdr, err := readAnyHeader(in)
	if err != nil {
		return err
	}

	// htslib injects the implicit ##FILTER=<ID=PASS> line (filter id 0)
	// when the header lacks an explicit PASS definition, so `head` must
	// emit it too for parity with genuine bcftools.
	ensurePASSFilter(hdr)

	if opts.SamplesOnly {
		for _, s := range hdr.Samples {
			if _, err := fmt.Fprintln(bw, s); err != nil {
				return err
			}
		}
		return nil
	}

	emitted := 0
	for _, m := range hdr.MetaInfo {
		if opts.NumLines > 0 && emitted >= opts.NumLines {
			return nil
		}
		if _, err := fmt.Fprintln(bw, m); err != nil {
			return err
		}
		emitted++
	}
	if opts.NumLines > 0 && emitted >= opts.NumLines {
		return nil
	}
	_, err = fmt.Fprintln(bw, buildChromLine(hdr))
	return err
}

// headWithRecords implements `head -n N`: it emits the full header followed
// by the first N variant records. NumLines is ignored in this mode (matches
// upstream — `-h` is a hard cap that applies independently of `-n`, but
// since our header pass always emits the complete header we approximate
// by ignoring NumLines when records are requested).
func headWithRecords(in io.Reader, bw *bufio.Writer, opts HeadOptions) error {
	hdr, recs, err := readAllVariants(in)
	if err != nil {
		return err
	}
	ensurePASSFilter(hdr)
	w := vcf.NewWriter(bw, hdr)
	if err := w.WriteHeader(); err != nil {
		return err
	}
	n := opts.NumRecords
	if n > len(recs) {
		n = len(recs)
	}
	for i := 0; i < n; i++ {
		if err := w.Write(recs[i]); err != nil {
			return err
		}
	}
	return w.Flush()
}

// readAnyHeader reads the header from either a VCF or BCF stream. It peeks at
// the magic bytes and dispatches to the appropriate parser. The returned
// *vcf.Header is the canonical header representation used throughout this
// package.
func readAnyHeader(in io.Reader) (*vcf.Header, error) {
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return nil, err
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		bh, err := bcf.NewReader(br)
		if err != nil {
			return nil, err
		}
		return bh.Header().VCF, nil
	}
	r := vcf.NewReader(br)
	return r.ReadHeader()
}

// buildChromLine reconstructs the `#CHROM\tPOS\tID\t...` line from the
// header, including the FORMAT and sample columns when samples are present.
func buildChromLine(hdr *vcf.Header) string {
	cols := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
	if len(hdr.Samples) > 0 {
		cols += "\tFORMAT\t" + strings.Join(hdr.Samples, "\t")
	}
	return cols
}
