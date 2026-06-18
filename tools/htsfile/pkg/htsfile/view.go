package htsfile

// View implements `htsfile -c` (--view): the textual, format-aware view of a
// file. Unlike a raw decompress, upstream routes each file through htslib's
// format reader and re-serialises it, so a VCF/BCF is emitted as canonical VCF
// text. This file implements the VCF/BCF path (the common case); other formats
// return ErrViewUnsupported so the CLI can report them rather than emit
// non-matching bytes.

import (
	"errors"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// ErrViewUnsupported is returned by View for formats whose canonical
// re-serialisation is not yet implemented.
var ErrViewUnsupported = errors.New("htsfile: --view not implemented for this format")

// View writes the canonical textual form of the file at path (already
// identified as f) to out, mirroring `htsfile -c`. path may be "-" for stdin.
func View(path string, f *Format, out io.Writer) error {
	switch f.Payload {
	case PayloadVCF, PayloadBCF:
		return viewVCF(path, out)
	default:
		return fmt.Errorf("%w: %s", ErrViewUnsupported, f.Payload)
	}
}

// viewVCF reads a (optionally gzip/BGZF-compressed) VCF and re-serialises it as
// canonical VCF text — the htslib round-trip htsfile -c performs. The implicit
// ##FILTER=<ID=PASS> line htslib injects is provenance-class boilerplate (the
// pipeline's comparison strips it), so it is not re-emitted here.
func viewVCF(path string, out io.Writer) error {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return fmt.Errorf("htsfile: open %s: %w", path, err)
	}
	defer in.Close()

	r := vcf.NewReader(in)
	hdr, err := r.ReadHeader()
	if err != nil {
		return fmt.Errorf("htsfile: read header: %w", err)
	}
	w := vcf.NewWriter(out, hdr)
	if err := w.WriteHeader(); err != nil {
		return err
	}
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("htsfile: read record: %w", err)
		}
		if err := w.Write(v); err != nil {
			return err
		}
	}
	return w.Flush()
}
