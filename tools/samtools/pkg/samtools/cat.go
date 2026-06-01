package samtools

import (
	"errors"
	"fmt"
	"io"
	"os"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// CatOptions configures the cat subcommand. Cat differs from merge in
// that it does NOT re-sort: input order is preserved record-by-record.
// All inputs must share the same @SQ table (we accept the union but
// reject contradictions).
type CatOptions struct {
	// HeaderOverride, when non-empty, replaces the header parsed from
	// the first input. The override file must be a SAM header text file
	// (matches upstream's `-h FILE`).
	HeaderOverride string
	// Threads is accepted for upstream-CLI compatibility; ignored.
	Threads int
}

// Cat concatenates the alignment blocks from each input BAM stream into out.
// It mirrors htslib's bam_cat byte-for-byte: a single output header is written
// (from the first input, or HeaderOverride when set), then each input's
// alignment BGZF blocks are copied verbatim (bgzf_raw_read / bgzf_raw_write)
// with that input's trailing EOF block stripped. A single canonical EOF block
// is appended at the end. Re-encoding is avoided so the alignment blocks stay
// byte-identical with the inputs.
func Cat(inputs []io.Reader, out io.Writer, opts CatOptions) error {
	if len(inputs) == 0 {
		return errors.New("samtools cat: no input files")
	}

	// Inflate each input's leading blocks far enough to recover its header
	// and the byte boundary of the first record, retaining the raw blocks so
	// the alignment data can be copied verbatim afterwards.
	parsed := make([]catInput, 0, len(inputs))
	for i, in := range inputs {
		ci, err := readCatInput(in)
		if err != nil {
			return fmt.Errorf("samtools cat: input %d: %w", i, err)
		}
		parsed = append(parsed, ci)
	}

	hdr := parsed[0].hdr
	if opts.HeaderOverride != "" {
		f, err := os.Open(opts.HeaderOverride)
		if err != nil {
			return err
		}
		defer f.Close()
		raw, err := io.ReadAll(f)
		if err != nil {
			return err
		}
		newHdr, err := sam.ParseHeaderText(string(raw))
		if err != nil {
			return err
		}
		hdr = newHdr
	}

	// Cross-check that every input shares the @SQ table with the chosen
	// header. Cat without `-h` requires this for the output to be
	// re-decodable; with `-h FILE` we let the user take responsibility.
	if opts.HeaderOverride == "" {
		for i := 1; i < len(parsed); i++ {
			if !sameRefTable(hdr.Refs, parsed[i].hdr.Refs) {
				return fmt.Errorf("samtools cat: input %d has a different @SQ table than input 0", i)
			}
		}
	}

	bw := sam.NewBAMWriter(out)
	if err := bw.WriteHeader(hdr); err != nil {
		return err
	}

	for _, ci := range parsed {
		if len(ci.leftover) > 0 {
			if err := bw.WriteUncompressed(ci.leftover); err != nil {
				return err
			}
			if err := bw.Flush(); err != nil {
				return err
			}
		}
		for _, rb := range ci.blocks {
			if err := bw.WriteRawBGZF(rb.Compressed); err != nil {
				return err
			}
		}
	}
	return bw.Close()
}

// catInput holds one input's parsed header plus the raw alignment blocks to be
// copied verbatim into the concatenated output.
type catInput struct {
	hdr      *sam.Header
	leftover []byte            // record bytes that shared the header's block
	blocks   []*bgzip.RawBlock // raw alignment blocks (EOF stripped)
}

// readCatInput inflates an input's leading BGZF blocks until the BAM header is
// fully available, then collects the remaining raw alignment blocks (excluding
// the trailing EOF marker) for verbatim copying.
func readCatInput(in io.Reader) (catInput, error) {
	var res catInput
	var decoded []byte
	var hdrLen int
	for {
		n, err := sam.BAMHeaderEncodedLen(decoded)
		if err == nil {
			hdrLen = n
			break
		}
		rb, rerr := bgzip.ReadRawBlock(in)
		if rerr != nil {
			if rerr == io.EOF {
				return res, io.ErrUnexpectedEOF
			}
			return res, rerr
		}
		if rb.IsEOF {
			return res, io.ErrUnexpectedEOF
		}
		decoded = append(decoded, rb.Uncompressed...)
	}
	hdr, err := sam.ParseEncodedBAMHeader(decoded[:hdrLen])
	if err != nil {
		return res, err
	}
	res.hdr = hdr
	res.leftover = decoded[hdrLen:]
	for {
		rb, rerr := bgzip.ReadRawBlock(in)
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return res, rerr
		}
		if rb.IsEOF {
			// Strip the input's EOF marker; a single canonical EOF is
			// appended after all inputs (matches htslib bam_cat).
			continue
		}
		res.blocks = append(res.blocks, rb)
	}
	return res, nil
}

// sameRefTable reports whether two @SQ tables are byte-for-byte
// equivalent (same length, same name+length pairs in the same order).
func sameRefTable(a, b []sam.Reference) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i].Name != b[i].Name || a[i].Length != b[i].Length {
			return false
		}
	}
	return true
}

// CatFiles is the high-level CLI entry point: opens each path, runs Cat,
// and closes the inputs. The output writer is owned by the caller.
func CatFiles(paths []string, out io.Writer, opts CatOptions) error {
	readers := make([]io.Reader, 0, len(paths))
	closers := make([]io.Closer, 0, len(paths))
	for _, p := range paths {
		f, err := os.Open(p)
		if err != nil {
			for _, c := range closers {
				_ = c.Close()
			}
			return err
		}
		readers = append(readers, f)
		closers = append(closers, f)
	}
	defer func() {
		for _, c := range closers {
			_ = c.Close()
		}
	}()
	return Cat(readers, out, opts)
}
