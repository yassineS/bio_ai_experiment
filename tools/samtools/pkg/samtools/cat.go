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

// Cat concatenates the records from each input BAM stream into out.
// The output's header is taken from the first input (or HeaderOverride
// when set); subsequent inputs must agree on the @SQ ordering.
func Cat(inputs []io.Reader, out io.Writer, opts CatOptions) error {
	if len(inputs) == 0 {
		return errors.New("samtools cat: no input files")
	}

	// Open all readers and parse headers up front so we fail fast on
	// header mismatches before writing any output.
	readers := make([]*sam.BAMReader, 0, len(inputs))
	for i, in := range inputs {
		br, err := sam.NewBAMReader(in)
		if err != nil {
			return fmt.Errorf("samtools cat: input %d: %w", i, err)
		}
		readers = append(readers, br)
	}

	hdr := readers[0].Header()
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
		for i := 1; i < len(readers); i++ {
			if !sameRefTable(hdr.Refs, readers[i].Header().Refs) {
				return fmt.Errorf("samtools cat: input %d has a different @SQ table than input 0", i)
			}
		}
	}

	// Fast path: mirror upstream bam_cat's raw BGZF-block passthrough. It does
	// zero (de)compression on records — after writing the merged header it
	// flushes each input's still-decompressed header-block tail, then copies the
	// remaining compressed blocks verbatim, trimming each input's trailing EOF
	// block, and finally emits exactly one EOF block. This is eligible only when
	// there is no header override (which may rewrite the @SQ table and, upstream,
	// is the user's responsibility) and every input is a real BGZF stream (so the
	// block-level raw copy is available). Anything else takes the decode fallback.
	fast := opts.HeaderOverride == ""
	if fast {
		for _, br := range readers {
			if br.BGZFReader() == nil {
				fast = false
				break
			}
		}
	}
	if fast {
		return catFast(readers, hdr, out)
	}
	return catDecode(readers, hdr, out)
}

// catFast concatenates inputs by copying their compressed BGZF blocks verbatim,
// never decoding or recompressing a record. It is the drop-in analogue of
// htslib's bam_cat: write the header, then for each input flush the leftover
// decompressed bytes of the block that also held the header (a tiny tail that
// is recompressed — harmless, and what htslib does too), copy the remaining
// compressed blocks minus their trailing EOF marker, and emit one final EOF
// block. Every input reader must be BGZF-backed (guaranteed by the caller).
func catFast(readers []*sam.BAMReader, hdr *sam.Header, out io.Writer) error {
	bgw := bgzip.NewWriter(out)
	headerBytes, err := sam.MarshalBAMHeader(hdr)
	if err != nil {
		return err
	}
	if _, err := bgw.Write(headerBytes); err != nil {
		return err
	}
	for i, br := range readers {
		bz := br.BGZFReader()
		// Re-emit the decompressed remainder of the block that held the header
		// (records packed into the same BGZF block as the header text). These
		// bytes are recompressed into the output stream.
		if rem := bz.DecompressedRemainder(); len(rem) > 0 {
			if _, err := bgw.Write(rem); err != nil {
				return err
			}
		}
		// Flush so the header/tail block(s) land in the output before the raw
		// compressed blocks that follow are appended directly to out.
		if err := bgw.Flush(); err != nil {
			return err
		}
		if _, err := bgzip.CopyRawTrimEOF(out, bz.RawRemaining()); err != nil {
			return fmt.Errorf("samtools cat: input %d: %w", i, err)
		}
	}
	// Exactly one EOF block terminates the concatenated stream.
	if _, err := out.Write(bgzip.EOFBlock); err != nil {
		return err
	}
	return nil
}

// catDecode is the record-by-record fallback used when the fast raw-block copy
// cannot apply (a header override, or a non-BGZF input). It decodes every record
// and re-encodes it against hdr, exactly as the original implementation did.
func catDecode(readers []*sam.BAMReader, hdr *sam.Header, out io.Writer) error {
	bw := sam.NewBAMWriter(out)
	if err := bw.WriteHeader(hdr); err != nil {
		return err
	}
	for i, br := range readers {
		for {
			rec, err := br.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return fmt.Errorf("samtools cat: input %d: %w", i, err)
			}
			if err := bw.Write(rec); err != nil {
				return err
			}
		}
	}
	return bw.Close()
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
