package samtools

import (
	"errors"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
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
