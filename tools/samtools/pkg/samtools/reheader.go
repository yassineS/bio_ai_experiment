package samtools

import (
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// ReheaderOptions configures Reheader. The new header source is either a
// SAM-text file path (HeaderText is empty), pre-loaded text (HeaderPath
// is empty), or — when Command is non-empty — the result of piping the
// existing header through Command (`-c CMD`).
type ReheaderOptions struct {
	// HeaderPath is the path to a SAM text file whose @-prefixed lines
	// replace the input BAM's header.
	HeaderPath string
	// HeaderText is an in-memory copy of the new header text. Takes
	// precedence over HeaderPath.
	HeaderText string
	// Command, when non-empty, is run via `sh -c` with the original
	// header on stdin; its stdout becomes the new header. Mirrors
	// upstream's `-c CMD` flag.
	Command string
	// InPlace requests an in-place rewrite (matches `-i FILE`); the
	// caller arranges the file swap. We simply emit the rewritten BAM
	// to the supplied writer either way.
	InPlace bool
	// NoPG suppresses the @PG record we would otherwise leave alone;
	// since v1 never injects @PG, setting this is a no-op kept for
	// flag-compat.
	NoPG bool
}

// Reheader emits a new BAM stream that has the original record bodies
// of the input but a replaced header. The operation is a streaming
// re-encode: read the entire input header (text + binary @SQ table),
// drop it, write the new header (re-derived @SQ table), then copy each
// record body across.
//
// The new header MUST contain an @SQ table whose order matches the
// original — record refIDs are stored as integer indices into that
// table and re-decoding requires the same ordering. Reheader fails
// loudly if the table size differs.
func Reheader(in io.Reader, out io.Writer, opts ReheaderOptions) error {
	br, err := sam.NewBAMReader(in)
	if err != nil {
		return err
	}
	defer br.Close()
	origHdr := br.Header()

	newHdrText := opts.HeaderText
	if newHdrText == "" && opts.HeaderPath != "" {
		f, err := os.Open(opts.HeaderPath)
		if err != nil {
			return err
		}
		raw, err := io.ReadAll(f)
		_ = f.Close()
		if err != nil {
			return err
		}
		newHdrText = string(raw)
	}
	if opts.Command != "" {
		piped, err := pipeCommand(opts.Command, origHdr.Text())
		if err != nil {
			return err
		}
		newHdrText = piped
	}
	if newHdrText == "" {
		return errors.New("samtools reheader: no header source (use -i FILE or -c CMD)")
	}
	newHdr, err := sam.ParseHeaderText(newHdrText)
	if err != nil {
		return err
	}
	if len(newHdr.Refs) != len(origHdr.Refs) {
		return fmt.Errorf("samtools reheader: new @SQ count %d != original %d (record refIDs would become invalid)",
			len(newHdr.Refs), len(origHdr.Refs))
	}

	bw := sam.NewBAMWriter(out)
	if err := bw.WriteHeader(newHdr); err != nil {
		return err
	}
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		// Refresh RName/RNext via the new header's order — the indices
		// we hold on disk don't change but the names might if the user
		// renamed an @SQ entry.
		if rec.RName != "" && rec.RName != "*" {
			idx := origHdr.RefIndex(rec.RName)
			if idx >= 0 && idx < len(newHdr.Refs) {
				rec.RName = newHdr.Refs[idx].Name
			}
		}
		if rec.RNext != "" && rec.RNext != "*" && rec.RNext != "=" {
			idx := origHdr.RefIndex(rec.RNext)
			if idx >= 0 && idx < len(newHdr.Refs) {
				rec.RNext = newHdr.Refs[idx].Name
			}
		}
		if err := bw.Write(rec); err != nil {
			return err
		}
	}
	return bw.Close()
}

// ReheaderFile is the high-level CLI entry point: opens inPath, runs
// Reheader, and writes to outPath (or, when InPlace is set, atomically
// swaps the rewritten file in over the input).
func ReheaderFile(inPath, outPath string, opts ReheaderOptions) error {
	in, err := os.Open(inPath)
	if err != nil {
		return err
	}
	defer in.Close()

	if opts.InPlace {
		tmp, err := os.CreateTemp(".", ".reheader.tmp.")
		if err != nil {
			return err
		}
		tmpName := tmp.Name()
		if err := Reheader(in, tmp, opts); err != nil {
			_ = tmp.Close()
			_ = os.Remove(tmpName)
			return err
		}
		if err := tmp.Close(); err != nil {
			_ = os.Remove(tmpName)
			return err
		}
		return os.Rename(tmpName, inPath)
	}

	out, err := os.Create(outPath)
	if err != nil {
		return err
	}
	defer out.Close()
	return Reheader(in, out, opts)
}

// pipeCommand runs `sh -c <cmd>` with stdinText on stdin and returns
// stdout. Used to implement `samtools reheader -c CMD`.
func pipeCommand(cmd, stdinText string) (string, error) {
	c := exec.Command("sh", "-c", cmd)
	c.Stdin = stringReader(stdinText)
	out, err := c.Output()
	if err != nil {
		return "", fmt.Errorf("samtools reheader: -c CMD failed: %w", err)
	}
	return string(out), nil
}

// stringReader returns an io.Reader over s without pulling in
// strings.NewReader at the top-level imports list (kept local so the
// reheader file is self-contained).
func stringReader(s string) io.Reader {
	return &strReader{s: s}
}

type strReader struct {
	s string
	i int
}

func (r *strReader) Read(p []byte) (int, error) {
	if r.i >= len(r.s) {
		return 0, io.EOF
	}
	n := copy(p, r.s[r.i:])
	r.i += n
	return n, nil
}
