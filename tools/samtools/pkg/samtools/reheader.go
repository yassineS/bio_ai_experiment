package samtools

import (
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
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
	// NoPG suppresses the @PG line injection. By default, reheader
	// appends a @PG line documenting the command via the shared
	// InjectPG helper.
	NoPG bool
	// PGCommand is the raw command-line stored under @PG:CL when NoPG
	// is false. The CLI populates this with os.Args. (Distinct from
	// `Command` above, which is the `-c CMD` shell pipeline applied
	// to the existing header text.)
	PGCommand string
}

// Reheader emits a new BAM stream that has the original alignment blocks of the
// input but a replaced header. It mirrors htslib's bam_reheader byte-for-byte:
// the new header is written into its own BGZF block(s), then the input's
// alignment BGZF blocks are copied verbatim (bgzf_raw_read / bgzf_raw_write) so
// they remain byte-identical with the input — only the header block is
// rewritten.
//
// The new header MUST contain an @SQ table whose order matches the original —
// record refIDs are integer indices into that table and the records are copied
// without rewriting, so the ordering must be preserved. Reheader fails loudly
// if the table size differs.
func Reheader(in io.Reader, out io.Writer, opts ReheaderOptions) error {
	// Read leading BGZF blocks, inflating until the full BAM header is
	// available, so we can locate the header/record boundary exactly as
	// htslib's bam_hdr_read does (which may leave leftover record bytes in
	// the header's final block).
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
				return io.ErrUnexpectedEOF
			}
			return rerr
		}
		if rb.IsEOF {
			return io.ErrUnexpectedEOF
		}
		decoded = append(decoded, rb.Uncompressed...)
	}

	origHdr, err := sam.ParseEncodedBAMHeader(decoded[:hdrLen])
	if err != nil {
		return err
	}

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

	newHdr = InjectPG(newHdr, "samtools", "samtools", "0.1.0", opts.PGCommand, opts.NoPG)
	bw := sam.NewBAMWriter(out)
	if err := bw.WriteHeader(newHdr); err != nil {
		return err
	}
	// Leftover record bytes that shared the header's final BGZF block must be
	// re-compressed; htslib does the same (bgzf_write + bgzf_flush) before the
	// raw-block copy begins.
	if leftover := decoded[hdrLen:]; len(leftover) > 0 {
		if err := bw.WriteUncompressed(leftover); err != nil {
			return err
		}
		if err := bw.Flush(); err != nil {
			return err
		}
	}
	// Copy the remaining input BGZF blocks verbatim. The input's own EOF block
	// is copied too (htslib copies all raw bytes through end of file); the
	// writer's Close then appends the canonical EOF, reproducing the upstream
	// double-EOF layout.
	for {
		rb, rerr := bgzip.ReadRawBlock(in)
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return rerr
		}
		if err := bw.WriteRawBGZF(rb.Compressed); err != nil {
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
