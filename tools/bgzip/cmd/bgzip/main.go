// Command bgzip is a pure-Go reimplementation of htslib's bgzip CLI.
//
// It compresses a file into a Blocked GNU Zip Format (BGZF) stream — the
// foundational on-disk codec used by .vcf.gz, BAM, BCF, and tabix indices —
// and provides decoding, in-place replacement, and the -b/-s/-r index
// utilities that downstream tools rely on.
package main

import (
	"compress/flate"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

const version = "1.0.0"

const usage = `bgzip - Block-gzip a file (BGZF) or decompress one.

Usage:
  bgzip [options] [file]
  bgzip -d [options] file.gz

Options:
  -c, --stdout              Write output to stdout, keep input intact.
  -d, --decompress          Decompress instead of compress.
  -f, --force               Overwrite existing output file.
  -k, --keep                Keep input file (do not delete it after success).
  -l, --compress-level N    Compression level 0-9 (default 6).
  -o, --output FILE         Write output to the named FILE. Use '-' for stdout.
                            This is the way to name the output when the input is
                            stdin (otherwise stdin output goes to stdout).
  -@, -t, --threads N       Number of compression threads (default 1).
                            Blocks are compressed in parallel; output and the
                            .gzi index stay correct regardless of thread count.
  -b, --offset N            Print uncompressed offset at compressed offset N.
  -s, --size                Print the decompressed size of the file.
  -r, --reindex             Write a .gzi index alongside file.gz.
      --test                Test the integrity of a BGZF file: decompress it
                            fully, write nothing, and exit non-zero on any
                            error. (Upstream spells this -t/--test; this port
                            binds -t to --threads, so only the long --test
                            form is offered — see the README "Deviations".)
  -h, --help                Show this help and exit.
  -v, --version             Show version and exit.

Use '-' as a filename to read from stdin or write to stdout.
`

func main() {
	os.Exit(run(os.Args[1:], os.Stdin, os.Stdout, os.Stderr))
}

func run(args []string, stdin io.Reader, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("bgzip", flag.ContinueOnError)
	fs.SetOutput(io.Discard) // we print parse errors and usage ourselves

	var (
		stdoutFlag  bool
		decompress  bool
		force       bool
		keep        bool
		level       int
		threads     int
		writeFname  string
		offset      int64
		offsetSet   bool
		showSize    bool
		reindex     bool
		test        bool
		showHelp    bool
		showVersion bool
	)

	cliflag.BoolVar(fs, &stdoutFlag, "c", "stdout", false, "Write to stdout")
	cliflag.BoolVar(fs, &decompress, "d", "decompress", false, "Decompress")
	cliflag.BoolVar(fs, &force, "f", "force", false, "Force overwrite")
	cliflag.BoolVar(fs, &keep, "k", "keep", false, "Keep input file")
	cliflag.IntVar(fs, &level, "l", "compress-level", bgzip.DefaultCompression, "Compression level (0-9)")
	cliflag.StringVar(fs, &writeFname, "o", "output", "", "Write output to the named file")
	cliflag.IntVar(fs, &threads, "t", "threads", 1, "Number of compression threads")
	// Upstream bgzip's canonical short flag for threads is -@; accept it as an
	// alias so existing command lines keep working.
	fs.IntVar(&threads, "@", 1, "")
	// -b takes an int64; cliflag does not expose Int64Var, so register both
	// names directly.
	fs.Func("b", "", func(s string) error { return parseInt64Flag(s, &offset, &offsetSet) })
	fs.Func("offset", "Print uncompressed offset at compressed offset N", func(s string) error {
		return parseInt64Flag(s, &offset, &offsetSet)
	})
	cliflag.BoolVar(fs, &showSize, "s", "size", false, "Print decompressed size")
	cliflag.BoolVar(fs, &reindex, "r", "reindex", false, "Write .gzi index")
	// Upstream's integrity check is -t/--test, but this port already binds -t
	// to --threads, so --test is registered as a long-only flag (documented
	// deviation). It decompresses the whole stream and discards the output.
	fs.BoolVar(&test, "test", false, "Test the integrity of a BGZF file")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "Show help")
	cliflag.BoolVar(fs, &showVersion, "v", "version", false, "Show version")
	// Upstream bgzip (bgzip.c getopt "cdh?fb:@:s:iI:l:grtko:") also accepts a
	// few flags this port does not act on yet. Register them as accepted
	// stubs so legacy command lines — and any bundled cluster that includes
	// them — still parse, matching upstream's surface:
	//   -i        write a .gzi index while compressing (we don't auto-index;
	//             use -r/--reindex). Accepted no-op.
	//   -I FILE   name of the index file. Accepted no-op.
	//   -g        rebgzip a raw stream using an existing index. Accepted no-op.
	// (-t is intentionally our --threads alias rather than upstream's "test"
	// switch; see the README "Deviations" section.)
	var (
		indexOnCompress bool
		indexFname      string
		rebgzip         bool
	)
	cliflag.BoolVar(fs, &indexOnCompress, "i", "index", false, "Ignored: write .gzi while compressing (use -r) (legacy)")
	cliflag.StringVar(fs, &indexFname, "I", "index-name", "", "Ignored: index file name (legacy)")
	cliflag.BoolVar(fs, &rebgzip, "g", "rebgzip", false, "Ignored: rebgzip a raw stream (legacy)")

	fs.Usage = func() { fmt.Fprint(stderr, usage) }

	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// (-cd == -c -d) and value concatenation (-l6 == -l 6) work the way
	// upstream bgzip's getopt parser accepts them.
	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(stderr, err)
		fmt.Fprint(stderr, usage)
		return 2
	}
	_ = indexOnCompress
	_ = indexFname
	_ = rebgzip
	if showHelp {
		fmt.Fprint(stdout, usage)
		return 0
	}
	if showVersion {
		fmt.Fprintf(stdout, "bgzip %s (Go implementation)\n", version)
		return 0
	}

	if threads < 1 {
		fmt.Fprintf(stderr, "bgzip: --threads must be >= 1\n")
		return 2
	}

	// Upstream bgzip treats `-o -` and `--output -` as a request for stdout,
	// equivalent to -c. Normalise that here so the output-naming logic below
	// only deals with real file paths.
	if writeFname == "-" {
		writeFname = ""
		stdoutFlag = true
	}
	if writeFname != "" && stdoutFlag {
		fmt.Fprintf(stderr, "bgzip: cannot write to %s and stdout at the same time\n", writeFname)
		return 2
	}

	rest := fs.Args()
	var input string
	if len(rest) == 0 {
		input = "-"
	} else if len(rest) == 1 {
		input = rest[0]
	} else {
		fmt.Fprintf(stderr, "bgzip: at most one input file may be given\n")
		return 2
	}

	// Mode dispatch — query flags first since they don't modify the input.
	switch {
	case test:
		return runTest(input, stdin, stderr)
	case offsetSet:
		return runOffset(input, offset, stdout, stderr)
	case showSize:
		return runSize(input, stdout, stderr)
	case reindex:
		return runReindex(input, stdout, stderr)
	case decompress:
		return runDecompress(input, writeFname, stdoutFlag, force, keep, stdin, stdout, stderr)
	default:
		return runCompress(input, writeFname, stdoutFlag, force, keep, level, threads, stdin, stdout, stderr)
	}
}

func parseInt64Flag(s string, dest *int64, set *bool) error {
	var v int64
	_, err := fmt.Sscanf(s, "%d", &v)
	if err != nil {
		return fmt.Errorf("invalid integer %q", s)
	}
	*dest = v
	*set = true
	return nil
}

func runCompress(input, writeFname string, useStdout, force, keep bool, level, threads int, stdin io.Reader, stdout, stderr io.Writer) int {
	if level < flate.HuffmanOnly || level > flate.BestCompression {
		fmt.Fprintf(stderr, "bgzip: invalid compression level %d\n", level)
		return 2
	}

	in, closeIn, err := openInputBinary(input, stdin)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer closeIn()

	// Decide where output goes. An explicit -o/--output FILE always wins; it is
	// how the output is named when the input is stdin. Otherwise stdin (or -c)
	// writes to stdout and a named input file produces "<input>.gz".
	var (
		outName   string
		out       io.Writer
		outCloser io.Closer
	)
	switch {
	case writeFname != "":
		outName = writeFname
		if !force {
			if _, err := os.Stat(outName); err == nil {
				fmt.Fprintf(stderr, "bgzip: %s already exists; use -f to overwrite\n", outName)
				return 1
			}
		}
		f, err := os.Create(outName)
		if err != nil {
			fmt.Fprintf(stderr, "bgzip: %v\n", err)
			return 1
		}
		out = f
		outCloser = f
	case useStdout || input == "-":
		out = stdout
		outCloser = io.NopCloser(nil)
	default:
		outName = input + ".gz"
		if !force {
			if _, err := os.Stat(outName); err == nil {
				fmt.Fprintf(stderr, "bgzip: %s already exists; use -f to overwrite\n", outName)
				return 1
			}
		}
		f, err := os.Create(outName)
		if err != nil {
			fmt.Fprintf(stderr, "bgzip: %v\n", err)
			return 1
		}
		out = f
		outCloser = f
	}

	bw, err := newCompressor(out, level, threads)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	// Ensure the compressor is always closed so the MultiWriter's worker
	// goroutines drain even when the copy below fails partway through. Close is
	// idempotent, so the explicit Close on the success path is harmless.
	bwClosed := false
	closeCompressor := func() error {
		if bwClosed {
			return nil
		}
		bwClosed = true
		return bw.Close()
	}
	defer closeCompressor()

	if _, err := io.Copy(bw, in); err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	if err := closeCompressor(); err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	if outCloser != nil {
		if err := outCloser.Close(); err != nil {
			fmt.Fprintf(stderr, "bgzip: %v\n", err)
			return 1
		}
	}

	// Only delete the input when we compressed a real named file in place and
	// the user did not ask to keep it or redirect output elsewhere.
	if input != "-" && !keep && !useStdout && writeFname == "" {
		if err := os.Remove(input); err != nil {
			fmt.Fprintf(stderr, "bgzip: %v\n", err)
			return 1
		}
	}
	return 0
}

// newCompressor returns a BGZF write-closer using the single-threaded Writer
// when threads <= 1 and the concurrent MultiWriter otherwise. Both produce a
// valid BGZF stream terminated by the EOF marker on Close.
func newCompressor(out io.Writer, level, threads int) (io.WriteCloser, error) {
	if threads <= 1 {
		return bgzip.NewWriterLevel(out, level)
	}
	return bgzip.NewMultiWriter(out, level, threads)
}

func runDecompress(input, writeFname string, useStdout, force, keep bool, stdin io.Reader, stdout, stderr io.Writer) int {
	in, closeIn, err := openInputBinary(input, stdin)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer closeIn()

	br, err := bgzip.NewReader(in)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer br.Close()

	var (
		outName   string
		out       io.Writer
		outCloser io.Closer
	)
	switch {
	case writeFname != "":
		outName = writeFname
		if !force {
			if _, err := os.Stat(outName); err == nil {
				fmt.Fprintf(stderr, "bgzip: %s already exists; use -f to overwrite\n", outName)
				return 1
			}
		}
		f, err := os.Create(outName)
		if err != nil {
			fmt.Fprintf(stderr, "bgzip: %v\n", err)
			return 1
		}
		out = f
		outCloser = f
	case useStdout || input == "-":
		out = stdout
		outCloser = io.NopCloser(nil)
	default:
		if !strings.HasSuffix(input, ".gz") {
			fmt.Fprintf(stderr, "bgzip: input %s does not end in .gz\n", input)
			return 1
		}
		outName = strings.TrimSuffix(input, ".gz")
		if !force {
			if _, err := os.Stat(outName); err == nil {
				fmt.Fprintf(stderr, "bgzip: %s already exists; use -f to overwrite\n", outName)
				return 1
			}
		}
		f, err := os.Create(outName)
		if err != nil {
			fmt.Fprintf(stderr, "bgzip: %v\n", err)
			return 1
		}
		out = f
		outCloser = f
	}
	if _, err := io.Copy(out, br); err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	if outCloser != nil {
		if err := outCloser.Close(); err != nil {
			fmt.Fprintf(stderr, "bgzip: %v\n", err)
			return 1
		}
	}
	if input != "-" && !keep && !useStdout && writeFname == "" {
		if err := os.Remove(input); err != nil {
			fmt.Fprintf(stderr, "bgzip: %v\n", err)
			return 1
		}
	}
	return 0
}

// runTest verifies the integrity of a BGZF stream the way upstream bgzip's
// -t/--test does: it decompresses the entire input, discarding the output,
// and reports any error. The input is never modified, nothing is written to
// stdout, and a decode failure (truncated stream, bad CRC, corrupt block)
// yields a non-zero exit code.
func runTest(input string, stdin io.Reader, stderr io.Writer) int {
	in, closeIn, err := openInputBinary(input, stdin)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer closeIn()

	br, err := bgzip.NewReader(in)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer br.Close()

	if _, err := io.Copy(io.Discard, br); err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	return 0
}

func runOffset(input string, off int64, stdout, stderr io.Writer) int {
	f, closeIn, err := openInputBinary(input, nil)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer closeIn()
	u, err := bgzip.UncompressedOffsetAt(f, off)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	fmt.Fprintf(stdout, "%d\n", u)
	return 0
}

func runSize(input string, stdout, stderr io.Writer) int {
	f, closeIn, err := openInputBinary(input, nil)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer closeIn()
	n, err := bgzip.DecompressedSize(f)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	fmt.Fprintf(stdout, "%d\n", n)
	return 0
}

func runReindex(input string, stdout, stderr io.Writer) int {
	if input == "-" {
		fmt.Fprintln(stderr, "bgzip: --reindex requires a real file path")
		return 2
	}
	f, err := os.Open(input)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer f.Close()
	offsets, err := bgzip.Scan(f)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	gziName := input + ".gzi"
	out, err := os.Create(gziName)
	if err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := bgzip.WriteGZI(out, offsets); err != nil {
		fmt.Fprintf(stderr, "bgzip: %v\n", err)
		return 1
	}
	return 0
}

// openInputBinary opens path for binary reading, honouring "-" as stdin. The
// returned close function is safe to call even when reading from stdin.
func openInputBinary(path string, stdin io.Reader) (io.Reader, func() error, error) {
	if path == "-" || path == "" {
		if stdin == nil {
			return nil, nil, fmt.Errorf("stdin not available")
		}
		return stdin, func() error { return nil }, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	return f, f.Close, nil
}
