package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/skewer/pkg/skewer"
)

// Upstream skewer's default Illumina adapter prefixes (parameter.cpp:48,50).
// They are applied when the user does not specify -x/-y explicitly, matching
// upstream's behaviour for `skewer <reads.fastq> [paired.fastq]`.
const (
	defaultAdapter3      = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	defaultAdapter3Pair2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
)

// runUpstreamCLI implements skewer 0.2.2's positional command line:
//
//	skewer [options] <reads.fastq> [paired-reads.fastq]
//
// Single-end vs paired-end is selected by the number of positional FASTQ
// inputs (upstream's nFileCnt), not by the -m/--mode flag — exactly as
// upstream does (parameter.cpp:902-911). The -o/--output value is the *base*
// name; outputs are written as <base>-trimmed.fastq for SE and
// <base>-trimmed-pair1.fastq / <base>-trimmed-pair2.fastq for PE, with the
// default base derived from the first input file (extension stripped).
func runUpstreamCLI(args []string) {
	fs := flag.NewFlagSet("skewer", flag.ExitOnError)

	var (
		adapter3   string
		adapter2   string
		mode       string
		errorRate  float64
		indelRate  float64
		minOverlap int
		endQuality int
		meanQual   int
		minLength  int
		maxLength  int
		output     string
		compress   bool
		toStdout   bool
		format     string
		threads    int
		quiet      bool
		showHelp   bool
		showVer    bool
	)

	cliflag.StringVar(fs, &adapter3, "x", "", "", "3' adapter sequence/file")
	cliflag.StringVar(fs, &adapter2, "y", "", "", "3' adapter sequence/file for the paired mate")
	cliflag.StringVar(fs, &mode, "m", "mode", "", "trimming mode: head|tail|any|pe")
	cliflag.Float64Var(fs, &errorRate, "r", "", 0.1, "maximum allowed error rate")
	cliflag.Float64Var(fs, &indelRate, "d", "", 0.03, "maximum allowed indel error rate")
	cliflag.IntVar(fs, &minOverlap, "k", "", 0, "minimum overlap length for adapter detection")
	cliflag.IntVar(fs, &endQuality, "q", "end-quality", 0, "trim 3' end until the given quality is reached")
	cliflag.IntVar(fs, &meanQual, "Q", "mean-quality", 0, "lowest mean quality allowed before trimming")
	cliflag.IntVar(fs, &minLength, "l", "min", 18, "minimum read length allowed after trimming")
	cliflag.IntVar(fs, &maxLength, "L", "max", 0, "maximum read length allowed after trimming (0 = no limit)")
	cliflag.StringVar(fs, &output, "o", "output", "", "base name of output file")
	cliflag.BoolVar(fs, &compress, "z", "compress", false, "compress output in GZIP format")
	cliflag.BoolVar(fs, &toStdout, "1", "stdout", false, "redirect output to STDOUT")
	cliflag.StringVar(fs, &format, "f", "format", "auto", "FASTQ quality format: sanger|solexa|auto")
	cliflag.IntVar(fs, &threads, "t", "threads", 1, "number of concurrent threads (accepted; single-threaded)")
	cliflag.BoolVar(fs, &quiet, "", "quiet", false, "no progress update")
	cliflag.BoolVar(fs, &showHelp, "h", "help", false, "show help")
	cliflag.BoolVar(fs, &showVer, "v", "version", false, "show version")

	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(2)
	}

	if showHelp {
		fmt.Print(usage)
		os.Exit(0)
	}
	if showVer {
		fmt.Println("skewer version 1.0.0 (Go implementation)")
		os.Exit(0)
	}

	inputs := fs.Args()
	if len(inputs) == 0 {
		fmt.Fprintln(os.Stderr, "Error: no input FASTQ file specified")
		fs.Usage()
		os.Exit(1)
	}
	if len(inputs) > 2 {
		fmt.Fprintf(os.Stderr, "Error: too many positional inputs (%d); expected 1 (single-end) or 2 (paired-end)\n", len(inputs))
		os.Exit(1)
	}
	paired := len(inputs) == 2

	// Resolve adapters. Upstream falls back to the default Illumina prefixes
	// when -x is not given, and uses the pair-2 default for -y when only -x is
	// specified (parameter.cpp). For 5'/head mode the supplied adapter is a 5'
	// adapter rather than a 3' one.
	x := adapter3
	if x == "" {
		x = defaultAdapter3
	}
	y := adapter2

	// Map upstream's mode to our TrimOptions.
	//   head -> 5' adapter; tail/any -> 3' adapter; pe -> matrix paired-end.
	var adapter5, threePrime string
	peMatrix := paired
	switch strings.ToLower(mode) {
	case "head":
		adapter5 = x
	case "any":
		// "any" trims the adapter wherever it occurs; we apply it as a 3'
		// adapter, which our aligner already scans across the read.
		threePrime = x
	case "pe", "":
		// Default: tail for SE, pe (matrix) for PE.
		threePrime = x
		if paired {
			peMatrix = true
			if y == "" {
				y = defaultAdapter3Pair2
			}
		}
	case "tail":
		threePrime = x
	case "mp", "ap":
		// Mate-pair / amplicon sub-modes are not yet modelled; fall back to the
		// default paired/3' behaviour so the CLI still runs end-to-end.
		threePrime = x
		if paired && y == "" {
			y = defaultAdapter3Pair2
		}
	default:
		fmt.Fprintf(os.Stderr, "Error: unknown -m/--mode value %q (expected head|tail|any|pe)\n", mode)
		os.Exit(1)
	}

	encoding := getUpstreamEncoding(format)

	// minOverlap default: upstream uses max(1, int(4-10*r)) for single-end when
	// -k is not given. Our trimmer's historical default is 3; honour an explicit
	// -k but otherwise keep the established default for output stability.
	overlap := 3
	if minOverlap > 0 {
		overlap = minOverlap
	}

	opts := skewer.TrimOptions{
		Adapter3:         threePrime,
		Adapter5:         adapter5,
		MinLength:        minLength,
		QualThreshold:    endQuality,
		MinOverlap:       overlap,
		ErrorRate:        errorRate,
		IndelRate:        indelRate,
		ProgressReport:   false,
		ProgressInterval: 100000,
		PEMatrixMode:     peMatrix,
	}

	if paired {
		runUpstreamPaired(inputs, output, toStdout, compress, encoding, opts, quiet)
	} else {
		runUpstreamSingle(inputs[0], output, toStdout, compress, encoding, opts, quiet)
	}
}

// runUpstreamSingle trims one input and writes <base>-trimmed.fastq (or stdout
// with -1), mirroring upstream's single-end output naming.
func runUpstreamSingle(input, outBase string, toStdout, compress bool, encoding fastq.QualityEncoding, opts skewer.TrimOptions, quiet bool) {
	in, err := iohelper.OpenReader(input)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input file: %v\n", err)
		os.Exit(1)
	}
	defer in.Close()

	var output io.WriteCloser
	if toStdout {
		output, err = iohelper.OpenWriter("-")
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening stdout: %v\n", err)
			os.Exit(1)
		}
	} else {
		// The ".gz" suffix added by trimmedName tells iohelper.OpenWriter to
		// gzip-compress the stream itself, so we must not wrap it again.
		name := trimmedName(outBase, input, "", compress)
		output, err = iohelper.OpenWriter(name)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating output file: %v\n", err)
			os.Exit(1)
		}
	}
	defer output.Close()

	stats, err := skewer.TrimSingleEnd(in, output, encoding, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error during trimming: %v\n", err)
		os.Exit(1)
	}
	if !quiet {
		printStats(stats, "SE")
	}
}

// runUpstreamPaired trims two inputs and writes <base>-trimmed-pair1.fastq and
// <base>-trimmed-pair2.fastq (or interleaved-to-stdout with -1), mirroring
// upstream's paired-end output naming.
func runUpstreamPaired(inputs []string, outBase string, toStdout, compress bool, encoding fastq.QualityEncoding, opts skewer.TrimOptions, quiet bool) {
	in1, err := iohelper.OpenReader(inputs[0])
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening first input file: %v\n", err)
		os.Exit(1)
	}
	defer in1.Close()

	in2, err := iohelper.OpenReader(inputs[1])
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening second input file: %v\n", err)
		os.Exit(1)
	}
	defer in2.Close()

	var out1, out2 io.WriteCloser
	if toStdout {
		// Upstream with -1 interleaves both mates onto STDOUT.
		w, oerr := iohelper.OpenWriter("-")
		if oerr != nil {
			fmt.Fprintf(os.Stderr, "Error opening stdout: %v\n", oerr)
			os.Exit(1)
		}
		out1 = nopCloser{w}
		out2 = nopCloser{w}
		defer w.Close()
	} else {
		// The ".gz" suffix added by trimmedName tells iohelper.OpenWriter to
		// gzip-compress each stream itself, so we must not wrap them again.
		name1 := trimmedName(outBase, inputs[0], "-pair1", compress)
		name2 := trimmedName(outBase, inputs[0], "-pair2", compress)
		out1, err = iohelper.OpenWriter(name1)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating first output file: %v\n", err)
			os.Exit(1)
		}
		defer out1.Close()

		out2, err = iohelper.OpenWriter(name2)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error creating second output file: %v\n", err)
			os.Exit(1)
		}
		defer out2.Close()
	}

	stats, err := skewer.TrimPairedEnd(in1, in2, out1, out2, nil, encoding, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error during trimming: %v\n", err)
		os.Exit(1)
	}
	if !quiet {
		printStats(stats, "PE")
	}
}

// trimmedName builds upstream skewer's output filename. The base is either the
// explicit -o value or, when empty, the first input's path with its last
// extension stripped (parameter.cpp:1166-1176). The "-trimmed" decoration plus
// the pairSuffix (""/"-pair1"/"-pair2") and a ".fastq" extension are appended,
// with ".gz" added when compression is enabled.
func trimmedName(outBase, input, pairSuffix string, compress bool) string {
	base := outBase
	if base == "" {
		base = stripLastExt(input)
	}
	name := base + "-trimmed" + pairSuffix + ".fastq"
	if compress {
		name += ".gz"
	}
	return name
}

// stripLastExt removes the final ".<ext>" from a path's basename component,
// matching upstream's occOfLastDot handling (the dot must follow the last path
// separator). A leading-dot name with no other dot is returned unchanged.
func stripLastExt(path string) string {
	dir := filepath.Dir(path)
	bn := filepath.Base(path)
	if i := strings.LastIndex(bn, "."); i > 0 {
		bn = bn[:i]
	}
	if dir == "." && !strings.HasPrefix(path, "./") {
		return bn
	}
	return filepath.Join(dir, bn)
}

// getUpstreamEncoding maps upstream's -f/--format value to a quality encoding.
// "auto" defaults to Sanger/Phred+33, which covers modern Illumina 1.8+ data.
func getUpstreamEncoding(format string) fastq.QualityEncoding {
	switch strings.ToLower(format) {
	case "solexa", "illumina", "phred64":
		return fastq.Phred64
	case "sanger", "phred33", "auto", "":
		return fastq.Phred33
	default:
		fmt.Fprintf(os.Stderr, "Warning: unknown -f/--format %q, using auto (sanger)\n", format)
		return fastq.Phred33
	}
}

// nopCloser adapts an io.Writer to io.WriteCloser without closing the wrapped
// writer; used so both interleaved -1 mates can share one STDOUT writer.
type nopCloser struct{ io.Writer }

func (nopCloser) Close() error { return nil }
