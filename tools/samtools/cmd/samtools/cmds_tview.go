package main

// cmds_tview.go wires the `samtools tview` subcommand. Only the
// non-interactive display modes are supported: -d T (text) and -d H (html),
// which render the alignment-viewer frame to stdout. The interactive ncurses
// mode (-d C, and the bare default which upstream treats as curses) is a
// deliberate non-goal — it needs a TTY UI library this project does not
// depend on — so it is rejected with a clear pointer to -d T / -d H.

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/samtools/pkg/samtools"
)

const tviewUsage = `samtools tview - text/HTML alignment viewer.

Usage:
  samtools tview [options] <in.bam|in.cram> [ref.fasta]

Options:
  -d, --display T|H|C   Output as (T)ext or (H)tml. (C)urses is not supported
                        in this port; use -d T or -d H in a pipeline.
  -p, --position REG    Start at this region (chr or chr:pos).
  -w, --width INT       Display width in columns (default 80).
  -s, --sample STR      Show only reads from this sample or read group.
  -T, --reference FA    Reference FASTA (also given positionally as ref.fasta).
  -i, --hide-inserts    Accepted for compatibility (insertion columns are not
                        expanded in this port).
  -h, --help            Show this help.
  -v, --version         Show version.

Notes:
  - -d T emits the plain character grid (ruler, reference, consensus, then one
    row per read, non-overlapping reads sharing a row). Matches are '.'/','
    (forward/reverse), mismatches the read base (UPPER/lower by strand),
    deletions '*', reference skips '>'/'<'.
  - -d H emits the same grid as a coloured HTML document.
  - The interactive ncurses viewer (-d C, and the bare default) is not built:
    this port targets pipelines, so it has no TTY UI dependency.
`

func runTview(args []string) int {
	fs := flag.NewFlagSet("samtools tview", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		display     string
		position    string
		width       int
		sample      string
		refFasta    string
		hideInserts bool
		showHelp    bool
		showVer     bool
	)
	cliflag.StringVar(fs, &display, "d", "display", "", "Display mode: T (text) or H (html)")
	cliflag.StringVar(fs, &position, "p", "position", "", "Start region (chr or chr:pos)")
	cliflag.IntVar(fs, &width, "w", "width", 0, "Display width in columns")
	cliflag.StringVar(fs, &sample, "s", "sample", "", "Show only reads from this sample/read group")
	cliflag.StringVar(fs, &refFasta, "T", "reference", "", "Reference FASTA")
	cliflag.BoolVar(fs, &hideInserts, "i", "hide-inserts", false, "Hide insertion columns (accepted)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream tview (bam_tview.c getopt "s:p:d:Xw:i") also accepts -X (an
	// explicit index-file argument). This port reads the alignment file
	// directly without a sibling index, so -X is an accepted no-op kept for
	// compatibility so legacy/bundled command lines parse.
	var customIdx bool
	fs.BoolVar(&customIdx, "X", false, "")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, tviewUsage)
		return 2
	}
	_ = customIdx
	if showHelp {
		fmt.Print(tviewUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprintln(os.Stderr, "samtools tview: missing input file")
		fmt.Fprint(os.Stderr, tviewUsage)
		return 2
	}

	// Resolve the display mode. Upstream maps the first char of -d:
	// H/h -> html, T/t -> text, C/c (and the default) -> curses. We support
	// only the two non-interactive modes; curses (and the bare default) is a
	// hard error pointing at -d T / -d H.
	var mode samtools.TviewMode
	switch {
	case len(display) > 0 && (display[0] == 'H' || display[0] == 'h'):
		mode = samtools.TviewHTML
	case len(display) > 0 && (display[0] == 'T' || display[0] == 't'):
		mode = samtools.TviewText
	default:
		fmt.Fprintln(os.Stderr,
			"samtools tview: interactive curses mode not supported; use -d T or -d H")
		return 2
	}

	// The reference may be supplied positionally as the second argument or
	// via -T/--reference; an explicit -T wins.
	input := fs.Arg(0)
	if refFasta == "" && fs.NArg() >= 2 {
		refFasta = fs.Arg(1)
	}

	opts := samtools.TviewOptions{
		Input:       input,
		Reference:   refFasta,
		Position:    position,
		Width:       width,
		Sample:      sample,
		Mode:        mode,
		HideInserts: hideInserts,
	}
	if err := samtools.Tview(opts, os.Stdout); err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err)
		return 1
	}
	return 0
}
