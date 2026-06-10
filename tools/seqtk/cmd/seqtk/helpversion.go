package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
)

// subFlags bundles the -h/--help and -v/--version switches that every seqtk
// subcommand registers, so each subcommand responds uniformly to them.
type subFlags struct {
	help    bool
	version bool
}

// register wires --help and -v/--version (plus -h when withShortH is true) onto
// fs. mergefa already binds -h to --haploid for upstream parity, so it passes
// withShortH=false and keeps -h as the haploid switch while still exposing
// --help for help.
func (s *subFlags) register(fs *flag.FlagSet, withShortH bool) {
	shortH := "h"
	if !withShortH {
		shortH = ""
	}
	cliflag.BoolVar(fs, &s.help, shortH, "help", false, "Show this help message and exit")
	cliflag.BoolVar(fs, &s.version, "v", "version", false, "Show version information and exit")
}

// handle prints usage (help) or the version banner and exits 0 when the
// corresponding switch was set. It must be called immediately after parsing so
// --help/--version short-circuit before any work happens.
func (s *subFlags) handle(fs *flag.FlagSet) {
	if s.help {
		fs.Usage()
		os.Exit(0)
	}
	if s.version {
		fmt.Printf("seqtk version %s\n", version)
		os.Exit(0)
	}
}
