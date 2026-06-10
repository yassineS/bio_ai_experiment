package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
)

// subFlags bundles the -v/--version switch that every prinseq subcommand
// registers so each responds uniformly to it. The -h/--help switch is left to
// the flag package's built-in ErrHelp handling (which prints usage and exits 0
// under flag.ExitOnError), matching prinseq's reliance on plain fs.Parse for
// its single-dash long-option compatibility.
type subFlags struct {
	version bool
}

// register wires -v/--version onto fs. No prinseq subcommand binds -v to
// another option, so the short form is available everywhere.
func (s *subFlags) register(fs *flag.FlagSet) {
	cliflag.BoolVar(fs, &s.version, "v", "version", false, "Show version information and exit")
}

// handle prints the version banner and exits 0 when -v/--version was set. Call
// it immediately after parsing.
func (s *subFlags) handle() {
	if s.version {
		fmt.Printf("prinseq version %s\n", version)
		os.Exit(0)
	}
}
