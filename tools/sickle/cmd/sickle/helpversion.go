package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
)

// sickleVersion is the version banner emitted by -v/--version at both the
// top level and every subcommand.
const sickleVersion = "sickle version 1.0.0 (Go implementation)"

// subFlags bundles the -h/--help and -v/--version switches that every sickle
// subcommand registers so each responds uniformly to them.
type subFlags struct {
	help    bool
	version bool
}

// register wires -h/--help and -v/--version onto fs. None of sickle's
// subcommands bind -h or -v to another option, so both short forms are
// available.
func (s *subFlags) register(fs *flag.FlagSet) {
	cliflag.BoolVar(fs, &s.help, "h", "help", false, "Show this help message and exit")
	cliflag.BoolVar(fs, &s.version, "v", "version", false, "Show version information and exit")
}

// handle prints usage (help) or the version banner and exits 0 when the
// corresponding switch was set. Call it immediately after parsing.
func (s *subFlags) handle(fs *flag.FlagSet) {
	if s.help {
		fs.Usage()
		os.Exit(0)
	}
	if s.version {
		fmt.Println(sickleVersion)
		os.Exit(0)
	}
}
