// Command conformance runs the spec-conformance and silent-corruption
// edge-case test batteries and prints a human-readable PASS/SKIP/FAIL summary.
//
// It is a thin driver over `go test` for the two battery packages:
//
//	pipeline/conformance  — the originals' own corpora (htslib test fixtures,
//	                        htscodecs compliance vectors) through our binaries.
//	pipeline/edgecases    — discrete named tests for each silent-corruption
//	                        class (CRAM reference handling, bcftools norm
//	                        re-indexing, index byte-identity, sort stability,
//	                        MD/NM tags, QUAL/PL ULP non-impact).
//
// Tests SKIP (rather than fail) when a prerequisite corpus or upstream binary
// is absent, so a partial environment still produces a meaningful report.
//
// Usage:
//
//	go run ./pipeline/cmd/conformance            # both batteries, verbose
//	go run ./pipeline/cmd/conformance -run CRAM  # filter by test-name regexp
//	go run ./pipeline/cmd/conformance -json       # machine-readable go test -json
package main

import (
	"flag"
	"fmt"
	"os"
	"os/exec"
	"strings"
)

func main() {
	var (
		runFlag  = flag.String("run", "", "test-name regexp passed through to `go test -run`")
		jsonFlag = flag.Bool("json", false, "emit raw `go test -json` (machine-readable)")
		pkgsFlag = flag.String("pkgs", "conformance,edgecases", "comma-separated battery packages to run")
	)
	flag.Parse()

	pkgMap := map[string]string{
		"conformance": "./pipeline/conformance/...",
		"edgecases":   "./pipeline/edgecases/...",
	}
	var pkgs []string
	for _, p := range strings.Split(*pkgsFlag, ",") {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		path, ok := pkgMap[p]
		if !ok {
			fmt.Fprintf(os.Stderr, "unknown battery package %q (want one of conformance, edgecases)\n", p)
			os.Exit(2)
		}
		pkgs = append(pkgs, path)
	}

	args := []string{"test"}
	if *jsonFlag {
		args = append(args, "-json")
	} else {
		args = append(args, "-v")
	}
	if *runFlag != "" {
		args = append(args, "-run", *runFlag)
	}
	args = append(args, pkgs...)

	fmt.Fprintf(os.Stderr, "# go %s\n", strings.Join(args, " "))
	cmd := exec.Command("go", args...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		// go test exits non-zero on FAIL; propagate so CI / callers notice.
		if ee, ok := err.(*exec.ExitError); ok {
			os.Exit(ee.ExitCode())
		}
		fmt.Fprintf(os.Stderr, "running go test: %v\n", err)
		os.Exit(1)
	}
}
