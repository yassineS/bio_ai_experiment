package main

import (
	"flag"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
)

// This file holds the bcftools-call CLI helpers salvaged from PR #219
// (originally split across subcmds.go and main.go). They live in their own
// file to keep the salvage self-contained and avoid clobbering main's
// independent subcmds.go/main.go changes.

// parseFlags is the shared parse entry point for every bcftools
// subcommand. It (1) ensures `--no-version` is accepted (see
// registerNoVersionIfAbsent) and (2) routes parsing through
// cliflag.Parse, which applies POSIX getopt normalization — short-flag
// bundling (`-hG` == `-h -G`), attached/value-concatenated short flags
// (`-Ob` == `-O b`, `-m-` == `-m -`, `-mboth` == `-m both`), the `--`
// end-of-options terminator, and bare `-` (stdin/stdout) — before
// delegating to fs.Parse. This gives every bcftools subcommand the same
// upstream getopt-compatible argument handling that samtools view
// already uses. On error the caller prints usage and exits with code 2,
// exactly as for a plain fs.Parse failure.
func parseFlags(fs *flag.FlagSet, args []string) error {
	registerNoVersionIfAbsent(fs)
	return cliflag.Parse(fs, args)
}

// registerNoVersionIfAbsent registers a no-op `--no-version` boolean flag
// on fs when one is not already present. Upstream bcftools accepts
// `--no-version` as a per-subcommand option that suppresses the
// `##bcftools_*Version`/`##bcftools_*Command` provenance header lines.
// Subcommands that already register `--no-version` (and wire it into
// their options) keep their own registration; this only fills the gap for
// subcommands such as `view` and `norm` which never emit a provenance
// line, so accepting the flag is a safe no-op.
func registerNoVersionIfAbsent(fs *flag.FlagSet) {
	if fs.Lookup("no-version") != nil {
		return
	}
	var noVersion bool
	fs.BoolVar(&noVersion, "no-version", false, "Accepted; no provenance line is emitted.")
}

func preprocessOptionalArg(args []string, flagName, defaultVal string) []string {
	out := make([]string, 0, len(args))
	for i := 0; i < len(args); i++ {
		a := args[i]
		if a == flagName {
			// Peek next arg: if it looks like a value (no leading
			// '-'), assume the user meant `--flag VALUE` and pass
			// through; otherwise expand to `--flag=defaultVal`.
			if i+1 < len(args) {
				next := args[i+1]
				if next == "" || next[0] != '-' {
					// Check if next looks like a positional file (e.g. has
					// '/' or ends in .vcf / .bcf / .gz) — heuristic.
					if looksLikePositional(next) {
						out = append(out, flagName+"="+defaultVal)
						continue
					}
				}
			} else {
				// Trailing bare flag.
				out = append(out, flagName+"="+defaultVal)
				continue
			}
		}
		out = append(out, a)
	}
	return out
}

// looksLikePositional reports whether s looks more like an input file
// path than an option value (heuristic: contains '/', '.', or ends in
// common suffixes).
func looksLikePositional(s string) bool {
	if strings.ContainsAny(s, "/") {
		return true
	}
	for _, suf := range []string{".vcf", ".bcf", ".gz", ".bgz", ".tbi", ".csi"} {
		if strings.HasSuffix(s, suf) {
			return true
		}
	}
	return false
}

// priorANFromSpec returns the first field of "AN,AC".
func priorANFromSpec(s string) string {
	if s == "" {
		return ""
	}
	if i := strings.IndexByte(s, ','); i >= 0 {
		return strings.TrimSpace(s[:i])
	}
	return strings.TrimSpace(s)
}

// priorACFromSpec returns the second field of "AN,AC".
func priorACFromSpec(s string) string {
	if i := strings.IndexByte(s, ','); i >= 0 {
		return strings.TrimSpace(s[i+1:])
	}
	return ""
}
