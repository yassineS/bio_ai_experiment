package main

import (
	"flag"
	"strings"
)

// This file holds the bcftools-call CLI helpers salvaged from PR #219
// (originally split across subcmds.go and main.go). They live in their own
// file to keep the salvage self-contained and avoid clobbering main's
// independent subcmds.go/main.go changes.

// boolFlag is the optional interface the standard library's flag package
// uses to recognise boolean flags (those that do not consume a following
// value). We use it to distinguish value-taking short flags from boolean
// ones when normalising getopt-style attached values.
type boolFlag interface {
	IsBoolFlag() bool
}

// valueTakingShortFlags inspects fs and returns the set of registered
// single-character flag names that consume a value (i.e. are NOT boolean
// flags). These are the only short flags for which an attached value
// (`-Xvalue`) is meaningful in upstream getopt semantics.
func valueTakingShortFlags(fs *flag.FlagSet) map[byte]bool {
	set := make(map[byte]bool)
	fs.VisitAll(func(f *flag.Flag) {
		if len(f.Name) != 1 {
			return
		}
		if bf, ok := f.Value.(boolFlag); ok && bf.IsBoolFlag() {
			return
		}
		set[f.Name[0]] = true
	})
	return set
}

// normalizeShortFlags rewrites getopt-style attached short-flag values
// into the two-token form that Go's flag package accepts. Upstream
// bcftools is getopt-based and accepts a value attached directly to a
// single-letter flag (e.g. `-Ob`, `norm -m-`, `-m+`, `-mboth`); Go's
// flag package only accepts `-X value` or `-X=value`. For each argument
// of the form `-X...` where X is a registered value-taking short flag
// and extra characters follow X, it splits the token into `-X` and the
// remainder. So `-Ob` -> `-O b`, `-m-` -> `-m -`, `-mboth` -> `-m both`.
//
// It deliberately leaves untouched:
//   - long flags (`--foo`, `--foo=bar`),
//   - boolean short flags (which take no value),
//   - the `-X=value` form (already valid; passed through),
//   - a bare `-` (stdin/stdout),
//   - everything after a bare `--` (end-of-options marker).
func normalizeShortFlags(fs *flag.FlagSet, args []string) []string {
	values := valueTakingShortFlags(fs)
	out := make([]string, 0, len(args)+2)
	for i, a := range args {
		if a == "--" {
			out = append(out, args[i:]...)
			break
		}
		// Candidate: a single-dash flag with at least one char after
		// the flag letter, and not a long flag (`--`) or bare `-`.
		if len(a) > 2 && a[0] == '-' && a[1] != '-' {
			if values[a[1]] && a[2] != '=' {
				out = append(out, a[:2], a[2:])
				continue
			}
		}
		out = append(out, a)
	}
	return out
}

// parseFlags is the shared parse entry point for every bcftools
// subcommand. It (1) ensures `--no-version` is accepted (see
// registerNoVersionIfAbsent) and (2) normalises getopt-style attached
// short-flag values (see normalizeShortFlags) so all subcommands accept
// the upstream attached short-flag forms (`-Ob`, `-m-`, ...) and
// `--no-version` uniformly, before delegating to fs.Parse.
func parseFlags(fs *flag.FlagSet, args []string) error {
	registerNoVersionIfAbsent(fs)
	return fs.Parse(normalizeShortFlags(fs, args))
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
