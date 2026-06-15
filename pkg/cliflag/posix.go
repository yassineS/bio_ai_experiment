package cliflag

import (
	"flag"
	"fmt"
)

// boolFlag is the interface Go's flag package uses internally to mark a
// flag.Value as boolean (its *flag.boolValue implements it, returning true).
// We probe for it so the normalizer can tell value-taking flags from boolean
// switches without depending on the unexported boolValue type.
type boolFlag interface {
	IsBoolFlag() bool
}

// isBoolFlag reports whether the flag registered under name on fs is a
// boolean switch (one that takes no value, like -b). It returns false when
// the flag is unknown or is value-taking.
func isBoolFlag(fs *flag.FlagSet, name string) bool {
	f := fs.Lookup(name)
	if f == nil {
		return false
	}
	if bf, ok := f.Value.(boolFlag); ok {
		return bf.IsBoolFlag()
	}
	return false
}

// Normalize expands POSIX getopt-style short-option clusters in args into the
// canonical one-flag-per-token form that Go's flag.FlagSet.Parse understands,
// using fs to introspect which short flags are boolean switches and which take
// a value. It returns the rewritten argument slice without mutating args or
// parsing anything.
//
// The expansion implements standard getopt semantics:
//
//   - A token "-XYZ..." (single dash, two-or-more characters, not "--...") is
//     a short-option cluster. Each character is looked up on fs. While the
//     character names a boolean flag it is emitted as its own "-X" token and
//     scanning continues. The first character that names a value-taking flag
//     consumes the REST of the cluster as its value: "-qYZ" becomes
//     "-q" "YZ". When that value-taking flag is the last character of the
//     cluster ("-bSq"), only "-q" is emitted so fs.Parse takes the value from
//     the following argument ("-bSq 20").
//   - Long options ("--foo", "--foo=bar") pass through untouched — Go's flag
//     package parses them directly.
//   - The "--" terminator stops normalization; it and everything after it are
//     passed through verbatim as positionals.
//   - A bare "-" (stdin/stdout) is passed through as a positional, never
//     treated as a flag cluster.
//   - When a value-taking short flag ends a cluster (or stands alone) so that
//     its value is the following argument, that following argument is passed
//     through verbatim — even if it begins with a dash ("-q -5" keeps "-5" as
//     the value), matching getopt's "the next argument is the option-argument"
//     rule and Go's flag package, which likewise consumes the next token.
//   - A character in a cluster that is not a registered short flag yields a
//     "flag provided but not defined" error, mirroring how upstream getopt
//     tools reject unknown options. The caller can treat this as exit code 2.
//   - Tokens that are already canonical (a lone "-x", a value-less single
//     short flag, positionals) are returned unchanged, so Normalize is
//     idempotent.
//
// Normalize is the foundation for cliflag.Parse; tools that only need the
// rewritten args (for example to inspect them before parsing) can call it
// directly.
func Normalize(fs *flag.FlagSet, args []string) ([]string, error) {
	out := make([]string, 0, len(args))
	for i := 0; i < len(args); i++ {
		arg := args[i]
		switch {
		case arg == "--":
			// Everything from the terminator onward is positional.
			out = append(out, args[i:]...)
			return out, nil
		case arg == "-" || arg == "":
			// Bare dash is stdin/stdout; empty stays empty. Pass through.
			out = append(out, arg)
		case len(arg) >= 2 && arg[0] == '-' && arg[1] == '-':
			// Long option ("--foo" / "--foo=bar"): hand straight to flag.
			out = append(out, arg)
		case arg[0] == '-':
			// Short-option cluster: "-XYZ...".
			expanded, wantsNext, err := expandCluster(fs, arg)
			if err != nil {
				return nil, err
			}
			out = append(out, expanded...)
			// The cluster ended on a value-taking flag with no inline value,
			// so its option-argument is the next token. Pass it through
			// verbatim (it may legitimately start with '-', e.g. "-q -5") and
			// skip re-normalizing it as a cluster of its own.
			if wantsNext && i+1 < len(args) {
				out = append(out, args[i+1])
				i++
			}
		default:
			// Positional argument.
			out = append(out, arg)
		}
	}
	return out, nil
}

// expandCluster rewrites a single short-option cluster token (e.g. "-bSq20")
// into canonical tokens. It assumes arg starts with a single '-' followed by
// at least one character. The returned wantsNext is true when the cluster
// ended on a value-taking flag that had no inline value, meaning the caller
// should treat the following argument as that flag's option-argument.
func expandCluster(fs *flag.FlagSet, arg string) (tokens []string, wantsNext bool, err error) {
	chars := arg[1:] // strip the leading '-'
	out := make([]string, 0, len(chars))
	for j := 0; j < len(chars); j++ {
		name := string(chars[j])
		if fs.Lookup(name) == nil {
			return nil, false, fmt.Errorf("flag provided but not defined: -%s", name)
		}
		if isBoolFlag(fs, name) {
			// Boolean switch: emit it and keep scanning the cluster.
			out = append(out, "-"+name)
			continue
		}
		// Value-taking flag terminates the cluster: the rest of the cluster
		// is its inline value, or — when the flag was the last character —
		// its value is the next argument.
		out = append(out, "-"+name)
		if rest := chars[j+1:]; rest != "" {
			out = append(out, rest)
			return out, false, nil
		}
		return out, true, nil
	}
	return out, false, nil
}

// flagName extracts the registered flag name from an option token and reports
// whether the token already carries an inline value. For "--chrom=MT" it
// returns ("chrom", true); for "--chrom" or "-c" it returns the bare name and
// false. Both single-dash short flags and double-dash long flags are handled.
// It assumes arg begins with at least one '-'.
func flagName(arg string) (name string, hasInline bool) {
	s := arg
	if len(s) >= 2 && s[0] == '-' && s[1] == '-' {
		s = s[2:]
	} else if len(s) >= 1 && s[0] == '-' {
		s = s[1:]
	}
	for i := 0; i < len(s); i++ {
		if s[i] == '=' {
			return s[:i], true
		}
	}
	return s, false
}

// Permute reorders an already-Normalize-d argument slice so that every option
// token (and the value it consumes) precedes the positional arguments,
// mirroring the GNU getopt / docopt permutation that the upstream tools rely
// on. Go's flag.FlagSet.Parse stops scanning options at the first non-flag
// argument, so without this pass an interspersed command line such as
// "PREFIX in.bam --chrom MT" would treat "--chrom" and "MT" as a third and
// fourth (rejected) positional. After permutation, fs.Parse sees the canonical
// flags-first form.
//
// Permute expects the output of Normalize: short-flag clusters already expanded
// and inline short-flag values already split out, so each token is one of a
// canonical option ("-c", "--chrom", "--chrom=MT"), an option value, the "--"
// terminator, a bare "-", or a positional. It uses fs to tell value-taking
// flags from boolean switches: a value-taking flag with no inline value
// consumes the following token as its option-argument and that token travels
// with the flag — even when it begins with '-', as a negative-number value such
// as "-c" "-5" does. The "--" terminator and everything after it are treated as
// positional verbatim (the terminator itself is kept with the options so
// fs.Parse still stops option scanning at the right place). The relative order
// of the positionals is preserved, so subcommands whose positional order is
// significant (for example "<in.fa> <reg.bed>") are unaffected. When the input
// is already in flags-first order the output equals the input, making Permute a
// no-op for non-interspersed command lines.
func Permute(fs *flag.FlagSet, args []string) []string {
	options := make([]string, 0, len(args))
	positionals := make([]string, 0, len(args))
	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "--" {
			// End of options: the terminator and every remaining token are
			// positional. Keep the terminator with the options so fs.Parse
			// stops option scanning at the right place.
			options = append(options, arg)
			positionals = append(positionals, args[i+1:]...)
			break
		}
		if len(arg) >= 2 && arg[0] == '-' && arg != "-" {
			// An option token. Determine whether it carries an inline value
			// ("--chrom=MT" or, post-Normalize, a value-concatenated short
			// flag) or consumes the following token as its value.
			options = append(options, arg)
			name, hasInline := flagName(arg)
			if !hasInline && !isBoolFlag(fs, name) && i+1 < len(args) {
				// Value-taking flag with no inline value: the next token is its
				// option-argument and must travel with it, even if it begins
				// with '-' (a negative-number value such as "-5").
				options = append(options, args[i+1])
				i++
			}
			continue
		}
		// A bare "-" (stdin/stdout) or any other token is a positional.
		positionals = append(positionals, arg)
	}
	return append(options, positionals...)
}

// Parse normalizes args with Normalize, permutes options ahead of positionals
// with Permute, and then calls fs.Parse on the result. It gives any tool wired
// through it POSIX getopt-compatible short-flag bundling ("-bS" == "-b -S"),
// value concatenation ("-q20" == "-q 20"), and — like upstream's getopt/docopt
// parsers — acceptance of flags interspersed among or after the positional
// arguments ("comp file.fa -l 60" == "comp -l 60 file.fa"), all on top of Go's
// native long-flag handling. Only options move; the relative order of the
// positionals is preserved, and "--" still terminates option scanning. It is a
// drop-in replacement for fs.Parse(args): on the error path the caller should
// print usage and exit with code 2, exactly as it would for a plain fs.Parse
// failure.
func Parse(fs *flag.FlagSet, args []string) error {
	normalized, err := Normalize(fs, args)
	if err != nil {
		return err
	}
	return fs.Parse(Permute(fs, normalized))
}
