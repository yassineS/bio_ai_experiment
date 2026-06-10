// Package cliflag provides utilities for CLI flag handling with both short and
// long options.
//
// Every helper registers the same destination under an optional POSIX-style
// short name (single dash, e.g. -i) and an optional GNU-style long name (double
// dash, e.g. --input) on a single *flag.FlagSet, so a tool can accept both
// forms without hand-rolling two registrations. Either name may be empty to
// register only one form. The long name carries the usage string; the short
// name is registered with an empty usage so it does not produce a duplicate
// line in flag.PrintDefaults output. See docs/CLI_CONVENTIONS.md for the
// project-wide flag conventions these helpers implement.
package cliflag

import (
	"flag"
	"fmt"
	"strings"
	"time"
)

// StringVar defines a string flag with both short and long names.
func StringVar(fs *flag.FlagSet, p *string, short, long, defaultValue, usage string) {
	if short != "" {
		fs.StringVar(p, short, defaultValue, "")
	}
	if long != "" {
		fs.StringVar(p, long, defaultValue, usage)
	}
}

// IntVar defines an int flag with both short and long names.
func IntVar(fs *flag.FlagSet, p *int, short, long string, defaultValue int, usage string) {
	if short != "" {
		fs.IntVar(p, short, defaultValue, "")
	}
	if long != "" {
		fs.IntVar(p, long, defaultValue, usage)
	}
}

// Int64Var defines an int64 flag with both short and long names. It is useful
// for flags holding large counts or random seeds that may exceed the int range
// on 32-bit platforms.
func Int64Var(fs *flag.FlagSet, p *int64, short, long string, defaultValue int64, usage string) {
	if short != "" {
		fs.Int64Var(p, short, defaultValue, "")
	}
	if long != "" {
		fs.Int64Var(p, long, defaultValue, usage)
	}
}

// Uint64Var defines a uint64 flag with both short and long names.
func Uint64Var(fs *flag.FlagSet, p *uint64, short, long string, defaultValue uint64, usage string) {
	if short != "" {
		fs.Uint64Var(p, short, defaultValue, "")
	}
	if long != "" {
		fs.Uint64Var(p, long, defaultValue, usage)
	}
}

// Float64Var defines a float64 flag with both short and long names.
func Float64Var(fs *flag.FlagSet, p *float64, short, long string, defaultValue float64, usage string) {
	if short != "" {
		fs.Float64Var(p, short, defaultValue, "")
	}
	if long != "" {
		fs.Float64Var(p, long, defaultValue, usage)
	}
}

// BoolVar defines a bool flag with both short and long names.
func BoolVar(fs *flag.FlagSet, p *bool, short, long string, defaultValue bool, usage string) {
	if short != "" {
		fs.BoolVar(p, short, defaultValue, "")
	}
	if long != "" {
		fs.BoolVar(p, long, defaultValue, usage)
	}
}

// DurationVar defines a time.Duration flag with both short and long names.
func DurationVar(fs *flag.FlagSet, p *time.Duration, short, long string, defaultValue time.Duration, usage string) {
	if short != "" {
		fs.DurationVar(p, short, defaultValue, "")
	}
	if long != "" {
		fs.DurationVar(p, long, defaultValue, usage)
	}
}

// Var defines a flag backed by a custom flag.Value under both short and long
// names. Use it for repeatable flags (where each occurrence accumulates a
// value) or any flag with bespoke parsing. Because both names target the same
// value, repeated occurrences under either form accumulate into the same
// destination.
func Var(fs *flag.FlagSet, value flag.Value, short, long, usage string) {
	if short != "" {
		fs.Var(value, short, "")
	}
	if long != "" {
		fs.Var(value, long, usage)
	}
}

// FormatUsage formats the usage message to show both short and long option names.
func FormatUsage(short, long, valueType, description string) string {
	var names []string
	if short != "" {
		if valueType != "" {
			names = append(names, fmt.Sprintf("-%s %s", short, valueType))
		} else {
			names = append(names, fmt.Sprintf("-%s", short))
		}
	}
	if long != "" {
		if valueType != "" {
			names = append(names, fmt.Sprintf("--%s %s", long, valueType))
		} else {
			names = append(names, fmt.Sprintf("--%s", long))
		}
	}
	return fmt.Sprintf("  %-30s %s", strings.Join(names, ", "), description)
}
