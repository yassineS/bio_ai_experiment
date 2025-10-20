// Package cliflag provides utilities for CLI flag handling with both short and long options.
package cliflag

import (
	"flag"
	"fmt"
	"strings"
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
