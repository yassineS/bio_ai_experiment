// Package samtools — shared SAM-FLAG spec parser.
//
// Upstream samtools accepts the SAM FLAG specification in two forms on
// every subcommand that exposes `--rf/--ff/--incl-flags/--excl-flags`:
//
//   - A non-negative integer (decimal, 0x hex, or 0 octal).
//   - A comma-separated list of mnemonic names from the SAM spec table:
//     PAIRED, PROPER_PAIR, UNMAP, MUNMAP, REVERSE, MREVERSE, READ1,
//     READ2, SECONDARY, QCFAIL, DUP, SUPPLEMENTARY.
//
// ParseFlagSpec mirrors `bam_str2flag` (sam.c) — case-insensitive name
// matching, both '|' and ',' as separators, returning a uint16 mask.
package samtools

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// flagNames maps each SAM-spec flag mnemonic to its bit mask. Aliases
// match upstream's bam_str2flag accept set (e.g. "PAIRED" and "P_PAIR"
// for the same bit, though htslib's table primarily uses the canonical
// names below). Names compare case-insensitively.
var flagNames = map[string]uint16{
	"PAIRED":        sam.FlagPaired,
	"PROPER_PAIR":   sam.FlagProperPair,
	"PROPER":        sam.FlagProperPair,
	"UNMAP":         sam.FlagUnmapped,
	"MUNMAP":        sam.FlagMateUnmapped,
	"REVERSE":       sam.FlagReverse,
	"MREVERSE":      sam.FlagMateReverse,
	"READ1":         sam.FlagRead1,
	"READ2":         sam.FlagRead2,
	"SECONDARY":     sam.FlagSecondary,
	"QCFAIL":        sam.FlagQCFail,
	"DUP":           sam.FlagDuplicate,
	"DUPLICATE":     sam.FlagDuplicate,
	"SUPPLEMENTARY": sam.FlagSupplementary,
}

// ParseFlagSpec parses the upstream `--rf/--ff` flag spec. Numeric
// inputs are accepted in decimal, hex (0x prefix) or octal (leading 0).
// Mnemonic inputs may be comma- or '|'-separated; matching is
// case-insensitive. Empty input returns 0.
func ParseFlagSpec(s string) (uint16, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return 0, nil
	}
	// Try integer first.
	if v, err := strconv.ParseUint(s, 0, 16); err == nil {
		return uint16(v), nil
	}
	var out uint16
	for _, tok := range splitFlagSpec(s) {
		t := strings.ToUpper(strings.TrimSpace(tok))
		if t == "" {
			continue
		}
		bit, ok := flagNames[t]
		if !ok {
			return 0, fmt.Errorf("unrecognised SAM FLAG name %q", tok)
		}
		out |= bit
	}
	return out, nil
}

func splitFlagSpec(s string) []string {
	// Split on commas and pipes.
	out := []string{}
	cur := strings.Builder{}
	for _, r := range s {
		if r == ',' || r == '|' {
			out = append(out, cur.String())
			cur.Reset()
			continue
		}
		cur.WriteRune(r)
	}
	if cur.Len() > 0 {
		out = append(out, cur.String())
	}
	return out
}
