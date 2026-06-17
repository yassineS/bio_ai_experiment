// Package bedmakewindows partitions intervals (either entire chromosomes from
// a genome-sizes file, or arbitrary BED intervals) into fixed-size or
// fixed-count windows. It mirrors `bedtools makewindows`.
//
// Two partition strategies are supported:
//
//   - FixedWidth (Width > 0): each window is up to Width bases wide, sliding
//     forward by Step bases (Step defaults to Width when zero). The final
//     window per interval is clipped to the interval end; with Step < Width
//     you get multiple short tail windows, matching upstream.
//
//   - FixedCount (Count > 0): the interval is split into Count windows of
//     as-equal length as possible. Following upstream, length is floor-divided
//     by Count, then the final window absorbs any remainder. Intervals shorter
//     than Count bases are skipped with the upstream warning
//     "WARNING: Interval CHR:START-END is smaller than the number of windows
//     requested. Skipping."
//
// Window naming is configurable via Naming:
//
//   - NoName        — no name column emitted; output is BED3.
//   - NameWinNum    — append the 1-based window index (per source interval).
//   - NameSrcWinNum — append "<src_name>_<index>" using the source name
//     (must come from a BED interval input).
//   - NameSrc       — append the source name unchanged.
//
// Reverse flips the index direction so the last window per interval carries
// index 1 (only meaningful for NameWinNum / NameSrcWinNum).
package bedmakewindows

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// Naming selects the per-window name annotation.
type Naming int

const (
	// NoName emits no fourth column (output is BED3).
	NoName Naming = iota
	// NameWinNum appends the 1-based window index.
	NameWinNum
	// NameSrcWinNum appends "<src_name>_<index>".
	NameSrcWinNum
	// NameSrc appends the unmodified source name.
	NameSrc
)

// ParseNaming parses the bedtools `-i` argument.
//
// Upstream `bedtools makewindows` accepts only "src", "winnum", and
// "srcwinnum"; when -i is omitted entirely the ID method defaults to
// ID_NONE (no name column, BED3 output). We map the empty string to that
// same default so the CLI's "-i" default sentinel produces BED3, and we
// additionally accept the literal "none" as an explicit spelling of the
// default. Accepting "none" is a documented fix-on-port superset: upstream
// errors on "-i none", whereas we treat it as the no-name default.
func ParseNaming(s string) (Naming, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "", "none":
		return NoName, nil
	case "winnum":
		return NameWinNum, nil
	case "srcwinnum":
		return NameSrcWinNum, nil
	case "src":
		return NameSrc, nil
	}
	return NoName, fmt.Errorf("unknown -i mode %q (want src, winnum, or srcwinnum)", s)
}

// Options bundles the configuration for MakeWindows.
type Options struct {
	// Width is the fixed window size in bases. Mutually exclusive with Count.
	Width int
	// Step is the slide between consecutive windows. Defaults to Width when 0.
	Step int
	// Count selects the FixedCount strategy: partition each interval into
	// Count windows of as-equal size as possible. Mutually exclusive with
	// Width / Step.
	Count int
	// Reverse flips per-interval window numbering so the last window is "1".
	Reverse bool
	// Naming controls the name column annotation.
	Naming Naming
}

// validate normalises Options and returns an error for inconsistent settings.
func (o *Options) validate() error {
	if o.Width > 0 && o.Count > 0 {
		return fmt.Errorf("-w and -n are mutually exclusive")
	}
	if o.Width <= 0 && o.Count <= 0 {
		return fmt.Errorf("must specify -w or -n")
	}
	if o.Width > 0 {
		if o.Step <= 0 {
			o.Step = o.Width
		}
	} else {
		if o.Step != 0 {
			return fmt.Errorf("-s is only meaningful with -w")
		}
	}
	return nil
}

// Interval is a chromosome region to partition.
type Interval struct {
	Chrom string
	Start int
	End   int
	Name  string // source name; the chromosome name for genome-file intervals
}

// FromGenome reads a chrom-sizes file ("chrom<TAB>size" lines) and returns one
// interval per chromosome, sorted by name for deterministic output.
func FromGenome(r io.Reader) ([]Interval, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []Interval
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("genome line %q must have at least 2 fields", line)
		}
		size, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid chromosome size %q: %v", fields[1], err)
		}
		if size <= 0 {
			continue
		}
		// Upstream constructs each genome interval as BED(chrom,0,size,chrom,...)
		// — the source name is the chromosome name, so `-i src` / `-i srcwinnum`
		// over a genome file annotate with the chrom (windowMaker.cpp:56).
		out = append(out, Interval{Chrom: fields[0], Start: 0, End: size, Name: fields[0]})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	sort.SliceStable(out, func(i, j int) bool { return out[i].Chrom < out[j].Chrom })
	return out, nil
}

// FromBED reads BED records and returns them as intervals. The 4th column (if
// present) is captured as the source name.
func FromBED(r io.Reader) ([]Interval, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []Interval
	for scanner.Scan() {
		raw := scanner.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return nil, fmt.Errorf("BED record must have at least 3 fields: %q", raw)
		}
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		iv := Interval{Chrom: fields[0], Start: start, End: end}
		if len(fields) >= 4 {
			iv.Name = fields[3]
		}
		out = append(out, iv)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// MakeWindows partitions each interval according to opts and writes the
// resulting windows to out. warn receives non-fatal warnings (e.g. intervals
// dropped under FixedCount). Returns the number of windows written.
func MakeWindows(intervals []Interval, out io.Writer, warn io.Writer, opts Options) (int, error) {
	if err := opts.validate(); err != nil {
		return 0, err
	}
	bw := bufio.NewWriter(out)
	defer bw.Flush()
	written := 0
	for _, iv := range intervals {
		var wins [][2]int
		if opts.Width > 0 {
			wins = fixedWidth(iv.Start, iv.End, opts.Width, opts.Step)
		} else {
			ws, skip := fixedCount(iv.Start, iv.End, opts.Count)
			if skip && warn != nil {
				fmt.Fprintf(warn, "WARNING: Interval %s:%d-%d is smaller than the number of windows requested. Skipping.\n",
					iv.Chrom, iv.Start, iv.End)
			}
			wins = ws
		}
		for i, w := range wins {
			idx := i + 1
			if opts.Reverse {
				idx = len(wins) - i
			}
			row := formatRow(iv, w[0], w[1], idx, opts.Naming)
			if _, err := bw.WriteString(row); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			written++
		}
	}
	return written, nil
}

// fixedWidth returns the per-interval windows for the FixedWidth strategy.
// The slide is `step`; the last window per interval is clipped to `end`.
//
// To match upstream `bedtools makewindows`, the loop continues advancing by
// `step` for as long as `s < end`. The final clipped windows can therefore be
// shorter than `width`, and with `step < width` you get multiple short tail
// windows (see makewindows.t03 — BBB:73000-90000 with -w 5000 -s 2000 emits
// 85000-90000, 87000-90000, and 89000-90000).
func fixedWidth(start, end, width, step int) [][2]int {
	var wins [][2]int
	if width <= 0 || step <= 0 {
		return wins
	}
	for s := start; s < end; s += step {
		e := s + width
		if e > end {
			e = end
		}
		if e <= s {
			break
		}
		wins = append(wins, [2]int{s, e})
	}
	return wins
}

// fixedCount returns the per-interval windows for the FixedCount strategy.
// When length < count the caller emits a warning and the boolean return is
// true (no windows produced).
func fixedCount(start, end, count int) ([][2]int, bool) {
	length := end - start
	if count <= 0 || length < count {
		return nil, length > 0
	}
	sz := length / count
	if sz <= 0 {
		return nil, true
	}
	wins := make([][2]int, 0, count)
	for i := 0; i < count; i++ {
		s := start + i*sz
		e := s + sz
		if i == count-1 {
			// Final window absorbs any remainder.
			e = end
		}
		wins = append(wins, [2]int{s, e})
	}
	return wins, false
}

// formatRow renders one BED row according to the naming mode.
func formatRow(iv Interval, start, end, idx int, naming Naming) string {
	base := fmt.Sprintf("%s\t%d\t%d", iv.Chrom, start, end)
	switch naming {
	case NoName:
		return base
	case NameWinNum:
		return fmt.Sprintf("%s\t%d", base, idx)
	case NameSrcWinNum:
		return fmt.Sprintf("%s\t%s_%d", base, iv.Name, idx)
	case NameSrc:
		return fmt.Sprintf("%s\t%s", base, iv.Name)
	}
	return base
}
