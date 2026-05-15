// Package bednuc implements `bedtools nuc`: for each interval in a BED file
// it fetches the corresponding sequence from a random-access FASTA reference
// and emits a per-interval nucleotide composition profile.
//
// Output (default): the original BED columns followed by:
//
//	%AT  %GC  #A  #C  #G  #T  #N  #oth  seq_len
//
// With `-seq`, the actual sequence is appended as a final column. With
// `-pattern STR`, the count of substring occurrences is appended (after
// the optional sequence column, matching upstream's column order).
//
// A one-line header starting with `#` describes the columns:
//
//	#1_usercol	2_usercol	...	N_pct_at	N+1_pct_gc	...
//
// where the leading `1..N_usercol` block has one column per column observed
// in the *first* BED record (matching upstream's behaviour, which derives
// the bedType from the first record).
//
// Strand-aware mode (`-s`) reverse-complements `-` strand intervals before
// counting (so A/T/G/C/N counts are computed against the transcript
// orientation). The `-C` flag toggles case-insensitive pattern matching
// (upstream default is case-sensitive — counter to the option's mnemonic).
//
// The `-fullHeader` flag tells the FASTA index to look up contigs by the
// raw header (including whitespace). Our `pkg/bioformats/fasta` index
// always uses the first whitespace-delimited token, so when `-fullHeader`
// is set we additionally consult a "full header" name map built lazily on
// the fly.
package bednuc

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
)

// Options configures Run.
type Options struct {
	// PrintSeq emits the raw (post-strand) sequence as a column.
	PrintSeq bool
	// Pattern, if non-empty, enables substring counting per interval.
	Pattern string
	// HasPattern is true when Pattern was explicitly supplied; an empty
	// pattern with HasPattern=true is treated as a CLI error upstream.
	HasPattern bool
	// IgnoreCase: case-insensitive pattern matching (`-C` in upstream).
	IgnoreCase bool
	// ForceStrand: reverse-complement '-' strand sequences before counting.
	ForceStrand bool
	// FullHeader: tolerate contig names that contain whitespace by matching
	// against the full FASTA header line. Off by default.
	FullHeader bool
}

// Counts holds the per-interval nucleotide counts.
type Counts struct {
	A, C, G, T, N, Other int
	PatternHits          int
	SeqLen               int
}

// Profile computes Counts for a single sequence (`-s` is applied by the
// caller before calling this). The pattern, if non-empty, is counted as
// fixed-substring occurrences with overlaps allowed (matching upstream).
func Profile(seq []byte, pattern string, ignoreCase bool) Counts {
	var c Counts
	c.SeqLen = len(seq)
	for _, b := range seq {
		switch b {
		case 'A', 'a':
			c.A++
		case 'C', 'c':
			c.C++
		case 'G', 'g':
			c.G++
		case 'T', 't', 'U', 'u':
			c.T++
		case 'N', 'n':
			c.N++
		default:
			c.Other++
		}
	}
	if pattern != "" {
		c.PatternHits = countPattern(seq, pattern, ignoreCase)
	}
	return c
}

// countPattern returns the number of (overlapping) occurrences of pattern
// in seq. Upstream's algorithm steps one character at a time and does a
// substring compare; we reproduce the same semantics, including matches
// at every position.
func countPattern(seq []byte, pattern string, ignoreCase bool) int {
	if pattern == "" || len(seq) < len(pattern) {
		return 0
	}
	hay := seq
	pat := []byte(pattern)
	if ignoreCase {
		hay = bytesToUpper(seq)
		pat = bytesToUpper([]byte(pattern))
	}
	count := 0
	plen := len(pat)
	last := len(hay) - plen
	for i := 0; i <= last; i++ {
		match := true
		for j := 0; j < plen; j++ {
			if hay[i+j] != pat[j] {
				match = false
				break
			}
		}
		if match {
			count++
		}
	}
	return count
}

func bytesToUpper(b []byte) []byte {
	out := make([]byte, len(b))
	for i, c := range b {
		if c >= 'a' && c <= 'z' {
			c -= 'a' - 'A'
		}
		out[i] = c
	}
	return out
}

// ReverseComplement returns a new slice that is the reverse complement of
// seq. Unknown characters pass through unchanged; IUPAC characters are
// complemented with case preservation. This mirrors the upstream `-s`
// behaviour for nuc (which uppercases at fetch time, so case isn't really
// observed, but we keep parity with -seq output).
func ReverseComplement(seq []byte) []byte {
	out := make([]byte, len(seq))
	for i, b := range seq {
		out[len(seq)-1-i] = complement(b)
	}
	return out
}

func complement(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'a':
		return 't'
	case 'C':
		return 'G'
	case 'c':
		return 'g'
	case 'G':
		return 'C'
	case 'g':
		return 'c'
	case 'T':
		return 'A'
	case 't':
		return 'a'
	case 'U':
		return 'A'
	case 'u':
		return 'a'
	case 'R':
		return 'Y'
	case 'r':
		return 'y'
	case 'Y':
		return 'R'
	case 'y':
		return 'r'
	case 'K':
		return 'M'
	case 'k':
		return 'm'
	case 'M':
		return 'K'
	case 'm':
		return 'k'
	case 'B':
		return 'V'
	case 'b':
		return 'v'
	case 'V':
		return 'B'
	case 'v':
		return 'b'
	case 'D':
		return 'H'
	case 'd':
		return 'h'
	case 'H':
		return 'D'
	case 'h':
		return 'd'
	case 'N':
		return 'N'
	case 'n':
		return 'n'
	case 'S':
		return 'S'
	case 's':
		return 's'
	case 'W':
		return 'W'
	case 'w':
		return 'w'
	}
	return b
}

// FormatHeader builds the `#`-prefixed column-header line for a record
// whose original column-count is bedType. Order:
//
//	1..bedType_usercol, then pct_at, pct_gc, num_A..num_oth, seq_len.
//
// Upstream appends `(bedType+10)_seq` when -seq is set and a pattern
// count column at the end when -pattern is set. We replicate that
// ordering verbatim.
func FormatHeader(bedType int, printSeq, hasPattern bool) string {
	var b strings.Builder
	b.WriteByte('#')
	for i := 1; i <= bedType; i++ {
		fmt.Fprintf(&b, "%d_usercol\t", i)
	}
	fmt.Fprintf(&b, "%d_pct_at\t", bedType+1)
	fmt.Fprintf(&b, "%d_pct_gc\t", bedType+2)
	fmt.Fprintf(&b, "%d_num_A\t", bedType+3)
	fmt.Fprintf(&b, "%d_num_C\t", bedType+4)
	fmt.Fprintf(&b, "%d_num_G\t", bedType+5)
	fmt.Fprintf(&b, "%d_num_T\t", bedType+6)
	fmt.Fprintf(&b, "%d_num_N\t", bedType+7)
	fmt.Fprintf(&b, "%d_num_oth\t", bedType+8)
	fmt.Fprintf(&b, "%d_seq_len", bedType+9)
	switch {
	case printSeq && hasPattern:
		fmt.Fprintf(&b, "\t%d_seq", bedType+10)
		fmt.Fprintf(&b, "\t%d_user_patt_count", bedType+11)
	case printSeq:
		fmt.Fprintf(&b, "\t%d_seq", bedType+10)
	case hasPattern:
		fmt.Fprintf(&b, "\t%d_user_patt_count", bedType+10)
	}
	b.WriteByte('\n')
	return b.String()
}

// FormatRow renders a result row given the original columns and the
// computed counts. printedSeq, when non-empty AND printSeq=true, is
// emitted as the seq column.
func FormatRow(cols []string, c Counts, printSeq bool, printedSeq []byte, hasPattern bool) string {
	var b strings.Builder
	b.WriteString(strings.Join(cols, "\t"))
	pctAT := 0.0
	pctGC := 0.0
	if c.SeqLen > 0 {
		pctAT = float64(c.A+c.T) / float64(c.SeqLen)
		pctGC = float64(c.C+c.G) / float64(c.SeqLen)
	}
	// Upstream uses `%f` which is 6-decimal-digit fixed format.
	fmt.Fprintf(&b, "\t%f\t%f", pctAT, pctGC)
	fmt.Fprintf(&b, "\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
		c.A, c.C, c.G, c.T, c.N, c.Other, c.SeqLen)
	if printSeq {
		b.WriteByte('\t')
		b.Write(printedSeq)
	}
	if hasPattern {
		fmt.Fprintf(&b, "\t%d", c.PatternHits)
	}
	b.WriteByte('\n')
	return b.String()
}

// Run streams BED records from bedR, looks each one up in the FASTA at
// fastaPath, and writes the profile to out. warn receives non-fatal
// warnings (unknown contig, zero-length feature, feature past chrom end).
// Returns the number of records emitted.
func Run(bedR io.Reader, fastaPath string, out io.Writer, warn io.Writer, opts Options) (int, error) {
	ra, err := fasta.OpenRandomAccess(fastaPath)
	if err != nil {
		return 0, err
	}
	defer ra.Close()

	// If FullHeader is requested, build a fallback map from raw-header
	// (everything after `>`) to the canonical first-token name. This lets
	// users pass `chr1 something` and still match a contig indexed as
	// `chr1` (upstream's behaviour pivots on whether the FASTA index used
	// the full header at indexing time; we deduplicate at lookup time so
	// the common case still works).
	var fullHeaderMap map[string]string
	if opts.FullHeader {
		fullHeaderMap = buildFullHeaderMap(fastaPath, ra)
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	sc := bufio.NewScanner(bedR)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	headerWritten := false
	written := 0
	lineNo := 0
	for sc.Scan() {
		lineNo++
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return written, fmt.Errorf("line %d: BED record needs >=3 columns: %q", lineNo, raw)
		}
		chrom := fields[0]
		start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
		if err != nil {
			return written, fmt.Errorf("line %d: invalid chromStart %q: %v", lineNo, fields[1], err)
		}
		end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
		if err != nil {
			return written, fmt.Errorf("line %d: invalid chromEnd %q: %v", lineNo, fields[2], err)
		}
		var strand string
		if len(fields) >= 6 {
			strand = fields[5]
		}

		if !headerWritten {
			fmt.Fprint(bw, FormatHeader(len(fields), opts.PrintSeq, opts.HasPattern))
			headerWritten = true
		}

		// Resolve contig (FullHeader fallback).
		resolved := chrom
		if opts.FullHeader && fullHeaderMap != nil {
			if alias, ok := fullHeaderMap[chrom]; ok {
				resolved = alias
			}
		}
		seqLength := ra.Length(resolved)
		if seqLength < 0 {
			if warn != nil {
				fmt.Fprintf(warn, "WARNING. chromosome (%s) was not found in the FASTA file. Skipping.\n", chrom)
			}
			continue
		}
		if end <= start {
			if warn != nil {
				fmt.Fprintf(warn, "Feature (%s:%d-%d) has length = 0, Skipping.\n", chrom, start+1, end-1)
			}
			continue
		}
		if int64(end) > seqLength {
			if warn != nil {
				fmt.Fprintf(warn, "Feature (%s:%d-%d) beyond the length of %s size (%d bp). Skipping.\n", chrom, start, end, chrom, seqLength)
			}
			continue
		}
		seq, err := ra.Fetch(resolved, int64(start), int64(end))
		if err != nil {
			return written, fmt.Errorf("line %d: %v", lineNo, err)
		}
		if opts.ForceStrand && strand == "-" {
			seq = ReverseComplement(seq)
		}
		counts := Profile(seq, opts.Pattern, opts.IgnoreCase)
		row := FormatRow(fields, counts, opts.PrintSeq, seq, opts.HasPattern)
		if _, err := bw.WriteString(row); err != nil {
			return written, err
		}
		written++
	}
	if err := sc.Err(); err != nil {
		return written, err
	}
	return written, nil
}

// buildFullHeaderMap walks the FASTA once, harvesting each `>...` header
// and mapping the full-header text (including spaces/tabs) to the
// first-token contig name that the index keys on. Returns nil on any I/O
// error — `-fullHeader` is best-effort.
func buildFullHeaderMap(path string, ra *fasta.RandomAccess) map[string]string {
	f, err := os.Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	known := make(map[string]bool)
	for _, n := range ra.Index().Names() {
		known[n] = true
	}
	out := make(map[string]string)
	br := bufio.NewScanner(f)
	br.Buffer(make([]byte, 64*1024), 16*1024*1024)
	for br.Scan() {
		line := br.Text()
		if len(line) == 0 || line[0] != '>' {
			continue
		}
		full := strings.TrimSpace(line[1:])
		first := full
		if i := strings.IndexAny(full, " \t"); i >= 0 {
			first = full[:i]
		}
		if !known[first] {
			continue
		}
		out[full] = first
	}
	return out
}
