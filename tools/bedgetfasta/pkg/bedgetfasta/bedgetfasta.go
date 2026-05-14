// Package bedgetfasta implements `bedtools getfasta`: it pulls FASTA
// subsequences for each interval in a BED file. The tool reads BED records
// from an io.Reader and looks up the sequence on disk via a FAI-indexed
// FASTA (the .fai is built on the fly when missing — same convenience as
// upstream).
//
// Supported options (mirror the upstream flags of the same name):
//
//   - Name   — use the BED name column as the FASTA header. Default
//     header is `chrom:start-end`; with Name the header becomes
//     `<name>::<chrom>:<start>-<end>` (matches upstream's modern
//     `-name` / `-name+` rendering).
//   - NameOnly — like Name but emits only `<name>` (no `::chrom:start-end`).
//   - Tab    — emit a TSV ("name<TAB>seq") instead of FASTA.
//   - Strand — when set and the BED record has a `-` strand, the emitted
//     sequence is reverse-complemented (case-preserving IUPAC,
//     same as upstream).
//   - Split  — if the BED record is BED12 with block columns, emit the
//     concatenation of its blocks (per-block stranded when both
//     Strand and Split are on, mirroring upstream).
//   - RNA    — emit the sequence with `T → U` (uppercase) / `t → u`
//     (lowercase). For -s on reverse strand, the complement is
//     computed first and then T→U applied.
//
// Implementation note: the package implements its own case-preserving
// Fetch on top of the FASTA index because the shared
// pkg/bioformats/fasta.RandomAccess.Fetch uppercases for downstream
// case-insensitive comparison; getfasta needs the original case.
package bedgetfasta

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
	Name     bool // -name / -name+ — header is `<name>::chrom:start-end`
	NameOnly bool // -nameOnly      — header is just `<name>`
	Tab      bool // -tab           — TSV output
	Strand   bool // -s             — reverse-complement '-' strand intervals
	Split    bool // -split         — concatenate BED12 blocks
	RNA      bool // -rna           — emit U instead of T (case-preserved)
}

// Run reads BED records from bedR, looks up sequence from the indexed FASTA
// at fastaPath, and writes the result to out. warn receives non-fatal
// warnings (e.g. unknown chrom). Returns the number of records emitted.
func Run(bedR io.Reader, fastaPath string, out io.Writer, warn io.Writer, opts Options) (int, error) {
	ra, err := openFasta(fastaPath)
	if err != nil {
		return 0, err
	}
	defer ra.Close()

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	sc := bufio.NewScanner(bedR)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
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
		var name string
		if len(fields) >= 4 {
			name = fields[3]
		}
		var strand string
		if len(fields) >= 6 {
			strand = fields[5]
		}

		// Fetch the sequence(s).
		var seq []byte
		if opts.Split && len(fields) >= 12 {
			blocks, err := parseBlocks(fields, start)
			if err != nil {
				return written, fmt.Errorf("line %d: %v", lineNo, err)
			}
			parts := make([][]byte, 0, len(blocks))
			for _, b := range blocks {
				s, err := ra.FetchPreserveCase(chrom, int64(b[0]), int64(b[1]))
				if err != nil {
					if isMissingChrom(err) {
						warnMissingChrom(warn, chrom)
						s = nil
						break
					}
					return written, fmt.Errorf("line %d: %v", lineNo, err)
				}
				if opts.Strand && strand == "-" {
					s = reverseComplement(s)
				}
				parts = append(parts, s)
			}
			if len(parts) == 0 && !isAllPresent(parts, blocks) {
				// chrom missing; skip
				continue
			}
			if opts.Strand && strand == "-" {
				// Upstream emits blocks in reverse order on negative strand
				// so the concatenation reads 5'->3' on the transcript.
				for i, j := 0, len(parts)-1; i < j; i, j = i+1, j-1 {
					parts[i], parts[j] = parts[j], parts[i]
				}
			}
			seq = bytesJoin(parts)
		} else {
			s, err := ra.FetchPreserveCase(chrom, int64(start), int64(end))
			if err != nil {
				if isMissingChrom(err) {
					warnMissingChrom(warn, chrom)
					continue
				}
				return written, fmt.Errorf("line %d: %v", lineNo, err)
			}
			if opts.Strand && strand == "-" {
				s = reverseComplement(s)
			}
			seq = s
		}

		if opts.RNA {
			seq = dnaToRNA(seq)
		}

		header := formatHeader(chrom, start, end, name, strand, opts)
		if opts.Tab {
			if _, err := bw.WriteString(header); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\t'); err != nil {
				return written, err
			}
			if _, err := bw.Write(seq); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
		} else {
			if _, err := bw.WriteString(">"); err != nil {
				return written, err
			}
			if _, err := bw.WriteString(header); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
			if _, err := bw.Write(seq); err != nil {
				return written, err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return written, err
			}
		}
		written++
	}
	if err := sc.Err(); err != nil {
		return written, err
	}
	return written, nil
}

// parseBlocks parses BED12 block columns 10-12 (blockCount, blockSizes,
// blockStarts) and returns absolute [start,end) ranges per block.
func parseBlocks(fields []string, recordStart int) ([][2]int, error) {
	count, err := strconv.Atoi(strings.TrimSpace(fields[9]))
	if err != nil {
		return nil, fmt.Errorf("invalid blockCount %q: %v", fields[9], err)
	}
	sizes := strings.Split(strings.TrimRight(strings.TrimSpace(fields[10]), ","), ",")
	starts := strings.Split(strings.TrimRight(strings.TrimSpace(fields[11]), ","), ",")
	if len(sizes) != count || len(starts) != count {
		return nil, fmt.Errorf("blockCount=%d but %d sizes / %d starts", count, len(sizes), len(starts))
	}
	out := make([][2]int, 0, count)
	for i := 0; i < count; i++ {
		sz, err := strconv.Atoi(strings.TrimSpace(sizes[i]))
		if err != nil {
			return nil, fmt.Errorf("invalid block size %q: %v", sizes[i], err)
		}
		st, err := strconv.Atoi(strings.TrimSpace(starts[i]))
		if err != nil {
			return nil, fmt.Errorf("invalid block start %q: %v", starts[i], err)
		}
		out = append(out, [2]int{recordStart + st, recordStart + st + sz})
	}
	return out, nil
}

// formatHeader builds the header string per the option-set, mirroring the
// upstream rendering rules:
//
//   - default          ->  `chrom:start-end` (+ `(+/-)` when Strand)
//   - Name / Name+     ->  `name::chrom:start-end` (+ `(+/-)` when Strand)
//   - NameOnly         ->  `name` (+ `(+/-)` when Strand)
//   - Tab + Strand     ->  same render, the `(+/-)` segment is still included
func formatHeader(chrom string, start, end int, name, strand string, opts Options) string {
	coord := fmt.Sprintf("%s:%d-%d", chrom, start, end)
	strandSuffix := ""
	if opts.Strand {
		s := strand
		if s != "-" {
			s = "+"
		}
		strandSuffix = "(" + s + ")"
	}
	switch {
	case opts.NameOnly:
		if name == "" {
			return coord + strandSuffix
		}
		return name + strandSuffix
	case opts.Name:
		if name == "" {
			return coord + strandSuffix
		}
		return name + "::" + coord + strandSuffix
	default:
		return coord + strandSuffix
	}
}

// FetchPreserveCase is an alternate, case-preserving fetch on a FASTA
// random-access handle. Reused from pkg/bioformats/fasta's geometry but
// emits the raw bytes (other than line terminators) without uppercasing.
type RandomAccess struct {
	idx     *fasta.Index
	r       readerAtCloser
	closeFn func() error
}

type readerAtCloser interface {
	io.ReaderAt
}

func openFasta(path string) (*RandomAccess, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	idx, err := fasta.LoadIndex(path + ".fai")
	if err != nil {
		if !os.IsNotExist(err) {
			f.Close()
			return nil, err
		}
		// Build on the fly.
		idx, err = fasta.BuildIndex(path)
		if err != nil {
			f.Close()
			return nil, err
		}
	}
	return &RandomAccess{idx: idx, r: f, closeFn: f.Close}, nil
}

// Close releases the underlying file.
func (ra *RandomAccess) Close() error {
	if ra.closeFn != nil {
		return ra.closeFn()
	}
	return nil
}

// FetchPreserveCase reads [start,end) on contig name, returning the raw
// bases (preserving case) with embedded newlines stripped.
func (ra *RandomAccess) FetchPreserveCase(name string, start, end int64) ([]byte, error) {
	entry, ok := ra.idx.Get(name)
	if !ok {
		return nil, errMissingChrom{name: name}
	}
	if start < 0 || end < start || end > entry.Length {
		return nil, fmt.Errorf("range %s:%d-%d out of bounds (length %d)", name, start, end, entry.Length)
	}
	if start == end {
		return []byte{}, nil
	}
	startLine := start / entry.LineBases
	startCol := start % entry.LineBases
	endLine := (end - 1) / entry.LineBases
	endCol := (end - 1) % entry.LineBases
	startByte := entry.Offset + startLine*entry.LineWidth + startCol
	endByte := entry.Offset + endLine*entry.LineWidth + endCol + 1
	if endByte <= startByte {
		return nil, fmt.Errorf("bad byte range %d-%d", startByte, endByte)
	}
	buf := make([]byte, endByte-startByte)
	if _, err := ra.r.ReadAt(buf, startByte); err != nil && err != io.EOF {
		return nil, err
	}
	out := buf[:0]
	for _, b := range buf {
		if b == '\n' || b == '\r' {
			continue
		}
		out = append(out, b)
	}
	if int64(len(out)) != end-start {
		return nil, fmt.Errorf("requested %d bases but read %d on %s", end-start, len(out), name)
	}
	return out, nil
}

// errMissingChrom signals an unknown contig in a typed way so the loop can
// emit the upstream-style warning and continue rather than aborting.
type errMissingChrom struct{ name string }

func (e errMissingChrom) Error() string {
	return fmt.Sprintf("chromosome %q not in index", e.name)
}

func isMissingChrom(err error) bool {
	_, ok := err.(errMissingChrom)
	return ok
}

func warnMissingChrom(w io.Writer, name string) {
	if w == nil {
		return
	}
	fmt.Fprintf(w, "WARNING. chromosome (%s) was not found in the FASTA file. Skipping.\n", name)
}

// isAllPresent returns true when parts has the same length as blocks. (Used
// to detect a partially-filled set of blocks after a missing-chrom warning.)
func isAllPresent(parts [][]byte, blocks [][2]int) bool {
	return len(parts) == len(blocks)
}

// bytesJoin concatenates parts without a separator.
func bytesJoin(parts [][]byte) []byte {
	total := 0
	for _, p := range parts {
		total += len(p)
	}
	out := make([]byte, 0, total)
	for _, p := range parts {
		out = append(out, p...)
	}
	return out
}

// dnaToRNA replaces T with U and t with u in-place; other bases are
// unaffected. Used by -rna.
func dnaToRNA(s []byte) []byte {
	for i, b := range s {
		switch b {
		case 'T':
			s[i] = 'U'
		case 't':
			s[i] = 'u'
		}
	}
	return s
}

// reverseComplement reverses s and complements each base, preserving case
// for IUPAC characters. Unknown characters map to themselves.
func reverseComplement(s []byte) []byte {
	out := make([]byte, len(s))
	for i, b := range s {
		out[len(s)-1-i] = complement(b)
	}
	return out
}

// complement maps a single base to its IUPAC complement, preserving case.
// Mapping (uppercase): A<->T, C<->G, U->A, R<->Y, S<->S, W<->W, K<->M,
// B<->V, D<->H, N->N, X->X. Anything else passes through.
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
	case 'S':
		return 'S'
	case 's':
		return 's'
	case 'W':
		return 'W'
	case 'w':
		return 'w'
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
	case 'X':
		return 'X'
	case 'x':
		return 'x'
	}
	return b
}
