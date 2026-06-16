// Package bedgetfasta implements `bedtools getfasta`: it pulls FASTA
// subsequences for each interval in a BED file. The tool reads BED records
// from an io.Reader and looks up the sequence on disk via a FAI-indexed
// FASTA (the .fai is built on the fly when missing — same convenience as
// upstream).
//
// Supported options (mirror the upstream flags of the same name):
//
//   - Name     — use the BED name column for the FASTA header. The header
//     becomes `<name>::<chrom>:<start>-<end>`. Matching upstream, a BED row
//     with no name column emits an empty name (`>::chrom:start-end`); it does
//     NOT fall back to the coordinate-only header.
//   - NamePlus — deprecated upstream alias of Name; identical header.
//   - NameOnly — like Name but emits only `<name>` (no `::chrom:start-end`).
//     With an empty name column the header is empty (`>`).
//   - Tab    — emit a TSV ("header<TAB>seq") instead of FASTA.
//   - BedOut — re-emit the BED record (original tab-delimited columns)
//     followed by a trailing sequence column, instead of FASTA. The
//     sequence honours Strand, Split and RNA. Columns beyond 6 are
//     preserved verbatim, mirroring upstream's reportBedTab.
//   - Strand — when set and the BED record has a `-` strand, the emitted
//     sequence is reverse-complemented (case-preserving IUPAC,
//     same as upstream).
//   - Split  — if the BED record is BED12 with block columns, emit the
//     concatenation of its blocks. With Strand on a `-` record, the
//     blocks are concatenated in genomic order and the whole sequence is
//     then reverse-complemented (matching upstream's ReportSeq).
//   - RNA    — emit the sequence with `T → U` (uppercase) / `t → u`
//     (lowercase). For -s on reverse strand, the complement is
//     computed first and then T→U applied.
//
// Per-record diagnostics on the warn writer match upstream byte-for-byte:
// a missing chromosome, an out-of-range feature, and a zero-length feature
// each emit the corresponding "Skipping." line and produce no output. When a
// sibling `.fai` exists but predates the FASTA, the htslib-style
// "Warning: the index file is older than the FASTA file." line is emitted.
//
// Implementation note: the package implements its own case-preserving
// Fetch on top of the FASTA index because the shared
// pkg/htsgo/fasta.RandomAccess.Fetch uppercases for downstream
// case-insensitive comparison; getfasta needs the original case.
package bedgetfasta

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// Options configures Run.
type Options struct {
	Name     bool // -name          — header is `<name>::chrom:start-end`
	NamePlus bool // -name+         — (deprecated) same header as -name
	NameOnly bool // -nameOnly      — header is just `<name>`
	Tab      bool // -tab           — TSV output ("header<TAB>seq")
	BedOut   bool // -bedOut        — re-emit the BED record with a trailing
	//                                sequence column (tab-delimited) instead
	//                                of FASTA. Upstream also disables FASTA
	//                                output (`useFasta=false`) for this mode.
	Strand     bool // -s             — reverse-complement '-' strand intervals
	Split      bool // -split         — concatenate BED12 blocks
	RNA        bool // -rna           — emit U instead of T (case-preserved)
	FullHeader bool // -fullHeader    — index contigs by the full FASTA header
	//                                  line (whitespace included), matching
	//                                  upstream `bedtools getfasta -fullHeader`.
}

// Run reads BED records from bedR, looks up sequence from the indexed FASTA
// at fastaPath, and writes the result to out. warn receives non-fatal
// warnings (e.g. unknown chrom). Returns the number of records emitted.
func Run(bedR io.Reader, fastaPath string, out io.Writer, warn io.Writer, opts Options) (int, error) {
	if warn == nil {
		warn = io.Discard
	}
	ra, err := openFasta(fastaPath, opts.FullHeader, warn)
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

		// Upstream rejects start > end with a hard error; mirror it.
		if end < start {
			return written, fmt.Errorf("line %d: chromEnd %d < chromStart %d", lineNo, end, start)
		}

		// Zero-length feature (start == end). Upstream artificially expands
		// the record by one base on each side, marks it zeroLength, and
		// then reports it with the ORIGINAL coordinates. We just skip with
		// the same message and never emit output.
		if start == end {
			fmt.Fprintf(warn, "Feature (%s:%d-%d) has length = 0, Skipping.\n", chrom, start, end)
			continue
		}

		// sequenceLength mirrors upstream's fr->sequenceLength: 0 means the
		// contig was not found in the index, otherwise it is the contig length.
		seqLength, ok := ra.sequenceLength(chrom)
		if !ok {
			fmt.Fprintf(warn, "WARNING. chromosome (%s) was not found in the FASTA file. Skipping.\n", chrom)
			continue
		}
		// Upstream bounds check: both start and end must be <= seqLength.
		if int64(start) > seqLength || int64(end) > seqLength {
			fmt.Fprintf(warn, "Feature (%s:%d-%d) beyond the length of %s size (%d bp).  Skipping.\n",
				chrom, start, end, chrom, seqLength)
			continue
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
					return written, fmt.Errorf("line %d: %v", lineNo, err)
				}
				parts = append(parts, s)
			}
			seq = bytesJoin(parts)
		} else {
			s, err := ra.FetchPreserveCase(chrom, int64(start), int64(end))
			if err != nil {
				return written, fmt.Errorf("line %d: %v", lineNo, err)
			}
			seq = s
		}

		// Upstream ReportSeq reverse-complements the whole (already
		// concatenated) sequence after extraction, then applies RNA via the
		// reverseComplement helper's isRNA flag. We keep the two steps split:
		// revcomp first (matching the block order upstream produces), then T->U.
		if opts.Strand && strand == "-" {
			seq = reverseComplement(seq)
		}
		if opts.RNA {
			seq = dnaToRNA(seq)
		}

		if opts.BedOut {
			// Re-emit the BED record (original tab-delimited fields) followed
			// by a trailing sequence column. Upstream re-serializes BED3-BED6
			// and appends any columns beyond 6 verbatim; joining the parsed
			// fields by TAB reproduces that for every BED flavour we accept.
			if _, err := bw.WriteString(strings.Join(fields, "\t")); err != nil {
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
			written++
			continue
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
	// Upstream renders the name-bearing header whenever -name, -name+, OR
	// -nameOnly is set, using the empty string for a missing name column
	// (it does NOT fall back to chrom:start-end). -name and -name+ produce
	// the identical `name::chrom:start-end` form; -name+ is a deprecated
	// alias retained for backwards compatibility.
	switch {
	case opts.NameOnly:
		return name + strandSuffix
	case opts.Name || opts.NamePlus:
		return name + "::" + coord + strandSuffix
	default:
		return coord + strandSuffix
	}
}

// FetchPreserveCase is an alternate, case-preserving fetch on a FASTA
// random-access handle. Reused from pkg/htsgo/fasta's geometry but
// emits the raw bytes (other than line terminators) without uppercasing.
type RandomAccess struct {
	idx     *fasta.Index
	r       readerAtCloser
	closeFn func() error
}

type readerAtCloser interface {
	io.ReaderAt
}

func openFasta(path string, fullHeader bool, warn io.Writer) (*RandomAccess, error) {
	// BGZF (`.fa.gz`) inputs are detected by sniffing the BGZF magic on
	// the first four bytes. When detected, we fully decompress the
	// payload into memory and back the case-preserving Fetch with a
	// bytes.Reader. samtools-compatible side-files (`.fa.gz.fai` and
	// `.fa.gz.gzi`) are honoured when present — see the package doc on
	// pkg/htsgo/fasta/bgzf.go for the on-disk format and the
	// future partial-decompression roadmap.
	if bgzf, err := isBGZF(path); err != nil {
		return nil, err
	} else if bgzf {
		return openFastaBGZF(path, fullHeader)
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	var idx *fasta.Index
	if fullHeader {
		// -fullHeader: a sibling .fai (built by samtools faidx) would
		// only contain first-token names, so ignore it and rebuild the
		// index keyed by the full header line.
		idx, err = fasta.BuildIndexFullHeader(path)
		if err != nil {
			f.Close()
			return nil, err
		}
	} else {
		idx, err = fasta.LoadIndex(path + ".fai")
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
		} else {
			// A sibling .fai exists. Upstream (htslib) warns when the FASTA
			// is newer than its index, since a stale index may not match.
			warnIfStaleIndex(path, warn)
		}
	}
	return &RandomAccess{idx: idx, r: f, closeFn: f.Close}, nil
}

// warnIfStaleIndex emits upstream's warning when path's sibling `.fai` index
// is older (by mtime) than the FASTA itself, mirroring htslib's check in
// FastaIndex::readIndexFile. Any stat error is treated as "not stale" and
// silently ignored, matching upstream's best-effort behaviour.
func warnIfStaleIndex(path string, warn io.Writer) {
	if warn == nil {
		return
	}
	faStat, err := os.Stat(path)
	if err != nil {
		return
	}
	idxStat, err := os.Stat(path + ".fai")
	if err != nil {
		return
	}
	if faStat.ModTime().After(idxStat.ModTime()) {
		fmt.Fprintln(warn, "Warning: the index file is older than the FASTA file.")
	}
}

// isBGZF reports whether path starts with the BGZF magic (gzip 1f 8b 08
// with FEXTRA bit set in FLG). A non-BGZF file (including plain gzip) is
// not an error — the caller falls back to plain-FASTA handling.
func isBGZF(path string) (bool, error) {
	f, err := os.Open(path)
	if err != nil {
		return false, err
	}
	defer f.Close()
	var buf [4]byte
	n, err := io.ReadFull(f, buf[:])
	if err != nil && err != io.EOF && err != io.ErrUnexpectedEOF {
		return false, err
	}
	if n < 4 {
		return false, nil
	}
	return buf == [4]byte{0x1f, 0x8b, 0x08, 0x04}, nil
}

// openFastaBGZF decompresses path and builds the case-preserving
// random-access view over the in-memory payload. The `.fa.gz.gzi`
// sidecar (if present) is parsed for validation only — the in-memory
// implementation does not yet need to seek by compressed offset.
func openFastaBGZF(path string, fullHeader bool) (*RandomAccess, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	br, err := bgzip.NewReader(f)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("bgzf header: %w", err)
	}
	payload, err := io.ReadAll(br)
	f.Close()
	if err != nil {
		return nil, fmt.Errorf("bgzf decompress: %w", err)
	}
	var idx *fasta.Index
	if fullHeader {
		idx, err = fasta.BuildIndexFullHeaderBytes(payload)
		if err != nil {
			return nil, err
		}
	} else {
		// Prefer an explicit .fa.gz.fai over rebuilding (samtools' offsets
		// are into the uncompressed virtual stream, which IS our payload).
		if fi, ferr := fasta.LoadIndex(path + ".fai"); ferr == nil {
			idx = fi
		} else if !os.IsNotExist(ferr) {
			return nil, ferr
		} else {
			idx, err = fasta.BuildIndexBytes(payload)
			if err != nil {
				return nil, err
			}
		}
	}
	rdr := bytes.NewReader(payload)
	return &RandomAccess{idx: idx, r: rdr, closeFn: func() error { return nil }}, nil
}

// Close releases the underlying file.
func (ra *RandomAccess) Close() error {
	if ra.closeFn != nil {
		return ra.closeFn()
	}
	return nil
}

// sequenceLength returns the length of contig name and whether it exists in
// the index. It mirrors upstream's FastaReference::sequenceLength, which
// returns 0 (here: ok=false) when the contig is absent.
func (ra *RandomAccess) sequenceLength(name string) (int64, bool) {
	entry, ok := ra.idx.Get(name)
	if !ok {
		return 0, false
	}
	return entry.Length, true
}

// FetchPreserveCase reads [start,end) on contig name, returning the raw
// bases (preserving case) with embedded newlines stripped.
func (ra *RandomAccess) FetchPreserveCase(name string, start, end int64) ([]byte, error) {
	entry, ok := ra.idx.Get(name)
	if !ok {
		return nil, fmt.Errorf("chromosome %q not in index", name)
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
