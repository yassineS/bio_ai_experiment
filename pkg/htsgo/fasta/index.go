// Index implements a samtools-compatible FASTA index (.fai) reader.
//
// The .fai format is a tab-separated plaintext file with one record per
// reference sequence and the following columns:
//
//	NAME    contig name (first whitespace-delimited token after '>')
//	LENGTH  total number of bases in the sequence
//	OFFSET  byte offset in the FASTA where the sequence's first base appears
//	LINEBLEN bases per line (excluding the line terminator)
//	LINELEN bytes per line (including the line terminator)
//
// The index lets us extract any sub-range of any contig in O(1) seeks plus
// O(span) bytes read — which is exactly what `bcftools norm` needs to
// left-align indels against the reference without slurping the whole genome
// into memory.
//
// Only the v1 format (5 columns) is supported here; the v2 format used for
// bgzipped FASTAs (`.gzi`) is out of scope for this slice.
package fasta

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// IndexEntry is one row of a .fai file.
type IndexEntry struct {
	// Name is the contig identifier (the part of the FASTA header up to
	// the first whitespace).
	Name string
	// Length is the total number of bases in the sequence.
	Length int64
	// Offset is the byte offset of the first base of the sequence inside
	// the underlying FASTA file.
	Offset int64
	// LineBases is the number of bases per line (excluding the newline
	// terminator) — typically 60 or 80 for well-formed FASTAs.
	LineBases int64
	// LineWidth is the number of bytes per line (including the newline
	// terminator).
	LineWidth int64
}

// Index is an in-memory map from contig name to its .fai entry, preserving
// the original on-disk order in IndexEntry slots.
type Index struct {
	entries []IndexEntry
	byName  map[string]int // contig name → index into entries
}

// LoadIndex reads a .fai file from disk.
func LoadIndex(path string) (*Index, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ReadIndex(f)
}

// ReadIndex parses a .fai file from an io.Reader.
func ReadIndex(r io.Reader) (*Index, error) {
	idx := &Index{byName: make(map[string]int)}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	line := 0
	for sc.Scan() {
		line++
		text := sc.Text()
		if text == "" {
			continue
		}
		fields := strings.Split(text, "\t")
		if len(fields) < 5 {
			return nil, fmt.Errorf("fasta: bad .fai line %d: want 5 columns, got %d", line, len(fields))
		}
		length, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("fasta: bad length on .fai line %d: %w", line, err)
		}
		offset, err := strconv.ParseInt(fields[2], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("fasta: bad offset on .fai line %d: %w", line, err)
		}
		linebases, err := strconv.ParseInt(fields[3], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("fasta: bad lineblen on .fai line %d: %w", line, err)
		}
		linewidth, err := strconv.ParseInt(fields[4], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("fasta: bad linelen on .fai line %d: %w", line, err)
		}
		if linebases <= 0 || linewidth < linebases {
			return nil, fmt.Errorf("fasta: invalid line geometry on .fai line %d (lineblen=%d linelen=%d)", line, linebases, linewidth)
		}
		idx.byName[fields[0]] = len(idx.entries)
		idx.entries = append(idx.entries, IndexEntry{
			Name:      fields[0],
			Length:    length,
			Offset:    offset,
			LineBases: linebases,
			LineWidth: linewidth,
		})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return idx, nil
}

// Entries returns all index entries in on-disk order.
func (i *Index) Entries() []IndexEntry {
	return i.entries
}

// Get returns the entry for a contig, or false if absent.
func (i *Index) Get(name string) (IndexEntry, bool) {
	pos, ok := i.byName[name]
	if !ok {
		return IndexEntry{}, false
	}
	return i.entries[pos], true
}

// Names returns all contig names in on-disk order.
func (i *Index) Names() []string {
	out := make([]string, 0, len(i.entries))
	for _, e := range i.entries {
		out = append(out, e.Name)
	}
	return out
}

// BuildIndex scans a plain (non-bgzipped) FASTA file and produces the
// equivalent of `samtools faidx` output. The function assumes uniform line
// length within each contig (the same constraint samtools enforces) and
// returns an error otherwise.
//
// This is convenient for tests where we want to build the index without
// shelling out to samtools, and for users that haven't run faidx ahead of
// time. For real-world genome files prefer the ahead-of-time `.fai` next to
// the FASTA.
//
// Contig names are taken to be the first whitespace-delimited token after
// '>' — the same convention samtools faidx uses. For `bedtools getfasta
// -fullHeader` semantics (the whole header line, whitespace included),
// use BuildIndexFullHeader instead.
func BuildIndex(path string) (*Index, error) {
	return buildIndex(path, false)
}

// BuildIndexFullHeader scans a plain (non-bgzipped) FASTA and indexes
// contigs by the full header line (everything after '>' up to but not
// including the terminating newline). It implements the `bedtools
// getfasta -fullHeader` lookup mode: a BED row keyed by
// `chr1 assembled by consortium X` resolves to the sequence whose header
// is exactly `>chr1 assembled by consortium X`. Names emitted into the
// .fai will contain spaces and are not interoperable with samtools faidx
// — only use this when the caller is going to consume the index in-memory.
func BuildIndexFullHeader(path string) (*Index, error) {
	return buildIndex(path, true)
}

func buildIndex(path string, fullHeader bool) (*Index, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return buildIndexFromReader(f, fullHeader)
}

// buildIndexFromBytes is the in-memory counterpart of buildIndex, used
// by the BGZF path which has already decompressed the FASTA payload.
func buildIndexFromBytes(data []byte, fullHeader bool) (*Index, error) {
	return buildIndexFromReader(bytes.NewReader(data), fullHeader)
}

// BuildIndexBytes builds a samtools-faidx-compatible index over an
// in-memory FASTA payload (the bytes the file would contain after any
// outer compression layer has been stripped). It is the in-memory
// counterpart of BuildIndex and is used by the BGZF wrapper in
// `tools/bedgetfasta` (and any future callers) to avoid an extra round
// trip through the filesystem.
func BuildIndexBytes(data []byte) (*Index, error) {
	return buildIndexFromBytes(data, false)
}

// BuildIndexFullHeaderBytes is the `-fullHeader` analogue of
// BuildIndexBytes: contigs are keyed by the full header line
// (whitespace included) instead of the first-token name. Use this when
// implementing `bedtools getfasta -fullHeader` against an in-memory
// FASTA payload.
func BuildIndexFullHeaderBytes(data []byte) (*Index, error) {
	return buildIndexFromBytes(data, true)
}

// buildIndexFromReader is the streaming, source-agnostic core of
// BuildIndex / BuildIndexFullHeader. It walks the FASTA line by line,
// recording header offsets and verifying that each contig has uniform
// line geometry (mirrors samtools faidx's invariant).
func buildIndexFromReader(r io.Reader, fullHeader bool) (*Index, error) {
	idx := &Index{byName: make(map[string]int)}
	br := bufio.NewReader(r)
	var (
		offset    int64
		curName   string
		curLength int64
		curOffset int64
		curLineB  int64
		curLineW  int64
		seenAny   bool
		mixed     bool // a contig had non-uniform line widths
	)
	flush := func() {
		if !seenAny {
			return
		}
		idx.byName[curName] = len(idx.entries)
		idx.entries = append(idx.entries, IndexEntry{
			Name:      curName,
			Length:    curLength,
			Offset:    curOffset,
			LineBases: curLineB,
			LineWidth: curLineW,
		})
	}
	for {
		line, err := br.ReadBytes('\n')
		lineLen := int64(len(line))
		if len(line) > 0 {
			if line[0] == '>' {
				flush()
				// Parse the contig name: the first whitespace-delimited
				// token after '>' for default-mode samtools/bedtools, or
				// the full header line (excluding the trailing CR/LF) for
				// -fullHeader mode.
				name := string(line[1:])
				name = strings.TrimRight(name, "\r\n")
				if !fullHeader {
					if i := strings.IndexAny(name, " \t"); i >= 0 {
						name = name[:i]
					}
				}
				curName = name
				curLength = 0
				curLineB = 0
				curLineW = 0
				seenAny = true
				mixed = false
				offset += lineLen
				curOffset = offset
			} else {
				// sequence line
				body := line
				// strip trailing newline + optional CR
				for len(body) > 0 && (body[len(body)-1] == '\n' || body[len(body)-1] == '\r') {
					body = body[:len(body)-1]
				}
				bases := int64(len(body))
				if curLineB == 0 {
					curLineB = bases
					curLineW = lineLen
				} else if mixed {
					// We already saw a short line; any further sequence line is illegal.
					return nil, fmt.Errorf("fasta: contig %q has non-uniform line width (cannot build .fai)", curName)
				} else if bases != curLineB {
					if bases < curLineB {
						mixed = true // first short line, treat as final-line terminator
					} else {
						return nil, fmt.Errorf("fasta: contig %q has non-uniform line width (cannot build .fai)", curName)
					}
				}
				curLength += bases
				offset += lineLen
			}
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
	}
	flush()
	return idx, nil
}

// WriteTo serialises the index in the standard .fai layout.
func (i *Index) WriteTo(w io.Writer) (int64, error) {
	bw := bufio.NewWriter(w)
	var n int64
	for _, e := range i.entries {
		written, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%d\t%d\n",
			e.Name, e.Length, e.Offset, e.LineBases, e.LineWidth)
		n += int64(written)
		if err != nil {
			return n, err
		}
	}
	return n, bw.Flush()
}

// Save writes the index to path (typically `<fasta>.fai`).
func (i *Index) Save(path string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	if _, err := i.WriteTo(f); err != nil {
		return err
	}
	return nil
}

// RandomAccess pairs an open FASTA file with its index for cheap range
// extraction. Callers should Close the file when finished.
type RandomAccess struct {
	r     io.ReaderAt
	idx   *Index
	close func() error
}

// OpenRandomAccess opens a FASTA file together with its sibling .fai index
// (or builds one on the fly when missing). The returned RandomAccess holds
// onto an open *os.File; call Close to release it.
//
// BGZF-compressed FASTAs (`.fa.gz` produced by `bgzip` or `samtools faidx
// --output-fmt fasta`) are detected automatically by sniffing the BGZF
// magic — callers do NOT need a separate code path for `.gz` inputs. See
// OpenRandomAccessBGZF for the underlying implementation; this thin
// wrapper exists so most callers can use a single entry point.
func OpenRandomAccess(path string) (*RandomAccess, error) {
	if bgzf, err := isBGZFFile(path); err != nil {
		return nil, err
	} else if bgzf {
		return OpenRandomAccessBGZF(path)
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	idx, err := LoadIndex(path + ".fai")
	if err != nil {
		if !os.IsNotExist(err) {
			f.Close()
			return nil, err
		}
		// Build the index on the fly when no sibling .fai exists.
		idx, err = BuildIndex(path)
		if err != nil {
			f.Close()
			return nil, err
		}
	}
	return &RandomAccess{r: f, idx: idx, close: f.Close}, nil
}

// OpenRandomAccessFullHeader opens a FASTA in `-fullHeader` mode: contigs
// are keyed by the full header line (whitespace included) instead of the
// first-token name. The on-disk `.fai`, which is built by samtools faidx
// and contains only the first-token names, is intentionally ignored —
// we always rebuild from the FASTA header lines so the in-memory index
// keys match the full header. Use this when implementing `bedtools
// getfasta -fullHeader` (or any consumer that expects a BED row keyed by
// a multi-word contig identifier to match the corresponding FASTA
// sequence).
func OpenRandomAccessFullHeader(path string) (*RandomAccess, error) {
	if bgzf, err := isBGZFFile(path); err != nil {
		return nil, err
	} else if bgzf {
		return OpenRandomAccessBGZFFullHeader(path)
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	idx, err := BuildIndexFullHeader(path)
	if err != nil {
		f.Close()
		return nil, err
	}
	return &RandomAccess{r: f, idx: idx, close: f.Close}, nil
}

// NewRandomAccess wraps an arbitrary ReaderAt plus an in-memory index.
// Useful in tests where the FASTA contents live in a bytes.Reader.
func NewRandomAccess(r io.ReaderAt, idx *Index) *RandomAccess {
	return &RandomAccess{r: r, idx: idx, close: func() error { return nil }}
}

// newRandomAccessWithCloser wraps a ReaderAt + index and installs a custom
// Close hook (used by the .gzi partial-decompression backend, which owns an
// open file handle that must be released on Close).
func newRandomAccessWithCloser(r io.ReaderAt, idx *Index, closeFn func() error) *RandomAccess {
	return &RandomAccess{r: r, idx: idx, close: closeFn}
}

// Close releases the underlying file (if any).
func (ra *RandomAccess) Close() error {
	if ra.close != nil {
		return ra.close()
	}
	return nil
}

// Index returns the parsed .fai index.
func (ra *RandomAccess) Index() *Index { return ra.idx }

// Fetch returns the sequence bases for the half-open range [start, end) on
// contig name. Coordinates are 0-based. Bases are returned as uppercase
// ASCII to keep downstream comparisons case-insensitive — left-alignment
// always works on a canonical case.
func (ra *RandomAccess) Fetch(name string, start, end int64) ([]byte, error) {
	pos, ok := ra.idx.byName[name]
	if !ok {
		return nil, fmt.Errorf("fasta: contig %q not in index", name)
	}
	entry := ra.idx.entries[pos]
	if start < 0 || end < start || end > entry.Length {
		return nil, fmt.Errorf("fasta: range %s:%d-%d out of bounds (length %d)", name, start, end, entry.Length)
	}
	if start == end {
		return []byte{}, nil
	}

	// Translate base coordinates into file byte coordinates using the
	// uniform line geometry.
	startLine := start / entry.LineBases
	startCol := start % entry.LineBases
	endLine := (end - 1) / entry.LineBases
	endCol := (end - 1) % entry.LineBases

	startByte := entry.Offset + startLine*entry.LineWidth + startCol
	endByte := entry.Offset + endLine*entry.LineWidth + endCol + 1
	if endByte <= startByte {
		return nil, fmt.Errorf("fasta: bad byte range %d-%d", startByte, endByte)
	}

	buf := make([]byte, endByte-startByte)
	if _, err := ra.r.ReadAt(buf, startByte); err != nil && err != io.EOF {
		return nil, err
	}
	// Strip embedded newlines (and CR) from the raw file slice.
	out := buf[:0]
	for _, b := range buf {
		if b == '\n' || b == '\r' {
			continue
		}
		if b >= 'a' && b <= 'z' {
			b -= 'a' - 'A'
		}
		out = append(out, b)
	}
	if int64(len(out)) != end-start {
		return nil, fmt.Errorf("fasta: requested %d bases but read %d on %s", end-start, len(out), name)
	}
	return out, nil
}

// Length returns the contig length, or -1 if the contig is unknown.
func (ra *RandomAccess) Length(name string) int64 {
	pos, ok := ra.idx.byName[name]
	if !ok {
		return -1
	}
	return ra.idx.entries[pos].Length
}
