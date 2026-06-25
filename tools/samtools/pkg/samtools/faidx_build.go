package samtools

// Index-build and FASTQ-quality support for faidx/fqidx.
//
// FaidxBuild writes the sibling .fai (and, for BGZF input, .gzi) for path,
// mirroring `samtools faidx <file>` / `samtools fqidx <file>` when called with
// no regions. The FASTA path reuses the shared pkg/htsgo/fasta index builder;
// the FASTQ path uses a dedicated 6-column builder that records the quality
// offset, matching htslib's fai_build_core.

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// FaidxBuild builds and writes the index sidecar(s) for path. For a plain
// FASTA/FASTQ it writes <path>.fai (or opts.FaiName when set). For a
// BGZF-compressed input it additionally writes <path>.gzi (or opts.GziName).
func FaidxBuild(path string, opts FaidxOptions) error {
	faiName := opts.FaiName
	if faiName == "" {
		faiName = path + ".fai"
	}

	bgzipped, err := isBGZFPath(path)
	if err != nil {
		return err
	}

	if opts.Format == FaidxFASTQ {
		return buildFastqIndex(path, faiName, opts, bgzipped)
	}

	// FASTA. Build the 5-column index over the (decompressed, if bgzipped)
	// payload using the shared builder so offsets land in the uncompressed
	// virtual stream — exactly what samtools faidx records. The bgzipped path
	// STREAMS the decompressed bytes through the index builder rather than
	// slurping the whole genome into memory: a genome-scale reference (e.g.
	// hs37d5, ~3 GB decompressed) would otherwise exhaust RAM, whereas upstream
	// builds the .fai in a few MB.
	var idx *fasta.Index
	if bgzipped {
		f, derr := os.Open(path)
		if derr != nil {
			return derr
		}
		br, derr := bgzf.NewReader(f)
		if derr != nil {
			f.Close()
			return derr
		}
		idx, err = fasta.BuildIndexReader(br)
		f.Close()
	} else {
		idx, err = fasta.BuildIndex(path)
	}
	if err != nil {
		return err
	}
	if err := idx.Save(faiName); err != nil {
		return err
	}
	if bgzipped {
		if err := writeGziSidecar(path, opts); err != nil {
			return err
		}
	}
	return nil
}

// writeGziSidecar scans the BGZF blocks of path and writes the .gzi block
// index (htslib's bgzf_index_dump format).
func writeGziSidecar(path string, opts FaidxOptions) error {
	gziName := opts.GziName
	if gziName == "" {
		gziName = path + ".gzi"
	}
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	offsets, err := bgzf.Scan(f)
	if err != nil {
		return err
	}
	out, err := os.Create(gziName)
	if err != nil {
		return err
	}
	defer out.Close()
	bw := bufio.NewWriter(out)
	if err := bgzf.WriteGZI(bw, offsets); err != nil {
		return err
	}
	return bw.Flush()
}

// openFastqReader returns a streaming reader over the (decompressed, if
// bgzipped) FASTQ payload of path, plus a closer for the underlying handle.
// It never buffers the whole payload: a BGZF input is inflated block by block
// through bgzf.NewReader, exactly like the FASTA build's BuildIndexReader path.
func openFastqReader(path string, bgzipped bool) (io.Reader, func() error, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	if !bgzipped {
		return f, f.Close, nil
	}
	br, err := bgzf.NewReader(f)
	if err != nil {
		f.Close()
		return nil, nil, err
	}
	return br, f.Close, nil
}

// fastqIndexEntry is one row of a 6-column FASTQ .fai (.fqi) file.
type fastqIndexEntry struct {
	name       string
	seqLen     int64
	seqOffset  int64
	lineBases  int64
	lineWidth  int64
	qualOffset int64
}

// buildFastqIndex builds the 6-column FASTQ index for path and writes it to
// faiName (plus a .gzi for BGZF input). It walks the decompressed payload with
// the same state machine semantics as htslib's fai_build_core.
func buildFastqIndex(path, faiName string, opts FaidxOptions, bgzipped bool) error {
	r, closeFn, err := openFastqReader(path, bgzipped)
	if err != nil {
		return err
	}
	defer closeFn()
	entries, err := scanFastqIndexReader(r)
	if err != nil {
		return err
	}
	out, err := os.Create(faiName)
	if err != nil {
		return err
	}
	defer out.Close()
	bw := bufio.NewWriter(out)
	for _, e := range entries {
		fmt.Fprintf(bw, "%s\t%d\t%d\t%d\t%d\t%d\n",
			e.name, e.seqLen, e.seqOffset, e.lineBases, e.lineWidth, e.qualOffset)
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	if bgzipped {
		if err := writeGziSidecar(path, opts); err != nil {
			return err
		}
	}
	return nil
}

// scanFastqIndex parses an in-memory FASTQ payload into index entries. It is a
// thin wrapper over scanFastqIndexReader retained for tests and in-memory
// callers; the build path streams instead (see buildFastqIndex).
func scanFastqIndex(data []byte) ([]fastqIndexEntry, error) {
	return scanFastqIndexReader(bytes.NewReader(data))
}

// scanFastqIndexReader parses a FASTQ payload from r into index entries with a
// single streaming pass, holding only O(1) state (one record's running counts)
// regardless of how large the FASTQ is. It records the byte offset (into the
// uncompressed stream) of each record's first sequence base and first quality
// base, the bases-per-line and the bytes-per-line, exactly like htslib's FASTQ
// indexer for uniform records — and produces byte-identical .fqi output to the
// previous slurp-then-scan implementation.
//
// A line-oriented scanner walks the stream; rather than materialising the whole
// payload, it tracks the byte offset of the start of each line so seqOffset and
// qualOffset land in the uncompressed stream just as before.
func scanFastqIndexReader(r io.Reader) ([]fastqIndexEntry, error) {
	var entries []fastqIndexEntry
	br := bufio.NewReader(r)

	var (
		offset  int64 // byte offset of the next line's first byte
		lineNum = 1
	)

	// readLine reads up to and including the next '\n' (or EOF). It returns the
	// line's bytes WITHOUT the trailing '\n' (line), the byte width INCLUDING
	// the consumed '\n' (width), the offset of the line's first byte (start),
	// and ok==false at end of input. offset is advanced past the line.
	readLine := func() (line []byte, width int64, start int64, ok bool, err error) {
		start = offset
		raw, rerr := br.ReadBytes('\n')
		if len(raw) == 0 && rerr != nil {
			if rerr == io.EOF {
				return nil, 0, start, false, nil
			}
			return nil, 0, start, false, rerr
		}
		width = int64(len(raw))
		offset += width
		// Strip the trailing '\n' (kept in width) for the returned line body.
		body := raw
		if len(body) > 0 && body[len(body)-1] == '\n' {
			body = body[:len(body)-1]
		}
		if rerr != nil && rerr != io.EOF {
			return body, width, start, true, rerr
		}
		return body, width, start, true, nil
	}

	for {
		// Peek the first byte of the next line to decide record boundary / blank.
		b, perr := br.ReadByte()
		if perr == io.EOF {
			break
		}
		if perr != nil {
			return nil, perr
		}
		if err := br.UnreadByte(); err != nil {
			return nil, err
		}
		// Skip blank lines between records.
		if b == '\n' {
			if _, _, _, _, err := readLine(); err != nil {
				return nil, err
			}
			lineNum++
			continue
		}
		if b != '@' {
			return nil, fmt.Errorf("fqidx: format error, unexpected %q at line %d", string(b), lineNum)
		}
		// Header line: name is the first whitespace-delimited token after '@'.
		hdr, _, _, _, err := readLine()
		if err != nil {
			return nil, err
		}
		lineNum++
		name := fastqName(hdr[1:])
		seqOffset := offset

		// Sequence lines until a line beginning with '+'.
		var seqLen, lineBases, lineWidth int64
		for {
			nb, nperr := br.ReadByte()
			if nperr == io.EOF {
				break
			}
			if nperr != nil {
				return nil, nperr
			}
			if err := br.UnreadByte(); err != nil {
				return nil, err
			}
			if nb == '+' {
				break
			}
			body, width, _, ok, lerr := readLine()
			if lerr != nil {
				return nil, lerr
			}
			if !ok {
				break
			}
			lineNum++
			bases := int64(len(body))
			if lineBases == 0 {
				lineBases = bases
				lineWidth = width
			}
			seqLen += bases
		}
		// '+' separator line.
		plus, perr := br.ReadByte()
		if perr != nil || plus != '+' {
			return nil, fmt.Errorf("fqidx: missing '+' for %q at line %d", name, lineNum)
		}
		if err := br.UnreadByte(); err != nil {
			return nil, err
		}
		if _, _, _, _, err := readLine(); err != nil { // consume the '+' line
			return nil, err
		}
		lineNum++
		qualOffset := offset

		// Quality lines until we have collected seqLen quality bytes.
		var qualLen int64
		for qualLen < seqLen {
			body, _, _, ok, lerr := readLine()
			if lerr != nil {
				return nil, lerr
			}
			if !ok {
				break
			}
			lineNum++
			qualLen += int64(len(body))
		}

		entries = append(entries, fastqIndexEntry{
			name:       name,
			seqLen:     seqLen,
			seqOffset:  seqOffset,
			lineBases:  lineBases,
			lineWidth:  lineWidth,
			qualOffset: qualOffset,
		})
	}
	return entries, nil
}

// fastqName extracts the first whitespace-delimited token from a FASTQ header
// (the bytes after '@', newline already stripped).
func fastqName(b []byte) string {
	b = bytes.TrimRight(b, "\r")
	for j := 0; j < len(b); j++ {
		if b[j] == ' ' || b[j] == '\t' {
			return string(b[:j])
		}
	}
	return string(b)
}

// faidxQualAccess provides quality-string fetches for FASTQ extraction. It is
// backed by the 6-column index (built on the fly when no .fqi sidecar exists)
// and reads the quality bytes from disk on demand via an io.ReaderAt over the
// (decompressed) payload — it never holds the whole FASTQ in memory, mirroring
// fasta.RandomAccess for the sequence path.
type faidxQualAccess struct {
	r       io.ReaderAt
	byName  map[string]fastqIndexEntry
	closeFn func() error
}

// openFaidxQual prepares quality access for a FASTQ input. It builds the 6
// column index by streaming the payload (bounded RSS), then opens a
// random-access reader over the (decompressed) payload so FetchQual can slice
// quality bytes directly from disk — without slurping the whole FASTQ.
func openFaidxQual(path string, opts FaidxOptions) (*faidxQualAccess, error) {
	bgzipped, err := isBGZFPath(path)
	if err != nil {
		return nil, err
	}
	// Build the index with a streaming pass (closes its own handle).
	r, closeFn, err := openFastqReader(path, bgzipped)
	if err != nil {
		return nil, err
	}
	entries, serr := scanFastqIndexReader(r)
	closeFn()
	if serr != nil {
		return nil, serr
	}
	byName := make(map[string]fastqIndexEntry, len(entries))
	for _, e := range entries {
		byName[e.name] = e
	}
	// Open a random-access reader over the decompressed payload for the actual
	// quality byte fetches (partial BGZF decompression when applicable).
	ra, raClose, err := fasta.OpenPayloadReaderAt(path)
	if err != nil {
		return nil, err
	}
	return &faidxQualAccess{r: ra, byName: byName, closeFn: raClose.Close}, nil
}

// FetchQual returns the quality bytes for [beg, end) (0-based half-open) on
// the named record, with embedded newlines stripped and case preserved. The
// quality bytes are read from the underlying payload via ReadAt — the same
// uniform line-geometry arithmetic fasta.RandomAccess.FetchRaw uses for
// sequence — so no in-memory copy of the FASTQ is required.
func (q *faidxQualAccess) FetchQual(name string, beg, end int64) ([]byte, error) {
	e, ok := q.byName[name]
	if !ok {
		return nil, fmt.Errorf("fqidx: record %q not found", name)
	}
	if beg < 0 || end < beg || end > e.seqLen {
		return nil, fmt.Errorf("fqidx: range %s:%d-%d out of bounds", name, beg, end)
	}
	if beg == end {
		return []byte{}, nil
	}
	lineBases := e.lineBases
	if lineBases <= 0 {
		lineBases = e.seqLen
	}
	startLine := beg / lineBases
	startCol := beg % lineBases
	endLine := (end - 1) / lineBases
	endCol := (end - 1) % lineBases
	startByte := e.qualOffset + startLine*e.lineWidth + startCol
	endByte := e.qualOffset + endLine*e.lineWidth + endCol + 1
	if startByte < 0 || endByte <= startByte {
		return nil, fmt.Errorf("fqidx: bad qual byte range %d-%d", startByte, endByte)
	}
	raw := make([]byte, endByte-startByte)
	if _, err := q.r.ReadAt(raw, startByte); err != nil && err != io.EOF {
		return nil, err
	}
	out := make([]byte, 0, end-beg)
	for _, b := range raw {
		if b == '\n' || b == '\r' {
			continue
		}
		out = append(out, b)
	}
	return out, nil
}

// Close releases any resources held by the quality accessor.
func (q *faidxQualAccess) Close() error {
	if q.closeFn != nil {
		return q.closeFn()
	}
	return nil
}
