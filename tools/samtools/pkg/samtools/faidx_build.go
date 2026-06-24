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

// decompressFile returns the fully decompressed payload of a BGZF file.
func decompressFile(path string) ([]byte, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	br, err := bgzf.NewReader(f)
	if err != nil {
		return nil, err
	}
	return io.ReadAll(br)
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
	var data []byte
	var err error
	if bgzipped {
		data, err = decompressFile(path)
	} else {
		data, err = os.ReadFile(path)
	}
	if err != nil {
		return err
	}
	entries, err := scanFastqIndex(data)
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

// scanFastqIndex parses an in-memory FASTQ payload into index entries. It
// records the byte offset (into the uncompressed stream) of each record's
// first sequence base and first quality base, the bases-per-line and the
// bytes-per-line, exactly like htslib's FASTQ indexer for uniform records.
func scanFastqIndex(data []byte) ([]fastqIndexEntry, error) {
	var entries []fastqIndexEntry
	i := 0
	n := len(data)
	lineNum := 1

	readLine := func() (start, end, next int) {
		start = i
		for i < n && data[i] != '\n' {
			i++
		}
		end = i // exclusive of '\n'
		if i < n {
			i++ // consume '\n'
		}
		return start, end, i
	}

	for i < n {
		// Skip blank lines between records.
		if data[i] == '\n' {
			i++
			lineNum++
			continue
		}
		if data[i] != '@' {
			return nil, fmt.Errorf("fqidx: format error, unexpected %q at line %d", string(data[i]), lineNum)
		}
		// Header line: name is the first whitespace-delimited token after '@'.
		hs, he, _ := readLine()
		lineNum++
		name := fastqName(data[hs+1 : he])
		seqOffset := int64(i)

		// Sequence lines until a line beginning with '+'.
		var seqLen, lineBases, lineWidth int64
		for i < n && data[i] != '+' {
			ls, le, nx := readLine()
			lineNum++
			bases := int64(le - ls)
			width := int64(nx - ls)
			if lineBases == 0 {
				lineBases = bases
				lineWidth = width
			}
			seqLen += bases
		}
		// '+' separator line.
		if i >= n || data[i] != '+' {
			return nil, fmt.Errorf("fqidx: missing '+' for %q at line %d", name, lineNum)
		}
		readLine() // consume the '+' line
		lineNum++
		qualOffset := int64(i)

		// Quality lines until we have collected seqLen quality bytes.
		var qualLen int64
		for i < n && qualLen < seqLen {
			ls, le, _ := readLine()
			lineNum++
			qualLen += int64(le - ls)
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
// and reads the quality bytes directly from the (decompressed) payload.
type faidxQualAccess struct {
	data    []byte
	byName  map[string]fastqIndexEntry
	closeFn func() error
}

// openFaidxQual prepares quality access for a FASTQ input. It loads the 6
// column index from the sibling .fai when present (else builds it), and holds
// the decompressed payload in memory for direct quality slicing.
func openFaidxQual(path string, opts FaidxOptions) (*faidxQualAccess, error) {
	bgzipped, err := isBGZFPath(path)
	if err != nil {
		return nil, err
	}
	var data []byte
	if bgzipped {
		data, err = decompressFile(path)
	} else {
		data, err = os.ReadFile(path)
	}
	if err != nil {
		return nil, err
	}
	entries, err := scanFastqIndex(data)
	if err != nil {
		return nil, err
	}
	byName := make(map[string]fastqIndexEntry, len(entries))
	for _, e := range entries {
		byName[e.name] = e
	}
	return &faidxQualAccess{data: data, byName: byName, closeFn: func() error { return nil }}, nil
}

// FetchQual returns the quality bytes for [beg, end) (0-based half-open) on
// the named record, with embedded newlines stripped and case preserved.
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
	if startByte < 0 || endByte > int64(len(q.data)) || endByte <= startByte {
		return nil, fmt.Errorf("fqidx: bad qual byte range %d-%d", startByte, endByte)
	}
	raw := q.data[startByte:endByte]
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
