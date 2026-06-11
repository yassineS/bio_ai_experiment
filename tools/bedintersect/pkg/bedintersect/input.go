// Input-format handling for bedintersect. Upstream `bedtools intersect`
// accepts BED, VCF, GFF, and BAM input on -a and -b. This file converts each
// of those into a common interval record that preserves the original input
// columns (so they can be echoed verbatim) and knows how to re-render its
// coordinates when clipped to an overlap, mirroring upstream's per-record-type
// print methods.
package bedintersect

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// recordFormat tags the detected input format of a single file.
type recordFormat int

const (
	fmtBED recordFormat = iota
	fmtVCF
	fmtGFF
	fmtBAM
)

// inRecord is the common interval model used by every intersect output path.
// It keeps the original input columns (line) so they round-trip verbatim, the
// 0-based half-open span used for overlap detection, the strand (for -s), and
// the source format so clipped output can be re-encoded the way upstream does.
type inRecord struct {
	chrom  string
	start  int // 0-based, half-open
	end    int
	strand string // "+", "-", ".", "*" or "" when absent
	line   string // original tab-joined columns, verbatim
	fields []string
	format recordFormat
}

// clippedLine renders this record clipped to the overlap span [s,e) (0-based),
// echoing every original column but with the coordinate fields adjusted, the
// same way upstream prints the default intersection output for each record
// type:
//
//   - BED: columns 2,3 become the clipped 0-based start/end.
//   - GFF: column 4 becomes start+1 (GFF is 1-based), column 5 becomes end.
//   - VCF: VCF records are never clipped; the full original line is echoed
//     (upstream's vcfRecord::print ignores the clip coordinates).
//   - BAM: handled separately (printed as a BED12 line); not reached here.
func (r *inRecord) clippedLine(s, e int) string {
	switch r.format {
	case fmtVCF:
		return r.line
	case fmtGFF:
		f := append([]string(nil), r.fields...)
		f[3] = strconv.Itoa(s + 1)
		f[4] = strconv.Itoa(e)
		return strings.Join(f, "\t")
	default: // fmtBED
		f := append([]string(nil), r.fields...)
		f[1] = strconv.Itoa(s)
		f[2] = strconv.Itoa(e)
		return strings.Join(f, "\t")
	}
}

// isInteger reports whether s is a base-10 integer (optionally signed),
// matching upstream's ParseTools::isInteger used by BedFile::parseLine.
func isInteger(s string) bool {
	if s == "" {
		return false
	}
	_, err := strconv.Atoi(s)
	return err == nil
}

// detectTextFormat classifies a tokenized data line as BED, VCF, or GFF using
// upstream BedFile::parseLine's precedence: BED first (cols 2,3 integer), then
// VCF (col 2 integer and >= 8 cols), then GFF (8 or 9 cols with cols 4,5
// integer). It returns ok=false when none match.
func detectTextFormat(fields []string) (recordFormat, bool) {
	n := len(fields)
	if n < 3 {
		return 0, false
	}
	if isInteger(fields[1]) && isInteger(fields[2]) {
		return fmtBED, true
	}
	if isInteger(fields[1]) && n >= 8 {
		return fmtVCF, true
	}
	if (n == 8 || n == 9) && isInteger(fields[3]) && isInteger(fields[4]) {
		return fmtGFF, true
	}
	return 0, false
}

// hasEnoughFields reports whether a tokenized data line has enough columns for
// the locked format's coordinate fields: BED needs >= 3 (chrom,start,end), VCF
// needs >= 8 (through INFO, with REF at col 4), and GFF needs >= 5 (through
// end). This guards parseTextRecord's field indexing on every line.
func hasEnoughFields(fields []string, format recordFormat) bool {
	switch format {
	case fmtVCF:
		return len(fields) >= 8
	case fmtGFF:
		return len(fields) >= 5
	default: // fmtBED
		return len(fields) >= 3
	}
}

// parseTextRecord converts a tokenized data line into an inRecord with a
// 0-based half-open span, applying the resolved format's coordinate
// conventions exactly as upstream does:
//
//   - BED: start = col2, end = col3 (already 0-based half-open).
//   - VCF: start = POS-1, end = start + len(REF) (col4).
//   - GFF: start = col4-1, end = col5 (1-based inclusive -> 0-based half-open).
func parseTextRecord(line string, fields []string, format recordFormat) (*inRecord, error) {
	rec := &inRecord{line: line, fields: fields, format: format}
	switch format {
	case fmtBED:
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid chromStart %q: %w", fields[1], err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("invalid chromEnd %q: %w", fields[2], err)
		}
		rec.chrom, rec.start, rec.end = fields[0], start, end
		if len(fields) > 5 {
			rec.strand = fields[5]
		}
	case fmtVCF:
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid VCF POS %q: %w", fields[1], err)
		}
		start := pos - 1
		end := start + len(fields[3])
		if start < 0 || end < start {
			return nil, fmt.Errorf("malformed VCF entry: start=%d end=%d", start, end)
		}
		rec.chrom, rec.start, rec.end = fields[0], start, end
		// VCF has no BED-style strand column; -s on VCF is unsupported upstream.
	case fmtGFF:
		gstart, err := strconv.Atoi(fields[3])
		if err != nil {
			return nil, fmt.Errorf("invalid GFF start %q: %w", fields[3], err)
		}
		gend, err := strconv.Atoi(fields[4])
		if err != nil {
			return nil, fmt.Errorf("invalid GFF end %q: %w", fields[4], err)
		}
		start := gstart - 1
		end := gend
		if start < 0 || end < start {
			return nil, fmt.Errorf("malformed GFF entry: start=%d end=%d", start, end)
		}
		rec.chrom, rec.start, rec.end = fields[0], start, end
		if len(fields) > 6 {
			rec.strand = fields[6]
		}
	default:
		return nil, fmt.Errorf("unknown text format")
	}
	return rec, nil
}

// readInRecords reads every interval from r. It autodetects whether the stream
// is BAM (by sniffing the leading bytes) or a text format (BED/VCF/GFF). The
// text format is locked on the first non-header data line, mirroring upstream's
// BedFile::parseLine. Header, track, browser and blank lines are skipped.
// `##fileformat=VCF` forces VCF detection.
func readInRecords(r io.Reader) ([]*inRecord, error) {
	br := bufio.NewReaderSize(r, 64*1024)
	magic, _ := br.Peek(4)
	if len(magic) >= 4 && magic[0] == 0x1f && magic[1] == 0x8b {
		// Gzip/BGZF magic: could be a BGZF-wrapped BAM. Probe the BAM reader.
		recs, ok, err := readBAMRecords(br)
		if err != nil {
			return nil, err
		}
		if ok {
			return recs, nil
		}
		// Not BAM (a plain-gzip text file would already have been
		// decompressed by iohelper.OpenReader before reaching here, so a
		// stream that still has gzip magic and is not BAM is unexpected).
		return nil, fmt.Errorf("gzip-compressed non-BAM input is not supported on -a/-b; decompress first")
	}
	if len(magic) >= 4 && string(magic) == "BAM\x01" {
		recs, ok, err := readBAMRecords(br)
		if err != nil {
			return nil, err
		}
		if ok {
			return recs, nil
		}
	}
	return readTextRecords(br)
}

// readTextRecords reads BED/VCF/GFF records from an already-peeked reader.
func readTextRecords(br *bufio.Reader) ([]*inRecord, error) {
	sc := bufio.NewScanner(br)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	var out []*inRecord
	var format recordFormat
	formatSet := false
	headerForcedVCF := false
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r")
		trimmed := strings.TrimSpace(line)
		if trimmed == "" {
			continue
		}
		if strings.HasPrefix(trimmed, "#") || strings.HasPrefix(trimmed, "track") ||
			strings.HasPrefix(trimmed, "browser") {
			if strings.HasPrefix(line, "##fileformat=VCF") {
				headerForcedVCF = true
			}
			continue
		}
		fields := strings.Split(line, "\t")
		if !formatSet {
			if headerForcedVCF {
				format = fmtVCF
			} else {
				f, ok := detectTextFormat(fields)
				if !ok {
					return nil, fmt.Errorf("unexpected file format: please use tab-delimited BED, GFF, or VCF (line: %q)", line)
				}
				format = f
			}
			formatSet = true
		}
		// Validate the field count against the locked format on EVERY line, not
		// just the first: upstream's type checker errors on a data line with the
		// wrong number of fields rather than indexing past the end.
		if !hasEnoughFields(fields, format) {
			return nil, fmt.Errorf("type checker found wrong number of fields while tokenizing data line (line: %q)", line)
		}
		rec, err := parseTextRecord(line, fields, format)
		if err != nil {
			return nil, err
		}
		out = append(out, rec)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// readBAMRecords decodes a BAM stream into inRecords, converting each primary,
// mapped alignment into a BED12 line exactly as upstream's BamRecord print path
// does (chrom start end name MAPQ strand thickStart thickEnd 0,0,0 numBlocks
// blockSizes, blockStarts; blocks derived from the CIGAR splitting on N ops).
// It returns ok=false (without consuming the reader past the header probe) when
// the stream is not BAM. Unmapped reads are skipped, matching `bedtools
// intersect -abam ... -bed` which only emits mapped alignments.
func readBAMRecords(br *bufio.Reader) ([]*inRecord, bool, error) {
	reader, err := sam.NewReader(br)
	if err != nil {
		if err == sam.ErrNotBAM {
			return nil, false, nil
		}
		return nil, false, err
	}
	var out []*inRecord
	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, true, fmt.Errorf("error reading BAM: %w", err)
		}
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" {
			continue
		}
		in := bamToBED12(rec)
		out = append(out, in)
	}
	return out, true, nil
}

// bamToBED12 converts a single mapped SAM/BAM alignment into a BED12 inRecord,
// mirroring upstream BamRecord::print + printRemainingBamFields. The blocks are
// the CIGAR's reference-consuming runs, broken on N ops (M/=/X advance the
// block, D extends it, N closes it). The colour is fixed at 0,0,0 and the
// block-size / block-start lists carry the trailing comma upstream emits.
func bamToBED12(rec *sam.Record) *inRecord {
	start := int(rec.Pos) - 1 // BAM Pos is 1-based; BED is 0-based.
	blocks := cigarBlocks(rec, start)
	end := start
	if len(blocks) > 0 {
		end = blocks[len(blocks)-1].end
	} else {
		end = int(rec.EndPosition()) - 1 // fallback; should not happen for mapped reads
	}
	strand := "+"
	if rec.Flag&sam.FlagReverse != 0 {
		strand = "-"
	}
	var sizes, starts strings.Builder
	for _, b := range blocks {
		sizes.WriteString(strconv.Itoa(b.end - b.start))
		sizes.WriteByte(',')
		starts.WriteString(strconv.Itoa(b.start - start))
		starts.WriteByte(',')
	}
	fields := []string{
		rec.RName,
		strconv.Itoa(start),
		strconv.Itoa(end),
		rec.QName,
		strconv.Itoa(int(rec.MapQ)),
		strand,
		strconv.Itoa(start),
		strconv.Itoa(end),
		"0,0,0",
		strconv.Itoa(len(blocks)),
		sizes.String(),
		starts.String(),
	}
	line := strings.Join(fields, "\t")
	return &inRecord{
		chrom:  rec.RName,
		start:  start,
		end:    end,
		strand: strand,
		line:   line,
		fields: fields,
		format: fmtBAM,
	}
}

// cigarBlocks returns the reference-consuming blocks of an alignment starting
// at 0-based refStart, breaking the alignment into separate blocks on each N
// (skip) op. M/=/X/D consume the reference; D extends the current block (so a
// deletion does not split it, matching upstream's breakOnDeletionOps=false); N
// closes the current block and starts a new one after the gap; I/S/H/P do not
// consume reference and are ignored.
func cigarBlocks(rec *sam.Record, refStart int) []block {
	var blocks []block
	cur := block{start: refStart, end: refStart}
	open := false
	pos := refStart
	for _, op := range rec.Cigar {
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch, sam.CigarDeletion:
			if !open {
				cur = block{start: pos, end: pos}
				open = true
			}
			pos += int(op.Length())
			cur.end = pos
		case sam.CigarSkipped:
			if open {
				blocks = append(blocks, cur)
				open = false
			}
			pos += int(op.Length())
		default:
			// I, S, H, P: no reference advance.
		}
	}
	if open {
		blocks = append(blocks, cur)
	}
	if len(blocks) == 0 {
		// No reference-consuming ops; fall back to a single zero-length block.
		blocks = append(blocks, block{start: refStart, end: refStart})
	}
	return blocks
}
