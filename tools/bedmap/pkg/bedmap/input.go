// Input-format handling for bedmap. Upstream `bedtools map` accepts BED, VCF,
// GFF, and BAM input on -b (the database file). This file converts each of
// those into a common interval record that preserves the original input
// columns (so a -c column extraction reads the literal source columns), the
// 0-based half-open span used for overlap detection, and the strand (for -s).
//
// The detection precedence mirrors upstream BedFile::parseLine: BED first
// (cols 2,3 integer), then VCF (col 2 integer and >= 8 cols), then GFF (8 or 9
// cols with cols 4,5 integer). BAM is sniffed by its leading magic bytes.
package bedmap

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
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
// end).
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

// parseTextRecord converts a tokenized data line into a rawRecord with a
// 0-based half-open span, applying the resolved format's coordinate
// conventions exactly as upstream does:
//
//   - BED: start = col2, end = col3 (already 0-based half-open).
//   - VCF: start = POS-1, end = start + len(REF) (col4).
//   - GFF: start = col4-1, end = col5 (1-based inclusive -> 0-based half-open).
//
// The original fields are preserved verbatim so a -c column extraction reads
// the literal source columns (matching upstream `bedtools map -b <vcf|gff>`).
func parseTextRecord(fields []string, format recordFormat) (rawRecord, error) {
	switch format {
	case fmtBED:
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return rawRecord{}, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return rawRecord{}, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
		}
		rec := &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}
		if len(fields) > 5 {
			rec.Strand = fields[5]
		}
		return rawRecord{rec: rec, fields: fields}, nil
	case fmtVCF:
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			return rawRecord{}, fmt.Errorf("invalid VCF POS %q: %v", fields[1], err)
		}
		start := pos - 1
		end := start + len(fields[3])
		if start < 0 || end < start {
			return rawRecord{}, fmt.Errorf("malformed VCF entry: start=%d end=%d", start, end)
		}
		// VCF has no BED-style strand column; -s on VCF is unsupported upstream.
		return rawRecord{rec: &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}, fields: fields}, nil
	case fmtGFF:
		gstart, err := strconv.Atoi(fields[3])
		if err != nil {
			return rawRecord{}, fmt.Errorf("invalid GFF start %q: %v", fields[3], err)
		}
		gend, err := strconv.Atoi(fields[4])
		if err != nil {
			return rawRecord{}, fmt.Errorf("invalid GFF end %q: %v", fields[4], err)
		}
		start := gstart - 1
		end := gend
		if start < 0 || end < start {
			return rawRecord{}, fmt.Errorf("malformed GFF entry: start=%d end=%d", start, end)
		}
		rec := &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}
		if len(fields) > 6 {
			rec.Strand = fields[6]
		}
		return rawRecord{rec: rec, fields: fields}, nil
	}
	return rawRecord{}, fmt.Errorf("unknown text format")
}

// readBRecords reads every B record from r, auto-detecting BAM (by sniffing the
// leading bytes) or a text format (BED/VCF/GFF). The text format is locked on
// the first non-header data line, mirroring upstream's BedFile::parseLine.
// Header, track, browser and blank lines are skipped. It returns the parsed
// records and the maximum number of fields seen on any record (used for the
// upstream "only has fields 1 - N" column-range error).
func readBRecords(r io.Reader) ([]rawRecord, int, error) {
	br := bufio.NewReaderSize(r, 64*1024)
	magic, _ := br.Peek(4)
	if (len(magic) >= 4 && magic[0] == 0x1f && magic[1] == 0x8b) ||
		(len(magic) >= 4 && string(magic) == "BAM\x01") {
		return readBAMBRecords(br)
	}
	return readTextBRecords(br)
}

// readTextBRecords reads BED/VCF/GFF records from an already-peeked reader.
func readTextBRecords(br *bufio.Reader) ([]rawRecord, int, error) {
	sc := bufio.NewScanner(br)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	var out []rawRecord
	var format recordFormat
	formatSet := false
	headerForcedVCF := false
	maxFields := 0
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
					second := ""
					if len(fields) > 1 {
						second = fields[1]
					}
					return nil, 0, fmt.Errorf("invalid record: column 2 %q is not a numeric BED start and the line is not a GFF feature", second)
				}
				format = f
			}
			formatSet = true
		}
		if !hasEnoughFields(fields, format) {
			return nil, 0, fmt.Errorf("type checker found wrong number of fields while tokenizing data line (line: %q)", line)
		}
		rr, err := parseTextRecord(fields, format)
		if err != nil {
			return nil, 0, err
		}
		if len(fields) > maxFields {
			maxFields = len(fields)
		}
		out = append(out, rr)
	}
	if err := sc.Err(); err != nil {
		return nil, 0, err
	}
	return out, maxFields, nil
}

// readBAMBRecords decodes a BAM stream into rawRecords, converting each primary,
// mapped alignment into a BED12 line exactly as upstream's BamRecord print path
// does (chrom start end name MAPQ strand thickStart thickEnd 0,0,0 numBlocks
// blockSizes, blockStarts). Unmapped reads are skipped, matching `bedtools map
// -b <bam>`. The maximum field count is the fixed 12 BED12 columns.
func readBAMBRecords(br *bufio.Reader) ([]rawRecord, int, error) {
	reader, err := sam.NewReader(br)
	if err != nil {
		if err == sam.ErrNotBAM {
			return nil, 0, fmt.Errorf("not a valid BAM file")
		}
		return nil, 0, err
	}
	var out []rawRecord
	maxFields := 0
	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, 0, fmt.Errorf("error reading BAM: %w", err)
		}
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" {
			continue
		}
		rr := bamToRaw(rec)
		if len(rr.fields) > maxFields {
			maxFields = len(rr.fields)
		}
		out = append(out, rr)
	}
	return out, maxFields, nil
}

// bamToRaw converts a single mapped SAM/BAM alignment into a BED12 rawRecord,
// mirroring upstream BamRecord::print. The blocks are the CIGAR's
// reference-consuming runs, broken on N ops (M/=/X advance the block, D extends
// it, N closes it). The colour is fixed at 0,0,0 and the block-size / block-
// start lists carry the trailing comma upstream emits.
func bamToRaw(rec *sam.Record) rawRecord {
	start := int(rec.Pos) - 1 // BAM Pos is 1-based; BED is 0-based.
	blocks := cigarBlocks(rec, start)
	end := start
	if len(blocks) > 0 {
		end = blocks[len(blocks)-1].end
	} else {
		end = int(rec.EndPosition()) - 1
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
	return rawRecord{
		rec:    &bed.Record{Chrom: rec.RName, ChromStart: start, ChromEnd: end, Strand: strand},
		fields: fields,
	}
}

// block is a single reference-consuming run of a BAM alignment.
type block struct {
	start int
	end   int
}

// cigarBlocks returns the reference-consuming blocks of an alignment starting
// at 0-based refStart, breaking the alignment into separate blocks on each N
// (skip) op. M/=/X/D consume the reference; D extends the current block; N
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
		blocks = append(blocks, block{start: refStart, end: refStart})
	}
	return blocks
}
