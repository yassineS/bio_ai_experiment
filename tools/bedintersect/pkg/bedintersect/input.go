// Input-format handling for bedintersect. Upstream `bedtools intersect`
// accepts BED, VCF, GFF, and BAM input on -a and -b. This file converts each
// of those into a common interval record that preserves the original input
// columns (so they can be echoed verbatim) and knows how to re-render its
// coordinates when clipped to an overlap, mirroring upstream's per-record-type
// print methods.
package bedintersect

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
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
	// order is the record's position within its chromosome's B slice. It lets
	// the interval-tree path restore B-file order (upstream echoes overlapping
	// B records in input order) after an out-of-order tree query.
	order int
	// dbID is the 0-based index of the B file this record came from. It is 0
	// for A records and for the single-B-file case; with multiple B files it
	// drives the DB-id column (-names/-filenames/numeric) and per-file -C
	// counts. A records ignore it.
	dbID int
	// unmapped marks an unmapped BAM alignment. Such a record is printed with
	// its stored RNAME (which may name the mate's chromosome) but never
	// participates in an overlap, so it is reported only by the modes that emit
	// non-overlapping A records (-v/-loj/-wao/-c), exactly like upstream's
	// printUnmapped path.
	unmapped bool
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

// vcfEnd computes the 0-based half-open end of a VCF record, mirroring upstream
// VcfRecord::initFromFile. A structural-variant ALT (one beginning with "<")
// derives its length from the INFO SVLEN/END tags (vcfSVlen), except a "<INS...>"
// insertion which is treated as zero-length (end == start). A non-SV record's
// end is the start plus the REF allele length.
func vcfEnd(start int, fields []string) int {
	alt := ""
	if len(fields) > 4 {
		alt = fields[4]
	}
	ref := ""
	if len(fields) > 3 {
		ref = fields[3]
	}
	if len(alt) > 0 && alt[0] == '<' {
		if strings.HasPrefix(alt, "<INS") {
			return start // insertions are zero-length
		}
		svlen := vcfSVlen(start, fields)
		if svlen == intMin {
			// No SVLEN/END found; fall back to the REF length like a plain record.
			return start + len(ref)
		}
		return start + svlen
	}
	return start + len(ref)
}

// intMin sentinels "SVLEN/END not found", matching upstream getVcfSVlen's
// INT_MIN return.
const intMin = -1 << 31

// vcfSVlen parses the INFO field (column 8) for a structural-variant length,
// mirroring upstream SingleLineDelimTextFileReader::getVcfSVlen: an SVLEN tag
// wins (the absolute value, or the absolute-max when comma-separated); otherwise
// an END tag yields END-POS+1. Returns intMin when neither is present. POS here
// is the 0-based start, so END-(start+1)+1 == END-start matches upstream's
// END-POS+1 on the 1-based POS.
func vcfSVlen(start int, fields []string) int {
	if len(fields) <= 7 {
		return intMin
	}
	for _, f := range strings.Split(fields[7], ";") {
		if f == "." {
			continue
		}
		eq := strings.IndexByte(f, '=')
		if eq < 0 {
			continue
		}
		key, val := f[:eq], f[eq+1:]
		switch key {
		case "SVLEN":
			best := 0
			found := false
			for _, part := range strings.Split(val, ",") {
				n, err := strconv.Atoi(strings.TrimSpace(part))
				if err != nil {
					continue
				}
				if n < 0 {
					n = -n
				}
				if !found || n > best {
					best = n
					found = true
				}
			}
			if found {
				return best
			}
		case "END":
			end, err := strconv.Atoi(strings.TrimSpace(val))
			if err != nil {
				continue
			}
			// Upstream: END - POS + 1, with POS the 1-based start (start+1).
			return end - (start + 1) + 1
		}
	}
	return intMin
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
		end := vcfEnd(start, fields)
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
// is BAM/CRAM (by sniffing the leading bytes) or a text format (BED/VCF/GFF).
// gzip/BGZF-compressed text (e.g. a piped `.bed.gz` on stdin, which iohelper
// does not transparently decompress) is gunzipped first. The text format is
// locked on the first non-header data line, mirroring upstream's
// BedFile::parseLine. Header, track, browser and blank lines are skipped;
// `##fileformat=VCF` forces VCF detection.
func readInRecords(r io.Reader) ([]*inRecord, error) {
	br := bufio.NewReaderSize(r, 64*1024)
	magic, _ := br.Peek(4)
	if len(magic) >= 4 && magic[0] == 0x1f && magic[1] == 0x8b {
		// Gzip/BGZF magic. BAM is BGZF-wrapped, so gunzip the leading block and
		// check for the BAM magic; a BAM stream is decoded by the BAM reader,
		// otherwise the decompressed bytes are (BGZF/gzip-compressed) text.
		buf, err := io.ReadAll(br)
		if err != nil {
			return nil, err
		}
		if isGzippedBAM(buf) {
			recs, _, berr := readBAMRecords(bufio.NewReader(bytes.NewReader(buf)))
			return recs, berr
		}
		gz, err := gzip.NewReader(bytes.NewReader(buf))
		if err != nil {
			return nil, fmt.Errorf("error opening gzip input: %w", err)
		}
		gz.Multistream(true)
		return readTextRecords(bufio.NewReaderSize(gz, 64*1024))
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
	if len(magic) >= 4 && string(magic) == "CRAM" {
		return readCRAMRecords(br)
	}
	return readTextRecords(br)
}

// isGzippedBAM reports whether a gzip/BGZF-compressed buffer decompresses to a
// BAM stream (its decompressed prefix is the "BAM\x01" magic). It reads only the
// first decompressed bytes, so a plain gzipped/BGZF text file is cheaply
// distinguished from a BAM without parsing the whole stream.
func isGzippedBAM(buf []byte) bool {
	gz, err := gzip.NewReader(bytes.NewReader(buf))
	if err != nil {
		return false
	}
	defer gz.Close()
	head := make([]byte, 4)
	n, _ := io.ReadFull(gz, head)
	return n == 4 && string(head) == "BAM\x01"
}

// readAllB reads every B file in order, tagging each record with its 0-based
// file id (dbID) so multi-database output can render the DB-id column and
// per-file -C counts. It also returns the first file's classification (record
// type + field count), which upstream uses to shape the null-B placeholder
// under -loj/-wao; the placeholder is determined by the first B file alone.
func readAllB(readersB []io.Reader) (recs []*inRecord, dbType dbRecordType, dbFields int, err error) {
	dbType, dbFields = dbBed3, 3
	for i, rb := range readersB {
		fileRecs, rerr := readInRecords(rb)
		if rerr != nil {
			return nil, dbType, dbFields, fmt.Errorf("error reading B intervals: %w", rerr)
		}
		if i == 0 {
			dbType, dbFields = classifyDB(fileRecs)
		}
		for _, r := range fileRecs {
			r.dbID = i
		}
		recs = append(recs, fileRecs...)
	}
	return recs, dbType, dbFields, nil
}

// readCRAMRecords decodes a CRAM stream into inRecords. Like the BAM path it
// only needs each mapped alignment's coordinate fields (RName, Pos, CIGAR,
// MAPQ, strand) — which CRAM stores directly — so no decode reference is
// required; the reader honours REF_CACHE/REF_PATH if set but the 'N'-base
// fallback is irrelevant to interval extraction. Unmapped reads are skipped,
// matching the BAM path and upstream `bedtools intersect` over CRAM.
func readCRAMRecords(br *bufio.Reader) ([]*inRecord, error) {
	reader, err := alnio.NewReaderWithReference(br, "")
	if err != nil {
		return nil, fmt.Errorf("error opening CRAM: %w", err)
	}
	var out []*inRecord
	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading CRAM: %w", err)
		}
		if rec.IsUnmapped() || rec.RName == "" || rec.RName == "*" {
			continue
		}
		out = append(out, bamToBED12(rec))
	}
	return out, nil
}

// readInRecordsWithHeader reads every interval from r and, when wantHeader is
// set, also returns the file's leading header text (the comment/track/browser
// lines before the first data record) so -header can echo it verbatim. The
// returned header string already ends with a newline when non-empty. Header
// capture applies to text inputs, including gzip/BGZF-compressed text;
// BAM/CRAM return an empty header.
func readInRecordsWithHeader(r io.Reader, wantHeader bool) ([]*inRecord, string, error) {
	br := bufio.NewReaderSize(r, 64*1024)
	magic, _ := br.Peek(4)
	isGzip := len(magic) >= 4 && magic[0] == 0x1f && magic[1] == 0x8b
	isBAMorCRAM := len(magic) >= 4 && (string(magic) == "BAM\x01" || string(magic) == "CRAM")
	if !wantHeader || isBAMorCRAM {
		recs, err := readInRecords(br)
		return recs, "", err
	}
	if isGzip {
		// Buffer, probe for BAM; if not BAM, gunzip and header-scan the text.
		buf, err := io.ReadAll(br)
		if err != nil {
			return nil, "", err
		}
		if isGzippedBAM(buf) {
			recs, _, berr := readBAMRecords(bufio.NewReader(bytes.NewReader(buf)))
			return recs, "", berr // BAM carries no text header to echo here
		}
		gz, err := gzip.NewReader(bytes.NewReader(buf))
		if err != nil {
			return nil, "", fmt.Errorf("error opening gzip input: %w", err)
		}
		gz.Multistream(true)
		return readTextRecordsHeader(bufio.NewReaderSize(gz, 64*1024))
	}
	return readTextRecordsHeader(br)
}

// readTextRecords reads BED/VCF/GFF records from an already-peeked reader.
func readTextRecords(br *bufio.Reader) ([]*inRecord, error) {
	recs, _, err := scanTextRecords(br, false)
	return recs, err
}

// readTextRecordsHeader is readTextRecords but also captures the leading header
// text (for -header).
func readTextRecordsHeader(br *bufio.Reader) ([]*inRecord, string, error) {
	return scanTextRecords(br, true)
}

// scanTextRecords parses BED/VCF/GFF records and, when wantHeader is set,
// collects the leading header lines (everything before the first data record)
// into the returned header string, newline-terminated.
func scanTextRecords(br *bufio.Reader, wantHeader bool) ([]*inRecord, string, error) {
	sc := bufio.NewScanner(br)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	var out []*inRecord
	var header strings.Builder
	var format recordFormat
	formatSet := false
	headerForcedVCF := false
	expectedFields := 0
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
			if wantHeader && !formatSet {
				header.WriteString(line)
				header.WriteByte('\n')
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
					return nil, "", fmt.Errorf("unexpected file format: please use tab-delimited BED, GFF, or VCF (line: %q)", line)
				}
				format = f
			}
			formatSet = true
			expectedFields = len(fields)
		}
		// Upstream's type checker locks the column count on the first data line
		// and errors on any later line with a different count (e.g. a stray
		// trailing tab adds an empty field). It also needs enough columns for the
		// locked format's coordinate fields.
		if len(fields) != expectedFields || !hasEnoughFields(fields, format) {
			return nil, "", &fieldCountError{}
		}
		rec, err := parseTextRecord(line, fields, format)
		if err != nil {
			return nil, "", err
		}
		out = append(out, rec)
	}
	if err := sc.Err(); err != nil {
		return nil, "", err
	}
	return out, header.String(), nil
}

// readBAMRecords decodes a BAM stream into inRecords. Each mapped alignment
// becomes a BED12 line exactly as upstream's BamRecord print path does (chrom
// start end name MAPQ strand thickStart thickEnd 0,0,0 numBlocks blockSizes,
// blockStarts; blocks derived from the CIGAR splitting on N ops). Unmapped reads
// are kept as null-placeholder records (chrom ".", start/end -1) so they flow
// through -v/-loj/-wao/-c exactly as upstream's printUnmapped does, rather than
// being silently dropped. It returns ok=false when the stream is not BAM.
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
			out = append(out, unmappedBAMRecord(rec))
			continue
		}
		out = append(out, bamToBED12(rec))
	}
	return out, true, nil
}

// bamReadName returns the alignment's query name with the upstream "/1" (first
// mate) or "/2" (second mate) suffix appended for paired reads, mirroring
// BamFileReader::getName.
func bamReadName(rec *sam.Record) string {
	name := rec.QName
	switch {
	case rec.Flag&sam.FlagRead1 != 0:
		name += "/1"
	case rec.Flag&sam.FlagRead2 != 0:
		name += "/2"
	}
	return name
}

// unmappedBAMRecord builds the null-placeholder inRecord upstream emits for an
// unmapped alignment (BamRecord::printUnmapped): chrom ".", start/end -1, the
// read name (with /1 or /2), MAPQ as the score, and a fixed empty BED12 tail.
// Its chrom is "." so it never overlaps any B record, which is exactly how an
// unmapped read behaves: absent in default mode, count 0 under -c, reported
// under -v, and paired with a null B under -loj/-wao.
func unmappedBAMRecord(rec *sam.Record) *inRecord {
	name := bamReadName(rec)
	if name == "" {
		name = "."
	}
	// Upstream prints the record's RNAME: "." when truly unaligned, or the
	// (mate's) chromosome htslib propagated onto the read. The record is still
	// flagged unmapped below so it never overlaps regardless of this name.
	chrom := rec.RName
	if chrom == "" || chrom == "*" {
		chrom = "."
	}
	line := chrom + "\t-1\t-1\t" + name + "\t" + strconv.Itoa(int(rec.MapQ)) +
		"\t.\t-1\t-1\t-1\t0,0,0\t0\t.\t."
	return &inRecord{
		chrom:    chrom,
		start:    -1,
		end:      -1,
		strand:   ".",
		line:     line,
		fields:   strings.Split(line, "\t"),
		format:   fmtBAM,
		unmapped: true,
	}
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
		bamReadName(rec),
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
