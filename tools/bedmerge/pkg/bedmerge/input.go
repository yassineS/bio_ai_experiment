// Input-format handling for bedmerge. Upstream `bedtools merge` accepts BED,
// GFF, VCF, and BAM input on -i (auto-detected) and merges them. This file
// converts each of those formats into a common interval record that preserves
// the original input columns (so -c/-o column operations can address them) and
// records the 0-based half-open span used for merging, mirroring upstream's
// per-record-type coordinate conventions.
package bedmerge

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// Header/track prefixes recognised on text input, as byte slices so the read
// loop can test the scanner buffer without allocating a string per line.
var (
	trackPrefix         = []byte("track")
	browserPrefix       = []byte("browser")
	fileformatVCFPrefix = []byte("##fileformat=VCF")
)

// inputFormat tags the detected input format of the -i stream.
type inputFormat int

const (
	fmtBED inputFormat = iota
	fmtGFF
	fmtVCF
	fmtBAM
)

// record is the common interval model used by every merge path. It keeps the
// original input columns (fields) so -c/-o can address them, the 0-based
// half-open span used for overlap detection, and the strand for -s/-S.
type record struct {
	chrom  string
	start  int // 0-based, half-open
	end    int
	strand string
	fields []string // columns addressable by -c (1-based index into this slice)
	isBAM  bool     // true for BAM-derived records (empty -c fields print as ".")
}

// detectFormat classifies the first tokenized data line. BED takes precedence
// (cols 2,3 integer), then VCF (col 2 integer with >= 8 cols), then GFF (8 or 9
// cols with cols 4,5 integer), matching upstream BedFile::parseLine.
func detectFormat(fields []string) (inputFormat, bool) {
	n := len(fields)
	if n < 3 {
		return 0, false
	}
	if isChrPos(fields[1]) && isChrPos(fields[2]) {
		return fmtBED, true
	}
	if isChrPos(fields[1]) && n >= 8 {
		return fmtVCF, true
	}
	if (n == 8 || n == 9) && isChrPos(fields[3]) && isChrPos(fields[4]) {
		return fmtGFF, true
	}
	return 0, false
}

// FieldCountError reports an input line whose tab-delimited field count differs
// from the number established by the first valid data line. It mirrors the two
// distinct messages upstream bedtools emits for this condition:
//
//   - TypeChecker == true: the inconsistency appears within the first four valid
//     data lines, so upstream's format auto-detector rejects the file with
//     "Error: Type checker found wrong number of fields while tokenizing data
//     line. / Perhaps you have extra TAB at the end of your line? ...". The
//     LineNum/Got/Want fields are unused in this case.
//   - TypeChecker == false: a later line disagrees, so upstream's per-line reader
//     reports "Error: line number N of file F has X fields, but Y were
//     expected." (LineNum=N, Got=X, Want=Y).
//
// The CLI (reportMergeError) formats the exact upstream wording, substituting
// the input file name (Filename).
type FieldCountError struct {
	TypeChecker bool
	Filename    string
	LineNum     int
	Got         int
	Want        int
}

func (e *FieldCountError) Error() string {
	if e.TypeChecker {
		return "type checker found wrong number of fields while tokenizing data line"
	}
	return fmt.Sprintf("line number %d of file %s has %d fields, but %d were expected",
		e.LineNum, e.Filename, e.Got, e.Want)
}

// SortOrderError reports input that is not coordinate-sorted, matching upstream
// bedtools merge's requirement that input be sorted by chromosome then start.
// Upstream prints "Error: Sorted input specified, but the file F has the
// following out of order record\n<record>" and exits 1; the CLI formats that
// wording from Filename and Line (the BED-reconstructed offending record).
type SortOrderError struct {
	Filename string
	Line     string
}

func (e *SortOrderError) Error() string {
	return fmt.Sprintf("sorted input specified, but the file %s has the following out of order record\n%s",
		e.Filename, e.Line)
}

// sortState tracks the running coordinate-sort check across the input stream,
// reproducing upstream FileRecordMgr::testInputSortOrder: input must be sorted
// by start within each chromosome, and a chromosome may not reappear once a
// different chromosome has been seen. The sort key is the record's original
// start (upstream's zero-length adjustment cancels out in the comparison).
type sortState struct {
	prevChrom   string
	prevStart   int
	havePrev    bool
	foundChroms map[string]struct{}
}

// check tests one record (in input order) against the running sort state and
// returns a SortOrderError when it is out of order, leaving recs unsorted at the
// point of failure. line is the BED-reconstructed offending-record text upstream
// would print.
func (s *sortState) check(chrom string, start int, line string, filename string) error {
	if s.foundChroms == nil {
		s.foundChroms = make(map[string]struct{})
	}
	if !s.havePrev || chrom != s.prevChrom {
		if _, seen := s.foundChroms[chrom]; seen {
			// A chromosome reappears after a different one: out of order.
			return &SortOrderError{Filename: filename, Line: line}
		}
		s.foundChroms[chrom] = struct{}{}
		s.prevChrom = chrom
		s.prevStart = start
		s.havePrev = true
		return nil
	}
	// Same chromosome as the previous record.
	if start < s.prevStart {
		return &SortOrderError{Filename: filename, Line: line}
	}
	s.prevStart = start
	return nil
}

// sortErrorLine renders the offending record exactly as upstream's Record
// operator<< does for the sort-order message: chrom, the (possibly
// zero-length-adjusted) start and end, then any original columns from index 4
// on. For GFF/VCF, which upstream prints in their native column layout, it falls
// back to the original line text.
func sortErrorLine(rec record, format inputFormat, rawLine string) string {
	if format != fmtBED {
		return rawLine
	}
	start, end := rec.start, rec.end
	if start == end { // zero-length: upstream adjusts start-1, end+1 for display.
		start--
		end++
	}
	var b strings.Builder
	b.WriteString(rec.chrom)
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(start))
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(end))
	// Append original columns beyond chrom/start/end when available.
	if len(rec.fields) > 3 {
		for _, f := range rec.fields[3:] {
			b.WriteByte('\t')
			b.WriteString(f)
		}
	} else {
		// Fast-path records keep no fields; recover trailing columns from the
		// raw line (everything after the third tab).
		if extra := trailingColumns(rawLine); extra != "" {
			b.WriteByte('\t')
			b.WriteString(extra)
		}
	}
	return b.String()
}

// trailingColumns returns the substring of a tab-delimited line after its third
// field (columns 4+), or "" when the line has three or fewer fields.
func trailingColumns(line string) string {
	tabs := 0
	for i := 0; i < len(line); i++ {
		if line[i] == '\t' {
			tabs++
			if tabs == 3 {
				return line[i+1:]
			}
		}
	}
	return ""
}

// fieldCount returns the number of tab-delimited fields in a line.
func fieldCount(line string) int {
	if line == "" {
		return 0
	}
	return strings.Count(line, "\t") + 1
}

// readText reads BED/GFF/VCF records from an already-peeked reader. It validates
// input the way upstream bedtools merge does: every data line must have the same
// field count as the first valid data line, and the records must be sorted by
// chromosome then start. Malformed or unsorted input yields a *FieldCountError
// or *SortOrderError, which the CLI renders with upstream-exact wording.
func readText(br *bufio.Reader, opts MergeOptions) ([]record, inputFormat, error) {
	src := newTextRecordSource(br, opts)
	var out []record
	for {
		rec, ok := src.next()
		if !ok {
			break
		}
		out = append(out, rec)
	}
	if src.err != nil {
		return nil, src.format, src.err
	}
	return out, src.format, nil
}

// textRecordSource yields BED/GFF/VCF records one at a time from a scanner,
// running upstream bedtools merge's per-line validation (field-count parity,
// format detection, coordinate sort order) as each record is read — exactly the
// loop readText used to run before collecting into a slice. next returns
// ok=false at end of input or when a validation error is hit (in s.err), so a
// caller can stream records into the merge without materialising them all.
type textRecordSource struct {
	sc              *bufio.Scanner
	opts            MergeOptions
	keepFields      bool
	ci              chromInterner
	sorts           sortState
	format          inputFormat
	formatSet       bool
	headerForcedVCF bool
	expectFields    int
	validData       int
	lineNum         int
	filename        string
	err             error
}

func newTextRecordSource(br *bufio.Reader, opts MergeOptions) *textRecordSource {
	sc := bufio.NewScanner(br)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	fn := opts.Filename
	if fn == "" {
		fn = "-"
	}
	return &textRecordSource{sc: sc, opts: opts, keepFields: opts.needsFields(), filename: fn}
}

func (s *textRecordSource) next() (record, bool) {
	if s.err != nil {
		return record{}, false
	}
	for s.sc.Scan() {
		s.lineNum++
		lineB := bytes.TrimRight(s.sc.Bytes(), "\r")
		trimmed := bytes.TrimSpace(lineB)
		if len(trimmed) == 0 {
			continue
		}
		if trimmed[0] == '#' || bytes.HasPrefix(trimmed, trackPrefix) ||
			bytes.HasPrefix(trimmed, browserPrefix) {
			if bytes.HasPrefix(lineB, fileformatVCFPrefix) {
				s.headerForcedVCF = true
			}
			continue
		}
		rawLine := string(lineB)
		nf := fieldCount(rawLine)
		s.validData++
		if s.expectFields == 0 {
			s.expectFields = nf
		} else if nf != s.expectFields {
			if s.validData <= 4 {
				s.err = &FieldCountError{TypeChecker: true, Filename: s.filename}
			} else {
				s.err = &FieldCountError{Filename: s.filename, LineNum: s.lineNum, Got: nf, Want: s.expectFields}
			}
			return record{}, false
		}
		if !s.formatSet {
			if s.headerForcedVCF {
				s.format = fmtVCF
			} else {
				f, ok := detectFormat(strings.Split(rawLine, "\t"))
				if !ok {
					s.err = fmt.Errorf("unexpected file format: please use tab-delimited BED, GFF, or VCF (line: %q)", rawLine)
					return record{}, false
				}
				s.format = f
			}
			if s.format == fmtVCF && s.opts.StrandSpec {
				s.err = errStrandedVCF
				return record{}, false
			}
			s.formatSet = true
		}
		var (
			rec record
			err error
		)
		if s.keepFields || s.format != fmtBED {
			rec, err = parseTextRecord(strings.Split(rawLine, "\t"), s.format)
		} else {
			rec, err = parseBEDFast(lineB, s.opts.StrandSpec || s.opts.StrandFilter != "", &s.ci)
		}
		if err != nil {
			s.err = err
			return record{}, false
		}
		if serr := s.sorts.check(rec.chrom, rec.start, "", s.filename); serr != nil {
			soErr := serr.(*SortOrderError)
			soErr.Line = sortErrorLine(rec, s.format, rawLine)
			s.err = soErr
			return record{}, false
		}
		return rec, true
	}
	if e := s.sc.Err(); e != nil {
		s.err = e
	}
	return record{}, false
}

// chromInterner caches the most recently interned chromosome name so a run of
// records sharing a chromosome (the norm for sorted BED input) allocates the
// name string once rather than per record. The `c.last == string(b)` compare is
// allocation-free: the compiler does not heap-allocate a []byte->string
// conversion used only as a comparison operand.
type chromInterner struct {
	last string
}

func (c *chromInterner) intern(b []byte) string {
	if c.last != "" && c.last == string(b) {
		return c.last
	}
	s := string(b)
	c.last = s
	return s
}

// internStrand maps the byte form of a strand column to a shared string,
// allocating nothing for the only values that occur in practice.
func internStrand(b []byte) string {
	switch string(b) {
	case "":
		return ""
	case "+":
		return "+"
	case "-":
		return "-"
	case ".":
		return "."
	}
	return string(b)
}

// parseBEDFast parses a BED data line directly from the scanner's byte buffer,
// copying only the chromosome (and the strand, when stranded merging is active)
// so the line buffer need not be retained. It mirrors parseTextRecord's BED case
// but skips building and keeping the full field slice, which is unused when no
// -c column operation or field-echoing output mode is requested. chrom is
// interned through ci to avoid an allocation per record on runs of one
// chromosome.
func parseBEDFast(line []byte, keepStrand bool, ci *chromInterner) (record, error) {
	var cols [6][]byte
	n := 0
	begin := 0
	for i := 0; i <= len(line); i++ {
		if i == len(line) || line[i] == '\t' {
			if n < len(cols) {
				cols[n] = line[begin:i]
			}
			n++
			begin = i + 1
		}
	}
	if n < 3 {
		return record{}, fmt.Errorf("BED record must have at least 3 fields, got %d", n)
	}
	start, err := parseChrPosBytes(cols[1])
	if err != nil {
		return record{}, fmt.Errorf("invalid chromStart %q: %w", cols[1], err)
	}
	end, err := parseChrPosBytes(cols[2])
	if err != nil {
		return record{}, fmt.Errorf("invalid chromEnd %q: %w", cols[2], err)
	}
	strand := ""
	if keepStrand && n > 5 {
		strand = internStrand(cols[5])
	}
	return record{chrom: ci.intern(cols[0]), start: start, end: end, strand: strand}, nil
}

// errStrandedVCF is the sentinel for the unsupported "-s with VCF" combination;
// the CLI prints the upstream-formatted message.
var errStrandedVCF = fmt.Errorf("stranded merge not supported for VCF file")

// parseTextRecord converts a tokenized data line into a record with a 0-based
// half-open span, applying the format's coordinate conventions exactly as
// upstream does:
//
//   - BED: start = col2, end = col3 (already 0-based half-open); strand = col6.
//   - GFF: start = col4-1, end = col5 (1-based inclusive -> 0-based); strand = col7.
//   - VCF: start = POS-1; end = start + SV length when the ALT is a symbolic
//     structural variant, else start + len(REF).
func parseTextRecord(fields []string, format inputFormat) (record, error) {
	switch format {
	case fmtBED:
		if len(fields) < 3 {
			return record{}, fmt.Errorf("BED record must have at least 3 fields, got %d", len(fields))
		}
		start, err := parseChrPos(fields[1])
		if err != nil {
			return record{}, fmt.Errorf("invalid chromStart %q: %w", fields[1], err)
		}
		end, err := parseChrPos(fields[2])
		if err != nil {
			return record{}, fmt.Errorf("invalid chromEnd %q: %w", fields[2], err)
		}
		strand := ""
		if len(fields) > 5 {
			strand = fields[5]
		}
		return record{chrom: fields[0], start: start, end: end, strand: strand, fields: fields}, nil
	case fmtGFF:
		if len(fields) < 5 {
			return record{}, fmt.Errorf("GFF record must have at least 5 fields, got %d", len(fields))
		}
		gstart, err := parseChrPos(fields[3])
		if err != nil {
			return record{}, fmt.Errorf("invalid GFF start %q: %w", fields[3], err)
		}
		gend, err := parseChrPos(fields[4])
		if err != nil {
			return record{}, fmt.Errorf("invalid GFF end %q: %w", fields[4], err)
		}
		strand := ""
		if len(fields) > 6 {
			strand = fields[6]
		}
		return record{chrom: fields[0], start: gstart - 1, end: gend, strand: strand, fields: fields}, nil
	case fmtVCF:
		if len(fields) < 8 {
			return record{}, fmt.Errorf("VCF record must have at least 8 fields, got %d", len(fields))
		}
		pos, err := parseChrPos(fields[1])
		if err != nil {
			return record{}, fmt.Errorf("invalid VCF POS %q: %w", fields[1], err)
		}
		start := pos - 1
		end := start + vcfEnd(fields)
		return record{chrom: fields[0], start: start, end: end, fields: fields}, nil
	default:
		return record{}, fmt.Errorf("unknown text format")
	}
}

// vcfEnd returns the reference length of a VCF record's interval, mirroring
// upstream VcfRecord::initFromFile and SingleLineDelimTextFileReader::getVcfSVlen.
// For a symbolic ALT (e.g. <DEL>, <DUP>) that is not an insertion it derives the
// span from the INFO SVLEN (abs of the largest-magnitude value) or, failing
// that, END-POS+1. Insertions (<INS...>) are zero length. A plain ALT uses
// len(REF).
func vcfEnd(fields []string) int {
	ref := fields[3]
	alt := fields[4]
	if len(alt) > 0 && alt[0] == '<' {
		if strings.HasPrefix(alt, "<INS") {
			return 0
		}
		if l, ok := vcfSVLen(fields); ok {
			return l
		}
		// Fall through to REF length if no SV length info present.
	}
	return len(ref)
}

// vcfSVLen parses the INFO column (field 8, 0-based index 7) for SVLEN or END,
// returning the reference span. SVLEN is preferred (abs of the entry with the
// largest magnitude across a comma list); END yields END-POS+1. Returns ok=false
// when neither tag is present.
func vcfSVLen(fields []string) (int, bool) {
	if len(fields) < 8 {
		return 0, false
	}
	pos, err := parseChrPos(fields[1])
	if err != nil {
		return 0, false
	}
	for _, kv := range strings.Split(fields[7], ";") {
		if kv == "." {
			continue
		}
		eq := strings.IndexByte(kv, '=')
		if eq < 0 {
			continue
		}
		key, val := kv[:eq], kv[eq+1:]
		switch key {
		case "SVLEN":
			best := 0
			have := false
			for _, s := range strings.Split(val, ",") {
				n, err := strconv.Atoi(strings.TrimSpace(s))
				if err != nil {
					continue
				}
				if a := abs(n); !have || a > best {
					best = a
					have = true
				}
			}
			if have {
				return best, true
			}
		case "END":
			end, err := parseChrPos(val)
			if err == nil {
				return end - pos + 1, true
			}
		}
	}
	return 0, false
}

func abs(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

// parseChrPos parses a coordinate the way upstream str2chrPos does: leading
// signed digits, with a fallback to floating-point parsing (truncated toward
// zero) when scientific notation like "8e02" is encountered.
func parseChrPos(s string) (int, error) {
	s = strings.TrimSpace(s)
	if n, err := strconv.Atoi(s); err == nil {
		return n, nil
	}
	// Scientific notation / decimal fallback: upstream casts strtod to integer.
	if f, err := strconv.ParseFloat(s, 64); err == nil {
		return int(math.Trunc(f)), nil
	}
	return 0, fmt.Errorf("illegal number %q", s)
}

// parseChrPosBytes parses a coordinate from a byte slice with the same
// semantics as parseChrPos, but without allocating a string on the common
// integer fast path (plain optional-signed digits). Anything else — embedded
// whitespace, scientific notation, decimals — falls back to parseChrPos.
func parseChrPosBytes(b []byte) (int, error) {
	if n := len(b); n > 0 {
		i := 0
		neg := false
		if b[0] == '+' || b[0] == '-' {
			neg = b[0] == '-'
			i++
		}
		if i < n {
			val := 0
			ok := true
			for ; i < n; i++ {
				c := b[i]
				if c < '0' || c > '9' {
					ok = false
					break
				}
				val = val*10 + int(c-'0')
			}
			if ok {
				if neg {
					val = -val
				}
				return val, nil
			}
		}
	}
	return parseChrPos(string(b))
}

// isChrPos reports whether s parses as a coordinate (integer or scientific).
func isChrPos(s string) bool {
	_, err := parseChrPos(s)
	return err == nil
}

// readBAM decodes a BAM stream into records, converting each mapped primary
// alignment into the SAM-field layout that upstream BamRecord::getField exposes
// to -c column operations (1=QNAME, 2=FLAG[unsupported], 3=RNAME, 4=POS[0-based],
// 5=MAPQ, 6=CIGAR[bedtools OpLen format], 7=RNEXT, 8=PNEXT[0-based], 9=TLEN,
// 10=SEQ, 11=QUAL). The merge span is the reference footprint of the alignment.
// It returns ok=false (without error) when the stream is not BAM. Unmapped reads
// are skipped, matching upstream merge over BAM.
func readBAM(br *bufio.Reader) ([]record, bool, error) {
	reader, err := sam.NewReader(br)
	if err != nil {
		if err == sam.ErrNotBAM {
			return nil, false, nil
		}
		return nil, false, err
	}
	var out []record
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
		out = append(out, bamToRecord(rec))
	}
	return out, true, nil
}

// bamToRecord converts a single mapped alignment into a record whose fields
// follow upstream BamRecord::getField's SAM-field ordering. POS and PNEXT are
// rendered 0-based (BamTools' internal convention, which upstream column ops
// read verbatim). The strand comes from the reverse flag.
func bamToRecord(rec *sam.Record) record {
	start := int(rec.Pos) - 1     // BAM Pos is 1-based; merge span is 0-based.
	end := int(rec.EndPosition()) // 0-based half-open end = Pos + refLen - 1.
	strand := "+"
	if rec.Flag&sam.FlagReverse != 0 {
		strand = "-"
	}
	// BamTools' getName appends /1 or /2 for paired reads.
	name := rec.QName
	if rec.Flag&sam.FlagRead1 != 0 {
		name += "/1"
	} else if rec.Flag&sam.FlagRead2 != 0 {
		name += "/2"
	}
	mateChr := ""
	switch rec.RNext {
	case "=":
		mateChr = rec.RName
	case "*", "":
	default:
		mateChr = rec.RNext
	}
	matePos := "-1"
	if rec.PNext > 0 {
		matePos = strconv.Itoa(int(rec.PNext) - 1) // 0-based, BamTools convention.
	}
	// FLAG (field 2) is an empty placeholder: upstream errors if a column op
	// requests it, which the CLI enforces separately.
	fields := []string{
		name,
		"",
		rec.RName,
		strconv.Itoa(start), // POS, 0-based
		strconv.Itoa(int(rec.MapQ)),
		bamCigarStr(rec),
		mateChr,
		matePos,
		strconv.Itoa(int(rec.TLen)),
		rec.Seq,
		bamQualStr(rec),
	}
	return record{chrom: rec.RName, start: start, end: end, strand: strand, fields: fields, isBAM: true}
}

// bamCigarStr renders the CIGAR in upstream BamRecord::buildCigarStr's format:
// each op is "<Op><Length>" (e.g. 100M -> "M100", 3S97M -> "S3M97"), with no
// separator. An empty CIGAR yields "".
func bamCigarStr(rec *sam.Record) string {
	var b strings.Builder
	for _, op := range rec.Cigar {
		b.WriteByte(op.Char())
		b.WriteString(strconv.FormatUint(uint64(op.Length()), 10))
	}
	return b.String()
}

// bamQualStr renders QUAL as ASCII Phred+33, matching upstream. A missing
// quality string (empty or all 0xff) is rendered as "*".
func bamQualStr(rec *sam.Record) string {
	if len(rec.Qual) == 0 {
		return "*"
	}
	allFF := true
	for _, q := range rec.Qual {
		if q != 0xff {
			allFF = false
			break
		}
	}
	if allFF {
		return "*"
	}
	b := make([]byte, len(rec.Qual))
	for i, q := range rec.Qual {
		b[i] = q + 33
	}
	return string(b)
}

// WriteHeader copies leading comment/track/browser lines from r to w, stopping
// at the first data (non-header) line. It mirrors upstream `bedtools merge
// -header`, which echoes the input file's header before the merged output. BAM
// input has no text header to echo, so a BAM stream yields nothing.
func WriteHeader(r io.Reader, w io.Writer) error {
	br := bufio.NewReader(r)
	if magic, _ := br.Peek(4); len(magic) >= 4 &&
		(magic[0] == 0x1f && magic[1] == 0x8b || string(magic) == "BAM\x01") {
		return nil
	}
	sc := bufio.NewScanner(br)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		trimmed := strings.TrimSpace(strings.TrimRight(line, "\r"))
		if trimmed != "" && !strings.HasPrefix(trimmed, "#") &&
			!strings.HasPrefix(trimmed, "track") && !strings.HasPrefix(trimmed, "browser") {
			break // first data line
		}
		if _, err := fmt.Fprintln(w, strings.TrimRight(line, "\r")); err != nil {
			return err
		}
	}
	return sc.Err()
}

// sortRecords sorts records by (chrom, start) only, the order merging requires,
// and preserves input order on equal (chrom, start) keys.
//
// Upstream `bedtools merge` consumes an already-sorted stream and collects the
// records of each merged interval in the order it reads them. Its sort key is
// (chrom, start); chromEnd is NOT a tie-break. So for records with an equal
// (chrom, start), the order-sensitive aggregations (-o collapse, distinct, …)
// emit values in INPUT order. Using a stable start-only sort here reproduces
// that — an earlier end-ascending tie-break reordered the collapsed/distinct
// values relative to upstream. The merge itself only needs ascending starts;
// each group's end is taken as the max over the group, so dropping the end
// tie-break does not affect the merged coordinates.
func sortRecords(recs []record) {
	sort.SliceStable(recs, func(i, j int) bool {
		if recs[i].chrom != recs[j].chrom {
			return recs[i].chrom < recs[j].chrom
		}
		return recs[i].start < recs[j].start
	})
}
