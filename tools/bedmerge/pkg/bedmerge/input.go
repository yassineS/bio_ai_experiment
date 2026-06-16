// Input-format handling for bedmerge. Upstream `bedtools merge` accepts BED,
// GFF, VCF, and BAM input on -i (auto-detected) and merges them. This file
// converts each of those formats into a common interval record that preserves
// the original input columns (so -c/-o column operations can address them) and
// records the 0-based half-open span used for merging, mirroring upstream's
// per-record-type coordinate conventions.
package bedmerge

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
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

// readInput reads every interval from r, auto-detecting BAM (by sniffing the
// leading magic bytes) or a text format (BED/GFF/VCF). It returns the parsed
// records and the detected format. Header/track/browser/blank lines are skipped
// for text input. The text format is locked on the first data line, matching
// upstream's per-file type detection.
func readInput(r io.Reader, opts MergeOptions) ([]record, inputFormat, error) {
	br := bufio.NewReaderSize(r, 64*1024)
	magic, _ := br.Peek(4)
	if len(magic) >= 4 && (magic[0] == 0x1f && magic[1] == 0x8b || string(magic) == "BAM\x01") {
		recs, ok, err := readBAM(br)
		if err != nil {
			return nil, fmtBAM, err
		}
		if ok {
			return recs, fmtBAM, nil
		}
		// Gzip magic but not BAM: iohelper.OpenReader already transparently
		// decompresses plain-gzip text, so a still-compressed non-BAM stream is
		// unexpected here. Fall through to the text reader.
	}
	recs, f, err := readText(br, opts)
	return recs, f, err
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

// readText reads BED/GFF/VCF records from an already-peeked reader.
func readText(br *bufio.Reader, opts MergeOptions) ([]record, inputFormat, error) {
	sc := bufio.NewScanner(br)
	sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	var out []record
	var format inputFormat
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
				f, ok := detectFormat(fields)
				if !ok {
					return nil, format, fmt.Errorf("unexpected file format: please use tab-delimited BED, GFF, or VCF (line: %q)", line)
				}
				format = f
			}
			// -s is not supported for VCF (matching upstream merge.t11).
			if format == fmtVCF && opts.StrandSpec {
				return nil, format, errStrandedVCF
			}
			formatSet = true
		}
		rec, err := parseTextRecord(fields, format)
		if err != nil {
			return nil, format, err
		}
		out = append(out, rec)
	}
	if err := sc.Err(); err != nil {
		return nil, format, err
	}
	return out, format, nil
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

// sortRecords sorts records by (chrom, start, end), the order merging requires.
func sortRecords(recs []record) {
	sort.SliceStable(recs, func(i, j int) bool {
		if recs[i].chrom != recs[j].chrom {
			return recs[i].chrom < recs[j].chrom
		}
		if recs[i].start != recs[j].start {
			return recs[i].start < recs[j].start
		}
		return recs[i].end < recs[j].end
	})
}
