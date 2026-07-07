// Package samtools — import subcommand (FASTQ → BAM/SAM).
//
// FastqImport reads one or more FASTQ files and emits BAM (or SAM text)
// alignment records with all reference fields blanked out, mirroring
// upstream's `samtools import`. The implementation follows
// `reference_code/samtools/bam_import.c::import_fastq`:
//
//   - Three positional shapes: a single file (single-ended or
//     interleaved /1 /2 records), a pair "R1 R2", or any combination of
//     `-0 unpaired`, `-1 R1 -2 R2`, and `-s singletons` (via OptionsPath
//     fields). The CLI driver maps short flags to these fields.
//   - Each FASTQ record becomes an unmapped (flag bit 0x4) record. R1 gets
//     0x1|0x40, R2 gets 0x1|0x80. If both R1 and R2 inputs are present, the
//     mate-unmapped bit 0x8 is set on both. Singletons get just 0x4.
//   - The FASTQ description text after the read ID is parsed for SAM aux
//     fields (TAG:TYPE:VALUE separated by whitespace) when AddTags is "*"
//     or a list of specific tags. Empty AddTags means no aux is parsed.
//   - The `/1` / `/2` suffix is stripped from the read name on the way
//     into BAM. With StripPairSuffix=false the suffix is preserved.
//   - When ReadGroup is non-empty an `@RG` header line is emitted and every
//     record gets a `RG:Z:<id>` aux tag.
//   - When OrderTag is non-empty (e.g. "oi"), each record gets a per-tag
//     incrementing counter aux ("oi:i:0", "oi:i:1", ...). The "TAG:LENGTH"
//     form requested by upstream is also recognised (`oi:5` → zero-padded
//     5-char string written as Z).
package samtools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// FastqImportOptions configures FastqImport.
type FastqImportOptions struct {
	// SinglePath is the `-s` positional shape: a single FASTQ that
	// contains both R1 and R2 (interleaved by /1 /2 suffix) or singletons.
	SinglePath string
	// UnpairedPath is the `-0` shape: a single unpaired FASTQ.
	UnpairedPath string
	// Read1Path / Read2Path are the `-1` / `-2` shapes; both must be set
	// to enable paired output.
	Read1Path string
	Read2Path string
	// Positional, when non-nil, is consumed before the named paths: a
	// single positional ends up as SinglePath, a pair "R1 R2" ends up as
	// Read1Path + Read2Path. Other lengths return an error.
	Positional []string
	// ReadGroup is the upstream `-R` shape: just the @RG ID. ReadGroupLine
	// is the `-r` shape: the full "@RG\t..." line. Either causes an @RG
	// header line to be emitted and an RG:Z:<id> aux to be appended to
	// every record.
	ReadGroup     string
	ReadGroupLine string
	// AuxTags governs FASTQ-comment → SAM-aux extraction. "*" means take
	// everything that looks like a TAG:TYPE:VALUE; a non-empty list
	// restricts to those tags. The empty value disables aux extraction.
	AuxTags string
	// OrderTag is the `--order` value. Passing "oi" appends an oi:i:<n>
	// counter to each record. The form "oi:N" (a non-zero N) emits a
	// zero-padded N-character oi:Z:<n> instead.
	OrderTag string
	// StripPairSuffix removes any trailing `/1` or `/2` from QNAME before
	// emission. The CLI default (matching upstream behaviour) is true; set
	// to false to match `-N` semantics, which keep the suffix.
	StripPairSuffix bool
	// OutputBAM (default) emits BGZF-wrapped BAM. Setting to false emits
	// text SAM.
	OutputBAM bool
	// Uncompressed implies OutputBAM and selects compression level 0 (stored
	// BGZF blocks) on the underlying writer (the -u flag). It takes precedence
	// over the Threads-driven level so -u output is level-0, and — because BGZF
	// blocks are independent gzip members — byte-identical for any worker count.
	Uncompressed bool
	// NoPG suppresses @PG injection (we never inject @PG so this is a
	// no-op).
	NoPG bool
	// Threads is upstream's -@/--threads worker count. When > 1 it spreads the
	// compressed-BAM OUTPUT's BGZF deflate across that many goroutines. Only
	// the output BGZF compression is parallelised: the FASTQ input decode
	// (plain-gzip DEFLATE, not block-parallel) and the record building stay
	// single-threaded, so the emitted bytes are identical for any worker count.
	// Ignored for SAM-text output, which has no BGZF framing to parallelise.
	Threads int
}

// FastqImport reads FASTQ records from the configured inputs and writes
// unmapped SAM/BAM records to out. The returned count is the number of
// records emitted.
func FastqImport(out io.Writer, opts FastqImportOptions) (int, error) {
	// Resolve positional arguments first.
	switch len(opts.Positional) {
	case 0:
		// nothing positional; named flags carry everything
	case 1:
		if opts.SinglePath == "" && opts.Read1Path == "" && opts.UnpairedPath == "" {
			opts.SinglePath = opts.Positional[0]
		}
	case 2:
		if opts.Read1Path == "" && opts.Read2Path == "" {
			opts.Read1Path = opts.Positional[0]
			opts.Read2Path = opts.Positional[1]
		}
	default:
		return 0, fmt.Errorf("samtools import: too many positional fastq files (%d); expected 0, 1, or 2", len(opts.Positional))
	}

	// Construct header: @HD with SO:unsorted, GO:query plus optional
	// @CO (matching upstream's "Reverse with: samtools fastq ..." line)
	// and optional @RG.
	hdr := &sam.Header{}
	hdr.Lines = append(hdr.Lines, sam.HeaderLine{Tag: "HD", Fields: []sam.HeaderField{
		{Tag: "VN", Value: "1.6"},
		{Tag: "SO", Value: "unsorted"},
		{Tag: "GO", Value: "query"},
	}})
	hdr.HDFields = hdr.Lines[len(hdr.Lines)-1].Fields
	if co := buildReverseComment(opts); co != "" {
		hdr.Lines = append(hdr.Lines, sam.HeaderLine{Tag: "CO", Fields: []sam.HeaderField{{Tag: "", Value: co}}})
		hdr.Comments = append(hdr.Comments, co)
	}
	rgID, rgLineText, rgErr := resolveRG(opts.ReadGroup, opts.ReadGroupLine)
	if rgErr != nil {
		return 0, rgErr
	}
	if rgLineText != "" {
		hl, err := parseRGHeaderLine(rgLineText)
		if err != nil {
			return 0, err
		}
		hdr.Lines = append(hdr.Lines, hl)
		rg := sam.ReadGroup{}
		for _, f := range hl.Fields {
			if f.Tag == "ID" {
				rg.ID = f.Value
			} else {
				rg.Extra = append(rg.Extra, f)
			}
		}
		hdr.ReadGroups = append(hdr.ReadGroups, rg)
	}

	var w sam.Writer
	if opts.OutputBAM || opts.Uncompressed {
		bw, err := sam.NewBAMWriterOptions(out, sam.BAMWriterOptions{Uncompressed: opts.Uncompressed, Threads: opts.Threads})
		if err != nil {
			return 0, err
		}
		w = bw
	} else {
		w = sam.NewSAMWriter(out)
	}
	if err := w.WriteHeader(hdr); err != nil {
		return 0, err
	}

	// Decide which reader path to walk.
	count := 0
	emitter := &recordEmitter{
		w:       w,
		opts:    &opts,
		rgID:    rgID,
		auxAll:  opts.AuxTags == "*",
		auxList: parseAuxTagList(opts.AuxTags),
		strip:   opts.StripPairSuffix,
	}
	switch {
	case opts.Read1Path != "" && opts.Read2Path != "":
		n, err := walkPaired(opts.Read1Path, opts.Read2Path, emitter)
		if err != nil {
			return n, err
		}
		count += n
	case opts.Read1Path != "" && opts.Read2Path == "":
		// Just R1 — treat as singletons that pre-claim FREAD1.
		n, err := walkSingle(opts.Read1Path, emitter, fqSingleR1)
		if err != nil {
			return n, err
		}
		count += n
	case opts.Read2Path != "":
		// Just R2 — pre-claim FREAD2.
		n, err := walkSingle(opts.Read2Path, emitter, fqSingleR2)
		if err != nil {
			return n, err
		}
		count += n
	case opts.UnpairedPath != "":
		n, err := walkSingle(opts.UnpairedPath, emitter, fqSingleR0)
		if err != nil {
			return n, err
		}
		count += n
	case opts.SinglePath != "":
		n, err := walkSingle(opts.SinglePath, emitter, fqSingleAuto)
		if err != nil {
			return n, err
		}
		count += n
	default:
		return 0, fmt.Errorf("samtools import: no input fastq files")
	}
	if err := w.Close(); err != nil {
		return count, err
	}
	return count, nil
}

// FastqImportFiles is the path-based entry point used by the CLI driver.
// It accepts opts already populated by the flag parser plus any leftover
// positional args.
func FastqImportFiles(positional []string, out io.Writer, opts FastqImportOptions) (int, error) {
	opts.Positional = positional
	return FastqImport(out, opts)
}

// fqSingleMode encodes how walkSingle should label a record from a
// single-file input.
type fqSingleMode int

const (
	fqSingleAuto fqSingleMode = iota // honour /1 /2 suffix in QNAME
	fqSingleR0                       // unpaired, just FUNMAP
	fqSingleR1                       // pre-claim FREAD1 (no mate-unmapped)
	fqSingleR2                       // pre-claim FREAD2 (no mate-unmapped)
)

// recordEmitter wraps a writer with the per-record decoration logic.
type recordEmitter struct {
	w       sam.Writer
	opts    *FastqImportOptions
	rgID    string
	auxAll  bool
	auxList map[string]struct{}
	strip   bool
	orderN  uint64
}

// emit writes one record after applying flag, name-suffix, aux-parsing, RG,
// and --order policies. The mate set is applied by the caller (paired vs
// singleton flag flips happen there).
func (e *recordEmitter) emit(name, desc string, seq, qual []byte, flag uint16) error {
	qname := name
	if e.strip {
		qname = stripPairSuffix(qname)
	}
	rec := &sam.Record{
		QName: qname,
		Flag:  flag,
		RName: "",
		Pos:   0,
		MapQ:  0,
		Cigar: nil,
		RNext: "",
		PNext: 0,
		TLen:  0,
		Seq:   string(seq),
		// SAM quality is stored as raw Phred; FASTQ uses ASCII+33.
		Qual: rawQualityFromFastq(qual),
	}
	if e.auxAll || len(e.auxList) > 0 {
		auxes := parseFastqAux(desc, e.auxAll, e.auxList)
		rec.Aux = append(rec.Aux, auxes...)
	}
	if e.rgID != "" {
		rec.Aux = append(rec.Aux, sam.Aux{Tag: "RG", Type: 'Z', Value: e.rgID})
	}
	if e.opts.OrderTag != "" {
		tag, width := parseOrderTag(e.opts.OrderTag)
		if width > 0 {
			val := fmt.Sprintf("%0*d", width, e.orderN)
			rec.Aux = append(rec.Aux, sam.Aux{Tag: tag, Type: 'Z', Value: val})
		} else {
			rec.Aux = append(rec.Aux, sam.Aux{Tag: tag, Type: 'i', Value: int64(e.orderN)})
		}
		e.orderN++
	}
	return e.w.Write(rec)
}

// walkSingle reads every FASTQ record from path and emits it via e according
// to the requested mode.
func walkSingle(path string, e *recordEmitter, mode fqSingleMode) (int, error) {
	rc, err := iohelper.OpenReader(path)
	if err != nil {
		return 0, fmt.Errorf("samtools import: open %s: %w", path, err)
	}
	defer rc.Close()
	rd := newFastqLineReader(rc)
	n := 0
	for {
		rec, ferr := rd.next()
		if ferr == io.EOF {
			break
		}
		if ferr != nil {
			return n, fmt.Errorf("samtools import: read %s: %w", path, ferr)
		}
		// Layer 1 — the htslib FASTQ reader flags (htslib sam.c
		// fastq_parse1): a read name ending in `/<digit>` is treated as
		// mate-of-a-pair, so it gets FUNMAP|FMUNMAP|FPAIRED plus FREAD1
		// (`/1`), FREAD2 (`/2`), or both (any other digit). A name without
		// the suffix is just FUNMAP. These flags are applied uniformly for
		// every input shape because upstream funnels all of -0/-s/-1/-2
		// through the same reader; the suffix is then stripped (see emit).
		flag := pairSuffixFlag(rec.name)
		// Layer 2 — the bam_import.c per-file-role flags
		// (bam_import.c:349). FQ_R0 (-0) and FQ_SINGLE (-s) add nothing
		// beyond the reader's suffix flags. FQ_R1 (-1) force-claims FREAD1
		// (when neither read bit is set yet) and FPAIRED; FQ_R2 (-2)
		// force-claims FREAD2 and FPAIRED. The neighbour-driven FMUNMAP that
		// upstream sets when an R1 is immediately followed by an R2 only
		// applies when both files are present, which is handled by
		// walkPaired; a lone -1/-2 here never gains it.
		switch mode {
		case fqSingleR0, fqSingleAuto:
			// reader flags only.
		case fqSingleR1:
			if flag&(sam.FlagRead1|sam.FlagRead2) == 0 {
				flag |= sam.FlagRead1
			}
			flag |= sam.FlagPaired
		case fqSingleR2:
			flag |= sam.FlagPaired | sam.FlagRead2
		}
		if err := e.emit(rec.name, rec.desc, rec.seq, rec.qual, flag); err != nil {
			return n, err
		}
		n++
	}
	return n, nil
}

// walkPaired reads two FASTQs in lock-step and emits R1 then R2 per pair.
// Both records get FPAIRED | FUNMAP | FMUNMAP (mate is unmapped too) plus
// FREAD1 or FREAD2 as appropriate.
func walkPaired(r1Path, r2Path string, e *recordEmitter) (int, error) {
	r1, err := iohelper.OpenReader(r1Path)
	if err != nil {
		return 0, fmt.Errorf("samtools import: open %s: %w", r1Path, err)
	}
	defer r1.Close()
	r2, err := iohelper.OpenReader(r2Path)
	if err != nil {
		return 0, fmt.Errorf("samtools import: open %s: %w", r2Path, err)
	}
	defer r2.Close()
	rd1 := newFastqLineReader(r1)
	rd2 := newFastqLineReader(r2)
	n := 0
	for {
		rec1, err1 := rd1.next()
		rec2, err2 := rd2.next()
		if err1 == io.EOF && err2 == io.EOF {
			break
		}
		if err1 != nil && err1 != io.EOF {
			return n, fmt.Errorf("samtools import: read %s: %w", r1Path, err1)
		}
		if err2 != nil && err2 != io.EOF {
			return n, fmt.Errorf("samtools import: read %s: %w", r2Path, err2)
		}
		if err1 != err2 {
			return n, fmt.Errorf("samtools import: input files with differing number of records (%s vs %s)", r1Path, r2Path)
		}
		flag1 := uint16(sam.FlagPaired | sam.FlagUnmapped | sam.FlagMateUnmapped | sam.FlagRead1)
		flag2 := uint16(sam.FlagPaired | sam.FlagUnmapped | sam.FlagMateUnmapped | sam.FlagRead2)
		if err := e.emit(rec1.name, rec1.desc, rec1.seq, rec1.qual, flag1); err != nil {
			return n, err
		}
		n++
		if err := e.emit(rec2.name, rec2.desc, rec2.seq, rec2.qual, flag2); err != nil {
			return n, err
		}
		n++
	}
	return n, nil
}

// fastqRecord is a lightweight FASTQ record. We don't use pkg/htsgo/fastq
// because we need access to the raw description tail (everything after the
// first whitespace) for SAM-aux parsing, including embedded tabs.
type fastqRecord struct {
	name string
	desc string
	seq  []byte
	qual []byte
}

// fastqLineReader wraps a bufio.Reader for record-at-a-time consumption.
type fastqLineReader struct {
	br *bufio.Reader
}

// newFastqLineReader wraps an arbitrary reader.
func newFastqLineReader(r io.Reader) *fastqLineReader {
	return &fastqLineReader{br: bufio.NewReaderSize(r, 1<<16)}
}

// next reads one record (4 lines). Trailing CRLF is stripped.
func (fr *fastqLineReader) next() (*fastqRecord, error) {
	hdr, err := fr.readLine()
	if err != nil {
		return nil, err
	}
	if hdr == "" {
		return nil, io.EOF
	}
	if hdr[0] != '@' {
		return nil, fmt.Errorf("expected '@' at start of FASTQ header, got %q", hdr)
	}
	rest := hdr[1:]
	name := rest
	desc := ""
	// Split on the first run of whitespace (space or tab).
	for i := 0; i < len(rest); i++ {
		if rest[i] == ' ' || rest[i] == '\t' {
			name = rest[:i]
			desc = rest[i+1:]
			break
		}
	}
	seqLine, err := fr.readLine()
	if err != nil {
		return nil, fmt.Errorf("samtools import: truncated record (no seq): %w", err)
	}
	sep, err := fr.readLine()
	if err != nil {
		return nil, fmt.Errorf("samtools import: truncated record (no sep): %w", err)
	}
	if len(sep) == 0 || sep[0] != '+' {
		return nil, fmt.Errorf("samtools import: bad separator line %q", sep)
	}
	qualLine, err := fr.readLine()
	if err != nil {
		return nil, fmt.Errorf("samtools import: truncated record (no qual): %w", err)
	}
	if len(seqLine) != len(qualLine) {
		return nil, fmt.Errorf("samtools import: seq/qual length mismatch (%d vs %d)", len(seqLine), len(qualLine))
	}
	return &fastqRecord{
		name: name,
		desc: desc,
		seq:  []byte(seqLine),
		qual: []byte(qualLine),
	}, nil
}

// readLine consumes one newline-terminated line. The trailing '\n' (and
// any '\r' before it) is stripped. An EOF on an empty buffer returns
// io.EOF; otherwise the last line without a trailing newline is returned.
func (fr *fastqLineReader) readLine() (string, error) {
	line, err := fr.br.ReadString('\n')
	if err == io.EOF && line == "" {
		return "", io.EOF
	}
	if err != nil && err != io.EOF {
		return "", err
	}
	line = strings.TrimRight(line, "\r\n")
	return line, nil
}

// stripPairSuffix drops a trailing "/<digit>" from a read name, matching the
// upstream htslib FASTQ reader which strips any `/0`–`/9` mate suffix. The
// name must be longer than two characters (htslib requires name.l > 2) so a
// bare "/1" is left untouched.
func stripPairSuffix(name string) string {
	if len(name) > 2 && name[len(name)-2] == '/' {
		c := name[len(name)-1]
		if c >= '0' && c <= '9' {
			return name[:len(name)-2]
		}
	}
	return name
}

// pairSuffixFlag returns the SAM flag bits the upstream htslib FASTQ reader
// (htslib sam.c fastq_parse1) assigns to a read based solely on its name
// suffix. A name longer than two characters ending in `/<digit>` is treated
// as one mate of a pair: it always gets FUNMAP|FMUNMAP|FPAIRED, plus FREAD1
// for `/1`, FREAD2 for `/2`, or both for any other digit. A name without the
// suffix gets just FUNMAP.
func pairSuffixFlag(name string) uint16 {
	flag := uint16(sam.FlagUnmapped)
	if len(name) > 2 && name[len(name)-2] == '/' {
		c := name[len(name)-1]
		if c >= '0' && c <= '9' {
			flag |= sam.FlagMateUnmapped | sam.FlagPaired
			switch c {
			case '1':
				flag |= sam.FlagRead1
			case '2':
				flag |= sam.FlagRead2
			default:
				flag |= sam.FlagRead1 | sam.FlagRead2
			}
		}
	}
	return flag
}

// rawQualityFromFastq converts an ASCII+33-encoded quality string to raw
// Phred bytes. A nil/empty input produces a nil slice (which the SAM writer
// then represents as "*").
func rawQualityFromFastq(qual []byte) []byte {
	if len(qual) == 0 {
		return nil
	}
	out := make([]byte, len(qual))
	for i, q := range qual {
		if q < 33 {
			out[i] = 0
		} else {
			out[i] = q - 33
		}
	}
	return out
}

// parseAuxTagList parses a comma-separated tag list into a set. The "*"
// (all tags) form is signalled separately by the caller and returns a nil
// map.
func parseAuxTagList(s string) map[string]struct{} {
	if s == "" || s == "*" {
		return nil
	}
	set := make(map[string]struct{})
	for _, t := range strings.Split(s, ",") {
		t = strings.TrimSpace(t)
		if len(t) == 2 {
			set[t] = struct{}{}
		}
	}
	return set
}

// parseFastqAux extracts SAM aux fields from a FASTQ description string. The
// description tail upstream allows is `TAG:TYPE:VALUE` records separated by
// whitespace (tab preferred, space tolerated). If allowAll is true every
// parseable field is kept; otherwise only those in allowList are.
func parseFastqAux(desc string, allowAll bool, allowList map[string]struct{}) []sam.Aux {
	if desc == "" {
		return nil
	}
	var out []sam.Aux
	// Tokenise on any run of whitespace; tab is the common separator but
	// space is also legal. We avoid strings.Fields because the upstream
	// parser also tolerates a multi-tab run.
	start := 0
	for i := 0; i <= len(desc); i++ {
		if i == len(desc) || desc[i] == '\t' || desc[i] == ' ' {
			if i > start {
				tok := desc[start:i]
				if a, ok := tryParseAuxField(tok); ok {
					if allowAll {
						out = append(out, a)
					} else if _, kept := allowList[a.Tag]; kept {
						out = append(out, a)
					}
				}
			}
			start = i + 1
		}
	}
	return out
}

// tryParseAuxField parses a single TAG:TYPE:VALUE token, returning false
// silently for malformed tokens (upstream behaviour: skip unknown comment
// junk rather than abort the import).
func tryParseAuxField(token string) (sam.Aux, bool) {
	if len(token) < 5 || token[2] != ':' || token[4] != ':' {
		return sam.Aux{}, false
	}
	a, err := sam.ParseAux(token)
	if err != nil {
		return sam.Aux{}, false
	}
	return a, true
}

// buildReverseComment constructs the "@CO Reverse with: samtools fastq ..."
// hint upstream emits. It is purely diagnostic — callers that want to roun
// d-trip can paste it into a `samtools fastq` command.
func buildReverseComment(opts FastqImportOptions) string {
	parts := []string{}
	if opts.UnpairedPath != "" {
		parts = append(parts, "-0 unpaired.fastq")
	}
	if opts.Read1Path != "" {
		parts = append(parts, "-1 R1.fastq")
	}
	if opts.Read2Path != "" {
		parts = append(parts, "-2 R2.fastq")
	}
	if opts.SinglePath != "" && opts.Read1Path == "" && opts.Read2Path == "" && opts.UnpairedPath == "" {
		parts = append(parts, "-N -o paired.fastq")
	}
	if len(opts.Positional) == 2 && opts.Read1Path == "" && opts.Read2Path == "" {
		parts = []string{"-1 R1.fastq", "-2 R2.fastq"}
	}
	if len(opts.Positional) == 1 && opts.SinglePath == "" && opts.Read1Path == "" && opts.UnpairedPath == "" {
		parts = []string{"-N -o paired.fastq"}
	}
	if len(parts) == 0 {
		return ""
	}
	// Upstream emits a trailing space before "\n"; preserve the quirk so
	// expected.sam fixtures diff byte-for-byte.
	return "Reverse with: samtools fastq " + strings.Join(parts, " ") + " "
}

// resolveRG handles the upstream -R / -r flag pair. -R is a bare ID;
// -r is a full @RG line ("ID:foo" or "@RG\tID:foo"). Returns the ID and
// the canonical header line text (already starting with "@RG\t").
func resolveRG(shortID, line string) (string, string, error) {
	if shortID == "" && line == "" {
		return "", "", nil
	}
	if line == "" {
		return shortID, "@RG\tID:" + shortID, nil
	}
	full := line
	if !strings.HasPrefix(full, "@") {
		full = "@RG\t" + full
	}
	// Extract ID.
	idx := strings.Index(full, "\tID:")
	if idx < 0 {
		return "", "", fmt.Errorf("samtools import: -r RG-line missing ID: field: %q", line)
	}
	rest := full[idx+4:]
	end := strings.IndexAny(rest, "\t")
	if end < 0 {
		end = len(rest)
	}
	id := rest[:end]
	return id, full, nil
}

// parseRGHeaderLine parses an @RG line into a HeaderLine, preserving the
// tag order.
func parseRGHeaderLine(line string) (sam.HeaderLine, error) {
	if !strings.HasPrefix(line, "@RG") {
		return sam.HeaderLine{}, fmt.Errorf("samtools import: expected @RG line, got %q", line)
	}
	rest := strings.TrimPrefix(line, "@RG")
	rest = strings.TrimPrefix(rest, "\t")
	hl := sam.HeaderLine{Tag: "RG"}
	for _, part := range strings.Split(rest, "\t") {
		if len(part) < 3 || part[2] != ':' {
			continue
		}
		hl.Fields = append(hl.Fields, sam.HeaderField{Tag: part[:2], Value: part[3:]})
	}
	return hl, nil
}

// parseOrderTag splits an --order spec into the 2-char tag and an optional
// zero-pad width (the upstream "TAG:N" form).
func parseOrderTag(s string) (string, int) {
	if len(s) < 2 {
		return s, 0
	}
	if len(s) == 2 {
		return s, 0
	}
	if s[2] != ':' {
		return s[:2], 0
	}
	n, err := strconv.Atoi(s[3:])
	if err != nil || n <= 0 {
		return s[:2], 0
	}
	return s[:2], n
}

// openImportStdout is a tiny helper used by the CLI driver to pick between
// stdout and a file path. Mirrors `openOut` in the cmd package without
// pulling the package dependency.
func openImportStdout(path string) (io.WriteCloser, error) {
	if path == "" || path == "-" {
		return nopWriteCloser{Writer: os.Stdout}, nil
	}
	return os.Create(path)
}

// nopWriteCloser wraps an io.Writer with a no-op Close.
type nopWriteCloser struct {
	io.Writer
}

// Close is a no-op so stdout isn't accidentally closed.
func (nopWriteCloser) Close() error { return nil }
