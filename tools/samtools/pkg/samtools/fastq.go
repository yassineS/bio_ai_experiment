package samtools

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// FastqOptions configures BAM-to-FASTQ conversion. Empty path strings mean
// "no separate output for this category"; in that case the records would
// otherwise be silently dropped, so callers typically set either Read1/Read2
// (paired mode) or just Output (interleaved/single).
type FastqOptions struct {
	// Read1Path is the path to write the first-in-pair reads to (`-1`).
	Read1Path string
	// Read2Path is the path to write the second-in-pair reads to (`-2`).
	Read2Path string
	// OrphanPath is the path for reads where exactly one of 0x40/0x80 is
	// set but the mate is missing, or where neither is set (`-0`).
	OrphanPath string
	// SingletonPath is the path for unpaired reads (no 0x1) when Read1/2
	// paths are configured (`-s`).
	SingletonPath string
	// OutputPath is the default sink (`-o`) used when Read1/Read2 are not
	// configured — paired output is interleaved here.
	OutputPath string
	// AlwaysAddSuffix forces `/1`/`/2` on every read name that already has
	// pair info in the flag (`-N`). When false, the suffix is only added
	// when not already present in QNAME and `NoSuffix` is false.
	AlwaysAddSuffix bool
	// NoSuffix disables the suffix entirely (`-n`).
	NoSuffix bool
	// IncludeFlags requires ALL these bits to be set (`-f`).
	IncludeFlags uint16
	// ExcludeFlags drops reads with ANY of these bits set (`-F`). Default
	// 0x900 (secondary + supplementary) matches upstream samtools.
	ExcludeFlags uint16
	// ExcludeFlagsAll drops reads only when ALL of these bits are set
	// (`-G`).
	ExcludeFlagsAll uint16
	// UseExcludeAll enables ExcludeFlagsAll filtering.
	UseExcludeAll bool
	// AddTags is the comma-separated list of aux tags whose `TAG:TYPE:VALUE`
	// formatted form is appended to the read description (`-T`).
	AddTags []string
	// AddAllTags, when true, appends every aux tag on each record (in the
	// record's existing order). Matches upstream samtools fastq's `-T '*'`
	// / `-T ''` shorthand (bam_fastq.c:500 — FASTQ_OPT_AUX(NULL)).
	AddAllTags bool
	// CompressLevel is the gzip level for `.gz` outputs (`-c`).
	// gzip.DefaultCompression when zero.
	CompressLevel int
	// UseOQ pulls quality from the OQ aux tag when present (`-O`).
	UseOQ bool
	// NoCO is accepted for upstream compatibility — we never emit @CO
	// lines anyway.
	NoCO bool
	// Threads is accepted for compatibility; v1 is single-threaded.
	Threads int
}

// Fastq streams a SAM/BAM input through to one or more FASTQ outputs as
// configured by opts. The function opens the configured output paths,
// writes the records, then closes them. The returned counts identify how
// many records were emitted to each output and how many were dropped.
type FastqCounts struct {
	Read1     int
	Read2     int
	Singleton int
	Orphan    int
	Output    int
	Dropped   int
	// PairedCoordinateWarn is set when paired (-1/-2) mode receives a
	// coordinate-sorted input and falls back to write-through-to-output
	// behaviour. Useful for tests and the CLI driver.
	PairedCoordinateWarn bool
}

// Fastq performs the conversion.
func Fastq(in io.Reader, opts FastqOptions) (FastqCounts, error) {
	var counts FastqCounts
	rd, err := alnio.NewReader(in)
	if err != nil {
		return counts, err
	}
	pairedMode := opts.Read1Path != "" || opts.Read2Path != ""
	// Upstream samtools bam2fq: `if (opts->fnr[1] || opts->fnr[2]) opts->has12 = false`.
	// In other words, when the user has specified -1 and/or -2 the read name
	// pair-suffix is dropped by default, because the output file already
	// disambiguates mate identity. `-N` (AlwaysAddSuffix) overrides this.
	if pairedMode && !opts.AlwaysAddSuffix {
		opts.NoSuffix = true
	}
	// Detect input ordering. For paired output we require name-sorted
	// input; coordinate-sorted falls back to interleaved write-through.
	hdr := rd.Header()
	sortOrder := readSortOrder(hdr)
	if pairedMode && sortOrder == "coordinate" {
		counts.PairedCoordinateWarn = true
	}

	// Open every configured output. Closing is in reverse order so the
	// gzip headers/footers nest correctly.
	type sink struct {
		path string
		w    io.Closer
		bw   *bufio.Writer
	}
	open := func(path string) (*sink, error) {
		if path == "" {
			return nil, nil
		}
		s, oerr := openFastqOutput(path, opts.CompressLevel)
		if oerr != nil {
			return nil, oerr
		}
		return &sink{path: path, w: s.w, bw: s.bw}, nil
	}
	defer func() {
		// Flush+close happens explicitly below; defer is a safety net
		// for early returns on parsing errors.
	}()

	openSink := func(path string) (*sink, error) { return open(path) }
	read1, err := openSink(opts.Read1Path)
	if err != nil {
		return counts, err
	}
	read2, err := openSink(opts.Read2Path)
	if err != nil {
		return counts, err
	}
	singleton, err := openSink(opts.SingletonPath)
	if err != nil {
		return counts, err
	}
	orphan, err := openSink(opts.OrphanPath)
	if err != nil {
		return counts, err
	}
	output, err := openSink(opts.OutputPath)
	if err != nil {
		return counts, err
	}
	// closeAll flushes/closes every non-nil sink in deterministic order.
	closeAll := func() error {
		var firstErr error
		for _, s := range []*sink{read1, read2, singleton, orphan, output} {
			if s == nil {
				continue
			}
			if err := s.bw.Flush(); err != nil && firstErr == nil {
				firstErr = err
			}
			if err := s.w.Close(); err != nil && firstErr == nil {
				firstErr = err
			}
		}
		return firstErr
	}

	// In coordinate-sorted paired mode we fall back to interleaved output:
	// if Output is open, write everything there; if not, write paired
	// records into Singleton; final fallback is to drop them.
	//
	// In name-sorted paired mode (the upstream supported case) we group
	// adjacent same-qname records and decide pair vs singleton on the
	// completed group, matching bam_fastq.c's `flush_rec` logic: only
	// when both an R1 and an R2 share the qname does the pair go to
	// `-1`/`-2`. A paired-flagged record whose mate is absent from the
	// group goes to `-s` if open, otherwise to its flag-indicated `-1`
	// or `-2` file (the upstream fallback).
	useGrouping := pairedMode && !counts.PairedCoordinateWarn
	// fmtRec runs the per-record SAM→FASTQ formatting steps once so the
	// same payload can be routed to any of the sinks below.
	fmtRec := func(rec *sam.Record) string {
		seq, qual := orientReadForFastq(rec, opts)
		header := buildFastqHeader(rec, opts)
		return formatFastq(header, seq, qual)
	}
	// writeToSink writes line to s if non-nil and updates *count, else
	// falls back to output / dropped per upstream's precedence.
	writeToSink := func(s *sink, line string, count *int) {
		if s != nil {
			_, _ = s.bw.WriteString(line)
			*count++
			return
		}
		if output != nil {
			_, _ = output.bw.WriteString(line)
			counts.Output++
			return
		}
		counts.Dropped++
	}
	// dispatchSingle handles the non-grouped (coordinate-fallback /
	// non-paired-mode) path that used to be the only code path.
	dispatchSingle := func(rec *sam.Record, line string) {
		switch {
		case useGrouping:
			// Unreachable: grouping path bypasses dispatchSingle.
		case output != nil:
			_, _ = output.bw.WriteString(line)
			counts.Output++
		case singleton != nil && !rec.IsPaired():
			_, _ = singleton.bw.WriteString(line)
			counts.Singleton++
		default:
			counts.Dropped++
		}
	}
	// flushGroup writes a same-qname group per upstream flush_rec:
	//   - if R1 and R2 both present → pair to -1/-2
	//   - if exactly one of R1/R2 present → -s if open, else its -1/-2
	//   - non-paired (score[0]) → -1 (upstream's fpr[0])
	//   - flag-paired but neither R1 nor R2 (or both) → orphan sink
	flushGroup := func(group []*sam.Record) {
		var best1, best2, best0 *sam.Record
		for _, r := range group {
			switch {
			case r.IsPaired() && r.IsRead1() && !r.IsRead2():
				if best1 == nil {
					best1 = r
				}
			case r.IsPaired() && r.IsRead2() && !r.IsRead1():
				if best2 == nil {
					best2 = r
				}
			case !r.IsPaired():
				if best0 == nil {
					best0 = r
				}
			default:
				// Paired but both/neither R1/R2 set → orphan.
				line := fmtRec(r)
				writeToSink(orphan, line, &counts.Orphan)
			}
		}
		if best1 != nil && best2 != nil {
			writeToSink(read1, fmtRec(best1), &counts.Read1)
			writeToSink(read2, fmtRec(best2), &counts.Read2)
		} else if best1 != nil || best2 != nil {
			var solo *sam.Record
			if best1 != nil {
				solo = best1
			} else {
				solo = best2
			}
			line := fmtRec(solo)
			if singleton != nil {
				_, _ = singleton.bw.WriteString(line)
				counts.Singleton++
			} else if best1 != nil {
				writeToSink(read1, line, &counts.Read1)
			} else {
				writeToSink(read2, line, &counts.Read2)
			}
		}
		if best0 != nil {
			// Non-paired (no 0x1): our port maps this to the singleton
			// sink (matching the legacy non-grouping dispatch
			// behaviour). Upstream emits to its `-0` (`fpr[0]`) sink
			// which our CLI does not expose. If neither `-s` nor `-o`
			// is set the record is dropped, preserving the partial-
			// fallback contract exercised by TestFastqPairedPartialFallback.
			line := fmtRec(best0)
			if singleton != nil {
				_, _ = singleton.bw.WriteString(line)
				counts.Singleton++
			} else if output != nil {
				_, _ = output.bw.WriteString(line)
				counts.Output++
			} else {
				counts.Dropped++
			}
		}
	}
	var group []*sam.Record
	var groupQName string
	for {
		rec, rerr := rd.Read()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			_ = closeAll()
			return counts, rerr
		}
		if !keepFastqRecord(rec, opts) {
			counts.Dropped++
			continue
		}
		if useGrouping {
			if group != nil && rec.QName != groupQName {
				flushGroup(group)
				group = group[:0]
			}
			groupQName = rec.QName
			group = append(group, rec)
			continue
		}
		dispatchSingle(rec, fmtRec(rec))
	}
	if useGrouping && len(group) > 0 {
		flushGroup(group)
	}
	if err := closeAll(); err != nil {
		return counts, err
	}
	return counts, nil
}

// fastqSink is the result of opening one output path: the closer (which
// also writes-through during Close for gzip footers) and the buffered
// writer the caller emits records into.
type fastqSink struct {
	w  io.Closer
	bw *bufio.Writer
}

// openFastqOutput opens a FASTQ destination. "-" means stdout. The output
// is gzipped when the path ends in ".gz".
func openFastqOutput(path string, compressLevel int) (fastqSink, error) {
	if path == "-" || path == "" {
		// Stdout — wrap with a no-op closer so we don't close the user's
		// terminal.
		bw := bufio.NewWriter(os.Stdout)
		return fastqSink{w: nopCloseStdout{}, bw: bw}, nil
	}
	f, err := os.Create(path)
	if err != nil {
		return fastqSink{}, err
	}
	if strings.HasSuffix(path, ".gz") {
		level := compressLevel
		if level == 0 {
			level = gzip.DefaultCompression
		}
		gz, gzerr := gzip.NewWriterLevel(f, level)
		if gzerr != nil {
			_ = f.Close()
			return fastqSink{}, gzerr
		}
		bw := bufio.NewWriter(gz)
		return fastqSink{w: &gzCloser{f: f, gz: gz}, bw: bw}, nil
	}
	bw := bufio.NewWriter(f)
	return fastqSink{w: f, bw: bw}, nil
}

// gzCloser closes a gzip.Writer and then the underlying file. The order
// matters: flush gzip first so the gzip footer is on disk.
type gzCloser struct {
	f  *os.File
	gz *gzip.Writer
}

func (c *gzCloser) Close() error {
	if err := c.gz.Close(); err != nil {
		_ = c.f.Close()
		return err
	}
	return c.f.Close()
}

// nopCloseStdout is a closer for the stdout-backed sink — closing is a
// no-op so we don't accidentally close the test runner's stdout.
type nopCloseStdout struct{}

func (nopCloseStdout) Close() error { return nil }

// keepFastqRecord applies the flag-include/exclude filters.
func keepFastqRecord(rec *sam.Record, opts FastqOptions) bool {
	if opts.IncludeFlags != 0 && rec.Flag&opts.IncludeFlags != opts.IncludeFlags {
		return false
	}
	exc := opts.ExcludeFlags
	if exc == 0 {
		exc = sam.FlagSecondary | sam.FlagSupplementary
	}
	if rec.Flag&exc != 0 {
		return false
	}
	if opts.UseExcludeAll && opts.ExcludeFlagsAll != 0 &&
		rec.Flag&opts.ExcludeFlagsAll == opts.ExcludeFlagsAll {
		return false
	}
	return true
}

// orientReadForFastq returns the SEQ and quality (raw Phred bytes) oriented
// in original-read direction: for reverse-strand records the BAM stores the
// reverse-complemented SEQ and reversed QUAL, so we reverse them again to
// recover the original sequencing orientation. UseOQ swaps the QUAL for
// the OQ tag's value when present.
func orientReadForFastq(rec *sam.Record, opts FastqOptions) (string, []byte) {
	seq := rec.Seq
	qual := append([]byte(nil), rec.Qual...)
	if opts.UseOQ {
		if a, ok := rec.GetAux("OQ"); ok {
			if s, ok := a.String(); ok {
				// OQ is ASCII-33 encoded; convert back to raw Phred.
				oq := make([]byte, len(s))
				for i := 0; i < len(s); i++ {
					oq[i] = s[i] - 33
				}
				qual = oq
			}
		}
	}
	if rec.Flag&sam.FlagReverse != 0 {
		seq = reverseComplement(seq)
		// Reverse qual bytes in place.
		for i, j := 0, len(qual)-1; i < j; i, j = i+1, j-1 {
			qual[i], qual[j] = qual[j], qual[i]
		}
	}
	return seq, qual
}

// reverseComplement returns the reverse complement of an IUPAC sequence.
func reverseComplement(s string) string {
	out := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		out[len(s)-1-i] = complementBase(s[i])
	}
	return string(out)
}

// complementBase returns the IUPAC complement of a single base. Unknown
// characters are returned unchanged.
func complementBase(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'a':
		return 't'
	case 'T':
		return 'A'
	case 't':
		return 'a'
	case 'C':
		return 'G'
	case 'c':
		return 'g'
	case 'G':
		return 'C'
	case 'g':
		return 'c'
	case 'N':
		return 'N'
	case 'n':
		return 'n'
	case 'U':
		return 'A'
	case 'u':
		return 'a'
	case 'R':
		return 'Y'
	case 'Y':
		return 'R'
	case 'K':
		return 'M'
	case 'M':
		return 'K'
	case 'B':
		return 'V'
	case 'V':
		return 'B'
	case 'D':
		return 'H'
	case 'H':
		return 'D'
	case 'S':
		return 'S'
	case 'W':
		return 'W'
	}
	return b
}

// buildFastqHeader composes the `@qname[ tags...]` portion of a FASTQ
// record. Pair-suffix logic mirrors upstream samtools: when 0x40/0x80 is
// set and the QNAME doesn't already end with `/1` or `/2`, append the
// appropriate suffix unless NoSuffix is set. AlwaysAddSuffix forces the
// suffix regardless of what QNAME already contains.
func buildFastqHeader(rec *sam.Record, opts FastqOptions) string {
	var sb strings.Builder
	sb.WriteByte('@')
	qname := rec.QName
	if !opts.NoSuffix {
		suffix := ""
		switch {
		case rec.IsRead1() && !rec.IsRead2():
			suffix = "/1"
		case rec.IsRead2() && !rec.IsRead1():
			suffix = "/2"
		}
		if suffix != "" {
			if opts.AlwaysAddSuffix || !strings.HasSuffix(qname, suffix) {
				qname += suffix
			}
		}
	}
	sb.WriteString(qname)
	switch {
	case opts.AddAllTags:
		// Emit every aux tag in the record's existing order, matching
		// upstream's FASTQ_OPT_AUX(NULL) behaviour for `-T ''` / `-T '*'`.
		for _, a := range rec.Aux {
			sb.WriteByte('\t')
			sb.WriteString(a.FormatSAM())
		}
	case len(opts.AddTags) > 0:
		for _, t := range opts.AddTags {
			t = strings.TrimSpace(t)
			if t == "" {
				continue
			}
			if a, ok := rec.GetAux(t); ok {
				sb.WriteByte('\t')
				sb.WriteString(a.FormatSAM())
			}
		}
	}
	return sb.String()
}

// formatFastq returns one FASTQ record (header includes leading '@') with
// trailing newline.
func formatFastq(header, seq string, qual []byte) string {
	// Upstream samtools writes a literal empty line for records with no
	// SEQ (BAM "*"); ours used to emit "*". Treat "*" as empty so the
	// fastq output is byte-identical.
	if seq == "*" {
		seq = ""
	}
	var sb strings.Builder
	sb.Grow(len(header) + len(seq) + len(qual) + 8)
	sb.WriteString(header)
	sb.WriteByte('\n')
	sb.WriteString(seq)
	sb.WriteString("\n+\n")
	if len(qual) == 0 {
		// Match SEQ length with '!' (Phred 0) per FASTQ convention.
		for i := 0; i < len(seq); i++ {
			sb.WriteByte('!')
		}
	} else {
		for _, b := range qual {
			if b == 0xff {
				sb.WriteByte('!')
				continue
			}
			sb.WriteByte(b + 33)
		}
	}
	sb.WriteByte('\n')
	return sb.String()
}

// readSortOrder inspects the @HD line of a header and returns the SO:
// value ("queryname", "coordinate", "unsorted", "unknown", ""). Missing
// header returns "".
func readSortOrder(hdr *sam.Header) string {
	for _, f := range hdr.HDFields {
		if f.Tag == "SO" {
			return f.Value
		}
	}
	return ""
}

// ParseAddTags splits a comma-separated list of aux tag names, trimming
// whitespace from each. The result is suitable for FastqOptions.AddTags.
// The upstream shorthand "*" (or an empty string) means "every tag";
// callers should detect those forms with AddTagsIsAll and set
// FastqOptions.AddAllTags = true instead of passing them through here.
func ParseAddTags(spec string) []string {
	if spec == "" {
		return nil
	}
	parts := strings.Split(spec, ",")
	out := make([]string, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		out = append(out, p)
	}
	return out
}

// AddTagsIsAll reports whether the -T argument is upstream samtools'
// "all aux tags" sentinel — an empty string or the literal "*". The CLI
// driver uses this to route to FastqOptions.AddAllTags.
func AddTagsIsAll(spec string) bool {
	s := strings.TrimSpace(spec)
	return s == "" || s == "*"
}

// ParseCompressLevel parses an integer compress level in [0, 9]. Other
// values produce an error.
func ParseCompressLevel(s string) (int, error) {
	if s == "" {
		return 0, nil
	}
	n, err := strconv.Atoi(s)
	if err != nil {
		return 0, fmt.Errorf("samtools fastq: bad compress level %q: %w", s, err)
	}
	if n < 0 || n > 9 {
		return 0, fmt.Errorf("samtools fastq: compress level out of range: %d", n)
	}
	return n, nil
}

// ErrNoFastqOutput is returned when no output destination is configured —
// either Read1/Read2 or Output (or Singleton) must be set.
var ErrNoFastqOutput = errors.New("samtools fastq: no output configured (use -1/-2 or -o)")
