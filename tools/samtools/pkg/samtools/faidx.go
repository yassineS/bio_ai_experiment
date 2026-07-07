package samtools

// faidx implements `samtools faidx` (and its FASTQ analogue `fqidx`):
// build a FASTA/FASTQ index (.fai, plus .gzi for bgzipped inputs) when no
// regions are supplied, or extract one or more regions to FASTA/FASTQ when
// they are. The behaviour mirrors upstream samtools/htslib byte-for-byte:
//
//   - The index build delegates to the shared pkg/htsgo/fasta index builder
//     (5-column NAME/LENGTH/OFFSET/LINEBASES/LINEWIDTH) and, for BGZF input,
//     scans the compressed blocks to emit the .gzi sidecar.
//   - Region extraction echoes the requested region string verbatim in the
//     '>'/'@' header, preserves the reference letter case, wraps the sequence
//     at the requested line length (default 60; 0 = single line), supports
//     reverse-complement (-i) with the configurable strand marker
//     (--mark-strand), and clamps past-the-end ranges exactly as htslib's
//     fai_get_val does. Region parsing reproduces htslib's hts_parse_region
//     coordinate conventions (1-based inclusive, chr:0 / chr: == whole
//     contig, chr:-N == chr:1-N, chr:N == chr:N-<end>).

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// FaidxFormat selects between the FASTA (faidx) and FASTQ (fqidx) flavours.
type FaidxFormat int

const (
	// FaidxFASTA is the `samtools faidx` flavour: '>' headers, no quality.
	FaidxFASTA FaidxFormat = iota
	// FaidxFASTQ is the `samtools fqidx` flavour: '@' headers plus a '+'
	// separator and a quality line per region.
	FaidxFASTQ
)

// faidxDefaultLineLen is the default sequence line wrap width sentinel (-60).
// Upstream encodes the default as a negative value meaning "same as the input
// data"; we resolve a negative LineLen to the contig's own line length at
// extraction time (fai_line_length), which equals the FASTA's wrap width.
const faidxDefaultLineLen = -60

// FaidxOptions configures the faidx/fqidx subcommand. Start from
// DefaultFaidxOptions so LineLen and the strand markers carry the upstream
// defaults.
type FaidxOptions struct {
	// Format selects FASTA (faidx) or FASTQ (fqidx) behaviour.
	Format FaidxFormat
	// Output is the destination FASTA/FASTQ path (-o/--output). Empty means
	// stdout (the CLI wrapper supplies the io.Writer).
	Output string
	// LineLen is the wrap width (-n/--length). Negative means "same as input"
	// (resolved to the contig's line length); 0 means a single unwrapped line.
	LineLen int
	// Continue, when true (-c/--continue), keeps processing after a region
	// that could not be retrieved instead of aborting with a non-zero exit.
	Continue bool
	// RegionFile names a file of regions, one per line (-r/--region-file).
	RegionFile string
	// ReverseComplement, when true (-i/--reverse-complement), reverse
	// complements each emitted sequence and appends the negative-strand mark.
	ReverseComplement bool
	// PosStrandName / NegStrandName are the strand markers appended to the
	// region name on the +/- strand (configured via --mark-strand). The
	// defaults are "" and "/rc".
	PosStrandName string
	NegStrandName string
	// FaiName overrides the .fai index path (--fai-idx).
	FaiName string
	// GziName overrides the .gzi index path (--gzi-idx).
	GziName string
	// Threads is upstream's -@/--threads worker count. When > 1 it drives
	// block-parallel BGZF inflate of the full-stream index-build pass over a
	// BGZF-compressed FASTA/FASTQ input (upstream attaches fai_thread_pool to
	// the input on the index-build side). Region EXTRACTION uses random-access
	// per-region seeks, not a stream, so it cannot benefit and stays serial;
	// there is no compressed faidx output to parallelise. The .fai/.gzi bytes
	// are identical for any worker count. Parallel decode is opt-in (0/1 stays
	// single-threaded) because each worker adds block buffers to peak RSS —
	// important for genome-scale references.
	Threads int
}

// DefaultFaidxOptions returns FaidxOptions populated with upstream defaults:
// the -60 line-length sentinel and the default "rc" strand markers.
func DefaultFaidxOptions(format FaidxFormat) FaidxOptions {
	return FaidxOptions{
		Format:        format,
		LineLen:       faidxDefaultLineLen,
		PosStrandName: "",
		NegStrandName: "/rc",
	}
}

// ParseMarkStrand applies a --mark-strand TYPE argument to opts, setting the
// positive/negative strand markers. It mirrors upstream's accepted TYPE
// values: "rc" (default), "no", "sign", and "custom,<pos>,<neg>". An
// unrecognised TYPE returns an error whose message matches upstream.
func ParseMarkStrand(opts *FaidxOptions, typ string) error {
	switch {
	case typ == "no":
		opts.PosStrandName = ""
		opts.NegStrandName = ""
	case typ == "sign":
		opts.PosStrandName = "(+)"
		opts.NegStrandName = "(-)"
	case typ == "rc":
		opts.PosStrandName = ""
		opts.NegStrandName = "/rc"
	case strings.HasPrefix(typ, "custom,"):
		rest := typ[len("custom,"):]
		// Split on the first comma: the part before is the +ve marker, the
		// part after (if any) is the -ve marker — matching faidx.c's malloc
		// split on the first comma after "custom,".
		comma := strings.IndexByte(rest, ',')
		if comma < 0 {
			opts.PosStrandName = rest
			opts.NegStrandName = ""
		} else {
			opts.PosStrandName = rest[:comma]
			opts.NegStrandName = rest[comma+1:]
		}
	default:
		return fmt.Errorf("[faidx] Unknown --mark-strand option \"%s\"", typ)
	}
	return nil
}

// Faidx is the high-level entry point: it indexes or extracts according to
// opts and the supplied regions, writing FASTA/FASTQ to out and diagnostics
// to errOut. It returns the process exit status (0 on success). path is the
// input FASTA/FASTQ. When regions is empty and no region file is set, Faidx
// builds the index instead of extracting.
func Faidx(path string, regions []string, opts FaidxOptions, out io.Writer, errOut io.Writer) int {
	// Index-build mode: exactly one positional (the FASTA) and no regions.
	if len(regions) == 0 && opts.RegionFile == "" {
		if err := FaidxBuild(path, opts); err != nil {
			faiName := opts.FaiName
			if faiName == "" {
				faiName = path + ".fai"
			}
			fmt.Fprintf(errOut, "[faidx] Could not build fai index %s\n", faiName)
			return 1
		}
		return 0
	}

	// Extraction mode.
	ra, err := openFaidxRandomAccess(path, opts)
	if err != nil {
		faiName := opts.FaiName
		if faiName == "" {
			faiName = path + ".fai"
		}
		fmt.Fprintf(errOut, "[faidx] Could not load fai index %s\n", faiName)
		return 1
	}
	defer ra.Close()

	// For FASTQ extraction we also need a quality fetcher backed by the
	// 6-column index; nil for FASTA.
	var qualFetch func(name string, beg, end int64) ([]byte, error)
	if opts.Format == FaidxFASTQ {
		qa, qerr := openFaidxQual(path, opts)
		if qerr == nil && qa != nil {
			defer qa.Close()
			qualFetch = qa.FetchQual
		}
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	exit := 0
	emit := func(region string) bool {
		if !writeFaidxRegion(ra, bw, errOut, region, opts, qualFetch) {
			if !opts.Continue {
				exit = 1
				return false
			}
		}
		return true
	}

	if opts.RegionFile != "" {
		if !faidxRegionsFromFile(opts.RegionFile, emit, errOut) {
			return 1
		}
		if exit != 0 {
			return exit
		}
	}
	for _, r := range regions {
		if !emit(r) {
			break
		}
	}
	return exit
}

// faidxRegionsFromFile reads a region file (one region per line) and calls
// emit for each. It returns false only when the file cannot be opened
// (upstream prints an error and aborts in that case).
func faidxRegionsFromFile(path string, emit func(string) bool, errOut io.Writer) bool {
	f, err := os.Open(path)
	if err != nil {
		fmt.Fprintf(errOut, "[faidx] Failed to open \"%s\" for reading.\n", path)
		return false
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r\n")
		if !emit(line) {
			return true // emit signalled stop; caller inspects exit code
		}
	}
	return true
}

// openFaidxRandomAccess opens path for random access. With no index override
// it defers to OpenRandomAccess (which finds the sibling .fai/.gzi). When a
// plain FASTA carries an explicit --fai-idx, the named index is loaded and
// paired with the input directly.
func openFaidxRandomAccess(path string, opts FaidxOptions) (*fasta.RandomAccess, error) {
	if opts.FaiName == "" {
		return fasta.OpenRandomAccess(path)
	}
	bgzipped, err := isBGZFPath(path)
	if err != nil {
		return nil, err
	}
	if bgzipped {
		// The partial-decompression backend expects the sidecars in their
		// default locations; defer to OpenRandomAccess for bgzipped input.
		return fasta.OpenRandomAccess(path)
	}
	idx, err := fasta.LoadIndex(opts.FaiName)
	if err != nil {
		return nil, err
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	return fasta.NewRandomAccessFile(f, idx), nil
}

// isBGZFPath reports whether path begins with the BGZF magic bytes.
func isBGZFPath(path string) (bool, error) {
	f, err := os.Open(path)
	if err != nil {
		return false, err
	}
	defer f.Close()
	var buf [4]byte
	n, _ := io.ReadFull(f, buf[:])
	if n < 4 {
		return false, nil
	}
	return buf[0] == 0x1f && buf[1] == 0x8b && buf[2] == 0x08 && buf[3] == 0x04, nil
}

// writeFaidxRegion extracts one region and writes it to w in FASTA/FASTQ
// form, echoing the requested region string in the header (with the strand
// marker) and preserving the reference case. It returns false when the region
// could not be retrieved (unknown contig / malformed region), printing the
// same diagnostics upstream does. Past-the-end ranges are clamped and emit a
// "Truncated sequence" warning but still count as success (return true).
//
// qualFetch supplies FASTQ quality bytes (nil for FASTA / when unavailable).
func writeFaidxRegion(ra *fasta.RandomAccess, w *bufio.Writer, errOut io.Writer, region string, opts FaidxOptions, qualFetch func(name string, beg, end int64) ([]byte, error)) bool {
	wrap := opts.LineLen

	pr, ok := parseFaidxRegion(ra.Index(), region)
	if !ok {
		// Region could not be parsed/resolved: upstream prints the warnings
		// from fai_get_val (once for the fai_line_length probe under the
		// default line length, plus once for the fetch) then "Failed to
		// fetch", and still writes the (empty) header to stdout.
		if wrap < 0 {
			fmt.Fprintf(errOut, "[W::fai_get_val] Reference %s not found in FASTA file, returning empty sequence\n", region)
		}
		fmt.Fprintf(errOut, "[W::fai_get_val] Reference %s not found in FASTA file, returning empty sequence\n", region)
		fmt.Fprintf(errOut, "[faidx] Failed to fetch sequence in %s\n", region)
		writeFaidxHeader(w, region, opts, false)
		return false
	}

	if wrap < 0 {
		wrap = int(pr.lineBases)
	}
	if wrap <= 0 {
		wrap = 1 << 62 // effectively a single line
	}

	seq, err := ra.FetchRaw(pr.name, pr.beg0, pr.end0)
	if err != nil {
		fmt.Fprintf(errOut, "[faidx] Failed to fetch sequence in %s\n", region)
		writeFaidxHeader(w, region, opts, false)
		return false
	}

	rev := opts.ReverseComplement && len(seq) > 0
	if rev {
		reverseComplementInPlace(seq)
	}

	writeFaidxHeader(w, region, opts, opts.ReverseComplement)

	// Zero-length and truncation diagnostics mirror upstream write_line: a
	// truncation warning fires when an explicit end was given (not the
	// open-ended sentinel) and the clamped length is short of the request.
	seqLen := int64(len(seq))
	if seqLen == 0 {
		fmt.Fprintf(errOut, "[faidx] Zero length sequence: %s\n", region)
	} else if pr.hasExplicitEnd && seqLen != pr.reqEnd0-pr.reqBeg0 {
		fmt.Fprintf(errOut, "[faidx] Truncated sequence: %s\n", region)
	}

	writeWrapped(w, seq, wrap)

	if opts.Format == FaidxFASTQ && qualFetch != nil {
		w.WriteString("+\n")
		qual, qerr := qualFetch(pr.name, pr.beg0, pr.end0)
		if qerr == nil {
			if rev {
				reverseInPlace(qual)
			}
			writeWrapped(w, qual, wrap)
		}
	}
	return true
}

// writeFaidxHeader writes the '>'/'@' header line for region, appending the
// +ve or -ve strand marker according to rev.
func writeFaidxHeader(w *bufio.Writer, region string, opts FaidxOptions, rev bool) {
	if opts.Format == FaidxFASTQ {
		w.WriteByte('@')
	} else {
		w.WriteByte('>')
	}
	w.WriteString(region)
	if rev {
		w.WriteString(opts.NegStrandName)
	} else {
		w.WriteString(opts.PosStrandName)
	}
	w.WriteByte('\n')
}

// writeWrapped writes seq to w wrapped at width bytes per line, each line
// terminated by '\n'. A zero-length seq writes nothing.
func writeWrapped(w *bufio.Writer, seq []byte, width int) {
	for i := 0; i < len(seq); i += width {
		end := i + width
		if end > len(seq) {
			end = len(seq)
		}
		w.Write(seq[i:end])
		w.WriteByte('\n')
	}
}

// faidxComplement is htslib's comp_base table restricted to the bases faidx
// reverse-complements; unmapped bytes map to themselves.
var faidxComplement = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = byte(i)
	}
	pairs := map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A',
		'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
		'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N',
		'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'u': 'a',
		'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm', 'm': 'k',
		'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n',
	}
	for k, v := range pairs {
		t[k] = v
	}
	return t
}()

// reverseComplementInPlace reverse-complements str in place, matching
// htslib's reverse_complement (comp_base lookup, midpoint swap including the
// centre element so it is complemented too).
func reverseComplementInPlace(str []byte) {
	i, j := 0, len(str)-1
	for i <= j {
		ci, cj := str[i], str[j]
		str[i] = faidxComplement[cj]
		str[j] = faidxComplement[ci]
		i++
		j--
	}
}

// reverseInPlace reverses str in place (used for FASTQ quality under -i).
func reverseInPlace(str []byte) {
	i, j := 0, len(str)-1
	for i < j {
		str[i], str[j] = str[j], str[i]
		i++
		j--
	}
}
