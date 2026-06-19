// Package bcftools — `bcftools isec` subcommand.
//
// `bcftools isec` performs set operations on N sorted VCF/BCF inputs. The
// canonical use is "find variants present in at least N of the inputs", but
// the flags also support exact membership counts and explicit inclusion
// bitmasks. Outputs can be either:
//
//   - the per-input projections written to `<prefix>/000<i>.vcf` (one file
//     per input, holding the records selected from that input — the
//     `-p`/`--prefix` mode), or
//   - one or more inputs written to stdout (the `-w`/`--write` mode).
//
// The two modes can be combined: `-p` always writes the per-input projections,
// and `-w` additionally streams the chosen inputs to stdout.
//
// Records collapse on (CHROM, POS, REF, ALT) by default. The `-c/--collapse`
// flag widens the comparison: `none` keeps the strict tuple, `snps` collapses
// any pair of SNP-records at the same site regardless of allele, `indels`
// likewise for indels, `both` is the union of `snps`+`indels`, `all`
// collapses everything at the same (CHROM, POS), `id` matches on the ID
// column, and `some` is upstream's hybrid (REF must match, any ALT in common
// is enough). The default is `none`.
package bcftools

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// CollapseMode controls how `isec` decides two records "match".
type CollapseMode int

const (
	// CollapseNone requires identical (CHROM, POS, REF, ALT). Default.
	CollapseNone CollapseMode = iota
	// CollapseSNPs collapses any two SNP records at the same site.
	CollapseSNPs
	// CollapseIndels collapses any two indel records at the same site.
	CollapseIndels
	// CollapseBoth is the union of CollapseSNPs and CollapseIndels.
	CollapseBoth
	// CollapseAll collapses every record at the same (CHROM, POS).
	CollapseAll
	// CollapseSome requires REF match and at least one ALT in common.
	CollapseSome
	// CollapseID matches on the ID column.
	CollapseID
)

// ParseCollapseMode parses the -c/--collapse flag.
func ParseCollapseMode(s string) (CollapseMode, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "", "none":
		return CollapseNone, nil
	case "snps":
		return CollapseSNPs, nil
	case "indels":
		return CollapseIndels, nil
	case "both":
		return CollapseBoth, nil
	case "all":
		return CollapseAll, nil
	case "some":
		return CollapseSome, nil
	case "id":
		return CollapseID, nil
	}
	return CollapseNone, fmt.Errorf("bcftools isec: unknown -c mode %q (expect none|snps|indels|both|all|some|id)", s)
}

// NfilesSpec encodes the `-n/--nfiles` argument.
type NfilesSpec struct {
	// Mode: '~' bitmask, '+' at-least-N, '=' exactly-N, 0 unset.
	Mode byte
	// N is the threshold for `+N` and `=N`. For `~bits` it is unused.
	N int
	// Bits is the inclusion mask for `~bits` mode (one bit per input,
	// bit 0 = first input).
	Bits []bool
}

// ParseNfilesSpec parses the -n/--nfiles argument. Empty s returns the
// "no constraint" spec (Mode == 0).
func ParseNfilesSpec(s string) (NfilesSpec, error) {
	out := NfilesSpec{}
	s = strings.TrimSpace(s)
	if s == "" {
		return out, nil
	}
	if len(s) >= 2 && (s[0] == '+' || s[0] == '=') {
		var n int
		_, err := fmt.Sscanf(s[1:], "%d", &n)
		if err != nil {
			return out, fmt.Errorf("bcftools isec: bad -n %q", s)
		}
		out.Mode = s[0]
		out.N = n
		return out, nil
	}
	if s[0] == '~' {
		bits := make([]bool, 0, len(s)-1)
		for _, r := range s[1:] {
			switch r {
			case '0':
				bits = append(bits, false)
			case '1':
				bits = append(bits, true)
			default:
				return out, fmt.Errorf("bcftools isec: bad -n bitmask %q (expect 0/1 chars)", s)
			}
		}
		out.Mode = '~'
		out.Bits = bits
		return out, nil
	}
	// Bare integer == "+N".
	var n int
	_, err := fmt.Sscanf(s, "%d", &n)
	if err != nil {
		return out, fmt.Errorf("bcftools isec: bad -n %q", s)
	}
	out.Mode = '+'
	out.N = n
	return out, nil
}

// IsecOptions controls the `bcftools isec` subcommand.
type IsecOptions struct {
	// Nfiles is the membership constraint. When unset (Mode == 0) all
	// records that appear in any input are emitted.
	Nfiles NfilesSpec
	// Collapse controls how records collapse. Defaults to CollapseNone.
	Collapse CollapseMode
	// Prefix is the per-input output directory: when non-empty, one file
	// `<prefix>/000<i>.vcf[.gz|.bcf]` is written per input.
	Prefix string
	// Write is the 1-based input indices to dump to stdout (the `-w` flag).
	// When empty and Prefix is also empty, every selected record is dumped
	// from input 1.
	Write []int
	// OutputFormat selects the encoding for both the prefix files and the
	// stdout stream.
	OutputFormat OutputFormat
	// CompressLevel is the gzip level for -O z output.
	CompressLevel int
	// Threads is upstream's -@/--threads. When greater than 1 it enables
	// parallel BGZF compression of -O z and -O b output via bgzf.MultiWriter;
	// the framed result decodes byte-identically regardless of thread count.
	Threads int
	// InputPaths are the input file names, used only for the README.txt legend
	// (which names each private/shared output file after its source). Set by
	// IsecFiles; may be empty when Isec is called directly.
	InputPaths []string
	// CmdArgv is the `isec ...` argument vector echoed into README.txt's "The
	// command line was:" line. Set by the CLI; may be empty.
	CmdArgv []string
}

// IsecFiles is the file-aware entry point. It opens each path through
// iohelper and runs the set operation.
//
// For the common case — every input is coordinate-sorted and chrom-contiguous
// with a consistent contig order — it takes a streaming fast-path that
// processes one contig at a time, so peak memory is bounded by a single
// contig's variants across all inputs rather than the full corpus. When any
// input violates that ordering (a contig reappears, POS decreases within a
// contig, or contigs are out of keyFor order), it falls back to reading every
// variant into memory and calling Isec, which is unchanged and remains correct
// for arbitrary/unsorted input.
func IsecFiles(paths []string, out io.Writer, opts IsecOptions) (int, error) {
	if len(paths) < 2 {
		return 0, fmt.Errorf("bcftools isec: need at least two input files (got %d)", len(paths))
	}
	if opts.InputPaths == nil {
		opts.InputPaths = paths
	}

	// Cheap pre-scan: read each file's header and (chrom, pos) stream, checking
	// for streamability and recording each file's contig sequence. No records
	// are retained, so this pass is memory-cheap.
	headers := make([]*vcf.Header, len(paths))
	fileContigs := make([][]string, len(paths))
	streamable := true
	for i, p := range paths {
		hdr, contigs, ok, err := isecPrescan(p)
		if err != nil {
			return 0, fmt.Errorf("bcftools isec: %s: %w", p, err)
		}
		headers[i] = hdr
		fileContigs[i] = contigs
		if !ok {
			streamable = false
		}
	}

	if streamable {
		// Build the merged contig order in keyFor order using header[0]'s
		// contig ranking (the same ordering Isec's keyOrder sort uses). Each
		// file's contigs already appear in keyFor-ascending order (verified in
		// the pre-scan); confirm the union is consistent across files, i.e. it
		// can be visited in a single ascending keyFor order. If not, fall back.
		primaryOrder := contigOrder(headers[0])
		mergedContigs, ok := isecMergeContigOrder(fileContigs, primaryOrder)
		if ok {
			return isecStreaming(paths, headers, mergedContigs, out, opts)
		}
		streamable = false
	}

	// Fallback: read all variants into memory and run Isec unchanged.
	groups := make([][]*vcf.Variant, len(paths))
	for i, p := range paths {
		in, err := iohelper.OpenReader(p)
		if err != nil {
			return 0, fmt.Errorf("bcftools isec: open %s: %w", p, err)
		}
		hdr, recs, err := readAllVariants(in)
		_ = in.Close()
		if err != nil {
			return 0, fmt.Errorf("bcftools isec: %s: %w", p, err)
		}
		headers[i] = hdr
		groups[i] = recs
	}
	return Isec(headers, groups, out, opts)
}

// isecStreamFile opens path through iohelper, reads its header, and invokes fn
// for each record in file order. It handles both VCF and BCF transparently
// (mirroring readAllVariants' format detection) and never retains records, so
// callers control retention. fn may return an error to stop early.
func isecStreamFile(path string, fn func(*vcf.Variant) error) (*vcf.Header, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("open %s: %w", path, err)
	}
	defer in.Close()
	br := bufio.NewReader(in)
	head, perr := br.Peek(5)
	if perr != nil && perr != io.EOF {
		return nil, perr
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		r, err := bcf.NewReader(br)
		if err != nil {
			return nil, err
		}
		hdr := r.Header()
		for {
			rec, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return hdr.VCF, err
			}
			if err := fn(rec.ToVariant(hdr)); err != nil {
				return hdr.VCF, err
			}
		}
		return hdr.VCF, nil
	}
	r := vcf.NewReader(br)
	hdr, err := r.ReadHeader()
	if err != nil {
		return nil, err
	}
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return hdr, err
		}
		if err := fn(v); err != nil {
			return hdr, err
		}
	}
	return hdr, nil
}

// isecPrescan streams path parsing only chrom/pos, returning its header, its
// contig sequence (in first-appearance order), and whether the file is
// internally streamable. A file is streamable when its records are
// non-decreasing in raw (contig-contiguous, pos non-decreasing within a
// contig) terms: a contig never reappears after a different contig, and POS
// never decreases within a contig run. Streamability here is independent of
// keyFor order; the cross-file keyFor consistency is checked later in
// isecMergeContigOrder.
func isecPrescan(path string) (*vcf.Header, []string, bool, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, nil, false, fmt.Errorf("open %s: %w", path, err)
	}
	defer in.Close()
	br := bufio.NewReader(in)
	head, perr := br.Peek(5)
	if perr != nil && perr != io.EOF {
		return nil, nil, false, perr
	}
	// BCF (and any malformed text) uses the full-parse scan; extracting chrom/pos
	// from BCF needs a binary decode, so it is not worth a bespoke fast path.
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		return isecPrescanFull(path)
	}

	sc := bufio.NewScanner(br)
	sc.Buffer(make([]byte, 0, 64*1024), 64*1024*1024)
	var hdrText []byte
	sawHeader := false
	for sc.Scan() {
		line := sc.Bytes()
		if len(line) == 0 {
			continue
		}
		if line[0] != '#' {
			return isecPrescanFull(path) // data before #CHROM — be safe
		}
		hdrText = append(hdrText, line...)
		hdrText = append(hdrText, '\n')
		if bytes.HasPrefix(line, []byte("#CHROM")) {
			sawHeader = true
			break
		}
	}
	if e := sc.Err(); e != nil {
		return nil, nil, false, e
	}
	if !sawHeader {
		return isecPrescanFull(path)
	}
	// Build the header by re-parsing the collected header text (small, one-time)
	// — identical to the header the full scan would return.
	hdr, herr := vcf.NewReader(bytes.NewReader(hdrText)).ReadHeader()
	if herr != nil {
		return nil, nil, false, herr
	}

	var contigs []string
	streamable := true
	seen := map[string]bool{}
	curChrom := ""
	curPos := 0
	have := false
	for sc.Scan() {
		line := sc.Bytes()
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		t1 := bytes.IndexByte(line, '\t')
		if t1 < 0 {
			return isecPrescanFull(path) // malformed — be safe
		}
		chromB := line[:t1]
		rest := line[t1+1:]
		posB := rest
		if t2 := bytes.IndexByte(rest, '\t'); t2 >= 0 {
			posB = rest[:t2]
		}
		pos, ok := isecParsePos(posB)
		if !ok {
			return isecPrescanFull(path) // non-integer POS — let the full parser report it
		}
		// string(chromB) in the comparison does not allocate; a fresh chrom
		// string is allocated only when a new contig is encountered.
		if !have || string(chromB) != curChrom {
			cs := string(chromB)
			if seen[cs] {
				streamable = false
			} else {
				seen[cs] = true
				contigs = append(contigs, cs)
			}
			curChrom = cs
			curPos = pos
			have = true
			continue
		}
		if pos < curPos {
			streamable = false
		}
		curPos = pos
	}
	if e := sc.Err(); e != nil {
		return nil, nil, false, e
	}
	return hdr, contigs, streamable, nil
}

// isecParsePos parses a VCF POS (a non-negative decimal integer) from bytes
// without allocating. It returns ok=false for any non-digit content so the
// caller can defer to the full parser (which reports the error).
func isecParsePos(b []byte) (int, bool) {
	if len(b) == 0 {
		return 0, false
	}
	n := 0
	for _, c := range b {
		if c < '0' || c > '9' {
			return 0, false
		}
		n = n*10 + int(c-'0')
	}
	return n, true
}

// isecPrescanFull is the full-parse pre-scan (used for BCF and any malformed
// text the fast path bails on). It mirrors isecPrescan's streamability logic.
func isecPrescanFull(path string) (*vcf.Header, []string, bool, error) {
	var contigs []string
	streamable := true
	seen := map[string]bool{}
	curChrom := ""
	curPos := 0
	have := false
	hdr, err := isecStreamFile(path, func(v *vcf.Variant) error {
		if !have || v.Chrom != curChrom {
			if seen[v.Chrom] {
				// Contig reappears after a different contig.
				streamable = false
			} else {
				seen[v.Chrom] = true
				contigs = append(contigs, v.Chrom)
			}
			curChrom = v.Chrom
			curPos = v.Pos
			have = true
			return nil
		}
		if v.Pos < curPos {
			streamable = false
		}
		curPos = v.Pos
		return nil
	})
	if err != nil {
		return nil, nil, false, err
	}
	return hdr, contigs, streamable, nil
}

// isecMergeContigOrder merges the per-file contig sequences into a single
// global order sorted by keyFor's contig ranking (using primaryOrder, the
// header[0] contig order). It returns the merged order and whether the union is
// consistent: each file's contigs must already be a subsequence of the
// keyFor-ascending global order (i.e. every file lists its contigs in
// keyFor-ascending order). If any file's contig sequence is out of keyFor order
// the result is not streamable and ok is false.
func isecMergeContigOrder(fileContigs [][]string, primaryOrder map[string]int) ([]string, bool) {
	rank := func(chrom string) int {
		if idx, ok := primaryOrder[chrom]; ok {
			return idx
		}
		return 1<<30 + sortFallback(chrom)
	}
	// Verify each file lists contigs in keyFor-ascending rank order, and that
	// no two distinct contigs share a rank in a way that would make the order
	// ambiguous across files. Collect the union.
	union := map[string]bool{}
	for _, contigs := range fileContigs {
		prev := -1
		first := true
		for _, c := range contigs {
			r := rank(c)
			if !first && r < prev {
				// File lists contigs out of keyFor order.
				return nil, false
			}
			if !first && r == prev {
				// Two distinct contigs collide on the same rank within a single
				// file; ambiguous ordering — fall back conservatively.
				return nil, false
			}
			prev = r
			first = false
			union[c] = true
		}
	}
	merged := make([]string, 0, len(union))
	for c := range union {
		merged = append(merged, c)
	}
	sort.SliceStable(merged, func(i, j int) bool {
		return rank(merged[i]) < rank(merged[j])
	})
	// Distinct contigs must not collide on the same rank globally either.
	for i := 1; i < len(merged); i++ {
		if rank(merged[i]) == rank(merged[i-1]) {
			return nil, false
		}
	}
	return merged, true
}

// isecStreaming runs the set operation one contig at a time. For each contig in
// keyFor order it re-streams every input, collecting only that contig's
// variants from each file into per-contig groups, runs them through isecCore
// against the shared open writers, then frees them. Peak memory is bounded by a
// single contig's variants across all inputs.
func isecStreaming(paths []string, headers []*vcf.Header, mergedContigs []string, out io.Writer, opts IsecOptions) (int, error) {
	state, err := isecSetup(headers, out, opts)
	if err != nil {
		return 0, err
	}
	n := len(paths)

	// Open one forward cursor per input and read each file exactly once: for
	// each contig (visited in keyFor order, which every input lists its contigs
	// in — verified by isecMergeContigOrder) pull that contig's contiguous run
	// from each cursor. A file lacking the contig contributes nothing; its
	// cursor simply stays parked on its next (later) contig.
	cursors := make([]*isecCursor, 0, n)
	closeCursors := func() {
		for _, c := range cursors {
			c.close()
		}
	}
	for _, p := range paths {
		c, err := openIsecCursor(p)
		if err != nil {
			closeCursors()
			_, _ = isecFinalize(state)
			return state.totalKept, fmt.Errorf("bcftools isec: %s: %w", p, err)
		}
		cursors = append(cursors, c)
	}
	defer closeCursors()

	for _, contig := range mergedContigs {
		groups := make([][]*vcf.Variant, n)
		for i, c := range cursors {
			for v := c.peek(); v != nil && v.Chrom == contig; v = c.peek() {
				groups[i] = append(groups[i], c.pop())
			}
			if c.err != nil && c.err != io.EOF {
				_, _ = isecFinalize(state)
				return state.totalKept, fmt.Errorf("bcftools isec: %s: %w", paths[i], c.err)
			}
		}
		if err := isecCore(state, groups); err != nil {
			_, _ = isecFinalize(state)
			return state.totalKept, err
		}
	}
	return isecFinalize(state)
}

// isecCursor is a forward, peekable pull reader over one input file (VCF or
// BCF, transparently). It yields owned *vcf.Variant values one at a time so the
// streaming merge can read each input exactly once.
type isecCursor struct {
	closer io.Closer
	next   func() (*vcf.Variant, error) // returns io.EOF when exhausted
	head   *vcf.Variant
	err    error
}

// openIsecCursor opens path and primes the first record. The reader uses Read
// (not ReadInto), so each yielded Variant is independently owned and safe to
// retain in a per-contig group.
func openIsecCursor(path string) (*isecCursor, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("open %s: %w", path, err)
	}
	br := bufio.NewReader(in)
	magic, perr := br.Peek(5)
	if perr != nil && perr != io.EOF {
		_ = in.Close()
		return nil, perr
	}
	c := &isecCursor{closer: in}
	if len(magic) >= 5 && magic[0] == 'B' && magic[1] == 'C' && magic[2] == 'F' {
		r, err := bcf.NewReader(br)
		if err != nil {
			_ = in.Close()
			return nil, err
		}
		bh := r.Header()
		c.next = func() (*vcf.Variant, error) {
			rec, err := r.Read()
			if err != nil {
				return nil, err
			}
			return rec.ToVariant(bh), nil
		}
	} else {
		r := vcf.NewReader(br)
		if _, err := r.ReadHeader(); err != nil {
			_ = in.Close()
			return nil, err
		}
		c.next = func() (*vcf.Variant, error) { return r.Read() }
	}
	c.advance()
	return c, nil
}

func (c *isecCursor) advance() {
	if c.err != nil {
		c.head = nil
		return
	}
	v, err := c.next()
	if err != nil {
		c.err = err
		c.head = nil
		return
	}
	c.head = v
}

func (c *isecCursor) peek() *vcf.Variant { return c.head }

func (c *isecCursor) pop() *vcf.Variant {
	v := c.head
	c.advance()
	return v
}

func (c *isecCursor) close() { _ = c.closer.Close() }

// isecState holds the open output writers and shared scratch for a single
// isec run, so the membership/emit core can be applied to either the whole
// corpus at once (Isec) or one contig at a time (isecStreaming). It is created
// by isecSetup and torn down by isecFinalize.
type isecState struct {
	opts    IsecOptions
	headers []*vcf.Header
	n       int

	// venn is true in the two-input no-constraint Venn decomposition mode.
	venn bool

	// Per-input (or four-way Venn) projection writers, opened when -p is set.
	perInputW     []variantWriter
	perInputClose []func()
	sitesF        *os.File

	// stdout stream writer, opened when -w is set or no -p/-w is given.
	openStdout      bool
	stdoutW         variantWriter
	stdoutFinish    func()
	writeFromInputs []int

	totalKept int
	bits      []byte // reused membership scratch, length n
}

// Isec performs the set operation on pre-loaded headers/groups, writing the
// per-input projections to `<prefix>/000<i>.vcf*` (when Prefix is set) and
// the stdout stream defined by `-w` (when Write is set, defaulting to input
// 1 when neither is set).
func Isec(headers []*vcf.Header, groups [][]*vcf.Variant, stdout io.Writer, opts IsecOptions) (int, error) {
	if len(headers) != len(groups) {
		return 0, fmt.Errorf("bcftools isec: header / variant count mismatch")
	}
	state, err := isecSetup(headers, stdout, opts)
	if err != nil {
		return 0, err
	}
	if err := isecCore(state, groups); err != nil {
		_, _ = isecFinalize(state)
		return state.totalKept, err
	}
	return isecFinalize(state)
}

// isecSetup validates the inputs, opens the per-input/venn projection writers,
// sites.txt, and the stdout stream writer, and returns the populated state.
// It contains all the output-opening logic; the membership build and emit loop
// live in isecCore so they can be driven incrementally.
func isecSetup(headers []*vcf.Header, stdout io.Writer, opts IsecOptions) (*isecState, error) {
	n := len(headers)
	if n == 0 {
		return nil, fmt.Errorf("bcftools isec: no inputs")
	}
	if opts.Nfiles.Mode == '~' && len(opts.Nfiles.Bits) != n {
		return nil, fmt.Errorf("bcftools isec: -n ~bits length %d != number of inputs %d", len(opts.Nfiles.Bits), n)
	}

	st := &isecState{
		opts:    opts,
		headers: headers,
		n:       n,
		bits:    make([]byte, n),
	}

	// VENN mode: with exactly two inputs and no membership constraint, upstream
	// (vcfisec.c: nfiles==2 && !isec_op -> OP_VENN) writes the four-way Venn
	// decomposition — 0000 records private to input 1, 0001 private to input 2,
	// 0002 the input-1 copy of shared records, 0003 the input-2 copy — rather
	// than one file per input.
	st.venn = n == 2 && opts.Nfiles.Mode == 0

	// Open per-input output files when -p is set. In VENN mode there are four
	// files; file i draws its header from input vennReader[i].
	vennReader := []int{0, 1, 0, 1}
	if opts.Prefix != "" {
		if err := os.MkdirAll(opts.Prefix, 0755); err != nil {
			return nil, fmt.Errorf("bcftools isec: mkdir %s: %w", opts.Prefix, err)
		}
		ext := outputExt(opts.OutputFormat)
		nOut := n
		hdrFor := func(i int) *vcf.Header { return headers[i] }
		if st.venn {
			nOut = 4
			hdrFor = func(i int) *vcf.Header { return headers[vennReader[i]] }
		}
		// README.txt names each output file after its source; build it with the
		// real output paths so the legend matches upstream.
		outNames := make([]string, nOut)
		for i := 0; i < nOut; i++ {
			outNames[i] = filepath.Join(opts.Prefix, fmt.Sprintf("%04d.vcf%s", i, ext))
		}
		readmePath := filepath.Join(opts.Prefix, "README.txt")
		_ = os.WriteFile(readmePath, []byte(isecReadme(opts, outNames, st.venn)), 0644)
		for i := 0; i < nOut; i++ {
			f, err := os.Create(outNames[i])
			if err != nil {
				return nil, fmt.Errorf("bcftools isec: create %s: %w", outNames[i], err)
			}
			w, finish, err := openOutput(f, ViewOptions{
				OutputFormat:  opts.OutputFormat,
				CompressLevel: opts.CompressLevel,
				Threads:       opts.Threads,
			}, hdrFor(i))
			if err != nil {
				_ = f.Close()
				return nil, err
			}
			closeF := f
			st.perInputW = append(st.perInputW, w)
			st.perInputClose = append(st.perInputClose, func() {
				finish()
				_ = closeF.Close()
			})
			if err := w.WriteHeader(); err != nil {
				return nil, err
			}
		}
		// sites.txt (one line per matched site, written in the per-occurrence
		// emit loop below) lists every site's CHROM POS REF ALT and membership
		// bits, matching upstream's `CHROM\tPOS\tREF\tALT\tBITS` shape.
		sitesPath := filepath.Join(opts.Prefix, "sites.txt")
		sf, err := os.Create(sitesPath)
		if err != nil {
			return nil, fmt.Errorf("bcftools isec: create %s: %w", sitesPath, err)
		}
		st.sitesF = sf
	}

	// Open the stdout writer if -w is set or both -p and -w are empty
	// (default = dump from the first input).
	st.writeFromInputs = opts.Write
	st.openStdout = opts.Prefix == "" || len(st.writeFromInputs) > 0
	if st.openStdout {
		// Header source: union (first input) when multiple write-targets,
		// or the chosen input when exactly one.
		hdrIdx := 0
		if len(st.writeFromInputs) == 1 {
			if st.writeFromInputs[0]-1 >= 0 && st.writeFromInputs[0]-1 < n {
				hdrIdx = st.writeFromInputs[0] - 1
			}
		}
		w, finish, err := openOutput(stdout, ViewOptions{
			OutputFormat:  opts.OutputFormat,
			CompressLevel: opts.CompressLevel,
			Threads:       opts.Threads,
		}, headers[hdrIdx])
		if err != nil {
			return nil, err
		}
		st.stdoutW = w
		st.stdoutFinish = finish
		if err := st.stdoutW.WriteHeader(); err != nil {
			return nil, err
		}
		if len(st.writeFromInputs) == 0 {
			st.writeFromInputs = []int{1}
		}
	}
	return st, nil
}

// isecCore builds the membership maps for the supplied groups, sorts the union
// key list by (contig order, POS), and runs the per-occurrence emit loop,
// writing to st's already-open writers and accumulating st.totalKept. It may be
// called repeatedly (once per contig in the streaming path); each call's
// keyOrder sort is self-contained because the global sort is contig-first, so
// concatenating per-contig outputs reproduces the single-pass output exactly.
func isecCore(st *isecState, groups [][]*vcf.Variant) error {
	n := st.n
	opts := st.opts
	if len(groups) != n {
		return fmt.Errorf("bcftools isec: header / variant count mismatch")
	}

	// Build per-input keyed maps and the union key list with membership.
	keyOf := func(v *vcf.Variant) string {
		return collapseKey(v, opts.Collapse)
	}
	type cell struct {
		variant   *vcf.Variant
		groupIdx  int
		variantIx int // index within its group
	}
	keyMembership := map[string]*[]bool{}
	keyOrder := []string{}
	keyVariants := map[string][]cell{}
	for gi, g := range groups {
		for vi, v := range g {
			k := keyOf(v)
			if _, ok := keyMembership[k]; !ok {
				m := make([]bool, n)
				keyMembership[k] = &m
				keyOrder = append(keyOrder, k)
			}
			(*keyMembership[k])[gi] = true
			keyVariants[k] = append(keyVariants[k], cell{variant: v, groupIdx: gi, variantIx: vi})
		}
	}

	// Sort the union key list by (contig order, POS) only, stably. Upstream's
	// synced reader presents records at the same position in file order, so the
	// tie-break must preserve first-appearance order (keyOrder is built in file
	// order) rather than re-sort by REF/ALT — otherwise two records at one POS
	// (e.g. A>T then A>C) would be emitted A>C-first.
	primaryOrder := contigOrder(st.headers[0])
	sort.SliceStable(keyOrder, func(i, j int) bool {
		ai, bi := keyVariants[keyOrder[i]][0].variant, keyVariants[keyOrder[j]][0].variant
		return keyFor(ai, primaryOrder).less(keyFor(bi, primaryOrder))
	})

	// The -n constraint is applied per matched site in the emit loop below (a
	// site is one paired occurrence across inputs), mirroring upstream's
	// per-record synced-reader test rather than a per-position one.

	// Walk the keys and emit one logical site per matched occurrence. When a
	// key holds several records from one input (intra-position duplicates),
	// upstream's synced reader pairs the k-th occurrence across inputs, so we
	// iterate occurrences rather than collapse them — each occurrence has its
	// own membership and its own sites.txt line.
	bits := st.bits
	for _, k := range keyOrder {
		// Group this key's records by input, preserving file order.
		byFile := make([][]*vcf.Variant, n)
		for _, c := range keyVariants[k] {
			byFile[c.groupIdx] = append(byFile[c.groupIdx], c.variant)
		}
		maxOcc := 0
		for gi := range byFile {
			if len(byFile[gi]) > maxOcc {
				maxOcc = len(byFile[gi])
			}
		}
		for occ := 0; occ < maxOcc; occ++ {
			present := make([]bool, n)
			firstGi := -1
			for gi := range byFile {
				if occ < len(byFile[gi]) {
					present[gi] = true
					if firstGi < 0 {
						firstGi = gi
					}
				}
			}
			if !nfilesPasses(present, opts.Nfiles) {
				continue
			}
			st.totalKept++
			rep := byFile[firstGi][occ]
			// sites.txt: CHROM POS REF ALT BITS (first present input's record).
			if st.sitesF != nil {
				for gi := 0; gi < n; gi++ {
					if present[gi] {
						bits[gi] = '1'
					} else {
						bits[gi] = '0'
					}
				}
				ref := rep.Ref
				if ref == "" {
					ref = "."
				}
				alt := "."
				if len(rep.Alt) > 0 {
					alt = strings.Join(rep.Alt, ",")
				}
				fmt.Fprintf(st.sitesF, "%s\t%d\t%s\t%s\t%s\n", rep.Chrom, rep.Pos, ref, alt, string(bits))
			}
			// Per-input projection.
			if opts.Prefix != "" {
				if st.venn && present[0] && present[1] {
					// Shared: input-1 copy -> 0002, input-2 copy -> 0003.
					if err := st.perInputW[2].Write(byFile[0][occ]); err != nil {
						return err
					}
					if err := st.perInputW[3].Write(byFile[1][occ]); err != nil {
						return err
					}
				} else {
					for gi := range byFile {
						if present[gi] {
							if err := st.perInputW[gi].Write(byFile[gi][occ]); err != nil {
								return err
							}
						}
					}
				}
			}
			// stdout dump: first requested input that is present at this site.
			if st.openStdout {
				for _, want := range st.writeFromInputs {
					wi := want - 1
					if wi < 0 || wi >= n || !present[wi] {
						continue
					}
					if err := st.stdoutW.Write(byFile[wi][occ]); err != nil {
						return err
					}
					break
				}
			}
		}
	}
	return nil
}

// isecFinalize flushes and closes every writer opened by isecSetup and returns
// the accumulated total kept count.
func isecFinalize(st *isecState) (int, error) {
	if st.openStdout {
		_ = st.stdoutW.Flush()
		st.stdoutFinish()
	}
	for _, w := range st.perInputW {
		_ = w.Flush()
	}
	for _, c := range st.perInputClose {
		c()
	}
	if st.sitesF != nil {
		_ = st.sitesF.Close()
	}
	return st.totalKept, nil
}

// collapseKey returns the canonical membership key for a variant under the
// given collapse mode.
func collapseKey(v *vcf.Variant, mode CollapseMode) string {
	switch mode {
	case CollapseNone:
		return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
	case CollapseAll:
		return fmt.Sprintf("%s\t%d", v.Chrom, v.Pos)
	case CollapseID:
		return v.ID
	case CollapseSNPs:
		if isSNPRecord(v) {
			return fmt.Sprintf("%s\t%d\tSNP", v.Chrom, v.Pos)
		}
		return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
	case CollapseIndels:
		if isIndelRecord(v) {
			return fmt.Sprintf("%s\t%d\tINDEL", v.Chrom, v.Pos)
		}
		return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
	case CollapseBoth:
		switch {
		case isSNPRecord(v):
			return fmt.Sprintf("%s\t%d\tSNP", v.Chrom, v.Pos)
		case isIndelRecord(v):
			return fmt.Sprintf("%s\t%d\tINDEL", v.Chrom, v.Pos)
		}
		return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
	case CollapseSome:
		return fmt.Sprintf("%s\t%d\t%s", v.Chrom, v.Pos, v.Ref)
	}
	return fmt.Sprintf("%s\t%d\t%s\t%s", v.Chrom, v.Pos, v.Ref, strings.Join(v.Alt, ","))
}

// nfilesPasses reports whether mem (the per-input membership bitmap) is
// admissible under spec.
func nfilesPasses(mem []bool, spec NfilesSpec) bool {
	if spec.Mode == 0 {
		// No constraint — keep every key that appears in at least one
		// input (which is always the case).
		return true
	}
	count := 0
	for _, b := range mem {
		if b {
			count++
		}
	}
	switch spec.Mode {
	case '+':
		return count >= spec.N
	case '=':
		return count == spec.N
	case '~':
		if len(spec.Bits) != len(mem) {
			return false
		}
		for i := range mem {
			if mem[i] != spec.Bits[i] {
				return false
			}
		}
		return true
	}
	return true
}

// outputExt returns the filename extension matching the chosen output
// format.
func outputExt(f OutputFormat) string {
	switch f {
	case OutputVCFGz:
		return ".gz"
	case OutputBCF:
		return "" // .vcf.bcf doesn't make sense — upstream uses .bcf, see below
	case OutputBCFUncompressed:
		return ""
	}
	return ""
}

// isecReadme renders the README.txt legend in upstream vcfisec's exact shape:
// a provenance preamble, then one "<outfile>\tfor records private to\t<input>"
// (and, in VENN mode, "<outfile>\tfor records from <input> shared by both\t..."
// ) line per output file. outNames are the full output paths; venn selects the
// four-file Venn legend over the one-file-per-input legend.
func isecReadme(opts IsecOptions, outNames []string, venn bool) string {
	var b strings.Builder
	b.WriteString("This file was produced by vcfisec.\n")
	b.WriteString("The command line was:\tbcftools")
	if len(opts.CmdArgv) > 0 {
		b.WriteString(" ")
		b.WriteString(opts.CmdArgv[0])
		b.WriteString(" ")
		for _, a := range opts.CmdArgv[1:] {
			b.WriteString(" ")
			b.WriteString(a)
		}
	}
	b.WriteString("\n\nUsing the following file names:\n")
	in := opts.InputPaths
	get := func(i int) string {
		if i < len(in) {
			return in[i]
		}
		return fmt.Sprintf("input%d", i+1)
	}
	if venn && len(outNames) == 4 {
		fmt.Fprintf(&b, "%s\tfor records private to\t%s\n", outNames[0], get(0))
		fmt.Fprintf(&b, "%s\tfor records private to\t%s\n", outNames[1], get(1))
		fmt.Fprintf(&b, "%s\tfor records from %s shared by both\t%s %s\n", outNames[2], get(0), get(0), get(1))
		fmt.Fprintf(&b, "%s\tfor records from %s shared by both\t%s %s\n", outNames[3], get(1), get(0), get(1))
	} else {
		for i, name := range outNames {
			fmt.Fprintf(&b, "%s\tfor stripped\t%s\n", name, get(i))
		}
	}
	return b.String()
}
