package samtools

// tview.go implements `samtools tview`'s NON-INTERACTIVE display backends:
// the text (`-d T`) and HTML (`-d H`) "alignment viewers". Both render the
// same in-memory character grid (the "screen") that upstream's
// bam_tview.c::base_draw_aln builds for a region, then serialise it — the
// text backend as a plain character grid, the HTML backend as a coloured
// <pre>/<span> markup. The interactive ncurses backend (`-d C`) is a
// deliberate non-goal here (it needs a TTY UI library, an external
// dependency this project avoids); the CLI emits a clear error pointing the
// user at -d T / -d H.
//
// # Faithfulness to upstream
//
// The layout is a direct port of bam_tview.c:
//
//   - The "screen" is a row-major grid of tixels (character + attribute),
//     allocated lazily one row at a time exactly as bam_tview_html.c's
//     html_mvaddch does. Cells default to a space with attribute 0; a write
//     past mcol columns is a no-op (clipped to the display width).
//   - Row 0 is the ruler: at every 0-based reference position that is a
//     multiple of the marker interval (10 below 1e9, 20 at/above it), the
//     1-based position is printed left-justified into consecutive ruler
//     cells, provided at least 10 columns remain.
//   - Row 1 is the reference line (FASTA base, or 'N' when no reference).
//   - Row 2 is the consensus line, computed by the same bam2bcf
//     bcf_call_glfgen + errmod genotype-likelihood call upstream uses
//     (tv_pl_func, bam_tview.c:191-217).
//   - Rows 3+ are the read rows, packed greedily so non-overlapping reads
//     share a row. Row assignment is the exact level-pool algorithm from
//     bam_lpileup.c::tview_func (TV_GAP=2: a freed level is reusable only
//     after a 2-column gap, lowest reusable level first). A base equal to
//     the reference renders as '.' (forward) / ',' (reverse); a mismatch as
//     the read base (UPPER forward / lower reverse); a deletion as '*'; a
//     reference skip as '>' (forward) / '<' (reverse).
//
// The text backend, when writing to a pipe (the only mode this port
// supports — there is no TTY here), emits just the characters: upstream's
// text_drawaln only inserts ANSI colour escapes when isatty(out) is true,
// so a piped `-d T` is a plain character grid. The HTML backend ports
// bam_tview_html.c::html_drawaln verbatim, including its CSS colour table
// and per-attribute <span> coalescing.
//
// # Insertion columns
//
// Upstream optionally expands insertion sub-columns (the max_ins loop in
// tv_pl_func) unless -i is given. This port does not expand insertion
// sub-columns; the reference/consensus/read grid is rendered at one column
// per reference position. -i is accepted for CLI compatibility. See
// docs/PARITY_ROADMAP.md.

import (
	"fmt"
	"io"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/alnio"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// DefaultTviewWidth is the display width (columns) used for the text and
// HTML tview backends when -w is not given. Upstream's `-d C` curses mode
// derives the width from the terminal; the non-interactive backends instead
// read the COLUMNS environment variable, falling back to 80
// (bam_tview_html.c::html_tv_init). We use 80 as the fixed default so the
// output is deterministic in a pipeline regardless of the environment.
const DefaultTviewWidth = 80

// tvMinAlnRow is upstream's TV_MIN_ALNROW: read rows begin at the grid row
// strictly greater than this, i.e. the first read (level 1) lands on row 3.
const tvMinAlnRow = 2

// tvGap is upstream's TV_GAP: the minimum number of empty columns required
// after a read ends before its row (level) may be reused by another read.
const tvGap = 2

// tvTenDigits is upstream's TEN_DIGITS threshold: below it ruler markers are
// spaced every 10 bp, at/above it every 20 bp (so 10-digit numbers do not
// collide).
const tvTenDigits = 1000000000

// tvErrModDepCorr is the depcorr passed to errmod_init for the consensus
// caller: upstream tview calls bcf_call_init(0.83, 13), which does
// errmod_init(1 - 0.83).
const tvErrModDepCorr = 1.0 - 0.83

// tvMinBaseQ is the consensus caller's per-base quality floor: upstream
// bcf_call_init(0.83, 13).
const tvMinBaseQ = 13

// tvUnderlineFlag is bam_tview_html.c's UNDERLINE_FLAG bit. It is set as an
// attribute on cells that should be underlined (anomalous-pair / secondary
// reads). The colour attributes are 1<<colorpair; the underline shares the
// same attribute word at bit 10.
const tvUnderlineFlag = 10

// seqNt16IntForTview maps a 4-bit IUPAC base code to a 2-bit base index
// (A=0,C=1,G=2,T=3, anything else 4), mirroring htslib's seq_nt16_int.
var seqNt16IntForTview = [16]int{4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4}

// asciiToNt16 maps an ASCII base to its 4-bit IUPAC code, mirroring the
// subset of htslib's seq_nt16_table used by tview's consensus caller. Any
// unmapped byte becomes 15 ('N').
func asciiToNt16(b byte) int {
	switch b {
	case 'A', 'a':
		return 1
	case 'C', 'c':
		return 2
	case 'M', 'm':
		return 3
	case 'G', 'g':
		return 4
	case 'R', 'r':
		return 5
	case 'S', 's':
		return 6
	case 'V', 'v':
		return 7
	case 'T', 't':
		return 8
	case 'W', 'w':
		return 9
	case 'Y', 'y':
		return 10
	case 'H', 'h':
		return 11
	case 'K', 'k':
		return 12
	case 'D', 'd':
		return 13
	case 'B', 'b':
		return 14
	default:
		return 15
	}
}

// TviewMode selects the non-interactive tview backend.
type TviewMode int

const (
	// TviewText is the `-d T` plain-text character-grid backend.
	TviewText TviewMode = iota
	// TviewHTML is the `-d H` coloured HTML backend.
	TviewHTML
)

// TviewOptions configures a tview render. Input/Reference are file paths;
// the region and width pin the displayed window.
type TviewOptions struct {
	// Input is the path to the alignment file (BAM/CRAM/SAM). "-" reads
	// from standard input.
	Input string
	// Reference is the optional reference FASTA path (positional ref.fa or
	// -T/--reference). When empty the reference line shows 'N' for every
	// column. A CRAM input also uses it as the decode reference.
	Reference string
	// Position is the -p region (chr[:pos]) the view starts at. When empty
	// the view starts at the first reference in the header (or, if a
	// reference FASTA is given, the first contig present in both).
	Position string
	// Width is the display width in columns (-w). Zero selects
	// DefaultTviewWidth.
	Width int
	// Sample restricts the displayed reads to those whose @RG ID or SM
	// matches this string (-s). Empty shows all reads.
	Sample string
	// Mode picks the text or HTML backend.
	Mode TviewMode
	// HideInserts suppresses insertion columns (-i). This port never expands
	// insertion sub-columns, so the flag is accepted for compatibility.
	HideInserts bool
}

// tixel is one screen cell: a character and an attribute word. The
// attribute word holds 1<<colorpair colour bits plus the underline bit
// (1<<tvUnderlineFlag), exactly as bam_tview_html.c's tixel_t.
type tixel struct {
	ch   byte
	attr uint32
}

// tvScreen is the lazily-grown row-major character grid. It mirrors the
// html_tview_t screen: rows are allocated on first write, every cell
// defaults to a space, and a write at or beyond mcol columns is clipped.
type tvScreen struct {
	mcol int
	rows [][]tixel
	// attr is the currently-active attribute (set via attron/attroff),
	// applied to every cell written by mvaddch — mirroring html_attron.
	attr uint32
}

// mvaddch writes ch at (y, x) with the current attribute, growing the grid
// to row y as needed. A write at x >= mcol is clipped (upstream
// html_mvaddch returns early).
func (s *tvScreen) mvaddch(y, x int, ch byte) {
	if x < 0 || x >= s.mcol {
		return
	}
	for len(s.rows) <= y {
		row := make([]tixel, s.mcol)
		for i := range row {
			row[i] = tixel{ch: ' ', attr: 0}
		}
		s.rows = append(s.rows, row)
	}
	s.rows[y][x] = tixel{ch: ch, attr: s.attr}
}

// mvprintw writes the decimal string of str starting at (y, x), one char per
// cell, mirroring html_mvprintw's per-character mvaddch loop. (Upstream's
// ruler uses the "%-PRIhts_pos" format, which has no field width and so is
// just the decimal digits.)
func (s *tvScreen) mvprintw(y, x int, str string) {
	for i := 0; i < len(str); i++ {
		s.mvaddch(y, x+i, str[i])
	}
}

// Tview renders the alignment viewer for opts to out. It is the single
// entry point for the text and HTML backends.
func Tview(opts TviewOptions, out io.Writer) error {
	mcol := opts.Width
	if mcol <= 0 {
		mcol = DefaultTviewWidth
	}

	// Open the alignment input (CRAM honours the reference, if given).
	rd, err := alnio.OpenReader(opts.Input, opts.Reference)
	if err != nil {
		return fmt.Errorf("samtools tview: %w", err)
	}
	defer rd.Close()
	hdr := rd.Header()

	// Optional reference FASTA for the reference line.
	var ref *fasta.RandomAccess
	if opts.Reference != "" {
		ref, err = fasta.OpenRandomAccess(opts.Reference)
		if err != nil {
			return fmt.Errorf("samtools tview: open reference %s: %w", opts.Reference, err)
		}
		defer ref.Close()
	}

	// Resolve the @RG ID set for -s (sample/read-group filtering).
	var rgSet map[string]struct{}
	if opts.Sample != "" {
		rgSet = tviewSampleReadGroups(hdr, opts.Sample)
		if len(rgSet) == 0 {
			return fmt.Errorf("samtools tview: the sample or read group %q is not present", opts.Sample)
		}
	}

	// Decide the displayed contig and left (0-based) position.
	chrom, leftPos0, err := tviewResolveStart(hdr, ref, opts.Position)
	if err != nil {
		return err
	}

	// Read every record on the contig that overlaps the display window
	// [leftPos0, leftPos0+mcol). alnio is a streaming reader (no index
	// iterator), so we scan and filter — the window is small (mcol columns),
	// so the retained set is bounded by coverage over those columns.
	recs, err := tviewCollectRecords(rd, chrom, leftPos0, leftPos0+mcol, rgSet)
	if err != nil {
		return err
	}

	// Fetch the reference slab for the window (clamped to the contig).
	var refSlab []byte
	if ref != nil {
		contigLen := ref.Length(chrom)
		end := int64(leftPos0 + mcol)
		if contigLen > 0 && end > contigLen {
			end = contigLen
		}
		if int64(leftPos0) < end {
			b, ferr := ref.Fetch(chrom, int64(leftPos0), end)
			if ferr != nil {
				return fmt.Errorf("samtools tview: fetch reference %s:%d-%d: %w", chrom, leftPos0, end, ferr)
			}
			refSlab = b
		}
	}

	em := errmod.Init(tvErrModDepCorr)
	screen := drawTviewScreen(recs, refSlab, leftPos0, mcol, em, opts.HideInserts)

	switch opts.Mode {
	case TviewHTML:
		return renderTviewHTML(out, screen, chrom, leftPos0)
	default:
		return renderTviewText(out, screen)
	}
}

// tviewResolveStart picks the displayed contig and 0-based left position.
// With a -p region it parses chrom:pos (1-based pos → 0-based left). Without
// one it starts at the first contig (or, with a reference, the first contig
// present in both the BAM header and the FASTA), position 0.
func tviewResolveStart(hdr *sam.Header, ref *fasta.RandomAccess, position string) (string, int, error) {
	if position != "" {
		reg, perr := region.ParseRegion(position)
		if perr != nil {
			return "", 0, fmt.Errorf("samtools tview: unknown reference or malformed region: %w", perr)
		}
		if hdr.RefIndex(reg.Chrom) < 0 {
			return "", 0, fmt.Errorf("samtools tview: unknown reference %q", reg.Chrom)
		}
		// region.Beg is 1-based inclusive; the 0-based left position is Beg-1.
		left := reg.Beg - 1
		if left < 0 {
			left = 0
		}
		return reg.Chrom, left, nil
	}
	if ref != nil {
		// First contig present in both the header and the reference index.
		for _, r := range hdr.Refs {
			if _, ok := ref.Index().Get(r.Name); ok {
				return r.Name, 0, nil
			}
		}
		return "", 0, fmt.Errorf("samtools tview: none of the BAM sequence names are present in the fasta file")
	}
	if len(hdr.Refs) == 0 {
		return "", 0, fmt.Errorf("samtools tview: no reference sequences in header")
	}
	return hdr.Refs[0].Name, 0, nil
}

// tviewSampleReadGroups builds the set of @RG IDs whose ID or SM matches
// sample, mirroring bam_tview.c::get_rg_sample. An empty result means the
// sample/read-group was not found.
func tviewSampleReadGroups(hdr *sam.Header, sample string) map[string]struct{} {
	set := map[string]struct{}{}
	for _, rg := range hdr.ReadGroups {
		if rg.ID == sample {
			set[rg.ID] = struct{}{}
			continue
		}
		for _, f := range rg.Extra {
			if f.Tag == "SM" && f.Value == sample {
				set[rg.ID] = struct{}{}
				break
			}
		}
	}
	return set
}

// tviewCollectRecords scans rd and returns the records on chrom whose
// aligned span overlaps the half-open 0-based window [beg0, end0), in
// coordinate (then input) order — the order htslib's pileup walker visits
// them. Records are filtered by the @RG set when rgSet is non-nil; unmapped
// records and records on other contigs are skipped.
func tviewCollectRecords(rd sam.Reader, chrom string, beg0, end0 int, rgSet map[string]struct{}) ([]*sam.Record, error) {
	var recs []*sam.Record
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("samtools tview: read: %w", err)
		}
		if rec.IsUnmapped() || rec.RName != chrom || rec.Pos <= 0 {
			continue
		}
		if rgSet != nil {
			rg, ok := rec.GetAux("RG")
			if !ok {
				continue
			}
			rgID, _ := rg.String()
			if _, keep := rgSet[rgID]; !keep {
				continue
			}
		}
		start0 := int(rec.Pos) - 1
		// EndPosition is the 1-based last reference position; its 0-based
		// exclusive end equals that value.
		endExcl := int(rec.EndPosition())
		if start0 >= end0 || endExcl <= beg0 {
			continue
		}
		recs = append(recs, rec)
	}
	// Stable sort by start position so the pileup column order matches
	// htslib's coordinate-ordered iterator.
	sort.SliceStable(recs, func(i, j int) bool {
		return recs[i].Pos < recs[j].Pos
	})
	return recs, nil
}

// tvColumnRead is one read's contribution to one pileup column.
type tvColumnRead struct {
	recIdx    int    // index into the record slice
	qpos      int    // 0-based query offset of the base at this column (valid when !isDel)
	isDel     bool   // D op covers this column
	isRefSkip bool   // N op covers this column
	isHead    bool   // first column this read appears in
	isTail    bool   // last column this read appears in
	insAfter  string // inserted bases immediately following this base (the bam_plp_insertion result), "" if none
	level     int    // row level assigned by the packing pass (1-based)
}

// tviewColumn is one pileup column: a 0-based reference position and the
// reads overlapping it, in coordinate order.
type tviewColumn struct {
	pos0  int
	reads []tvColumnRead
}

// drawTviewScreen builds the full screen grid: ruler, reference, consensus,
// and the packed read rows, for the window [leftPos0, leftPos0+mcol). When
// hideInserts is false (the upstream default) insertion sub-columns are
// expanded; -i sets it true to suppress them.
func drawTviewScreen(recs []*sam.Record, refSlab []byte, leftPos0, mcol int, em *errmod.Errmod, hideInserts bool) *tvScreen {
	s := &tvScreen{mcol: mcol}

	// Build per-column read lists. Upstream's region iterator returns every
	// read overlapping [leftPos0, leftPos0+mcol) and pushes it through the
	// pileup IN FULL, so the level-assignment walk runs over each read's
	// entire span — from the earliest read start (which may sit BEFORE the
	// display window) to the latest read end (which may sit AFTER it). The
	// out-of-window columns are never drawn (tv_pl_func's `pos < left_pos`
	// and `ccol > mcol` guards) but they shape is_head/is_tail and thus the
	// levels of the reads that ARE visible. We therefore build columns over
	// the full [minStart0, maxEnd0) span, assign levels over all of them, and
	// draw only the columns inside [leftPos0, leftPos0+mcol).
	buildBeg0 := leftPos0
	buildEnd0 := leftPos0 + mcol
	for _, rec := range recs {
		if s0 := int(rec.Pos) - 1; s0 < buildBeg0 {
			buildBeg0 = s0
		}
		if e0 := int(rec.EndPosition()); e0 > buildEnd0 {
			buildEnd0 = e0
		}
	}
	columns := buildTviewColumns(recs, buildBeg0, buildEnd0)

	// Assign row levels with the greedy free-slot pool (bam_lpileup.c).
	assignTviewLevels(columns)

	// last_pos tracks the previous emitted 0-based position so the ruler
	// fill between coverage gaps is reproduced (base_draw_aln). It starts at
	// leftPos0-1.
	lastPos := leftPos0 - 1
	ccol := 0

	emitRuler := func(pos0 int) {
		interval := 10
		if pos0 >= tvTenDigits {
			interval = 20
		}
		if pos0%interval == 0 && mcol-ccol >= 10 {
			s.mvprintw(0, ccol, strconv.Itoa(pos0+1))
		}
	}
	refBaseAt := func(pos0 int) byte {
		off := pos0 - leftPos0
		if refSlab != nil && off >= 0 && off < len(refSlab) {
			return refSlab[off]
		}
		return 'N'
	}

	for _, col := range columns {
		pos0 := col.pos0
		if pos0 < leftPos0 || ccol > mcol {
			continue
		}
		// Fill the gap of uncovered columns up to this position (the
		// `for cp = last_pos+1; cp < pos` loop in tv_pl_func).
		for cp := lastPos + 1; cp < pos0; cp++ {
			emitRuler(cp)
			s.attr = 0
			s.mvaddch(1, ccol, refBaseAt(cp))
			ccol++
		}
		// Ruler marker at this covered position.
		emitRuler(pos0)

		rb := refBaseAt(pos0)
		// Consensus call for this column. It is written ONCE at the j==0
		// (leftmost) column position, before the j-loop advances ccol
		// (tv_pl_func writes the consensus to row 2 at tv->ccol prior to the
		// insertion loop).
		consChar, consAttr := tviewConsensusCall(col.reads, recs, rb, em)
		s.attr = consAttr
		s.mvaddch(2, ccol, consChar)
		s.attr = 0

		// Maximum insertion length across the column's reads (unless -i).
		maxIns := 0
		if !hideInserts {
			for _, r := range col.reads {
				if len(r.insAfter) > maxIns {
					maxIns = len(r.insAfter)
				}
			}
		}

		// Core loop: j==0 is the base column; j>0 are inserted sub-columns.
		for j := 0; j <= maxIns; j++ {
			for _, r := range col.reads {
				row := tvMinAlnRow + r.level
				if row <= tvMinAlnRow {
					continue
				}
				rec := recs[r.recIdx]
				isRev := rec.Flag&sam.FlagReverse != 0
				var c byte
				if j == 0 {
					switch {
					case !r.isDel && !r.isRefSkip:
						var base byte = 'N'
						if r.qpos >= 0 && r.qpos < len(rec.Seq) {
							base = rec.Seq[r.qpos]
						}
						if upper(base) == upper(rb) {
							if isRev {
								c = ','
							} else {
								c = '.'
							}
						} else {
							c = base
						}
					case r.isRefSkip:
						if isRev {
							c = '<'
						} else {
							c = '>'
						}
					default:
						c = '*'
					}
				} else {
					// Inserted sub-column. A read with a shorter (or no)
					// insertion pads with '*'; otherwise it shows its j-th
					// inserted base verbatim (never dotted — tv_pl_func's dot
					// rule is gated on j==0, which never holds here).
					if j > len(r.insAfter) {
						c = '*'
					} else {
						c = r.insAfter[j-1]
					}
				}
				// Attribute: underline for anomalous-pair / secondary reads,
				// colour by MAPQ (the upstream default TV_COLOR_MAPQ). The
				// colour uses the base column's MAPQ for every sub-column.
				var attr uint32
				if (rec.Flag&sam.FlagPaired != 0 && rec.Flag&sam.FlagProperPair == 0) ||
					rec.Flag&sam.FlagSecondary != 0 {
					attr |= 1 << tvUnderlineFlag
				}
				x := int(rec.MapQ)/10 + 1
				if x > 4 {
					x = 4
				}
				attr |= 1 << uint(x)
				s.attr = attr
				if isRev {
					s.mvaddch(row, ccol, lower(c))
				} else {
					s.mvaddch(row, ccol, upper(c))
				}
				s.attr = 0
			}
			// Reference row: the real base at j==0, '*' (red, colorpair 8)
			// for each inserted sub-column.
			if j == 0 {
				s.attr = 0
				s.mvaddch(1, ccol, rb)
			} else {
				s.attr = 1 << 8
				s.mvaddch(1, ccol, '*')
				s.attr = 0
			}
			ccol++
		}
		lastPos = pos0
	}

	// Trailing reference / ruler fill out to mcol (base_draw_aln's final
	// while loop).
	for ccol < mcol {
		pos0 := lastPos + 1
		emitRuler(pos0)
		s.attr = 0
		s.mvaddch(1, ccol, refBaseAt(pos0))
		ccol++
		lastPos++
	}
	return s
}

// buildTviewColumns produces the per-column read lists for the window,
// reproducing htslib's pileup: a read occupies every reference position from
// its start through its last reference-consuming op (D and N gaps included),
// with is_head/is_tail at the first/last such positions and is_del/is_refskip
// flagged on gap columns. Columns are returned in ascending position order;
// only columns with at least one read are present.
func buildTviewColumns(recs []*sam.Record, beg0, end0 int) []tviewColumn {
	byPos := map[int][]tvColumnRead{}
	for idx, rec := range recs {
		appendTviewReadColumns(rec, idx, beg0, end0, byPos)
	}
	if len(byPos) == 0 {
		return nil
	}
	positions := make([]int, 0, len(byPos))
	for p := range byPos {
		positions = append(positions, p)
	}
	sort.Ints(positions)
	cols := make([]tviewColumn, 0, len(positions))
	for _, p := range positions {
		reads := byPos[p]
		// Within a column, preserve coordinate (then input) order so the
		// level-assignment pass sees reads in the same order htslib does.
		sort.SliceStable(reads, func(i, j int) bool {
			ri, rj := recs[reads[i].recIdx], recs[reads[j].recIdx]
			if ri.Pos != rj.Pos {
				return ri.Pos < rj.Pos
			}
			return reads[i].recIdx < reads[j].recIdx
		})
		cols = append(cols, tviewColumn{pos0: p, reads: reads})
	}
	return cols
}

// appendTviewReadColumns walks rec's CIGAR and appends a tvColumnRead to
// byPos[pos0] for every reference position the read covers within
// [beg0, end0). is_head / is_tail are set on the read's first / last covered
// reference position.
func appendTviewReadColumns(rec *sam.Record, recIdx, beg0, end0 int, byPos map[int][]tvColumnRead) {
	refPos := int(rec.Pos) - 1
	queryPos := 0

	// First pass: collect every covered (pos0, qpos, isDel, isRefSkip) so we
	// know the read's overall first/last reference position for head/tail.
	// An insertion (I op) attaches its bases to the preceding covered base as
	// insAfter, reproducing bam_plp_insertion (which reports the insertion
	// that follows the current pileup base).
	type cov struct {
		pos0      int
		qpos      int
		isDel     bool
		isRefSkip bool
		insAfter  string
	}
	var covs []cov
	for _, op := range rec.Cigar {
		l := int(op.Length())
		switch op.Op() {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			for k := 0; k < l; k++ {
				covs = append(covs, cov{pos0: refPos + k, qpos: queryPos + k})
			}
			refPos += l
			queryPos += l
		case sam.CigarDeletion:
			for k := 0; k < l; k++ {
				covs = append(covs, cov{pos0: refPos + k, qpos: queryPos, isDel: true})
			}
			refPos += l
		case sam.CigarSkipped:
			for k := 0; k < l; k++ {
				covs = append(covs, cov{pos0: refPos + k, qpos: queryPos, isRefSkip: true})
			}
			refPos += l
		case sam.CigarInsertion:
			// Attach the inserted bases to the preceding covered base.
			if len(covs) > 0 && queryPos+l <= len(rec.Seq) {
				covs[len(covs)-1].insAfter = rec.Seq[queryPos : queryPos+l]
			}
			queryPos += l
		case sam.CigarSoftClip:
			queryPos += l
		case sam.CigarHardClip, sam.CigarPadding:
			// Consume neither reference nor (for our accounting) query.
		}
	}
	if len(covs) == 0 {
		return
	}
	// is_head / is_tail are computed over the WINDOW-CLIPPED covered
	// positions, not the read's full span: upstream resets the pileup buffer
	// at left_pos and starts the region iterator there, so a read overlapping
	// the left edge gets is_head at the first column inside the window (and
	// is_tail at the last column inside it). This matches sam_itr_queryi +
	// bam_lplbuf_reset in base_draw_aln.
	firstIn, lastIn := -1, -1
	for i, c := range covs {
		if c.pos0 < beg0 || c.pos0 >= end0 {
			continue
		}
		if firstIn == -1 {
			firstIn = i
		}
		lastIn = i
	}
	if firstIn == -1 {
		return
	}
	for i, c := range covs {
		if c.pos0 < beg0 || c.pos0 >= end0 {
			continue
		}
		byPos[c.pos0] = append(byPos[c.pos0], tvColumnRead{
			recIdx:    recIdx,
			qpos:      c.qpos,
			isDel:     c.isDel,
			isRefSkip: c.isRefSkip,
			isHead:    i == firstIn,
			isTail:    i == lastIn,
			insAfter:  c.insAfter,
		})
	}
}

// freenode mirrors bam_lpileup.c's freenode_t: a node in the singly-linked
// list of freed row levels. level is the row index; cnt is the reuse
// countdown. The list is terminated by a sentinel node (sentinel == true),
// and is kept ordered by (cnt, level) ascending — the splaysort order in
// tview_func. Crucially, a freed level inherits the cnt of the *sentinel*
// node that held it: a sentinel freshly allocated by mp_alloc (calloc) has
// cnt 0, so the level it later carries is reusable with no gap; a sentinel
// recycled through mp_free has cnt = TV_GAP, enforcing the 2-column gap.
type freenode struct {
	level    int
	cnt      int
	next     *freenode
	sentinel bool
}

// tvLevelMempool mirrors bam_lpileup.c's mempool_t: it recycles freenodes.
// alloc returns a pool node (cnt = TV_GAP, set by free) or a fresh node
// (cnt = 0); free returns a node to the pool with cnt = TV_GAP.
type tvLevelMempool struct {
	buf []*freenode
}

func (mp *tvLevelMempool) alloc() *freenode {
	if len(mp.buf) == 0 {
		return &freenode{} // calloc: cnt 0
	}
	n := mp.buf[len(mp.buf)-1]
	mp.buf = mp.buf[:len(mp.buf)-1]
	return n
}

func (mp *tvLevelMempool) free(p *freenode) {
	p.next = nil
	p.cnt = tvGap
	p.sentinel = false
	mp.buf = append(mp.buf, p)
}

// assignTviewLevels assigns a 1-based row level to every read occurrence in
// every column, a faithful port of bam_lpileup.c::tview_func including its
// linked-list + mempool mechanics. Per column it:
//
//  1. decrements every non-sentinel node's reuse countdown toward 0;
//  2. for each read in order: an is_head read reuses the list head's level
//     when the head is a real node with cnt 0 (freeing that node), else
//     allocates a fresh ++maxLevel; a continuing read inherits its level; an
//     is_tail read parks its level on the current sentinel node and allocates
//     a new sentinel (so the freed level carries the old sentinel's cnt —
//     0 for the very first end, TV_GAP after a recycle);
//  3. discards nodes whose level exceeds the column maximum, resets maxLevel
//     to that maximum, and re-sorts the list by (cnt, level).
func assignTviewLevels(columns []tviewColumn) {
	mp := &tvLevelMempool{}
	sentinel := mp.alloc()
	sentinel.sentinel = true
	head := sentinel
	tail := sentinel
	maxLevel := 0
	// levelOf maps a still-open read (by recIdx) to its assigned level so
	// non-head occurrences inherit it.
	levelOf := map[int]int{}

	for ci := range columns {
		col := &columns[ci]
		// (1) age every non-sentinel node.
		for p := head; p.next != nil; p = p.next {
			if p.cnt > 0 {
				p.cnt--
			}
		}
		colMax := 0
		for ri := range col.reads {
			r := &col.reads[ri]
			if r.isHead {
				var lvl int
				if head.next != nil && head.cnt == 0 {
					// Reuse the head's level and free the node.
					lvl = head.level
					next := head.next
					mp.free(head)
					head = next
				} else {
					maxLevel++
					lvl = maxLevel
				}
				r.level = lvl
				levelOf[r.recIdx] = lvl
			} else {
				r.level = levelOf[r.recIdx]
			}
			if r.level > colMax {
				colMax = r.level
			}
			if r.isTail {
				// Park this level on the current sentinel; allocate a new
				// sentinel. The parked node keeps the sentinel's cnt.
				tail.level = r.level
				tail.sentinel = false
				newSentinel := mp.alloc()
				newSentinel.sentinel = true
				newSentinel.next = nil
				tail.next = newSentinel
				tail = newSentinel
				delete(levelOf, r.recIdx)
			}
		}
		// (3) collect real (non-sentinel) nodes, discard level > colMax,
		// re-sort by (cnt, level), relink with the sentinel at the end.
		var nodes []*freenode
		for p := head; p.next != nil; {
			next := p.next
			if p.level > colMax {
				mp.free(p)
			} else {
				nodes = append(nodes, p)
			}
			p = next
		}
		sort.SliceStable(nodes, func(i, j int) bool {
			if nodes[i].cnt != nodes[j].cnt {
				return nodes[i].cnt < nodes[j].cnt
			}
			return nodes[i].level < nodes[j].level
		})
		// Relink: nodes... -> tail (sentinel).
		head = tail
		for i := len(nodes) - 1; i >= 0; i-- {
			nodes[i].next = head
			head = nodes[i]
		}
		maxLevel = colMax
	}
}

// tviewConsensusCall computes the consensus character and its colour
// attribute for one column, porting bam_tview.c::tv_pl_func's consensus block
// (the bcf_call_glfgen + errmod genotype-likelihood call). reads are the
// column's per-read occurrences; rb is the reference base (ASCII).
func tviewConsensusCall(reads []tvColumnRead, recs []*sam.Record, rb byte, em *errmod.Errmod) (byte, uint32) {
	refNt16 := asciiToNt16(rb)

	// Build the packed base array glfgen feeds to errmod and the qsum[4]
	// per-base quality sums (bam2bcf.c::bcf_call_glfgen, SNP branch).
	bases := make([]uint16, 0, len(reads))
	var qsum [4]float64
	for _, r := range reads {
		if r.isDel || r.isRefSkip {
			continue
		}
		rec := recs[r.recIdx]
		if rec.Flag&sam.FlagUnmapped != 0 {
			continue
		}
		mapQ := int(rec.MapQ)
		if rec.MapQ == 255 {
			mapQ = 20 // DEF_MAPQ
		}
		q := 0
		if r.qpos >= 0 && r.qpos < len(rec.Qual) {
			q = int(rec.Qual[r.qpos])
		}
		if q < tvMinBaseQ {
			continue
		}
		// seqQ = 99 for SNPs, so q = min(q, seqQ) is a no-op. Cap by mapQ
		// (capQ = 60) then clamp to [4, 63].
		if mapQ > 60 {
			mapQ = 60
		}
		if q > mapQ {
			q = mapQ
		}
		if q > 63 {
			q = 63
		}
		if q < 4 {
			q = 4
		}
		var b int
		if r.qpos >= 0 && r.qpos < len(rec.Seq) {
			code := asciiToNt16(rec.Seq[r.qpos])
			if code == 0 {
				code = refNt16
			}
			b = seqNt16IntForTview[code&0xf]
		} else {
			b = 4 // N
		}
		var strand uint16
		if rec.Flag&sam.FlagReverse != 0 {
			strand = 1
		}
		bases = append(bases, uint16(q)<<5|strand<<4|uint16(b))
		if b < 4 {
			qsum[b] += float64(q)
		}
	}

	// errmod computes the 5x5 phred-scaled genotype-likelihood matrix p.
	p := make([]float32, 25)
	em.Cal(bases, 5, p)

	// Pick the call exactly as tv_pl_func: sort qsum descending, score the
	// top two alleles with the het prior, take the most likely genotype.
	var qs [4]int
	for i := 0; i < 4; i++ {
		qs[i] = int(qsum[i])<<2 | i
	}
	for i := 1; i < 4; i++ {
		for j := i; j > 0 && qs[j] > qs[j-1]; j-- {
			qs[j], qs[j-1] = qs[j-1], qs[j]
		}
	}
	a1 := qs[0] & 3
	a2 := qs[1] & 3
	const prior = 30.0
	pp := [3]float64{
		float64(p[a1*5+a1]),
		float64(p[a1*5+a2]) + prior,
		float64(p[a2*5+a2]),
	}
	rbUpper := upper(rb)
	if "ACGT"[a1] != rbUpper {
		pp[0] += prior + 3
	}
	if "ACGT"[a2] != rbUpper {
		pp[2] += prior + 3
	}
	var call uint32
	switch {
	case pp[0] < pp[1] && pp[0] < pp[2]:
		min := pp[1]
		if pp[2] < min {
			min = pp[2]
		}
		call = uint32(1<<uint(a1))<<16 | uint32(int(min-pp[0]+0.499))
	case pp[2] < pp[1] && pp[2] < pp[0]:
		min := pp[0]
		if pp[1] < min {
			min = pp[1]
		}
		call = uint32(1<<uint(a2))<<16 | uint32(int(min-pp[2]+0.499))
	default:
		min := pp[0]
		if pp[2] < min {
			min = pp[2]
		}
		call = uint32(1<<uint(a1)|1<<uint(a2))<<16 | uint32(int(min-pp[1]+0.499))
	}

	c := ",ACMGRSVTWYHKDBN"[call>>16&0xf]
	i := int(call&0xffff)/10 + 1
	if i > 4 {
		i = 4
	}
	// tv_pl_func always underlines the consensus character
	// (attr = my_underline(tv); attr |= my_colorpair(tv, i)).
	attr := uint32(1)<<tvUnderlineFlag | uint32(1)<<uint(i)
	if c == rbUpper {
		c = '.'
	}
	return c, attr
}

// renderTviewText writes the screen as a plain character grid: each row's
// mcol cells followed by a newline. This matches upstream text_drawaln when
// the output is NOT a terminal (no ANSI escapes), which is the only mode this
// pipeline-oriented port supports.
func renderTviewText(out io.Writer, s *tvScreen) error {
	buf := make([]byte, 0, (s.mcol+1)*len(s.rows))
	for _, row := range s.rows {
		for x := 0; x < s.mcol; x++ {
			buf = append(buf, row[x].ch)
		}
		buf = append(buf, '\n')
	}
	_, err := out.Write(buf)
	return err
}

// tviewCSSColors is bam_tview_html.c's colour table, indexed by colorpair
// id 0..9.
var tviewCSSColors = [10]string{
	"black", "blue", "green", "yellow", "black",
	"green", "cyan", "yellow", "red", "blue",
}

// renderTviewHTML writes the screen as the coloured HTML document
// bam_tview_html.c::html_drawaln emits: a styled <pre> block where each run
// of equal-attribute cells is wrapped in a <span class='tviewc[u]N'>.
func renderTviewHTML(out io.Writer, s *tvScreen, chrom string, leftPos0 int) error {
	w := &htmlWriter{out: out}
	pos1 := leftPos0 + 1
	w.s("<html><head>")
	w.s("<title>")
	w.s(chrom)
	w.s(":")
	w.s(strconv.Itoa(pos1))
	w.s("</title>")
	w.s("<style type='text/css'>\n")
	w.s(".tviewbody { margin:5px; background-color:white;text-align:center;}\n")
	w.s(".tviewtitle {text-align:center;}\n")
	w.s(".tviewpre { margin:5px; background-color:white;}\n")
	for id, col := range tviewCSSColors {
		w.s(".tviewc")
		w.s(strconv.Itoa(id))
		w.s(" {color:")
		w.s(col)
		w.s(";}\n.tviewcu")
		w.s(strconv.Itoa(id))
		w.s(" {color:")
		w.s(col)
		w.s(";text-decoration:underline;}\n")
	}
	w.s("</style>")
	w.s("</head><body>")
	w.s("<div class='tviewbody'><div class='tviewtitle'>")
	w.s(chrom)
	w.s(":")
	w.s(strconv.Itoa(pos1))
	w.s("</div>")
	w.s("<pre class='tviewpre'>")

	for y, row := range s.rows {
		for x := 0; x < s.mcol; x++ {
			cell := row[x]
			if x == 0 || cell.attr != row[x-1].attr {
				w.s("<span")
				// The lowest set colour bit selects the css class; the
				// underline bit adds the 'u' infix.
				css := tviewLowestBit(cell.attr)
				if css >= 0 && css < 32 {
					w.s(" class='tviewc")
					if cell.attr&(1<<tvUnderlineFlag) != 0 {
						w.s("u")
					}
					w.s(strconv.Itoa(css))
					w.s("'")
				}
				w.s(">")
			}
			switch cell.ch {
			case '<':
				w.s("&lt;")
			case '>':
				w.s("&gt;")
			case '&':
				w.s("&amp;")
			default:
				w.b(cell.ch)
			}
			if x+1 == s.mcol || cell.attr != row[x+1].attr {
				w.s("</span>")
			}
		}
		if y+1 < len(s.rows) {
			w.s("<br/>")
		}
	}
	w.s("</pre></div></body></html>")
	return w.err
}

// tviewLowestBit returns the index of the lowest set bit of attr (matching
// the html_drawaln `while(css<32)` scan), or -1 when no bit is set.
func tviewLowestBit(attr uint32) int {
	for css := 0; css < 32; css++ {
		if attr&(1<<uint(css)) != 0 {
			return css
		}
	}
	return -1
}

// htmlWriter is a small error-collecting writer so the HTML emitter stays
// terse.
type htmlWriter struct {
	out io.Writer
	err error
}

func (w *htmlWriter) s(str string) {
	if w.err != nil {
		return
	}
	_, w.err = io.WriteString(w.out, str)
}

func (w *htmlWriter) b(c byte) {
	if w.err != nil {
		return
	}
	_, w.err = w.out.Write([]byte{c})
}
