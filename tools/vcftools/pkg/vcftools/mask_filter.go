// FASTA-like mask filter for vcftools.
//
// `--mask FILE` and `--invert-mask FILE` accept a FASTA-style file where each
// record holds a chromosome name on a `>NAME` header line followed by one or
// more lines of digit characters (one per reference base, 1-based). Each digit
// is interpreted as the integer `char - '0'`; a site at (CHROM, POS) is kept
// when its mask digit is `<= --mask-min INT` (default 0). `--invert-mask`
// flips the keep/drop decision. `--mask-min` must satisfy 0 <= V <= 9 (upstream
// errors when V > 9; we mirror that and additionally validate the lower bound
// for clarity since negative values silently drop every site upstream).
//
// Mirrors upstream `entry::filter_sites_by_mask`
// (reference_code/vcftools/src/cpp/entry_filters.cpp:674-752). Important
// upstream behaviours we preserve verbatim:
//
//  1. The mask file is read in a single forward pass. The mask state advances
//     only forwards; if a VCF record names a chromosome that has already been
//     consumed (e.g. VCF records arrive out of chromosomal order relative to
//     the mask), the record is dropped. This matches upstream's stateful
//     `ifstream` walk and is observable: a VCF reordered to put chr2 before
//     chr1 against a mask listing chr1 then chr2 will lose the chr1 sites.
//  2. When the mask advances onto a new chromosome header BEFORE reaching the
//     requested POS, the current site is dropped (upstream's `passed_filters
//     = false` at entry_filters.cpp:729).
//  3. When the requested POS is past the end of the mask sequence for the
//     current chromosome AND EOF is hit, the site is dropped.
//  4. Header parsing tokenises on the first space/tab after `>`, matching
//     upstream's `line.substr(1, line.find_first_of(" \t")-1)` — comments
//     after the chromosome name are ignored.
//
// We diverge from the upstream `ifstream`-based incremental reader for
// implementation simplicity: we load the entire mask file into memory as a
// per-chromosome list of (offset, line) slabs and walk a cursor over it. The
// observable behaviour (including ordering quirks (1)/(2)/(3)) is preserved
// because the cursor only ever moves forward and is keyed by (chromIdx,
// slabIdx) — exactly the upstream `mask` state machine.
package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// maskSlab is one continuous line of mask digits for a chromosome along with
// the 1-based reference position of the first character.
type maskSlab struct {
	startPos int    // 1-based start position of slab[0] on the chromosome
	line     string // raw mask characters (digits 0-9)
}

// maskChromosome is the ordered list of slabs for a single chromosome plus the
// chromosome name (kept for parity-error reporting).
type maskChromosome struct {
	name  string
	slabs []maskSlab
}

// maskFilter is the streaming state for the mask filter. It is invoked from
// the hot loop with one VCF site at a time and advances its cursor forward
// only. The cursor mirrors upstream's `(mask_chr, mask_line, mask_pos)`
// triplet: chromIdx is the index into `chroms`, slabIdx is the index into the
// current chromosome's slabs, and exhausted is true once we've walked past
// every chromosome (matching upstream's `mask.eof()`).
type maskFilter struct {
	chroms    []maskChromosome
	chromIdx  int
	slabIdx   int
	exhausted bool
	invert    bool
	minKept   int
}

// loadMaskFilter parses a FASTA-style mask file. Lines that aren't preceded by
// a `>CHROM` header are silently dropped (upstream similarly never reads them
// once `mask_chr == ""` — they don't match any VCF site). The chromosome list
// preserves input order so the cursor can walk forward.
func loadMaskFilter(filename string, invert bool, minKept int) (*maskFilter, error) {
	if minKept < 0 || minKept > 9 {
		return nil, fmt.Errorf("--mask-min must be between 0 and 9 (got %d)", minKept)
	}
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, fmt.Errorf("opening mask file %s: %w", filename, err)
	}
	defer f.Close()

	chroms, err := parseMaskFile(f)
	if err != nil {
		return nil, fmt.Errorf("parsing mask file %s: %w", filename, err)
	}

	return &maskFilter{
		chroms:    chroms,
		chromIdx:  -1, // not yet positioned on any chromosome (matches upstream's mask_chr == "" start)
		slabIdx:   0,
		exhausted: len(chroms) == 0,
		invert:    invert,
		minKept:   minKept,
	}, nil
}

// parseMaskFile reads a FASTA-style mask file into per-chromosome slabs in
// input order. Trailing whitespace on each line is trimmed (upstream:
// `line.erase(line.find_last_not_of(" \t") + 1)`).
func parseMaskFile(r io.Reader) ([]maskChromosome, error) {
	scanner := bufio.NewScanner(r)
	// Allow long mask lines (a chromosome on a single line could easily exceed
	// the default 64KiB buffer for short references that emit one big slab).
	scanner.Buffer(make([]byte, 0, 1<<20), 1<<24)

	var chroms []maskChromosome
	var cur *maskChromosome
	curPos := 1
	for scanner.Scan() {
		line := strings.TrimRight(scanner.Text(), " \t\r")
		if line == "" {
			continue
		}
		if line[0] == '>' {
			// Header — split on first whitespace, like upstream.
			rest := line[1:]
			if idx := strings.IndexAny(rest, " \t"); idx >= 0 {
				rest = rest[:idx]
			}
			chroms = append(chroms, maskChromosome{name: rest})
			cur = &chroms[len(chroms)-1]
			curPos = 1
			continue
		}
		if cur == nil {
			// Mask data before any `>CHROM` — upstream never matches such
			// lines to any site (mask_chr stays ""), so we drop silently.
			continue
		}
		cur.slabs = append(cur.slabs, maskSlab{startPos: curPos, line: line})
		curPos += len(line)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return chroms, nil
}

// passes reports whether the site at (chrom, pos) is kept by the mask filter.
// Advances the internal cursor forward. After the cursor has walked past the
// end of the chromosome list (exhausted == true), every subsequent site is
// dropped, mirroring upstream's behaviour after `mask.eof()`.
//
// pos is the 1-based VCF POS. chrom is the VCF CHROM string.
func (m *maskFilter) passes(chrom string, pos int) bool {
	if m == nil {
		return true
	}
	if m.exhausted {
		return false
	}

	// Phase 1: advance chromIdx until either we hit a matching chromosome
	// header or run out of headers.
	for m.chromIdx == -1 || m.chroms[m.chromIdx].name != chrom {
		m.chromIdx++
		m.slabIdx = 0
		if m.chromIdx >= len(m.chroms) {
			m.exhausted = true
			return false
		}
	}

	// Phase 2: within the current chromosome, advance slabIdx until the
	// current slab's range covers pos, or we run off the end of this
	// chromosome's slabs. Mirroring upstream: if we walk to a `>` header for a
	// DIFFERENT chromosome before finding pos, drop. If we exhaust the file
	// entirely, drop.
	chromSlabs := m.chroms[m.chromIdx].slabs
	for m.slabIdx < len(chromSlabs) {
		slab := chromSlabs[m.slabIdx]
		slabEnd := slab.startPos + len(slab.line) - 1 // inclusive
		if pos > slabEnd {
			m.slabIdx++
			continue
		}
		if pos < slab.startPos {
			// VCF position is before the current slab — upstream's cursor
			// never moves backwards, so this only happens when records are
			// out-of-order on the same chromosome. Upstream would have moved
			// past this position already on a previous call and returns false.
			return false
		}
		// pos is within slab.
		digit := int(slab.line[pos-slab.startPos]) - '0'
		keep := digit <= m.minKept
		if m.invert {
			keep = !keep
		}
		return keep
	}

	// Ran off the end of this chromosome's slabs without finding pos. Upstream
	// would either advance to the next `>CHROM` header (and set
	// passed_filters = false because the chromosome changed before pos was
	// covered) or hit EOF and also drop. We advance the cursor onto the next
	// chromosome — but never reset slabIdx onto a new chromosome for a future
	// "same chrom" query (matches upstream forward-only state).
	m.chromIdx++
	m.slabIdx = 0
	if m.chromIdx >= len(m.chroms) {
		m.exhausted = true
	}
	return false
}
