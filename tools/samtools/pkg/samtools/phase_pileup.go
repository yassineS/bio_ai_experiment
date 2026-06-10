package samtools

import (
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// phasePileupColumn is a single pileup column produced by the
// streaming pileup driver below. It mirrors the subset of
// `bam_pileup1_t` upstream actually inspects in phase.c's main loop:
// the read, qpos (0-based query position), and the is_del / is_refskip
// flags (always false in this minimal port — phase.c skips them, so
// they would never reach the emit loop anyway).
type phasePileupColumn struct {
	pos   int32 // 1-based reference position
	tid   int32 // reference index (0 for our single-contig fixtures)
	plp   []phasePlpRead
	rname string
}

// phasePlpRead is one read's contribution at a pileup column.
// Matches the bam_pileup1_t fields phase.c actually reads:
//
//	b      → the read record
//	qpos   → 0-based index into b.Seq for the aligned base
//	isDel  → CIGAR D at this column
//	isRefSkip → CIGAR N at this column
type phasePlpRead struct {
	b         *sam.Record
	qpos      int32
	isDel     bool
	isRefSkip bool
}

// phaseStreamPileup is a per-reference streaming pileup driver
// mirroring htslib's bam_plp_auto semantics for the subset of CIGAR ops
// phase.c examines. It walks recs (assumed pre-sorted by Pos), emits
// one phasePileupColumn per 1-based reference position covered by any
// read, in coordinate order; within each column, reads appear in the
// order they entered the active buffer (i.e. sorted by Pos, then by
// the order in which equal-Pos reads appeared in the input — which is
// the BAM record order for coordinate-sorted input). Reads contribute
// to a column iff the CIGAR op at that column is M/=/X (a single base
// callable from the query); D/N columns are produced with is_del /
// is_refskip set so the downstream main loop can skip them as upstream
// does.
//
// The implementation pre-computes per-read CIGAR walks (since the
// records are in memory and the fixture is small) and emits columns
// from a min-heap of (pos, reading-pointer) — equivalent to upstream's
// active-read buffer.
type phaseStreamPileup struct {
	recs  []*sam.Record
	rname string
	// per-read state:
	cigarWalk []cigarColumn // sorted by (pos, readIdx)
	cur       int
}

// cigarColumn is one (refPos, readIdx, qpos, kind) tuple produced by
// expanding a single record's CIGAR. They are stably sorted by refPos
// and then by readIdx so the column-emit order matches upstream.
type cigarColumn struct {
	pos       int32
	readIdx   int32
	qpos      int32
	isDel     bool
	isRefSkip bool
}

// newPhaseStreamPileup expands every record into per-position
// cigarColumn entries and sorts them. The set is small for phase
// fixtures (read count × read length) so this is cheap.
func newPhaseStreamPileup(recs []*sam.Record, rname string) *phaseStreamPileup {
	var cols []cigarColumn
	for ri, rec := range recs {
		if rec.Pos <= 0 || len(rec.Cigar) == 0 {
			continue
		}
		refPos := rec.Pos // 1-based
		qpos := int32(0)
		for _, op := range rec.Cigar {
			o := op.Op()
			n := int(op.Length())
			switch o {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
				for k := 0; k < n; k++ {
					cols = append(cols, cigarColumn{
						pos:     refPos + int32(k),
						readIdx: int32(ri),
						qpos:    qpos + int32(k),
					})
				}
				refPos += int32(n)
				qpos += int32(n)
			case sam.CigarInsertion, sam.CigarSoftClip:
				qpos += int32(n)
			case sam.CigarDeletion:
				for k := 0; k < n; k++ {
					cols = append(cols, cigarColumn{
						pos:     refPos + int32(k),
						readIdx: int32(ri),
						qpos:    qpos,
						isDel:   true,
					})
				}
				refPos += int32(n)
			case sam.CigarSkipped:
				for k := 0; k < n; k++ {
					cols = append(cols, cigarColumn{
						pos:       refPos + int32(k),
						readIdx:   int32(ri),
						qpos:      qpos,
						isRefSkip: true,
					})
				}
				refPos += int32(n)
			case sam.CigarHardClip, sam.CigarPadding:
				// no axis movement
			}
		}
	}
	sort.SliceStable(cols, func(i, j int) bool {
		if cols[i].pos != cols[j].pos {
			return cols[i].pos < cols[j].pos
		}
		return cols[i].readIdx < cols[j].readIdx
	})
	return &phaseStreamPileup{recs: recs, rname: rname, cigarWalk: cols}
}

// next returns the next phasePileupColumn or nil when the stream is
// exhausted. Columns are emitted in coordinate order.
func (p *phaseStreamPileup) next() *phasePileupColumn {
	if p.cur >= len(p.cigarWalk) {
		return nil
	}
	pos := p.cigarWalk[p.cur].pos
	col := &phasePileupColumn{pos: pos, rname: p.rname}
	for p.cur < len(p.cigarWalk) && p.cigarWalk[p.cur].pos == pos {
		cc := p.cigarWalk[p.cur]
		col.plp = append(col.plp, phasePlpRead{
			b:         p.recs[cc.readIdx],
			qpos:      cc.qpos,
			isDel:     cc.isDel,
			isRefSkip: cc.isRefSkip,
		})
		p.cur++
	}
	return col
}
