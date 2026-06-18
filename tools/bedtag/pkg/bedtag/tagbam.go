// This file implements the upstream `bedtools tag` (aka tagBam) model: annotate
// the alignments of a BAM with an aux tag built from overlaps against one or
// more BED/GFF/VCF annotation files.
//
// Upstream reference: reference_code/bedtools/src/tagBam/tagBam.cpp. For each
// MAPPED alignment it forms the interval [POS, end) and, per annotation file,
// either records the file's label when there is ANY overlap (the default), or
// joins the overlapping records' name/score fields. The per-file fragments are
// joined with ';' and written as a single aux tag (default "YB", type Z). An
// alignment with no overlap (and every unmapped alignment) is emitted
// unchanged.

package bedtag

import (
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TagBAMMode selects what populates the tag for each overlapping annotation
// file.
type TagBAMMode int

const (
	// TagModeLabels writes the file's label when the alignment overlaps any
	// record in that file (upstream default).
	TagModeLabels TagBAMMode = iota
	// TagModeNames joins the overlapping records' name fields.
	TagModeNames
	// TagModeScores joins the overlapping records' score fields.
	TagModeScores
)

// TagBAMOptions configures TagBAM.
type TagBAMOptions struct {
	// Tag is the two-character aux tag to write (default "YB").
	Tag string
	// Mode selects labels/names/scores.
	Mode TagBAMMode
	// Labels are the per-file labels used in TagModeLabels (one per file).
	Labels []string
	// SameStrand (-s) only counts overlaps on the same strand; OppositeStrand
	// (-S) only the opposite. They are mutually exclusive.
	SameStrand     bool
	OppositeStrand bool
	// MinFraction (-f) is the minimum overlap as a fraction of the alignment
	// length. Default 1e-9 (effectively 1bp).
	MinFraction float64
}

// TagBAM reads BAM alignments from in, tags each mapped alignment per opts
// using the per-file annotation records in annoFiles (one slice of records per
// annotation file, in the same order as opts.Labels), and writes BAM to out.
// It returns the number of alignments written.
func TagBAM(in io.Reader, annoFiles [][]*bed.Record, out io.Writer, opts TagBAMOptions) (int, error) {
	if opts.Tag == "" {
		opts.Tag = "YB"
	}
	if len(opts.Tag) != 2 {
		return 0, fmt.Errorf("bedtag: tag must be exactly two characters, got %q", opts.Tag)
	}
	if opts.SameStrand && opts.OppositeStrand {
		return 0, fmt.Errorf("bedtag: -s and -S are mutually exclusive")
	}
	if opts.MinFraction <= 0 {
		opts.MinFraction = 1e-9
	}

	// Build a per-file, per-chrom interval index (with each record's original
	// file order for the UCSC-bin tie-break that allHits uses).
	annos := make([]annoIndex, len(annoFiles))
	for fi, recs := range annoFiles {
		byChrom := map[string][]*bed.Record{}
		idx := make(map[*bed.Record]int, len(recs))
		for i, r := range recs {
			byChrom[r.Chrom] = append(byChrom[r.Chrom], r)
			idx[r] = i
		}
		trees := make(map[string]*bed.IntervalTree, len(byChrom))
		for chrom, rs := range byChrom {
			sort.SliceStable(rs, func(i, j int) bool { return rs[i].ChromStart < rs[j].ChromStart })
			trees[chrom] = bed.NewIntervalTree(rs)
		}
		annos[fi] = annoIndex{trees: trees, fileIdx: idx}
	}

	br, err := sam.NewBAMReader(in)
	if err != nil {
		return 0, fmt.Errorf("bedtag: open BAM: %w", err)
	}
	bw := sam.NewBAMWriter(out)
	if err := bw.WriteHeader(br.Header()); err != nil {
		return 0, err
	}

	count := 0
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return count, err
		}
		if rec.IsMapped() {
			if tagVal := buildTag(rec, annos, opts); tagVal != "" {
				setAux(rec, opts.Tag, tagVal)
			}
		}
		if err := bw.Write(rec); err != nil {
			return count, err
		}
		count++
	}
	if err := bw.Close(); err != nil {
		return count, err
	}
	return count, nil
}

// annoIndex is a per-annotation-file overlap index: per-chrom interval trees
// plus each record's original file order for the allHits bin-order tie-break.
type annoIndex struct {
	trees   map[string]*bed.IntervalTree
	fileIdx map[*bed.Record]int
}

// buildTag returns the aux value for an alignment, or "" when nothing overlaps.
func buildTag(rec *sam.Record, annos []annoIndex, opts TagBAMOptions) string {
	aStart := int(rec.Pos) - 1 // 0-based start
	aEnd := int(rec.EndPosition())
	if aEnd <= aStart {
		return ""
	}
	aStrand := "+"
	if rec.Flag&sam.FlagReverse != 0 {
		aStrand = "-"
	}
	var fragments []string
	for fi := range annos {
		hits := overlapHits(rec.RName, aStart, aEnd, aStrand, annos[fi].trees, annos[fi].fileIdx, &opts)
		if len(hits) == 0 {
			continue
		}
		switch opts.Mode {
		case TagModeLabels:
			label := ""
			if fi < len(opts.Labels) {
				label = opts.Labels[fi]
			}
			fragments = append(fragments, label)
		case TagModeNames:
			parts := make([]string, len(hits))
			for i, h := range hits {
				parts[i] = h.Name
			}
			fragments = append(fragments, strings.Join(parts, ","))
		case TagModeScores:
			parts := make([]string, len(hits))
			for i, h := range hits {
				parts[i] = strconv.Itoa(h.Score)
			}
			fragments = append(fragments, strings.Join(parts, ","))
		}
	}
	return strings.Join(fragments, ";")
}

// overlapHits returns the annotation records overlapping [aStart,aEnd) by at
// least opts.MinFraction of the alignment length, in upstream allHits order
// (UCSC bin levels finest-first, then bin number, then file order). For the
// labels mode the caller only checks len()>0, so the order is harmless there.
func overlapHits(chrom string, aStart, aEnd int, aStrand string, trees map[string]*bed.IntervalTree, fileIdx map[*bed.Record]int, opts *TagBAMOptions) []*bed.Record {
	tree, ok := trees[chrom]
	if !ok {
		return nil
	}
	cands := tree.Query(&bed.Record{Chrom: chrom, ChromStart: aStart, ChromEnd: aEnd})
	if len(cands) == 0 {
		return nil
	}
	alignLen := aEnd - aStart
	out := make([]*bed.Record, 0, len(cands))
	for _, c := range cands {
		if !tagStrandOK(aStrand, c.Strand, opts) {
			continue
		}
		ov := overlapLen(aStart, aEnd, c.ChromStart, c.ChromEnd)
		if ov <= 0 {
			continue
		}
		if float64(ov)/float64(alignLen) < opts.MinFraction {
			continue
		}
		out = append(out, c)
	}
	sort.SliceStable(out, func(i, j int) bool {
		li, bi := binLevelIdx(out[i].ChromStart, out[i].ChromEnd)
		lj, bj := binLevelIdx(out[j].ChromStart, out[j].ChromEnd)
		switch {
		case li != lj:
			return li < lj
		case bi != bj:
			return bi < bj
		default:
			return fileIdx[out[i]] < fileIdx[out[j]]
		}
	})
	return out
}

// binLevelIdx returns the UCSC bin LEVEL (0 = finest) and bin index within that
// level for an interval, mirroring bedtools' getBin (bedFile.h).
func binLevelIdx(start, end int) (level, idx int) {
	s := start >> binFirstShift
	e := (end - 1) >> binFirstShift
	for i := 0; i < binLevels; i++ {
		if s == e {
			return i, s
		}
		s >>= binNextShift
		e >>= binNextShift
	}
	return binLevels - 1, s
}

const (
	binFirstShift = 14
	binNextShift  = 3
	binLevels     = 8
)

func tagStrandOK(aStrand, bStrand string, opts *TagBAMOptions) bool {
	if opts.SameStrand {
		return aStrand == bStrand
	}
	if opts.OppositeStrand {
		return bStrand != "" && bStrand != "." && aStrand != bStrand
	}
	return true
}

func overlapLen(aS, aE, bS, bE int) int {
	if aS < bS {
		aS = bS
	}
	if aE > bE {
		aE = bE
	}
	if aE <= aS {
		return 0
	}
	return aE - aS
}

// setAux replaces or appends a Z-type aux tag on a record.
func setAux(rec *sam.Record, tag, value string) {
	for i := range rec.Aux {
		if rec.Aux[i].Tag == tag {
			rec.Aux[i].Type = 'Z'
			rec.Aux[i].Value = value
			return
		}
	}
	rec.Aux = append(rec.Aux, sam.Aux{Tag: tag, Type: 'Z', Value: value})
}
