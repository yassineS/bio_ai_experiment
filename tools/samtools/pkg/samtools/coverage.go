package samtools

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// CoverageOptions configures Coverage. Defaults match `samtools coverage`'s
// tabular mode: one row per reference (or per region when Regions is set).
type CoverageOptions struct {
	// Regions restricts output to the given chr[:start-end] specs. When
	// empty, every @SQ entry is reported.
	Regions []string
	// MinMAPQ skips records with MAPQ < MinMAPQ.
	MinMAPQ uint8
	// MinBaseQ skips bases with quality < MinBaseQ from depth/baseq stats.
	MinBaseQ uint8
	// IncludeFlags requires every flag bit in the mask.
	IncludeFlags uint16
	// ExcludeFlags drops records with any of these flag bits. Defaults to
	// FlagUnmapped|FlagSecondary|FlagSupplementary|FlagDuplicate when zero.
	ExcludeFlags uint16
	// Histogram requests the ASCII-histogram output mode (`-A`). Falls
	// back to a stub message in v1 — the tabular form is the heavily
	// used variant.
	Histogram bool
	// HistogramBins is the column width used when Histogram is true.
	HistogramBins int
	// HeaderOff (-H) suppresses the column-header line.
	HeaderOff bool
}

// CoverageRow is one line of `samtools coverage`'s tabular output.
type CoverageRow struct {
	RName     string
	StartPos  int32
	EndPos    int32
	NumReads  uint64
	CovBases  uint64
	Coverage  float64
	MeanDepth float64
	MeanBaseQ float64
	MeanMapQ  float64
}

// coverageRefState is the per-reference accumulator used during the
// streaming pass.
type coverageRefState struct {
	length   int32
	numReads uint64
	// posDelta maps 0-based reference position -> delta of in-flight depth.
	posDelta map[int]int
	baseQSum uint64
	baseQCnt uint64
	mapQSum  uint64
	mapQCnt  uint64
}

// Coverage walks records from in (a single BAM/SAM stream — multi-file
// support layered by the CLI) and emits one row per @SQ entry (or per
// region) summarising depth statistics.
func Coverage(in io.Reader, w io.Writer, opts CoverageOptions) error {
	r, err := sam.NewReader(in)
	if err != nil {
		return err
	}
	hdr := r.Header()

	excl := opts.ExcludeFlags
	if excl == 0 {
		excl = sam.FlagUnmapped | sam.FlagSecondary | sam.FlagSupplementary | sam.FlagDuplicate
	}

	states := make([]*coverageRefState, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		states[i] = &coverageRefState{length: ref.Length, posDelta: map[int]int{}}
	}

	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if rec.Flag&excl != 0 {
			continue
		}
		if opts.IncludeFlags != 0 && rec.Flag&opts.IncludeFlags != opts.IncludeFlags {
			continue
		}
		if rec.MapQ < opts.MinMAPQ {
			continue
		}
		idx := -1
		if rec.RName != "" && rec.RName != "*" {
			idx = hdr.RefIndex(rec.RName)
		}
		if idx < 0 {
			continue
		}
		st := states[idx]
		st.numReads++
		st.mapQSum += uint64(rec.MapQ)
		st.mapQCnt++

		// Walk CIGAR, accumulating per-base depth and baseQ stats. Bases
		// below MinBaseQ contribute neither to depth nor baseq sums
		// (matches upstream).
		refPos := int(rec.Pos) - 1
		queryPos := 0
		for _, op := range rec.Cigar {
			length := int(op.Length())
			switch op.Op() {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
				for i := 0; i < length; i++ {
					var q uint8 = 0xff
					if queryPos+i < len(rec.Qual) {
						q = rec.Qual[queryPos+i]
					}
					if opts.MinBaseQ > 0 && q < opts.MinBaseQ {
						continue
					}
					if q != 0xff {
						st.baseQSum += uint64(q)
						st.baseQCnt++
					}
					pos := refPos + i
					st.posDelta[pos]++
					st.posDelta[pos+1]--
				}
				refPos += length
				queryPos += length
			case sam.CigarInsertion, sam.CigarSoftClip:
				queryPos += length
			case sam.CigarDeletion, sam.CigarSkipped:
				refPos += length
			case sam.CigarHardClip, sam.CigarPadding:
				// no movement on query or ref
			}
		}
	}

	// Resolve regions list — when empty, every @SQ counts as the full
	// reference.
	type spec struct {
		name     string
		idx      int
		startPos int32
		endPos   int32
	}
	var specs []spec
	if len(opts.Regions) == 0 {
		for i, ref := range hdr.Refs {
			specs = append(specs, spec{name: ref.Name, idx: i, startPos: 1, endPos: ref.Length})
		}
	} else {
		for _, r := range opts.Regions {
			rg, perr := region.ParseRegion(r)
			if perr != nil {
				return perr
			}
			i := hdr.RefIndex(rg.Chrom)
			if i < 0 {
				return fmt.Errorf("samtools coverage: unknown ref %q", rg.Chrom)
			}
			startPos := int32(1)
			endPos := hdr.Refs[i].Length
			if rg.Beg > 0 {
				startPos = int32(rg.Beg)
			}
			if rg.End > 0 {
				endPos = int32(rg.End)
			}
			specs = append(specs, spec{name: rg.Chrom, idx: i, startPos: startPos, endPos: endPos})
		}
	}

	if !opts.HeaderOff {
		fmt.Fprintln(w, "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tbaseq\tmapq")
	}
	for _, s := range specs {
		row := summariseCoverage(s.name, s.startPos, s.endPos, states[s.idx])
		fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.4f\t%.4f\n",
			row.RName, row.StartPos, row.EndPos, row.NumReads,
			row.CovBases, row.Coverage, row.MeanDepth,
			row.MeanBaseQ, row.MeanMapQ)
	}
	return nil
}

// summariseCoverage computes the per-region row from a refState by
// walking the delta map in sorted order.
func summariseCoverage(name string, start, end int32, st *coverageRefState) CoverageRow {
	row := CoverageRow{
		RName:    name,
		StartPos: start,
		EndPos:   end,
		NumReads: st.numReads,
	}
	if st.baseQCnt > 0 {
		row.MeanBaseQ = float64(st.baseQSum) / float64(st.baseQCnt)
	}
	if st.mapQCnt > 0 {
		row.MeanMapQ = float64(st.mapQSum) / float64(st.mapQCnt)
	}
	if len(st.posDelta) == 0 {
		return row
	}
	keys := make([]int, 0, len(st.posDelta))
	for k := range st.posDelta {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	depth := 0
	prev := keys[0]
	var covBases, depthSum uint64
	for _, k := range keys {
		if depth > 0 && k > prev {
			lo := prev
			hi := k
			if int(start)-1 > lo {
				lo = int(start) - 1
			}
			if int(end) < hi {
				hi = int(end)
			}
			if hi > lo {
				covBases += uint64(hi - lo)
				depthSum += uint64(depth) * uint64(hi-lo)
			}
		}
		depth += st.posDelta[k]
		prev = k
	}
	row.CovBases = covBases
	regionLen := uint64(end - start + 1)
	if regionLen > 0 {
		row.Coverage = 100.0 * float64(covBases) / float64(regionLen)
		row.MeanDepth = float64(depthSum) / float64(regionLen)
	}
	return row
}
