package samtools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// samFastPathEligible reports whether the BAM->SAM text output path can use the
// allocation-light direct serialiser (sam.BAMReader.WriteSAMBody) instead of
// decoding each record into a full sam.Record and formatting it through
// SAMWriter.
//
// The fast path serialises straight from the raw BAM bytes and decodes only the
// fixed-prefix fields (flag, MAPQ, RName, Pos, CIGAR reference span). It is
// therefore eligible only when:
//   - output is plain text SAM (not BAM/CRAM, not -c/Count), and
//   - every active per-record filter is decidable from those fixed-prefix
//     fields — i.e. flag include/exclude, MAPQ, region overlap and BED overlap.
//
// Filters that need the read name or an aux tag (read-group, -d/-D tag filters,
// -N qname sets, -s subsampling) fall back to the full-decode path, which still
// has the decoded Record those predicates require. The output bytes are
// identical either way; only the per-record work differs.
func samFastPathEligible(opts *ViewOptions) bool {
	if opts.OutputBAM || opts.OutputCRAM {
		return false
	}
	if opts.ReadGroup != "" || len(opts.ReadGroupSet) > 0 {
		return false
	}
	if len(opts.TagFilters) > 0 || len(opts.QNameSet) > 0 {
		return false
	}
	if opts.Subsample > 0 && opts.Subsample < 1 {
		return false
	}
	return true
}

// rawAuxBAMSinkEligible reports whether the CRAM→BAM decode may use the
// memory-lean raw-aux passthrough: the CRAM decoder builds each record's aux as
// a raw on-disk BAM aux byte block (sam.Record.RawAux) that the BAM writer emits
// verbatim, so the heavier []sam.Aux is never materialised.
//
// It is eligible only when the sink is a BAM writer AND nothing on the path
// reads a record's aux fields. The aux-reading touch points are the read-group
// filter (-r/-R, which calls GetAux("RG")) and the -d/-D tag filters (which read
// the tag value); every other view filter (flag/MAPQ/region/BED/-N qname/-s
// subsample) is decided from non-aux fields. A SAM text or CRAM sink reads the
// decoded aux to re-serialise it, so those are excluded. Count mode writes no
// records and gains nothing, so it stays on the eager path too. When unsure the
// predicate errs OFF, leaving today's eager behaviour and its byte output
// unchanged. The mode is a no-op for non-CRAM input (the BAM/SAM readers never
// set RawAux), so this only ever engages on the CRAM→BAM path.
func rawAuxBAMSinkEligible(opts *ViewOptions) bool {
	if opts.OutputCRAM {
		return false
	}
	if !opts.OutputBAM && !opts.Uncompressed {
		return false // SAM text or any non-BAM sink reads the decoded aux.
	}
	if opts.Count {
		return false // counting writes no records; no passthrough benefit.
	}
	if opts.ReadGroup != "" || len(opts.ReadGroupSet) > 0 {
		return false // -r/-R reads RG via GetAux.
	}
	if len(opts.TagFilters) > 0 {
		return false // -d/-D reads aux tag values.
	}
	return true
}

// fastRegionPredicate builds a region-overlap test over the cheaply decoded
// FastFields, keyed on the raw refID so no header name lookup is needed. It
// returns nil when no regions are configured (meaning "keep all"); otherwise a
// predicate that keeps a record when its [pos0, pos0+refSpan) range overlaps any
// region on the same reference. A record on a reference no region targets is
// dropped, matching buildRegionFilter.
func fastRegionPredicate(regions []region.ResolvedRegion) func(*sam.FastFields) bool {
	if len(regions) == 0 {
		return nil
	}
	return func(ff *sam.FastFields) bool {
		if ff.RefID < 0 {
			return false
		}
		rid := int(ff.RefID)
		pos0 := int(ff.Pos) - 1
		if pos0 < 0 {
			pos0 = 0
		}
		refLen := ff.RefSpan
		if refLen <= 0 {
			refLen = 1
		}
		recEnd := pos0 + refLen
		for i := range regions {
			r := &regions[i]
			if r.RefID != rid {
				continue
			}
			if pos0 < r.End0 && recEnd > r.Beg0 {
				return true
			}
		}
		return false
	}
}

// fastFlagFilter applies the flag/MAPQ predicates of keepRecord using only the
// FastFields (which carry Flag and MapQ). It is a pure subset of keepRecord:
// the read-group/aux/qname/subsample branches are excluded by
// samFastPathEligible, so this covers every filter that can still be active on
// the fast path.
func fastFlagFilter(ff *sam.FastFields, opts *ViewOptions) bool {
	if opts.IncludeFlags != 0 && ff.Flag&opts.IncludeFlags != opts.IncludeFlags {
		return false
	}
	if opts.ExcludeFlags != 0 && ff.Flag&opts.ExcludeFlags != 0 {
		return false
	}
	if opts.UseExcludeAll && ff.Flag&opts.ExcludeFlagsAll == opts.ExcludeFlagsAll {
		return false
	}
	if opts.MinMAPQ > 0 && ff.MapQ < opts.MinMAPQ {
		return false
	}
	return true
}

// fastBedFilter builds a BED-overlap predicate over FastFields, mirroring
// loadBedFilter's per-record test (a record is kept when its
// [pos0, pos0+refSpan) range intersects any BED interval on its RName). It
// returns (nil, nil) when path is empty. The per-chromosome interval trees and
// the half-open / zero-span semantics match loadBedFilter exactly so the kept
// set is identical to the full-decode path.
func fastBedFilter(path string) (func(*sam.FastFields) bool, error) {
	if path == "" {
		return nil, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("samtools view: open BED %q: %w", path, err)
	}
	defer f.Close()
	byChrom := map[string][]*bed.Record{}
	rd := bed.NewReader(f)
	for {
		rec, rerr := rd.Read()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return nil, fmt.Errorf("samtools view: read BED %q: %w", path, rerr)
		}
		byChrom[rec.Chrom] = append(byChrom[rec.Chrom], rec)
	}
	if len(byChrom) == 0 {
		return func(*sam.FastFields) bool { return false }, nil
	}
	trees := make(map[string]*bed.IntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		sort.Slice(recs, func(i, j int) bool {
			return recs[i].ChromStart < recs[j].ChromStart
		})
		trees[chrom] = bed.NewIntervalTree(recs)
	}
	return func(ff *sam.FastFields) bool {
		if ff.Flag&sam.FlagUnmapped != 0 || ff.RName == "" || ff.RName == "*" {
			return false
		}
		t, ok := trees[ff.RName]
		if !ok {
			return false
		}
		pos0 := int(ff.Pos) - 1
		if pos0 < 0 {
			pos0 = 0
		}
		refLen := ff.RefSpan
		if refLen <= 0 {
			return false
		}
		q := &bed.Record{Chrom: ff.RName, ChromStart: pos0, ChromEnd: pos0 + refLen}
		return len(t.Query(q)) > 0
	}, nil
}

// viewStreamFast is the streaming (linear-scan) counterpart of
// viewIndexedChunksFast: it serialises a whole-stream BAM->SAM conversion
// straight from the raw BAM bytes, applying only the fixed-prefix-decidable
// filters. It is used by View when the underlying reader is a *sam.BAMReader
// and samFastPathEligible(opts) holds. Output is byte-identical to the decode
// loop. When opts.Count is set it counts survivors and writes the count rather
// than the records.
func viewStreamFast(br *sam.BAMReader, out io.Writer, opts *ViewOptions, hdr *sam.Header, resolved []region.ResolvedRegion) (int, error) {
	regionPred := fastRegionPredicate(resolved)
	if regionPred == nil && len(opts.Regions) > 0 {
		regionPred = func(*sam.FastFields) bool { return false }
	}
	bedPred, berr := fastBedFilter(opts.BedPath)
	if berr != nil {
		return 0, berr
	}

	var bw *bufio.Writer
	if !opts.Count {
		bw = bufio.NewWriter(out)
		// SAM is text: emit the header only when -h/-H requested (matching
		// SAMWriter via openViewWriter's emitHeader rule for non-binary output).
		if opts.HeaderOnly || opts.WithHeader {
			if _, err := hdr.WriteTo(bw); err != nil {
				return 0, err
			}
		}
		if opts.HeaderOnly {
			return 0, bw.Flush()
		}
	}

	matched, err := fastSAMScan(br, bw, opts, regionPred, bedPred)
	if err != nil {
		return matched, err
	}
	if opts.Count {
		fmt.Fprintln(out, matched)
		return matched, nil
	}
	if err := bw.Flush(); err != nil {
		return matched, err
	}
	return matched, nil
}

// fastSAMScan reads BAM records from br, applies the fast-path filters and
// serialises every survivor straight to bw via br.WriteSAMBody. It is the
// allocation-light replacement for the decode-into-Record / SAMWriter.Write
// inner loop, used when samFastPathEligible reports the configuration permits
// it. When bw is nil (the -c/Count case) survivors are counted but not
// serialised, so the cheap ReadSAMInto decode covers counting too. It returns
// the number of records that passed the filters; bw, when non-nil, is flushed
// by the caller.
func fastSAMScan(br *sam.BAMReader, bw *bufio.Writer, opts *ViewOptions,
	regionPred func(*sam.FastFields) bool, bedPred func(*sam.FastFields) bool) (int, error) {
	matched := 0
	var ff sam.FastFields
	for {
		err := br.ReadSAMInto(&ff)
		if err == io.EOF {
			break
		}
		if err != nil {
			return matched, err
		}
		if !fastFlagFilter(&ff, opts) {
			continue
		}
		if regionPred != nil && !regionPred(&ff) {
			continue
		}
		if bedPred != nil && !bedPred(&ff) {
			continue
		}
		matched++
		if bw == nil {
			continue
		}
		if err := br.WriteSAMBody(bw, &ff); err != nil {
			return matched, err
		}
	}
	return matched, nil
}
