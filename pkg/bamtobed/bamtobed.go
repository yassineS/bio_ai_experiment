// Package bamtobed provides streaming converters that read BAM, SAM, GFF
// or VCF text/binary and emit BED-style records suitable for the various
// `bedtools` ports in this repository.
//
// The converters are designed for the parity-testing path: most BED tools
// in this repo are BED-only at their library boundary, but the upstream
// `bedtools` programs accept BAM / SAM / GFF / VCF inputs via flags like
// `-ibam` / `-i x.vcf` / `-b x.gff`. Rather than teach every tool's core
// to dispatch on input format, we provide thin streaming adapters here
// and the CLI wrappers / parity tests wrap their input reader with one.
//
// All converters return a plain `io.Reader` emitting newline-terminated
// BED text. They run their decoding in a background goroutine and pipe
// the result through an `io.Pipe`, so consumers can read at their own
// pace without buffering the whole conversion in memory.
package bamtobed

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// BAMRef is a (chromosome, length) pair pulled from a BAM header. Used by
// bedgenomecov to seed its per-chromosome depth arrays when the user passes
// `-ibam` without a separate `-g` file (matching upstream's behaviour of
// using the BAM SQ headers as the genome).
type BAMRef struct {
	Name   string
	Length int
}

// ReadBAMHeader parses just the BAM header from r and returns the SQ
// reference list in declared order. The returned io.Reader yields the
// remainder of the BAM stream (header consumed), ready to be passed to a
// fresh sam.NewBAMReader via NewBAMBodyReader, or to one of the FromBAM*
// converters here.
//
// Note: this helper is for callers that want to extract the genome AND
// the per-record BED stream from the same BAM. It reads the full BAM into
// memory because sam.BAMReader doesn't currently expose the half-consumed
// underlying byte stream. For small BAM fixtures (which is the only place
// we use this) the simplicity is worth the extra copy.
func ReadBAMHeader(r io.Reader) (*sam.Header, []BAMRef, error) {
	br, err := sam.NewBAMReader(r)
	if err != nil {
		return nil, nil, err
	}
	hdr := br.Header()
	refs := make([]BAMRef, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refs[i] = BAMRef{Name: ref.Name, Length: int(ref.Length)}
	}
	return hdr, refs, nil
}

// FromBAM wraps a BGZF-wrapped BAM byte stream and emits BED6 text lines:
// `chrom\tstart\tend\tname\tmapq\tstrand`. Unmapped, secondary,
// supplementary, duplicate and QC-fail alignments are dropped (matching the
// upstream `BamRecordMgr` default filter). The reference span is computed
// from the CIGAR (`Pos-1`, `Pos-1+ReferenceLength()`); strand is derived
// from the BAM FLAG's reverse bit. Records with zero reference length are
// skipped.
func FromBAM(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		br, err := sam.NewBAMReader(r)
		if err != nil {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: open BAM: %w", err))
			return
		}
		for {
			rec, err := br.Read()
			if err == io.EOF {
				return
			}
			if err != nil {
				_ = pw.CloseWithError(fmt.Errorf("bamtobed: read BAM: %w", err))
				return
			}
			if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
				rec.IsDuplicate() || rec.IsQCFail() {
				continue
			}
			refLen := rec.Cigar.ReferenceLength()
			if refLen <= 0 {
				continue
			}
			start := int(rec.Pos) - 1
			if start < 0 {
				continue
			}
			end := start + refLen
			strand := "+"
			if rec.Flag&sam.FlagReverse != 0 {
				strand = "-"
			}
			name := rec.QName
			if name == "" {
				name = "."
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, start, end, name, rec.MapQ, strand); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
	}()
	return pr
}

// DecodeBAMToBED reads the entire BAM byte stream from r and returns
// (bedText, refs, nil). bedText is the BED6 textual rendering of every
// surviving primary alignment (same filter / formatting as FromBAM). refs
// is the SQ list parsed from the BAM header.
//
// This single-pass helper exists so callers that need BOTH the per-record
// BED stream AND the genome (e.g. bedgenomecov's `-ibam` mode) don't have
// to decode the BAM twice. Buffers the whole BAM in memory; intended for
// fixture-scale inputs.
func DecodeBAMToBED(r io.Reader) ([]byte, []BAMRef, error) {
	return decodeBAM(r, false)
}

// DecodeBAMSplitToBED is like DecodeBAMToBED but emits one BED record per
// CIGAR M-run (N breaks blocks, D extends; matches upstream `-split`).
func DecodeBAMSplitToBED(r io.Reader) ([]byte, []BAMRef, error) {
	return decodeBAM(r, true)
}

func decodeBAM(r io.Reader, split bool) ([]byte, []BAMRef, error) {
	br, err := sam.NewBAMReader(r)
	if err != nil {
		return nil, nil, fmt.Errorf("bamtobed: open BAM: %w", err)
	}
	hdr := br.Header()
	refs := make([]BAMRef, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refs[i] = BAMRef{Name: ref.Name, Length: int(ref.Length)}
	}
	var out strings.Builder
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("bamtobed: read BAM: %w", err)
		}
		if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
			rec.IsDuplicate() || rec.IsQCFail() {
			continue
		}
		cur := int(rec.Pos) - 1
		if cur < 0 {
			continue
		}
		strand := "+"
		if rec.Flag&sam.FlagReverse != 0 {
			strand = "-"
		}
		name := rec.QName
		if name == "" {
			name = "."
		}
		if !split {
			refLen := rec.Cigar.ReferenceLength()
			if refLen <= 0 {
				continue
			}
			fmt.Fprintf(&out, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, cur, cur+refLen, name, rec.MapQ, strand)
			continue
		}
		// split mode
		blockStart := cur
		blockLen := 0
		flush := func() {
			if blockLen <= 0 {
				return
			}
			fmt.Fprintf(&out, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, blockStart, blockStart+blockLen, name, rec.MapQ, strand)
		}
		for _, op := range rec.Cigar {
			oc := op.Op()
			l := int(op.Length())
			switch oc {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch, sam.CigarDeletion:
				if blockLen == 0 {
					blockStart = cur
				}
				blockLen += l
				cur += l
			case sam.CigarSkipped:
				flush()
				blockLen = 0
				cur += l
				blockStart = cur
			}
		}
		flush()
	}
	return []byte(out.String()), refs, nil
}

// FromBAMSplit is like FromBAM but expands each alignment's CIGAR `N` (and
// optionally `D`) gaps so each contiguous reference-consuming block is
// emitted as its own BED record. Mirrors upstream `-split`: M/=/X extend
// a block, N breaks a block, D extends (matches multicov / coverage
// `breakOnDeletionOps=false`), I/S/H/P consume no reference.
func FromBAMSplit(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		br, err := sam.NewBAMReader(r)
		if err != nil {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: open BAM: %w", err))
			return
		}
		for {
			rec, err := br.Read()
			if err == io.EOF {
				return
			}
			if err != nil {
				_ = pw.CloseWithError(fmt.Errorf("bamtobed: read BAM: %w", err))
				return
			}
			if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
				rec.IsDuplicate() || rec.IsQCFail() {
				continue
			}
			cur := int(rec.Pos) - 1
			if cur < 0 {
				continue
			}
			strand := "+"
			if rec.Flag&sam.FlagReverse != 0 {
				strand = "-"
			}
			name := rec.QName
			if name == "" {
				name = "."
			}
			blockStart := cur
			blockLen := 0
			flushBlock := func() error {
				if blockLen <= 0 {
					return nil
				}
				_, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%s\t%d\t%s\n",
					rec.RName, blockStart, blockStart+blockLen, name, rec.MapQ, strand)
				return err
			}
			for _, op := range rec.Cigar {
				oc := op.Op()
				l := int(op.Length())
				switch oc {
				case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch, sam.CigarDeletion:
					if blockLen == 0 {
						blockStart = cur
					}
					blockLen += l
					cur += l
				case sam.CigarSkipped:
					if err := flushBlock(); err != nil {
						_ = pw.CloseWithError(err)
						return
					}
					blockLen = 0
					cur += l
					blockStart = cur
				default:
					// I, S, H, P don't advance reference.
				}
			}
			if err := flushBlock(); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
	}()
	return pr
}

// DecodeOpts configures the per-CIGAR-op behaviour of DecodeBAM/SAM*Opts
// variants. SplitOnN breaks blocks on CIGAR `N` (default for `-split`).
// SplitOnDel additionally breaks blocks on CIGAR `D` — upstream genomecov
// uses this combination (see ContextGenomeCoverage.cpp), while multicov
// and coverage do not (they keep `breakOnDeletionOps=false`).
type DecodeOpts struct {
	SplitOnN   bool
	SplitOnDel bool
}

// DecodeBAMOptsToBED is the option-flag variant of DecodeBAMToBED.
func DecodeBAMOptsToBED(r io.Reader, opts DecodeOpts) ([]byte, []BAMRef, error) {
	return decodeBAMOpts(r, opts)
}

// DecodeSAMOptsToBED is the option-flag variant of DecodeSAMToBED.
func DecodeSAMOptsToBED(r io.Reader, opts DecodeOpts) ([]byte, []BAMRef, error) {
	return decodeSAMOpts(r, opts)
}

// DecodeSAMSplitOptsToBED is a convenience alias kept for naming symmetry
// with DecodeSAMSplitToBED.
func DecodeSAMSplitOptsToBED(r io.Reader, opts DecodeOpts) ([]byte, []BAMRef, error) {
	return decodeSAMOpts(r, opts)
}

// DecodeSAMToBED is the SAM-text counterpart of DecodeBAMToBED: it reads
// a SAM text stream and returns (bedText, refs, nil). Same filter and
// formatting as DecodeBAMToBED. Used by parity tests that only have a SAM
// fixture (no precomputed BAM) when we can't call out to `htsutil samtobam`.
func DecodeSAMToBED(r io.Reader) ([]byte, []BAMRef, error) {
	return decodeSAM(r, false)
}

// DecodeSAMSplitToBED is the SAM-text counterpart of DecodeBAMSplitToBED.
func DecodeSAMSplitToBED(r io.Reader) ([]byte, []BAMRef, error) {
	return decodeSAM(r, true)
}

func decodeBAMOpts(r io.Reader, opts DecodeOpts) ([]byte, []BAMRef, error) {
	br, err := sam.NewBAMReader(r)
	if err != nil {
		return nil, nil, fmt.Errorf("bamtobed: open BAM: %w", err)
	}
	hdr := br.Header()
	refs := make([]BAMRef, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refs[i] = BAMRef{Name: ref.Name, Length: int(ref.Length)}
	}
	var out strings.Builder
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("bamtobed: read BAM: %w", err)
		}
		writeRecOpts(&out, rec, opts)
	}
	return []byte(out.String()), refs, nil
}

func decodeSAMOpts(r io.Reader, opts DecodeOpts) ([]byte, []BAMRef, error) {
	sr, err := sam.NewSAMReader(r)
	if err != nil {
		return nil, nil, fmt.Errorf("bamtobed: open SAM: %w", err)
	}
	hdr := sr.Header()
	refs := make([]BAMRef, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refs[i] = BAMRef{Name: ref.Name, Length: int(ref.Length)}
	}
	var out strings.Builder
	for {
		rec, err := sr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("bamtobed: read SAM: %w", err)
		}
		writeRecOpts(&out, rec, opts)
	}
	return []byte(out.String()), refs, nil
}

// writeRecOpts emits BED text for one alignment using the configured op
// split-mode. Filter rules and field layout match decodeBAM/decodeSAM.
func writeRecOpts(out *strings.Builder, rec *sam.Record, opts DecodeOpts) {
	if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
		rec.IsDuplicate() || rec.IsQCFail() {
		return
	}
	cur := int(rec.Pos) - 1
	if cur < 0 {
		return
	}
	strand := "+"
	if rec.Flag&sam.FlagReverse != 0 {
		strand = "-"
	}
	name := rec.QName
	if name == "" {
		name = "."
	}
	if !opts.SplitOnN && !opts.SplitOnDel {
		refLen := rec.Cigar.ReferenceLength()
		if refLen <= 0 {
			return
		}
		fmt.Fprintf(out, "%s\t%d\t%d\t%s\t%d\t%s\n",
			rec.RName, cur, cur+refLen, name, rec.MapQ, strand)
		return
	}
	blockStart := cur
	blockLen := 0
	flush := func() {
		if blockLen <= 0 {
			return
		}
		fmt.Fprintf(out, "%s\t%d\t%d\t%s\t%d\t%s\n",
			rec.RName, blockStart, blockStart+blockLen, name, rec.MapQ, strand)
	}
	for _, op := range rec.Cigar {
		oc := op.Op()
		l := int(op.Length())
		switch oc {
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			if blockLen == 0 {
				blockStart = cur
			}
			blockLen += l
			cur += l
		case sam.CigarDeletion:
			if opts.SplitOnDel {
				flush()
				blockLen = 0
				cur += l
				blockStart = cur
			} else {
				if blockLen == 0 {
					blockStart = cur
				}
				blockLen += l
				cur += l
			}
		case sam.CigarSkipped:
			if opts.SplitOnN {
				flush()
				blockLen = 0
				cur += l
				blockStart = cur
			} else {
				if blockLen == 0 {
					blockStart = cur
				}
				blockLen += l
				cur += l
			}
		}
	}
	flush()
}

func decodeSAM(r io.Reader, split bool) ([]byte, []BAMRef, error) {
	sr, err := sam.NewSAMReader(r)
	if err != nil {
		return nil, nil, fmt.Errorf("bamtobed: open SAM: %w", err)
	}
	hdr := sr.Header()
	refs := make([]BAMRef, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refs[i] = BAMRef{Name: ref.Name, Length: int(ref.Length)}
	}
	var out strings.Builder
	for {
		rec, err := sr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("bamtobed: read SAM: %w", err)
		}
		if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
			rec.IsDuplicate() || rec.IsQCFail() {
			continue
		}
		cur := int(rec.Pos) - 1
		if cur < 0 {
			continue
		}
		strand := "+"
		if rec.Flag&sam.FlagReverse != 0 {
			strand = "-"
		}
		name := rec.QName
		if name == "" {
			name = "."
		}
		if !split {
			refLen := rec.Cigar.ReferenceLength()
			if refLen <= 0 {
				continue
			}
			fmt.Fprintf(&out, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, cur, cur+refLen, name, rec.MapQ, strand)
			continue
		}
		blockStart := cur
		blockLen := 0
		flush := func() {
			if blockLen <= 0 {
				return
			}
			fmt.Fprintf(&out, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, blockStart, blockStart+blockLen, name, rec.MapQ, strand)
		}
		for _, op := range rec.Cigar {
			oc := op.Op()
			l := int(op.Length())
			switch oc {
			case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch, sam.CigarDeletion:
				if blockLen == 0 {
					blockStart = cur
				}
				blockLen += l
				cur += l
			case sam.CigarSkipped:
				flush()
				blockLen = 0
				cur += l
				blockStart = cur
			}
		}
		flush()
	}
	return []byte(out.String()), refs, nil
}

// FromBAMPaired wraps a BGZF-wrapped BAM byte stream and emits BED3 lines
// covering the full fragment span [POS-1, MatePos-1+something) for
// properly-paired reads, mirroring `bedtools genomecov -ibam -pc`
// (paired-end coverage).
//
// Upstream logic: skip if not a proper pair (unpaired records are also
// skipped — `-pc` requires both mates). For each surviving record:
//
//   - IsFirstMate && IsReverseStrand    -> emit [MatePos-1, EndPosition)
//     (this read is to the right of its mate; cover from mate's left edge
//     to this read's right edge)
//   - IsFirstMate && IsMateReverseStrand -> emit [POS-1, POS-1+|ISIZE|)
//     (this read is to the left of its mate; cover from POS to ISIZE
//     downstream)
//
// All other records (second mates and oddly-oriented pairs) contribute
// nothing — the first mate's pair already covered the fragment, and
// mixed-orientation pairs are dropped wholesale (matches upstream's
// silent skip). Output is BED6 (`chrom\tstart\tend\tQNAME\tMAPQ\tstrand`)
// because consumers downstream still want a sortable name/mapq column.
func FromBAMPaired(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		br, err := sam.NewBAMReader(r)
		if err != nil {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: open BAM: %w", err))
			return
		}
		for {
			rec, err := br.Read()
			if err == io.EOF {
				return
			}
			if err != nil {
				_ = pw.CloseWithError(fmt.Errorf("bamtobed: read BAM: %w", err))
				return
			}
			if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
				rec.IsDuplicate() || rec.IsQCFail() {
				continue
			}
			if !rec.IsPaired() || !rec.IsProperPair() || rec.IsMateUnmapped() {
				continue
			}
			start, end, ok := pairedCoverageInterval(rec)
			if !ok {
				continue
			}
			strand := "+"
			if rec.Flag&sam.FlagReverse != 0 {
				strand = "-"
			}
			name := rec.QName
			if name == "" {
				name = "."
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, start, end, name, rec.MapQ, strand); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
	}()
	return pr
}

// pairedCoverageInterval returns the [start, end) reference span for one
// alignment under bedtools' `-pc` rules, plus an ok flag. Only the
// leftmost first-mate of each pair contributes; everyone else returns
// ok=false. Upstream code: see genomeCoverageBed.cpp pair_chip branch.
func pairedCoverageInterval(rec *sam.Record) (int, int, bool) {
	if !rec.IsRead1() {
		return 0, 0, false
	}
	// Skip pairs where the orientation is wrong (reverse mate to the
	// left of, or forward mate to the right of, its partner).
	pos := int(rec.Pos) - 1
	matePos := int(rec.PNext) - 1
	if (pos < matePos && rec.Flag&sam.FlagReverse != 0) ||
		(matePos < pos && rec.Flag&sam.FlagMateReverse != 0) {
		return 0, 0, false
	}
	switch {
	case rec.Flag&sam.FlagReverse != 0:
		// Right mate of a forward+reverse pair: cover from mate's left
		// edge to this read's right edge.
		refLen := rec.Cigar.ReferenceLength()
		if refLen <= 0 {
			return 0, 0, false
		}
		end := pos + refLen
		if matePos < 0 || matePos >= end {
			return 0, 0, false
		}
		return matePos, end, true
	case rec.Flag&sam.FlagMateReverse != 0:
		// Left mate: cover [POS-1, POS-1+|ISIZE|).
		isize := int(rec.TLen)
		if isize < 0 {
			isize = -isize
		}
		if isize <= 0 {
			return 0, 0, false
		}
		return pos, pos + isize, true
	}
	return 0, 0, false
}

// FromSAMPaired is the SAM-text counterpart of FromBAMPaired (see
// FromBAMPaired for semantics).
func FromSAMPaired(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		sr, err := sam.NewSAMReader(r)
		if err != nil {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: open SAM: %w", err))
			return
		}
		for {
			rec, err := sr.Read()
			if err == io.EOF {
				return
			}
			if err != nil {
				_ = pw.CloseWithError(fmt.Errorf("bamtobed: read SAM: %w", err))
				return
			}
			if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
				rec.IsDuplicate() || rec.IsQCFail() {
				continue
			}
			if !rec.IsPaired() || !rec.IsProperPair() || rec.IsMateUnmapped() {
				continue
			}
			start, end, ok := pairedCoverageInterval(rec)
			if !ok {
				continue
			}
			strand := "+"
			if rec.Flag&sam.FlagReverse != 0 {
				strand = "-"
			}
			name := rec.QName
			if name == "" {
				name = "."
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, start, end, name, rec.MapQ, strand); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
	}()
	return pr
}

// FromSAMExtended is the SAM-text counterpart of FromBAMExtended.
func FromSAMExtended(r io.Reader, fragSize int) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		if fragSize <= 0 {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: -fs fragment size must be > 0"))
			return
		}
		sr, err := sam.NewSAMReader(r)
		if err != nil {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: open SAM: %w", err))
			return
		}
		for {
			rec, err := sr.Read()
			if err == io.EOF {
				return
			}
			if err != nil {
				_ = pw.CloseWithError(fmt.Errorf("bamtobed: read SAM: %w", err))
				return
			}
			if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
				rec.IsDuplicate() || rec.IsQCFail() {
				continue
			}
			refLen := rec.Cigar.ReferenceLength()
			if refLen <= 0 {
				continue
			}
			pos := int(rec.Pos) - 1
			if pos < 0 {
				continue
			}
			strand := "+"
			var start, end int
			if rec.Flag&sam.FlagReverse != 0 {
				strand = "-"
				end = pos + refLen
				start = end - fragSize
				if start < 0 {
					start = 0
				}
			} else {
				start = pos
				end = pos + fragSize
			}
			if end <= start {
				continue
			}
			name := rec.QName
			if name == "" {
				name = "."
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, start, end, name, rec.MapQ, strand); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
	}()
	return pr
}

// ReadSAMHeaderRefs parses the SAM header from r and returns the @SQ
// reference list, plus the remainder of the SAM body buffered into a
// fresh reader. Used by parity tests that only have a SAM fixture and
// need both the @SQ-derived genome and the body stream.
func ReadSAMHeaderRefs(r io.Reader) ([]BAMRef, io.Reader, error) {
	// Drain the whole input once — fixtures are small.
	buf, err := io.ReadAll(r)
	if err != nil {
		return nil, nil, err
	}
	sr, err := sam.NewSAMReader(strings.NewReader(string(buf)))
	if err != nil {
		return nil, nil, fmt.Errorf("bamtobed: open SAM: %w", err)
	}
	hdr := sr.Header()
	refs := make([]BAMRef, len(hdr.Refs))
	for i, ref := range hdr.Refs {
		refs[i] = BAMRef{Name: ref.Name, Length: int(ref.Length)}
	}
	// Hand back the full SAM (header + body) so callers can use any of
	// the FromSAM* converters on it.
	return refs, strings.NewReader(string(buf)), nil
}

// FromBAMExtended wraps a BGZF-wrapped BAM byte stream and emits BED6
// lines with each alignment extended downstream-or-upstream-of-the-5'-end
// to a fixed fragment length, mirroring `bedtools genomecov -ibam -fs N`.
//
// Forward-strand records become [POS-1, POS-1+fragSize); reverse-strand
// records become [POS-1+ReferenceLength-fragSize, POS-1+ReferenceLength).
// Reverse-strand extensions that would underflow are clamped to start at
// 0 (matching upstream's "if(end<fragSize) AddCoverage(0,end)" branch).
//
// Filtering matches FromBAM: unmapped/secondary/supplementary/duplicate/
// QC-fail records are dropped.
func FromBAMExtended(r io.Reader, fragSize int) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		if fragSize <= 0 {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: -fs fragment size must be > 0"))
			return
		}
		br, err := sam.NewBAMReader(r)
		if err != nil {
			_ = pw.CloseWithError(fmt.Errorf("bamtobed: open BAM: %w", err))
			return
		}
		for {
			rec, err := br.Read()
			if err == io.EOF {
				return
			}
			if err != nil {
				_ = pw.CloseWithError(fmt.Errorf("bamtobed: read BAM: %w", err))
				return
			}
			if rec.IsUnmapped() || rec.IsSecondary() || rec.IsSupplementary() ||
				rec.IsDuplicate() || rec.IsQCFail() {
				continue
			}
			refLen := rec.Cigar.ReferenceLength()
			if refLen <= 0 {
				continue
			}
			pos := int(rec.Pos) - 1
			if pos < 0 {
				continue
			}
			strand := "+"
			var start, end int
			if rec.Flag&sam.FlagReverse != 0 {
				strand = "-"
				end = pos + refLen
				start = end - fragSize
				if start < 0 {
					start = 0
				}
			} else {
				start = pos
				end = pos + fragSize
			}
			if end <= start {
				continue
			}
			name := rec.QName
			if name == "" {
				name = "."
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%s\t%d\t%s\n",
				rec.RName, start, end, name, rec.MapQ, strand); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
	}()
	return pr
}

// FromGFF wraps a GFF text reader and emits BED3 lines (chrom,
// col4-1, col5). Header / comment / blank lines and rows with fewer than
// 5 columns or unparseable coordinates are silently dropped.
func FromGFF(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		sc := bufio.NewScanner(r)
		sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
		for sc.Scan() {
			line := sc.Text()
			trim := strings.TrimSpace(line)
			if trim == "" || strings.HasPrefix(trim, "#") ||
				strings.HasPrefix(trim, "browser") || strings.HasPrefix(trim, "track") {
				continue
			}
			fields := strings.Split(line, "\t")
			if len(fields) < 5 {
				continue
			}
			s, err := strconv.Atoi(fields[3])
			if err != nil || s < 1 {
				continue
			}
			e, err := strconv.Atoi(fields[4])
			if err != nil || e < s {
				continue
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\n", fields[0], s-1, e); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
		if err := sc.Err(); err != nil {
			_ = pw.CloseWithError(err)
		}
	}()
	return pr
}

// FromVCF wraps a VCF text reader and emits BED3 lines (CHROM,
// POS-1, POS-1+len(REF)). All `#`-prefixed header lines and rows with
// fewer than 4 columns are silently dropped.
func FromVCF(r io.Reader) io.Reader {
	pr, pw := io.Pipe()
	go func() {
		bw := bufio.NewWriter(pw)
		defer func() {
			_ = bw.Flush()
			_ = pw.Close()
		}()
		sc := bufio.NewScanner(r)
		sc.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
		for sc.Scan() {
			line := sc.Text()
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			fields := strings.Split(line, "\t")
			if len(fields) < 4 {
				continue
			}
			pos, err := strconv.Atoi(fields[1])
			if err != nil || pos < 1 {
				continue
			}
			ref := fields[3]
			refLen := len(ref)
			if refLen == 0 {
				refLen = 1
			}
			if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\n", fields[0], pos-1, pos-1+refLen); err != nil {
				_ = pw.CloseWithError(err)
				return
			}
		}
		if err := sc.Err(); err != nil {
			_ = pw.CloseWithError(err)
		}
	}()
	return pr
}
