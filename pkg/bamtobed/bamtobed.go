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
