// `bcftools index` builds a CSI (or .tbi for VCF.gz) index for a sorted
// BCF/VCF file. CSI is the default for BCF and for VCFs with chromosomes
// longer than the BAI scheme can address. The .tbi path delegates to the
// shared tabix builder; this file owns the BCF→CSI traversal.
package bcftools

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bcf"
	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
	"github.com/yassineS/bio_ai_experiment/tools/tabix/pkg/tabix"
)

// IndexFormat selects the on-disk index flavour.
type IndexFormat int

const (
	// IndexCSI emits a `.csi` index (the bcftools default for BCF).
	IndexCSI IndexFormat = iota
	// IndexTBI emits a `.tbi` index (.vcf.gz only; refuses BCF input).
	IndexTBI
)

// IndexOptions controls `bcftools index`.
type IndexOptions struct {
	Format     IndexFormat
	MinShift   int32 // CSI only; 0 means the htslib default (14)
	OutputPath string
	Force      bool
}

// BuildIndex inspects the input file and dispatches to the right indexer.
// VCFs (gzipped) can be indexed as CSI or TBI; BCFs are always CSI.
func BuildIndex(inputPath string, opts IndexOptions) (string, error) {
	outPath := opts.OutputPath
	if outPath == "" {
		switch opts.Format {
		case IndexTBI:
			outPath = inputPath + ".tbi"
		default:
			outPath = inputPath + ".csi"
		}
	}
	if !opts.Force {
		if _, err := os.Stat(outPath); err == nil {
			return "", fmt.Errorf("bcftools index: %s already exists (use -f to overwrite)", outPath)
		}
	}
	isBCF, err := looksLikeBCF(inputPath)
	if err != nil {
		return "", err
	}
	switch {
	case isBCF && opts.Format == IndexTBI:
		return "", errors.New("bcftools index: --tbi is not valid for BCF input")
	case isBCF:
		csi, err := buildCSIForBCF(inputPath, opts.MinShift)
		if err != nil {
			return "", err
		}
		return outPath, csi.WriteFile(outPath)
	case opts.Format == IndexTBI:
		cfg, err := tabix.PresetConfig(tabix.PresetVCF)
		if err != nil {
			return "", err
		}
		idx, err := tabix.Build(inputPath, cfg)
		if err != nil {
			return "", err
		}
		return outPath, idx.WriteFile(outPath)
	default:
		csi, err := buildCSIForVCFGz(inputPath, opts.MinShift)
		if err != nil {
			return "", err
		}
		return outPath, csi.WriteFile(outPath)
	}
}

// looksLikeBCF returns true if the first decoded bytes of path are the BCF
// magic. We sniff through bgzip if the file is BGZF.
func looksLikeBCF(path string) (bool, error) {
	f, err := os.Open(path)
	if err != nil {
		return false, err
	}
	defer f.Close()
	br := bufio.NewReader(f)
	head, _ := br.Peek(16)
	if isBGZF(head) {
		rr, err := bgzip.NewReader(br)
		if err != nil {
			return false, err
		}
		defer rr.Close()
		var sig [5]byte
		if _, err := io.ReadFull(rr, sig[:]); err != nil {
			return false, nil
		}
		return sig[0] == 'B' && sig[1] == 'C' && sig[2] == 'F', nil
	}
	return len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F', nil
}

// isBGZF returns true if b starts with a BGZF block header (gzip magic + the
// `BC` extra-subfield used by htslib).
func isBGZF(b []byte) bool {
	if len(b) < 16 {
		return false
	}
	if b[0] != 0x1f || b[1] != 0x8b {
		return false
	}
	return b[12] == 'B' && b[13] == 'C'
}

// buildCSIForBCF traverses a BGZF-wrapped BCF file, recording per-record
// (chrom, pos, rlen, virtual offset) tuples and folding them into a CSI.
func buildCSIForBCF(path string, minShift int32) (*tabix.CSI, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	// Scan the BGZF blocks first; we need block offsets so we can resolve
	// each record's uncompressed start to a virtual offset.
	offsets, err := bgzip.Scan(f)
	if err != nil && !errors.Is(err, bgzip.ErrTruncated) {
		return nil, err
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		return nil, err
	}
	br, err := bgzip.NewReader(f)
	if err != nil {
		return nil, err
	}
	defer br.Close()
	data, err := io.ReadAll(br)
	if err != nil {
		return nil, err
	}

	// Helper: convert an uncompressed byte position into a virtual offset.
	uoffToV := func(pos int64) tabix.VOffset {
		lo, hi := 0, len(offsets)
		for lo < hi {
			mid := (lo + hi) / 2
			if int64(offsets[mid].UncompressedOffset) <= pos {
				lo = mid + 1
			} else {
				hi = mid
			}
		}
		i := lo - 1
		if i < 0 {
			i = 0
		}
		blk := offsets[i]
		uoff := int(pos - blk.UncompressedOffset)
		return tabix.MakeVOffset(blk.CompressedOffset, uoff)
	}

	// Parse the BCF header so we know where records start.
	rdr := bytes.NewReader(data)
	hdr, err := bcf.ReadHeader(rdr)
	if err != nil {
		return nil, err
	}
	bodyStart := int64(len(data)) - int64(rdr.Len())
	pos := bodyStart

	csi := tabix.NewCSI(minShift, 5)
	// Encode the contig names + tabix-style preset block so consumers can
	// look up refIDs by name (mirrors htslib's CSI-for-BCF auxiliary block).
	names := make([]string, len(hdr.Contigs))
	for i, c := range hdr.Contigs {
		names[i] = c.ID
	}
	csi.SetAuxFromTabix(tabix.Config{Format: tabix.FormatVCF, ColSeq: 1, ColBeg: 2, ColEnd: 0, Meta: '#', Skip: 0}, names)

	br2 := bcf.NewReaderWithHeader(hdr)
	for {
		recStart := pos
		rec, err := br2.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("bcftools index: read record: %w", err)
		}
		// Advance pos by the (l_shared + l_indiv + 8) bytes the reader just
		// consumed. We can recompute that from the record's encoded size.
		// Each record on the wire is l_shared:u32 + l_indiv:u32 + payload.
		// Recompute by recording the buffered position before/after the
		// read — but the bcf reader hides its bufio. Simpler: derive
		// payload bytes from the on-the-fly encode.
		shared, indiv, err := encodeRecordSize(rec)
		if err != nil {
			return nil, err
		}
		pos += 8 + int64(shared+indiv)
		if rec.ChromID < 0 {
			csi.NoCoor++
			continue
		}
		beg := int64(rec.Pos)
		end := beg + int64(rec.Rlen)
		if end <= beg {
			end = beg + 1
		}
		v := uoffToV(recStart)
		vEnd := uoffToV(pos)
		csi.AddRecord(int(rec.ChromID), beg, end, v, vEnd)
	}
	return csi, nil
}

// encodeRecordSize returns the (l_shared, l_indiv) of rec without writing
// any bytes. We do it the easy way — re-encode into a buffer and discard
// the result — because the reader doesn't expose the original byte run.
func encodeRecordSize(rec *bcf.Record) (int, int, error) {
	var buf bytes.Buffer
	w := bcf.NewWriter(&buf, &bcf.Header{}) // header is irrelevant for WriteRecord
	if err := w.WriteRecord(rec); err != nil {
		return 0, 0, err
	}
	if err := w.Flush(); err != nil {
		return 0, 0, err
	}
	all := buf.Bytes()
	// The Writer prepends the magic + length-prefixed text header on its
	// first call. We skip past that to find the (l_shared, l_indiv, body)
	// frame.
	off := 0
	if len(all) >= 5 && all[0] == 'B' && all[1] == 'C' && all[2] == 'F' {
		off = 5
		var lText uint32
		if off+4 > len(all) {
			return 0, 0, errors.New("bcftools index: encoded record too short")
		}
		lText = binary.LittleEndian.Uint32(all[off : off+4])
		off += 4 + int(lText)
	}
	if off+8 > len(all) {
		return 0, 0, errors.New("bcftools index: encoded record too short for frame")
	}
	lShared := binary.LittleEndian.Uint32(all[off : off+4])
	lIndiv := binary.LittleEndian.Uint32(all[off+4 : off+8])
	return int(lShared), int(lIndiv), nil
}

// buildCSIForVCFGz scans a bgzipped VCF and builds a CSI index using the
// same (chrom, pos, end) extraction the tabix builder uses.
func buildCSIForVCFGz(path string, minShift int32) (*tabix.CSI, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	offsets, err := bgzip.Scan(f)
	if err != nil && !errors.Is(err, bgzip.ErrTruncated) {
		return nil, err
	}
	if _, err := f.Seek(0, io.SeekStart); err != nil {
		return nil, err
	}
	br, err := bgzip.NewReader(f)
	if err != nil {
		return nil, err
	}
	defer br.Close()
	data, err := io.ReadAll(br)
	if err != nil {
		return nil, err
	}

	uoffToV := func(pos int64) tabix.VOffset {
		lo, hi := 0, len(offsets)
		for lo < hi {
			mid := (lo + hi) / 2
			if int64(offsets[mid].UncompressedOffset) <= pos {
				lo = mid + 1
			} else {
				hi = mid
			}
		}
		i := lo - 1
		if i < 0 {
			i = 0
		}
		blk := offsets[i]
		uoff := int(pos - blk.UncompressedOffset)
		return tabix.MakeVOffset(blk.CompressedOffset, uoff)
	}

	csi := tabix.NewCSI(minShift, 5)
	csi.SetAuxFromTabix(tabix.Config{Format: tabix.FormatVCF, ColSeq: 1, ColBeg: 2, ColEnd: 0, Meta: '#', Skip: 0}, nil)
	// Build the name list as we encounter chroms.
	nameID := map[string]int{}

	sc := bufio.NewScanner(bytes.NewReader(data))
	sc.Buffer(make([]byte, 0, 1<<16), 1<<24)
	var pos int64
	for sc.Scan() {
		raw := sc.Bytes()
		lineStart := pos
		pos += int64(len(raw)) + 1
		if len(raw) == 0 || raw[0] == '#' {
			continue
		}
		// Split on tabs to grab the first 4 columns (CHROM, POS, ID, REF).
		fields := bytes.SplitN(raw, []byte("\t"), 5)
		if len(fields) < 4 {
			continue
		}
		chrom := string(fields[0])
		var posVal int64
		for _, c := range fields[1] {
			if c < '0' || c > '9' {
				posVal = -1
				break
			}
			posVal = posVal*10 + int64(c-'0')
		}
		if posVal <= 0 {
			continue
		}
		refLen := int64(len(fields[3]))
		if refLen <= 0 {
			refLen = 1
		}
		id, ok := nameID[chrom]
		if !ok {
			id = len(nameID)
			nameID[chrom] = id
		}
		// VCF is 1-based; CSI bin scheme uses 0-based half-open.
		beg := posVal - 1
		end := beg + refLen
		v := uoffToV(lineStart)
		vEnd := uoffToV(pos)
		csi.AddRecord(id, beg, end, v, vEnd)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	// Rebuild the name slice in name-id order.
	names := make([]string, len(nameID))
	for n, i := range nameID {
		names[i] = n
	}
	csi.SetAuxFromTabix(tabix.Config{Format: tabix.FormatVCF, ColSeq: 1, ColBeg: 2, ColEnd: 0, Meta: '#', Skip: 0}, names)
	return csi, nil
}

// ChromIDInCSI returns the dictionary index of name in c.Names, or -1 if
// the name is unknown.
func ChromIDInCSI(c *tabix.CSI, name string) int {
	for i, n := range c.Names {
		if n == name {
			return i
		}
	}
	return -1
}
