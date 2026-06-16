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
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// tbxMaxShift mirrors htslib's TBX_MAX_SHIFT (tbx.h). CSI for VCF/BCF derives
// its starting number of bin levels from this and the requested min-shift.
const tbxMaxShift = 31

// htsBinMaxPos returns the (0-based, exclusive) maximum position a binning
// index addresses for the given (minShift, nLvls), matching htslib's
// hts_bin_maxpos (hts.h): 1 << (min_shift + n_lvls*3).
func htsBinMaxPos(minShift, nLvls int32) int64 {
	return int64(1) << uint32(minShift+nLvls*3)
}

// htsAdjustCSISettings ports htslib's hts_adjust_csi_settings (hts.c). Given
// the longest reference length, it grows the number of bin levels (and, only
// when the level ceiling is hit, the min-shift) so the index addresses every
// coordinate. It returns the adjusted (minShift, nLvls). The +256 slack and
// the max_n_lvls=9 ceiling match htslib exactly.
func htsAdjustCSISettings(maxLenIn int64, minShift, nLvls int32) (int32, int32) {
	const maxNLvls = 9
	maxLen := maxLenIn + 256
	if maxLen <= htsBinMaxPos(minShift, maxNLvls) {
		maxpos := htsBinMaxPos(minShift, nLvls)
		for maxLen > maxpos {
			nLvls++
			maxpos *= 8
		}
		return minShift, nLvls
	}
	nLvls = maxNLvls
	maxpos := htsBinMaxPos(minShift, nLvls)
	for maxLen > maxpos {
		minShift++
		maxpos *= 2
	}
	return minShift, nLvls
}

// maxContigLenFromMeta returns the largest `##contig=<...,length=N>` value in a
// VCF-style header's meta lines, mirroring htslib's idx_calc_n_lvls_ids, which
// scans the contig dictionary's stored length (info[0]). When no contig length
// is found it returns 0 so callers can apply htslib's broken-header fallback of
// 2^31-1.
func maxContigLenFromMeta(meta []string) int64 {
	var max int64
	for _, m := range meta {
		if !strings.HasPrefix(m, "##contig=") {
			continue
		}
		idx := strings.Index(m, "length")
		if idx < 0 {
			continue
		}
		rest := m[idx+len("length"):]
		// Skip the run of spaces / '=' that separates the key from the value.
		i := 0
		for i < len(rest) && (rest[i] == ' ' || rest[i] == '=') {
			i++
		}
		j := i
		for j < len(rest) && rest[j] >= '0' && rest[j] <= '9' {
			j++
		}
		if j == i {
			continue
		}
		if n, err := strconv.ParseInt(rest[i:j], 10, 64); err == nil && n > max {
			max = n
		}
	}
	return max
}

// countContigsInMeta returns the number of `##contig=` lines in a VCF-style
// header. htslib's BCF CSI emits one (possibly empty) reference slot per
// header contig, so the on-disk n_ref equals this count.
func countContigsInMeta(meta []string) int {
	n := 0
	for _, m := range meta {
		if strings.HasPrefix(m, "##contig=") {
			n++
		}
	}
	return n
}

// csiNLvlsBCF returns the (minShift, nLvls) htslib uses for a BCF CSI: a
// starting level count of 0 grown by the longest header contig
// (idx_calc_n_lvls_ids with starting_n_lvls==0 in bcf_idx_init / bcf_index).
func csiNLvlsBCF(minShift int32, meta []string) (int32, int32) {
	maxLen := maxContigLenFromMeta(meta)
	if maxLen == 0 {
		maxLen = (int64(1) << 31) - 1 // broken-header fallback, matching htslib
	}
	return htsAdjustCSISettings(maxLen, minShift, 0)
}

// csiNLvlsVCF returns the (minShift, nLvls) htslib uses for a VCF.gz CSI: a
// starting level count of (TBX_MAX_SHIFT-min_shift+2)/3 grown by the longest
// header contig (vcf_idx_init / idx_calc_n_lvls_ids).
func csiNLvlsVCF(minShift int32, meta []string) (int32, int32) {
	starting := (int32(tbxMaxShift) - minShift + 2) / 3
	if starting < 0 {
		starting = 0
	}
	maxLen := maxContigLenFromMeta(meta)
	if maxLen == 0 {
		maxLen = (int64(1) << 31) - 1
	}
	return htsAdjustCSISettings(maxLen, minShift, starting)
}

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
// magic. We sniff through bgzip if the file is BGZF. The open is routed through
// openSeekable so a remote URL (http(s)/s3/gs) is sniffed via a ranged read of
// just its leading bytes rather than a full download.
func looksLikeBCF(path string) (bool, error) {
	f, err := openSeekable(path)
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

// csiRefMeta accumulates the per-reference statistics htslib stores in the
// CSI metadata pseudo-bin: the virtual-offset span of the reference's records
// and the mapped/unmapped record counts.
type csiRefMeta struct {
	offBeg   tabix.VOffset
	offEnd   tabix.VOffset
	nMapped  uint64
	hasBeg   bool
	nUnmappd uint64
}

// observe folds one record's virtual-offset span into the per-ref meta.
func (m *csiRefMeta) observe(vBeg, vEnd tabix.VOffset) {
	if !m.hasBeg || vBeg < m.offBeg {
		m.offBeg = vBeg
		m.hasBeg = true
	}
	if vEnd > m.offEnd {
		m.offEnd = vEnd
	}
	m.nMapped++
}

// finishCSI mirrors htslib's hts_idx_finish: it (1) extends the very last
// record's data-bin chunk to the file's final virtual offset, (2) appends the
// per-reference metadata pseudo-bin (id BinLimit()+1, htslib's idx->n_bins+1)
// carrying the reference's virtual-offset span and the (n_mapped, n_unmapped)
// record counts, with the last reference's span likewise extended to the final
// offset. lastRef/lastBin identify the last pushed record's reference and bin;
// they are <0 for an empty index. final is bgzf_tell at EOF — the virtual
// offset of the BGZF EOF block. Matching this makes the emitted CSI
// byte-identical to `bcftools index` and the -W writer.
func finishCSI(csi *tabix.CSI, metas []csiRefMeta, lastRef, lastBin int, final tabix.VOffset) {
	// (1) The last record's chunk end is bumped to the final offset.
	if lastRef >= 0 && lastRef < len(csi.Refs) {
		for i := range csi.Refs[lastRef].Bins {
			if int(csi.Refs[lastRef].Bins[i].ID) == lastBin {
				ch := csi.Refs[lastRef].Bins[i].Chunks
				if n := len(ch); n > 0 && ch[n-1].End < final {
					ch[n-1].End = final
				}
				break
			}
		}
	}
	// (2) The metadata pseudo-bin, one per reference that has records.
	metaBinID := csi.BinLimit() + 1
	for r := range metas {
		m := metas[r]
		if !m.hasBeg || r >= len(csi.Refs) {
			continue
		}
		offEnd := m.offEnd
		if r == lastRef && final > offEnd {
			offEnd = final
		}
		csi.Refs[r].Bins = append(csi.Refs[r].Bins, tabix.CSIBin{
			ID:      metaBinID,
			LOffset: 0,
			Chunks: []tabix.CSIChunk{
				{Beg: m.offBeg, End: offEnd},
				{Beg: tabix.VOffset(m.nMapped), End: tabix.VOffset(m.nUnmappd)},
			},
		})
	}
}

// eofVOffset returns the virtual offset of the BGZF EOF block (bgzf_tell at
// end-of-file): the compressed offset just past the final data block, with a
// zero in-block offset. offsets is the block list from bgzip.Scan (which omits
// the EOF marker itself but lets us derive where it begins).
func eofVOffset(offsets []bgzip.BlockOffset) tabix.VOffset {
	if len(offsets) == 0 {
		return 0
	}
	last := offsets[len(offsets)-1]
	return tabix.MakeVOffset(last.CompressedOffset+int64(last.CompressedSize), 0)
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

	// Parse the BCF header so we know where records start. We must NOT derive
	// the body start from the reader's remaining length: bcf.ReadHeader buffers
	// ahead via an internal bufio, draining the bytes.Reader past the header and
	// into the record stream. Instead compute the on-wire header length directly
	// from the BCF framing — magic(5) + l_text(4) + l_text — which is exactly
	// where the first record begins (htslib flushes a BGZF block boundary there,
	// so this start coincides with a block boundary and matches the virtual
	// offsets `bcftools index` records).
	rdr := bytes.NewReader(data)
	hdr, err := bcf.ReadHeader(rdr)
	if err != nil {
		return nil, err
	}
	if len(data) < 9 {
		return nil, errors.New("bcftools index: BCF too short for header")
	}
	lText := binary.LittleEndian.Uint32(data[5:9])
	bodyStart := int64(9) + int64(lText)
	pos := bodyStart

	// Derive the bin-level count the way htslib's bcf_idx_init does: start at 0
	// and grow by the longest header contig. BCF CSI carries NO tabix-style
	// auxiliary block (l_aux == 0), unlike the VCF.gz CSI; chrom→refID lookups
	// for BCF go through the BCF header dictionary, not the index.
	effMinShift := minShift
	if effMinShift <= 0 {
		effMinShift = 14
	}
	ms, depth := csiNLvlsBCF(effMinShift, hdr.VCF.MetaInfo)
	// Construct directly rather than via NewCSI: a BCF CSI legitimately has
	// depth 0 (short contigs), which NewCSI would clamp up to the BAI default 5.
	csi := &tabix.CSI{MinShift: ms, Depth: depth}
	// htslib writes one (possibly empty) reference slot per header contig, so
	// pre-size the reference list to the full contig count; contigs without any
	// records are serialised as empty refs (n_bin == 0).
	nContigs := countContigsInMeta(hdr.VCF.MetaInfo)
	if nContigs > len(csi.Refs) {
		csi.Refs = make([]tabix.CSIRef, nContigs)
	}
	metas := make([]csiRefMeta, nContigs)
	lastRef, lastBin := -1, -1

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
		if rid := int(rec.ChromID); rid >= 0 && rid < len(metas) {
			metas[rid].observe(v, vEnd)
			lastRef = rid
			lastBin = int(csi.Reg2bin(beg, end))
		}
	}
	finishCSI(csi, metas, lastRef, lastBin, eofVOffset(offsets))
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

	// First pass: collect the ##contig header lines so the bin-level count can
	// be derived from the longest contig (vcf_idx_init -> idx_calc_n_lvls_ids).
	var contigMeta []string
	{
		hsc := bufio.NewScanner(bytes.NewReader(data))
		hsc.Buffer(make([]byte, 0, 1<<16), 1<<24)
		for hsc.Scan() {
			line := hsc.Text()
			if len(line) == 0 || line[0] != '#' {
				break
			}
			if strings.HasPrefix(line, "##contig=") {
				contigMeta = append(contigMeta, line)
			}
		}
	}
	effMinShift := minShift
	if effMinShift <= 0 {
		effMinShift = 14
	}
	ms, depth := csiNLvlsVCF(effMinShift, contigMeta)
	csi := &tabix.CSI{MinShift: ms, Depth: depth}
	csi.SetAuxFromTabix(tabix.Config{Format: tabix.FormatVCF, ColSeq: 1, ColBeg: 2, ColEnd: 0, Meta: '#', Skip: 0}, nil)
	// Build the name list as we encounter chroms.
	nameID := map[string]int{}
	var metas []csiRefMeta
	lastRef, lastBin := -1, -1

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
		for id >= len(metas) {
			metas = append(metas, csiRefMeta{})
		}
		metas[id].observe(v, vEnd)
		lastRef = id
		lastBin = int(csi.Reg2bin(beg, end))
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	finishCSI(csi, metas, lastRef, lastBin, eofVOffset(offsets))
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
