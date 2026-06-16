// Haplotype-aware consequence engine for bcftools csq.
//
// This file ports the haplotype engine from upstream csq.c: the
// hap_node_t haplotype tree, hap_init / hap_finalize / hap_add_csq,
// cds_translate, the splice_csq family with set_refalt (splice_build_hap),
// the vbuf / pos2vbuf variant buffer, csq_push / csq_stage, and the
// -p/--phase and -n/--ncsq handling.
//
// Variants overlapping a coding region are phased onto per-sample
// haplotypes and walked together, so compound consequences (multiple
// variants in the same codon) are called jointly and the reference vs
// haplotype-altered CDS are translated and diffed to produce the
// amino-acid change string and the true frameshift / inframe /
// elongation / truncation calls. UTR / splice / intron / non-coding
// consequences are routed through the same vbuf so the emitted
// INFO/BCSQ matches upstream byte for byte.

package bcftools

import (
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// nRefPad mirrors csq.c's N_REF_PAD: the number of bases padded onto
// both ends of the cached transcript reference to avoid boundary
// effects.
const nRefPad = 10

// Haplotype node types, mirroring HAP_CDS / HAP_ROOT / HAP_SSS.
const (
	hapCDS  = 0
	hapRoot = 1
	hapSSS  = 2 // start/stop/splice node, no codon prediction
)

// Phasing modes, mirroring PHASE_* in csq.c.
const (
	phaseRequire = 0 // -p r
	phaseMerge   = 1 // -p m
	phaseAsIs    = 2 // -p a
	phaseSkip    = 3 // -p s
	phaseNonRef  = 4 // -p R
	phaseDropGT  = 5 // no samples
)

// csqStartStop and csqPrintedUpstream / csqUpstreamStop / csqIncompleteCDS
// complete the bit set used by the engine (csq_classify.go defines the
// rest).
const (
	csqPrintedUpstream = 1 << 0
	csqUpstreamStop    = 1 << 19
	csqIncompleteCDS   = 1 << 20
)

// csqStartStop is the CSQ_START_STOP macro.
const csqStartStop = csqStopLost | csqStopGained | csqStopRetained | csqStartLost | csqStartRetained

// csqCompoundFull is the complete CSQ_COMPOUND macro from csq.c
// (csq_classify.go's csqCompound omits the engine-only bits).
const csqCompoundFull = csqSynonymous | csqMissense | csqStopLost | csqStopGained |
	csqInframeDel | csqInframeIns | csqFrameshift |
	csqStartLost | csqStopRetained | csqInframeAlter | csqIncompleteCDS |
	csqUpstreamStop | csqStartRetained | csqElongation | csqTruncation

// gencode is a single NCBI genetic-code table, the Go analogue of
// csq.c's gencode_t. Code holds the amino acid and Stop holds the
// start/stop annotation ('*', 'M' or '-') for each of the 64 codons,
// indexed as idx = (nt4[c0]<<4)|(nt4[c1]<<2)|nt4[c2] with
// A=0,C=1,G=2,T=3.
type gencode struct {
	// ID is the NCBI translation-table number selected by
	// -C/--genetic-code (with 0 being bcftools' "standard simplified"
	// variant of table 1).
	ID int
	// Name is the human-readable table name.
	Name string
	// Code maps each of the 64 codons to its amino-acid letter.
	Code string
	// Stop marks start ('M') and stop ('*') codons; '-' otherwise.
	Stop string
}

// gencodeTables transcribes csq.c's gencode_tables (generated upstream
// by misc/gencode-tables from the NCBI translation tables). Table 0 is
// bcftools' "standard simplified" variant of NCBI table 1. New tables
// can be added by appending the corresponding NCBI entry.
var gencodeTables = []gencode{
	{ID: 0, Name: "Standard simplified",
		Code: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
		Stop: "--------------M---------------------------------*-*-----*-------"},
	{ID: 1, Name: "Standard",
		Code: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
		Stop: "--------------M---------------M-----------------*-*-----*-----M-"},
	{ID: 2, Name: "Vertebrate Mitochondrial",
		Code: "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		Stop: "--------*-*-MMMM------------------------------M-*-*-------------"},
	{ID: 3, Name: "Yeast Mitochondrial",
		Code: "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		Stop: "------------M-M-------------------------------M-*-*-------------"},
	{ID: 5, Name: "Invertebrate Mitochondrial",
		Code: "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
		Stop: "------------MMMM------------------------------M-*-*-----------M-"},
}

// gencodeByID returns the genetic-code table for the given NCBI id (0
// is bcftools' standard simplified table). It reports false when no
// table with that id has been transcribed.
func gencodeByID(id int) (*gencode, bool) {
	for i := range gencodeTables {
		if gencodeTables[i].ID == id {
			return &gencodeTables[i], true
		}
	}
	return nil, false
}

// standardGencode is the default table (NCBI 0 / standard simplified),
// used by the engine when -C/--genetic-code is unset.
var standardGencode = &gencodeTables[0]

// GeneticCodeKnown reports whether the given NCBI translation-table id
// is supported by `bcftools csq -C/--genetic-code`.
func GeneticCodeKnown(id int) bool {
	_, ok := gencodeByID(id)
	return ok
}

// GeneticCodeIDs returns the supported -C/--genetic-code table ids as a
// comma-separated string (for error messages), e.g. "0, 1, 2, 3, 5".
func GeneticCodeIDs() string {
	var sb strings.Builder
	for i, gc := range gencodeTables {
		if i > 0 {
			sb.WriteString(", ")
		}
		sb.WriteString(strconv.Itoa(gc.ID))
	}
	return sb.String()
}

// GeneticCodeListing renders the supported genetic-code tables in the
// same id/name/code/stop layout upstream prints for `-C l`.
func GeneticCodeListing() string {
	var sb strings.Builder
	for _, gc := range gencodeTables {
		sb.WriteString(strconv.Itoa(gc.ID))
		sb.WriteByte('\t')
		sb.WriteString(gc.Name)
		sb.WriteString("\n\t")
		sb.WriteString(gc.Code)
		sb.WriteString("\n\t")
		sb.WriteString(gc.Stop)
		sb.WriteString("\n\n")
	}
	return sb.String()
}

// nt4 maps a DNA base to 0..3 (A,C,G,T); 4 for anything else.
var nt4 = func() [256]uint8 {
	var t [256]uint8
	for i := range t {
		t[i] = 4
	}
	t['A'], t['a'] = 0, 0
	t['C'], t['c'] = 1, 1
	t['G'], t['g'] = 2, 2
	t['T'], t['t'] = 3, 3
	return t
}()

// cnt4 maps a DNA base to the index of its complement (A<->T, C<->G).
var cnt4 = func() [256]uint8 {
	var t [256]uint8
	for i := range t {
		t[i] = 4
	}
	t['A'], t['a'] = 3, 3
	t['C'], t['c'] = 2, 2
	t['G'], t['g'] = 1, 1
	t['T'], t['t'] = 0, 0
	return t
}()

// dna2aa translates a forward-strand codon to its amino acid using the
// selected genetic-code table.
func (gc *gencode) dna2aa(c0, c1, c2 byte) byte {
	idx := int(nt4[c0])<<4 | int(nt4[c1])<<2 | int(nt4[c2])
	if idx > 63 {
		return 'X'
	}
	return gc.Code[idx]
}

// dna2stop returns the stop/start annotation ('*', 'M' or '-') of a
// forward-strand codon under the selected genetic-code table.
func (gc *gencode) dna2stop(c0, c1, c2 byte) byte {
	idx := int(nt4[c0])<<4 | int(nt4[c1])<<2 | int(nt4[c2])
	if idx > 63 {
		return 0
	}
	return gc.Stop[idx]
}

// cdna2aa translates a reverse-strand codon (read on the forward
// strand) to its amino acid using the selected genetic-code table.
func (gc *gencode) cdna2aa(c0, c1, c2 byte) byte {
	idx := int(cnt4[c2])<<4 | int(cnt4[c1])<<2 | int(cnt4[c0])
	if idx > 63 {
		return 'X'
	}
	return gc.Code[idx]
}

// cdna2stop returns the stop/start annotation of a reverse-strand
// codon read on the forward strand under the selected genetic-code
// table.
func (gc *gencode) cdna2stop(c0, c1, c2 byte) byte {
	idx := int(cnt4[c2])<<4 | int(cnt4[c1])<<2 | int(cnt4[c0])
	if idx > 63 {
		return 0
	}
	return gc.Stop[idx]
}

// strandCode converts a gff.Strand to csq.c's STRAND_FWD / STRAND_REV.
func strandCode(s gff.Strand) int {
	if s == gff.StrandReverse {
		return strandRev
	}
	return strandFwd
}

const (
	strandRev = 0
	strandFwd = 1
)

// hapTranscript is the engine's per-transcript working state, the Go
// analogue of csq.c's tscript_t. It is created lazily when the first
// coding variant overlapping the transcript is seen.
type hapTranscript struct {
	tr    *CSQTranscript
	beg0  int    // transcript start, 0-based
	end0  int    // transcript end, 0-based inclusive
	ref   []byte // padded reference: nRefPad + [beg..end] + nRefPad
	sref  []byte // spliced (CDS-only) reference, padded
	root  *hapNode
	hap   []*hapNode // per-haplotype leaves; len == 1 (no GT) or 2*nsmpl
	nsref int        // len(sref)
}

// hapNode is one node of the haplotype tree, the Go analogue of
// hap_node_t. The root represents the unaltered transcript; each child
// applies one more variant.
type hapNode struct {
	seq      string // CDS segment [parent_node, this_node)
	variant  string // "ref>alt"
	typ      int    // hapRoot / hapCDS / hapSSS
	csq      uint32 // this node's per-record consequence bits
	dlen     int    // alt length minus ref length
	rbeg     int    // variant VCF position, 0-based inclusive
	rlen     int    // variant rlen on the spliced reference
	sbeg     int    // position on the spliced reference, 0-based
	icds     int    // index of the overlapped CDS exon
	child    []*hapNode
	prev     *hapNode
	curRec   *vcf.Variant
	rec      *vcf.Variant
	vcfIal   int
	nend     int
	curChild []int // allele -> active child index, reset per record
	csqList  []*csqEntry
}

// csqEntry is the engine's analogue of csq_t: a top-level consequence
// tied to a haplotype node and a VCF record.
type csqEntry struct {
	pos  int // VCF position, 0-based
	vrec *vrecBuf
	idx  int
	typ  vcsq
}

// vcsq mirrors vcsq_t: everything needed to render one BCSQ entry.
type vcsq struct {
	strand  int
	typ     uint32
	trid    string
	vcfIal  int
	biotype string
	gene    string
	ref     *vcf.Variant // for CSQ_PRINTED_UPSTREAM, the @pos back-reference
	vstr    string       // variant string, eg "|2V>2I|103G>A"
	hasVstr bool
}

// txtCsq mirrors txt_csq_t: one staged (consequence-index, sample,
// haplotype) tuple for the -O t streaming-text output. ismpl is -1 when
// genotypes are dropped (the no-GT, sample-agnostic path); ihap is the
// 1-based haplotype number, or 0 when no haplotype applies.
type txtCsq struct {
	idx   int
	ismpl int
	ihap  int
}

// vrecBuf mirrors vrec_t: a single VCF record plus the consequences
// staged against it.
type vrecBuf struct {
	rec  *vcf.Variant
	vcsq []vcsq
	// fmtBM is the per-sample FORMAT/BCSQ bitmask block, length
	// nsamples*nfmtBcsq when populated. Bit (2*idx+ihap) at offset
	// nfmtBcsq*ismpl marks "INFO/BCSQ entry idx applies to haplotype
	// ihap of sample ismpl". Mirrors vrec_t.fmt_bm.
	fmtBM []uint32
	// nfmt tracks the highest int32 index touched by any bit, so the
	// emitted FORMAT/BCSQ can be trimmed to the minimum width.
	nfmt int
	// txt holds the staged (idx, ismpl, ihap) tuples for -O t text
	// output, in staging order. Mirrors vrec_t.txt. Only populated when
	// the engine runs in textMode.
	txt []txtCsq
}

// vbuf mirrors vbuf_t: VCF records sharing a position.
type vbuf struct {
	vrec      []*vrecBuf
	keepUntil int // maximum transcript end position seen
}

// hapEngine drives the whole haplotype-aware pipeline for one VCF
// stream. It owns the variant buffer and the active-transcript set.
type hapEngine struct {
	idx     *CSQIndex
	opts    CSQOptions
	hdr     *vcf.Header
	phase   int
	samples []int // header indices of samples to process; nil when phaseDropGT
	ncsq2   int
	// gencode is the selected NCBI translation table (-C/--genetic-code),
	// defaulting to the standard simplified table when unset.
	gencode *gencode
	// brief is the -b/--brief-predictions / -B/--trim-protein-seq value:
	// 0 means full amino-acid predictions, N>0 truncates each prediction
	// to the first N residues. Mirrors upstream's args->brief_predictions.
	brief int
	// nfmtBcsq is the per-sample width of FORMAT/BCSQ (number of
	// int32 ints needed to hold ncsq2 effective bits, 30 per int).
	// Mirrors upstream's args->nfmt_bcsq.
	nfmtBcsq int

	vcfBuf   []*vbuf       // round buffer of buffered VCF lines, ordered by pos
	pos2vbuf map[int]*vbuf // position -> vbuf
	hapTr    map[string]*hapTranscript
	activeTr []*hapTranscript // transcripts still receiving variants
	rid      string           // current contig

	out []*vcf.Variant // finalised, ready-to-write records in order

	// textMode selects upstream's -O t streaming-text output
	// (FT_TAB_TEXT). When set, the engine stages per-(sample,haplotype)
	// consequence tuples (vrecBuf.txt) instead of FORMAT/BCSQ bitmasks
	// and renders them into outText at flush time.
	textMode bool
	// outText accumulates the finalised "CSQ\t..." text lines (no
	// header), in emission order. Only populated when textMode is set.
	outText []string
}

// newHapEngine constructs an engine for the given index and options.
func newHapEngine(idx *CSQIndex, opts CSQOptions, hdr *vcf.Header) *hapEngine {
	e := &hapEngine{
		idx:      idx,
		opts:     opts,
		hdr:      hdr,
		pos2vbuf: map[int]*vbuf{},
		hapTr:    map[string]*hapTranscript{},
	}
	// Select the genetic-code table. Unknown ids fall back to the
	// standard table; the CLI validates -C/--genetic-code up front so
	// the engine never sees an unknown id in practice.
	if gc, ok := gencodeByID(opts.GeneticCode); ok {
		e.gencode = gc
	} else {
		e.gencode = standardGencode
	}
	// brief_predictions: -b sets 1, -B sets N. The CLI maps both onto
	// opts.TrimProteinSeq (1 for -b), matching upstream where -b and -B
	// share args->brief_predictions.
	e.brief = opts.TrimProteinSeq
	e.phase = phaseByteToMode(opts.Phase)
	if hdr == nil || len(hdr.Samples) == 0 {
		e.phase = phaseDropGT
	}
	if e.phase != phaseDropGT {
		// All samples are processed (subsetting is a later slice).
		e.samples = make([]int, len(hdr.Samples))
		for i := range e.samples {
			e.samples[i] = i
		}
	}
	// ncsq2 mirrors upstream's args->ncsq2 (2*--ncsq): the per-haplotype
	// consequence cap. slice 4: consumed by FORMAT/BCSQ emission.
	e.ncsq2 = opts.NCSQ
	if e.ncsq2 <= 0 {
		e.ncsq2 = 16 // upstream default --ncsq
	}
	e.ncsq2 *= 2
	e.nfmtBcsq = ncsq2ToNfmt(e.ncsq2)
	e.textMode = opts.TextOutput
	return e
}

// ncsq2ToNfmt mirrors upstream's ncsq2_to_nfmt: the number of 32-bit
// ints needed to hold ncsq2 effective bits, with 30 bits per int (the
// top two bits are reserved to avoid BCF missing/end values).
func ncsq2ToNfmt(ncsq2 int) int {
	if ncsq2 <= 0 {
		return 1
	}
	return 1 + (ncsq2-1)/30
}

// icsq2ToBit maps an icsq2 (2*idx+ihap) to its int32 index and bit
// position. Mirrors upstream's icsq2_to_bit.
func icsq2ToBit(icsq2 int) (ival, ibit int) {
	return icsq2 / 30, icsq2 % 30
}

// phaseByteToMode maps the -p/--phase byte (a|m|r|R|s) to a phase*
// constant. An unset (zero) byte defaults to require, as upstream does.
func phaseByteToMode(b byte) int {
	switch b {
	case 'a':
		return phaseAsIs
	case 'm':
		return phaseMerge
	case 'R':
		return phaseNonRef
	case 's':
		return phaseSkip
	default: // 'r' or unset
		return phaseRequire
	}
}

// csqGTAlleles parses a GT string into allele indices and a phased
// flag. Missing alleles are reported as -1.
func csqGTAlleles(gt string) (alleles []int, phased bool) {
	if gt == "" || gt == "." {
		return nil, false
	}
	phased = strings.ContainsRune(gt, '|')
	sep := func(r rune) bool { return r == '|' || r == '/' }
	for _, f := range strings.FieldsFunc(gt, sep) {
		if f == "." || f == "" {
			alleles = append(alleles, -1)
			continue
		}
		n, err := strconv.Atoi(f)
		if err != nil {
			alleles = append(alleles, -1)
			continue
		}
		alleles = append(alleles, n)
	}
	return alleles, phased
}

// getTranscript lazily builds the per-transcript engine state and its
// padded reference. Returns nil if the contig sequence is unavailable.
func (e *hapEngine) getTranscript(t *CSQTranscript) *hapTranscript {
	if ht, ok := e.hapTr[t.ID]; ok {
		return ht
	}
	refSeq := e.idx.Refs[t.Chrom]
	if len(refSeq) == 0 {
		e.hapTr[t.ID] = nil
		return nil
	}
	beg0 := t.Beg - 1
	end0 := t.End - 1
	padBeg := nRefPad
	if beg0 < padBeg {
		padBeg = beg0
	}
	total := (end0 - beg0 + 1) + 2*nRefPad
	ref := make([]byte, total)
	for i := range ref {
		ref[i] = 'N'
	}
	// Copy [beg0-padBeg .. end0+nRefPad] of the contig into ref so that
	// ref[nRefPad] aligns with the transcript's first base.
	srcStart := beg0 - padBeg
	for i := 0; i < total; i++ {
		dst := i + (nRefPad - padBeg)
		if dst < 0 || dst >= total {
			continue
		}
		src := srcStart + i
		if src >= 0 && src < len(refSeq) {
			ref[dst] = upper(refSeq[src])
		}
	}
	ht := &hapTranscript{
		tr:   t,
		beg0: beg0,
		end0: end0,
		ref:  ref,
	}
	ht.buildSplicedRef()
	ht.root = &hapNode{typ: hapRoot}
	nhap := 1
	if e.phase != phaseDropGT {
		nhap = 2 * len(e.samples)
	}
	ht.hap = make([]*hapNode, nhap)
	ht.root.nend = nhap
	e.hapTr[t.ID] = ht
	e.activeTr = append(e.activeTr, ht)
	return ht
}

// buildSplicedRef concatenates the CDS exons of the transcript into the
// spliced reference, padding nRefPad bases on each end. It mirrors
// tscript_splice_ref.
func (ht *hapTranscript) buildSplicedRef() {
	if len(ht.tr.CDSExons) == 0 {
		// Non-coding transcript: there is no spliced CDS to build. The
		// splice/intron/non-coding tests run directly on ht.ref.
		ht.nsref = 2 * nRefPad
		ht.sref = make([]byte, ht.nsref)
		return
	}
	var length int
	for _, c := range ht.tr.CDSExons {
		length += c.End - c.Start + 1
	}
	ht.nsref = length + 2*nRefPad
	sref := make([]byte, ht.nsref)
	cds := ht.tr.CDSExons
	// Left padding: nRefPad bases preceding the first CDS exon.
	copy(sref[:nRefPad], ht.ref[cds[0].Start-1-ht.beg0:cds[0].Start-1-ht.beg0+nRefPad])
	pos := nRefPad
	for _, c := range cds {
		clen := c.End - c.Start + 1
		copy(sref[pos:pos+clen], ht.ref[nRefPad+c.Start-1-ht.beg0:nRefPad+c.Start-1-ht.beg0+clen])
		pos += clen
	}
	last := cds[len(cds)-1]
	copy(sref[pos:pos+nRefPad], ht.ref[nRefPad+last.End-ht.beg0:nRefPad+last.End-ht.beg0+nRefPad])
	ht.sref = sref
}

// cdsPos returns the offset of a CDS exon's first base within the
// spliced reference (0-based, excluding nRefPad padding).
func (ht *hapTranscript) cdsPos(icds int) int {
	pos := 0
	for i := 0; i < icds; i++ {
		pos += ht.tr.CDSExons[i].End - ht.tr.CDSExons[i].Start + 1
	}
	return pos
}
