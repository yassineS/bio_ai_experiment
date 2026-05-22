// Shared consequence-type constants and the SO-term string table for
// bcftools csq.
//
// These mirror the CSQ_* macros and csq_strings[] from upstream csq.c.
// They are consumed by the haplotype engine (csq_splice.go /
// csq_process.go / csq_engine.go / csq_hap.go), which is the sole
// production consequence path. The former v1 per-record classifier
// (splice_csq family, test_cds/utr/splice/tscript dispatch) that also
// lived in this file has been removed — its splice_csq port now exists
// exactly once, in csq_splice.go.

package bcftools

// Consequence-type bits, mirroring the CSQ_* macros in csq.c. The bit
// positions match upstream so the precedence walk in kput_vcsq can
// iterate them in the same order.
const (
	csqSynonymous     = 1 << 1
	csqMissense       = 1 << 2
	csqStopLost       = 1 << 3
	csqStopGained     = 1 << 4
	csqInframeDel     = 1 << 5
	csqInframeIns     = 1 << 6
	csqFrameshift     = 1 << 7
	csqSpliceAcceptor = 1 << 8
	csqSpliceDonor    = 1 << 9
	csqStartLost      = 1 << 10
	csqSpliceRegion   = 1 << 11
	csqStopRetained   = 1 << 12
	csqUTR5           = 1 << 13
	csqUTR3           = 1 << 14
	csqNonCoding      = 1 << 15
	csqIntron         = 1 << 16
	csqInframeAlter   = 1 << 18
	csqCodingSeq      = 1 << 21
	csqElongation     = 1 << 22
	csqTruncation     = 1 << 23
	csqStartRetained  = 1 << 24
)

// csqStrings maps a consequence bit index to its SO-term string. The
// slot ordering is verbatim from csq.c's csq_strings[]; the precedence
// walk in kput_vcsq relies on it.
var csqStrings = [...]string{
	1:  "synonymous",
	2:  "missense",
	3:  "stop_lost",
	4:  "stop_gained",
	5:  "inframe_deletion",
	6:  "inframe_insertion",
	7:  "frameshift",
	8:  "splice_acceptor",
	9:  "splice_donor",
	10: "start_lost",
	11: "splice_region",
	12: "stop_retained",
	13: "5_prime_utr",
	14: "3_prime_utr",
	15: "non_coding",
	16: "intron",
	18: "inframe_altering",
	21: "coding_sequence",
	22: "feature_elongation",
	23: "feature_truncation",
	24: "start_retained",
}

// splice region/donor window sizes, from gff.h.
const (
	nSpliceDonor        = 2 // 2bp at the intron edge -> donor/acceptor
	nSpliceRegionIntron = 8 // up to 8bp into the intron -> splice_region
	nSpliceRegionExon   = 3 // up to 3bp into the exon -> splice_region
)

// splice-csq return codes, mirroring the SPLICE_* macros.
const (
	spliceVarRef  = 0 // ref==alt, not a variant
	spliceOutside = 1 // csq set, does not overlap the coding region
	spliceInside  = 2 // overlaps the coding region, prediction needed
	spliceOverlap = 3 // indel overlaps the region boundary
)
