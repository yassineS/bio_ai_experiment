// Package edgecases holds the silent-corruption edge-case test battery: one
// focused, named Go test per class of failure where a bioinformatics tool can
// produce *wrong* output without crashing or warning. Each test name maps to
// the manuscript's ranked list of silent-corruption risks:
//
//  1. TestCRAMReferenceHandling      — reference-relative CRAM base corruption.
//  2. TestBCFToolsNormReindex        — multiallelic split/join Number=A/R/G re-indexing.
//  3. TestIndexByteIdentity          — .bai/.csi/.tbi byte-identity vs upstream.
//  4. TestSortStabilityStrnumCmp     — coordinate & queryname (strnum) sort tie-breaks.
//  5. TestCalmdMDNMTags              — MD/NM recomputation across =/X/N CIGAR ops.
//  6. TestQualPLULPNonImpact         — QUAL/PL last-ULP differences never flip GT/FILTER.
//
// Each test builds our binaries via pipeline/internal/upstream.OurBinary,
// resolves the upstream binary via upstream.Binary, and SKIPs gracefully when a
// prerequisite (an upstream binary, a reference, vendored data) is absent.
// Generated fixtures are tiny and written into t.TempDir().
package edgecases
