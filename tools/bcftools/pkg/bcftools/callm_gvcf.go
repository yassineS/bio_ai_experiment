// gVCF block emission for `bcftools call -m --gvcf`. Wires the
// post-mcall record stream through the upstream gvcf.c block-flush
// algorithm — the same logic mpileup_gvcf.go ports for `bcftools
// mpileup --gvcf`, with the adaptations required by call's
// post-mcall record shape:
//
//   - Ref-only sites are recognised by ALT in {empty, "."}, not the
//     mpileup-side <*> placeholder.
//   - The emitted block carries no PL (mcall removed PL via
//     bcf_update_format_int32(... "PL" NULL 0) — mcall.c:1580) and no
//     QS (mcall.c:1529 removes QS), so the block's per-sample data
//     reduces to GT:DP. INFO is END (only when the block spans at
//     least 2 sites) + MIN_DP.
//   - The block's allele set is REF only: ALT renders as ".".
//
// Upstream call dispatch (vcfcall.c:1249-1254):
//
//	if ( (args.aux.flag & CALL_VARONLY) && ret==0 && !args.gvcf ) continue;
//	if ( args.gvcf )
//	    bcf_rec = gvcf_write(args.gvcf, args.out_fh, args.aux.hdr, bcf_rec, ret==1?1:0);
//	if ( bcf_rec && bcf_write1(args.out_fh, args.aux.hdr, bcf_rec)!=0 ) ...
//	at EOF: if ( args.gvcf ) gvcf_write(args.gvcf, args.out_fh, args.aux.hdr, NULL, 0);
//
// In our streaming pipeline the equivalent dispatch is baked into
// callGVCFBlocker.Write/Flush — the inner writer sees a sequence of
// records that already includes the collapsed block rows.

package bcftools

import (
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// augmentGVCFHeader appends the END / MIN_DP INFO declarations
// produced by upstream gvcf_update_header (gvcf.c:40-44). The lines
// are appended idempotently — if the user's input header already
// declares them (rare but possible), the originals are kept.
func augmentGVCFHeader(hdr *vcf.Header) *vcf.Header {
	if hdr == nil {
		return hdr
	}
	out := &vcf.Header{
		MetaInfo: append([]string(nil), hdr.MetaInfo...),
		Samples:  append([]string(nil), hdr.Samples...),
	}
	declarations := []struct {
		marker string
		line   string
	}{
		{`##INFO=<ID=END,`, `##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">`},
		{`##INFO=<ID=MIN_DP,`, `##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum per-sample depth in this gVCF block">`},
	}
	for _, d := range declarations {
		found := false
		for _, m := range out.MetaInfo {
			if strings.HasPrefix(m, d.marker) {
				found = true
				break
			}
		}
		if !found {
			out.MetaInfo = append(out.MetaInfo, d.line)
		}
	}
	return out
}

// callIsRefOnly reports whether v is a post-mcall ref-only record.
// Such records carry ALT=[] or ALT=["."] (mcall.c:1577-1581 — the
// als_new==1 branch that runs bcf_update_format_int32(... "PL" NULL 0)
// and leaves only REF on the record). The mpileup-side isRefOnly
// recognises the <*> placeholder instead.
func callIsRefOnly(v *vcf.Variant) bool {
	if len(v.Alt) == 0 {
		return true
	}
	if len(v.Alt) == 1 && (v.Alt[0] == "." || v.Alt[0] == "") {
		return true
	}
	return false
}

// callGVCFBlocker wraps an underlying variantWriter and bands
// consecutive post-mcall REF-only rows into INFO/END+MIN_DP blocks
// per the upstream gvcf.c::gvcf_write state machine.
type callGVCFBlocker struct {
	inner   variantWriter
	dpRange []int

	// Active block state. valid == false means no block in flight.
	valid     bool
	rid       string // chrom
	start     int    // 1-based start
	end       int    // 1-based last collapsed pos (inclusive)
	prevRange int    // dp_range bin (>=1) of the active block

	// Carried from the block's first record so the emitted block has
	// the right REF base and per-sample column names.
	refAllele string
	nSamples  int
	sampleNm  []string

	// Per-sample minima across the block. minDP[i] is min DP for
	// sample i. blockMinDP is the global min across all samples and
	// all sites — the value emitted as INFO/MIN_DP.
	minDP      []int
	blockMinDP int
}

// newCallGVCFBlocker returns a wrapper over inner that bands
// post-mcall REF-only records by --gvcf thresholds. dpRange must be
// non-empty and sorted (parseGVCFRanges already sorts).
func newCallGVCFBlocker(inner variantWriter, dpRange []int) *callGVCFBlocker {
	return &callGVCFBlocker{inner: inner, dpRange: dpRange}
}

// WriteHeader forwards to the inner writer; the header rewrite was
// already applied in callStreaming via augmentGVCFHeader.
func (g *callGVCFBlocker) WriteHeader() error { return g.inner.WriteHeader() }

// Flush emits any in-flight block then forwards to the inner writer.
// Mirror of upstream's "EOF: gvcf_write(..., NULL, 0)" from
// vcfcall.c:1254.
func (g *callGVCFBlocker) Flush() error {
	if g.valid {
		if err := g.emitBlock(); err != nil {
			return err
		}
	}
	return g.inner.Flush()
}

// Write is the per-record dispatcher. It mirrors gvcf_write
// (gvcf.c:88) including the SNP-row → indel-row position de-dup and
// the dp_range==0 pass-through (upstream sets can_collapse=0 and
// attaches MIN_DP to the record on output).
func (g *callGVCFBlocker) Write(v *vcf.Variant) error {
	canCollapse := callIsRefOnly(v)
	var minDP int
	var dpBin int
	var dpVec []int

	if canCollapse {
		var ok bool
		dpVec, ok = perSampleDP(v)
		if !ok {
			// Per-sample DP missing — upstream's gvcf.c:127 treats
			// this as a block-breaker.
			canCollapse = false
		} else {
			minDP = dpVec[0]
			for _, d := range dpVec[1:] {
				if d < minDP {
					minDP = d
				}
			}
			dpBin = bandFor(minDP, g.dpRange)
			if dpBin == 0 {
				// "DP too small" — upstream sets can_collapse=0 and
				// needs_flush=1, then attaches MIN_DP to the
				// passed-through record below.
				canCollapse = false
			}
		}
	}

	needsFlush := !canCollapse
	if g.valid && g.prevRange != dpBin {
		needsFlush = true
	}
	if g.valid && (g.rid != v.Chrom || v.Pos > g.end+1) {
		needsFlush = true
	}

	if g.valid && needsFlush {
		// SNP + indel at the same position: drop the trailing SNP
		// position from the block end (gvcf.c:139).
		if v.Chrom == g.rid && v.Pos == g.end {
			g.end--
		}
		if err := g.emitBlock(); err != nil {
			return err
		}
	}

	if canCollapse {
		if !g.valid {
			g.valid = true
			// v.Chrom / v.Ref / sample names may alias the reader's reused
			// line buffer (records sourced via ReadInto are consume-and-discard:
			// their strings are only valid until the next read). The block is
			// held in flight across many subsequent reads, so clone the strings
			// it retains to give the block its own stable backing.
			g.rid = strings.Clone(v.Chrom)
			g.start = v.Pos
			g.end = v.Pos
			g.prevRange = dpBin
			g.refAllele = strings.Clone(v.Ref)
			g.nSamples = len(v.Samples)
			g.sampleNm = make([]string, g.nSamples)
			for i, s := range v.Samples {
				g.sampleNm[i] = strings.Clone(s.Name)
			}
			g.minDP = append([]int{}, dpVec...)
			g.blockMinDP = minDP
		} else {
			g.end = v.Pos
			if minDP < g.blockMinDP {
				g.blockMinDP = minDP
			}
			for i, d := range dpVec {
				if d < g.minDP[i] {
					g.minDP[i] = d
				}
			}
		}
		return nil
	}

	// Non-collapsible record: pass through unchanged. Upstream
	// (gvcf.c:221-222) injects MIN_DP for the is_ref-but-DP-too-low
	// case.
	if dpBin == 0 && callIsRefOnly(v) && len(dpVec) > 0 {
		if v.Info == nil {
			v.Info = map[string]string{}
		}
		if _, present := v.Info["MIN_DP"]; !present {
			v.Info["MIN_DP"] = strconv.Itoa(minDP)
			v.InfoOrder = append(v.InfoOrder, "MIN_DP")
		}
	}
	return g.inner.Write(v)
}

// emitBlock builds and writes the collapsed gVCF record for the
// currently-buffered block, then resets the block state. The record
// shape mirrors gvcf.c:143-158: REF allele only (ALT="."), INFO=END
// (only when start<end) + MIN_DP, FORMAT=GT:DP per sample, GT 0/0,
// DP=per-sample-min.
func (g *callGVCFBlocker) emitBlock() error {
	rec := &vcf.Variant{
		Chrom:   g.rid,
		Pos:     g.start,
		ID:      ".",
		Ref:     g.refAllele,
		Alt:     []string{"."},
		Qual:    -1, // upstream emits "." for QUAL on gVCF blocks
		Filter:  nil,
		Info:    map[string]string{},
		Format:  []string{"GT", "DP"},
		Samples: make([]vcf.Sample, g.nSamples),
	}
	if g.start < g.end {
		rec.Info["END"] = strconv.Itoa(g.end)
		rec.InfoOrder = append(rec.InfoOrder, "END")
	}
	rec.Info["MIN_DP"] = strconv.Itoa(g.blockMinDP)
	rec.InfoOrder = append(rec.InfoOrder, "MIN_DP")
	for i := 0; i < g.nSamples; i++ {
		// Post-mcall ref-only records always carry GT=0/0 (mcall.c
		// mcall_set_ref_genotypes runs for the als_new==1 branch).
		// Block emits the same flat 0/0 for every sample.
		data := map[string]string{
			"GT": "0/0",
			"DP": strconv.Itoa(g.minDP[i]),
		}
		rec.Samples[i] = vcf.Sample{Name: g.sampleNm[i], Data: data}
	}
	if err := g.inner.Write(rec); err != nil {
		return err
	}
	g.valid = false
	g.rid = ""
	g.start = 0
	g.end = 0
	g.prevRange = 0
	g.minDP = nil
	g.sampleNm = nil
	g.nSamples = 0
	g.blockMinDP = 0
	return nil
}
