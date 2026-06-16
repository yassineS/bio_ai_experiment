// Native port of the `id`/`--use-id` mode of the upstream fixref plugin
// (plugins/fixref.c, MODE_USE_ID). This mode determines the correct REF allele
// from a separate dbSNP VCF/BCF keyed by the ID (rsID) column rather than from
// strand convention: for each input record it looks the rsID up in the dbSNP
// file and swaps REF/ALT (and the genotypes) when the input ALT matches the
// dbSNP REF base.
//
// Upstream consults the dbSNP file through a synced BCF reader restricted to the
// current record's chromosome (bcf_sr_set_regions(sr, chr, 0)); it builds a
// per-chromosome hash map keyed by the dbSNP ID string -> {pos, ref-base},
// rebuilt whenever the input chromosome changes. We reproduce the exact same
// map by streaming the dbSNP VCF once per chromosome and keeping only records on
// that chromosome. Because the result is fully determined by that map (an
// ID->{pos,ref} lookup), streaming with iohelper's transparent BGZF/gzip
// decoding yields byte-for-byte identical output without needing the .tbi/.csi
// index that upstream's synced reader requires.
package bcftools

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// fixrefMarker mirrors the C marker_t: the 0-based position and the 0/1/2/3
// reference-base index recorded for a dbSNP rsID.
type fixrefMarker struct {
	pos int // 0-based position of the dbSNP record
	ref int // 0/1/2/3 ref-base index (A/C/G/T)
}

// dbsnpBuildMap streams the dbSNP VCF and returns an rsID -> marker map holding
// only the records on chromosome chr. It mirrors fixref.c dbsnp_init: non-SNP
// (REF or ALT not a single base), non-[ACGT] REF, and missing/"." IDs are
// skipped, and a duplicate ID keeps the FIRST occurrence (the ambiguous-id
// skip). It is a pure function of the dbSNP file and the chromosome name.
func dbsnpBuildMap(dbsnpFname, chr string) (map[string]fixrefMarker, error) {
	rc, err := iohelper.OpenReader(dbsnpFname)
	if err != nil {
		return nil, fmt.Errorf("fixref: failed to open %s: %w", dbsnpFname, err)
	}
	defer rc.Close()

	rd := vcf.NewReader(rc)
	if _, err := rd.ReadHeader(); err != nil {
		return nil, fmt.Errorf("fixref: failed to read header of %s: %w", dbsnpFname, err)
	}

	m := make(map[string]fixrefMarker)
	for {
		rec, err := rd.Read()
		if err != nil || rec == nil {
			break
		}
		if rec.Chrom != chr {
			continue // synced reader is region-restricted to chr
		}
		dbsnpAddRecord(m, rec)
	}
	return m, nil
}

// dbsnpAddRecord inserts one dbSNP record into the rsID map following the
// upstream filtering and first-wins rules. It is split out as a pure helper for
// unit testing without any I/O.
func dbsnpAddRecord(m map[string]fixrefMarker, rec *vcf.Variant) {
	// Skip non-SNPs: upstream requires both REF and ALT to be a single base
	// (allele[0][1]==0 && allele[1][1]==0). A record with no ALT is skipped too.
	if len(rec.Ref) != 1 {
		return
	}
	if len(rec.Alt) < 1 || len(rec.Alt[0]) != 1 {
		return
	}
	ref := fixrefNt2int(rec.Ref[0])
	if ref < 0 {
		return // non-[ACGT] base
	}
	if rec.ID == "" || rec.ID == "." {
		return
	}
	if _, ok := m[rec.ID]; ok {
		return // ambiguous id: keep the first occurrence
	}
	m[rec.ID] = fixrefMarker{pos: rec.Pos - 1, ref: ref}
}

// applyUseID is MODE_USE_ID: it looks the record's ID up in the per-chromosome
// dbSNP map (rebuilding the map on a chromosome change), decides the orientation
// from the dbSNP REF base, and swaps REF/ALT (and the genotypes) when needed.
// It returns whether to keep the record. ia/ib are the 0/1/2/3 indices of the
// input REF/ALT bases; ir is the fetched forward reference-base index.
func (p *fixrefPlugin) applyUseID(v *vcf.Variant, ir, ia, ib int) (bool, error) {
	if p.dbsnpMap == nil || p.dbsnpRID != v.Chrom {
		p.dbsnpRID = v.Chrom
		p.dbsnpPrevPos = 0 // mirrors `args.pos = 0` on a chromosome change
		m, err := dbsnpBuildMap(p.dbsnpFname, v.Chrom)
		if err != nil {
			return false, err
		}
		p.dbsnpMap = m
	}

	keep, err := p.dbsnpCheck(v, ir, ia, ib)
	if err != nil {
		return false, err
	}

	// Upstream emits a one-shot warning when a corrected position makes the
	// output unsorted (the previous output position is greater than this one).
	if !p.dbsnpUnsorted && p.dbsnpPrevPos > v.Pos-1 {
		if p.stderr != nil {
			fmt.Fprintf(p.stderr,
				"Warning: corrected position(s) results in unsorted VCF, for example %s:%d comes after %s:%d\n"+
					"         The command `bcftools sort` can be used to fix the order.\n",
				v.Chrom, v.Pos, v.Chrom, p.dbsnpPrevPos+1)
		}
		p.dbsnpUnsorted = true
	}
	p.dbsnpPrevPos = v.Pos - 1
	return keep, nil
}

// dbsnpCheck mirrors fixref.c dbsnp_check. On an ID hit it corrects the position
// if the dbSNP record sits elsewhere (re-fetching the forward REF and counting a
// fixed pos), verifies the dbSNP REF base matches the fetched REF (a fatal error
// otherwise, exactly like upstream), and then: leaves the record unchanged when
// the input REF already equals the dbSNP REF (FIX_NONE), or swaps REF/ALT and
// the genotypes when the input ALT equals the dbSNP REF (FIX_SWAP). A missing
// ID, an unknown ID, or neither allele matching is left unresolved (annotated
// "skip", counted as unresolved, and dropped when --discard is set).
func (p *fixrefPlugin) dbsnpCheck(v *vcf.Variant, ir, ia, ib int) (bool, error) {
	if v.ID == "" || v.ID == "." {
		return p.dbsnpNoInfo(), nil
	}
	mk, ok := p.dbsnpMap[v.ID]
	if !ok {
		return p.dbsnpNoInfo(), nil
	}

	if mk.pos != v.Pos-1 {
		v.Pos = mk.pos + 1
		ir = p.fetchRefBase(v)
		if ir == -3 {
			return false, fmt.Errorf("fixref: faidx fetch failed at %s:%d", v.Chrom, v.Pos)
		}
		p.npesErr++
	}

	if mk.ref != ir {
		return false, fmt.Errorf("fixref: Reference base mismatch at %s:%d .. %c vs %c",
			v.Chrom, v.Pos, fixrefInt2nt(mk.ref), fixrefInt2nt(ir))
	}

	if ia == mk.ref {
		p.dirty = fixFixNone
		return true, nil
	}
	if ib == mk.ref {
		p.dirty = fixFixSwap
		p.nswap++
		p.setRefAlt(v, fixrefInt2nt(ib), fixrefInt2nt(ia), true)
		return true, nil
	}
	return p.dbsnpNoInfo(), nil
}

// dbsnpNoInfo handles the unresolved (no-info) branch shared by missing/unknown
// IDs and the "neither allele matches" case: it counts an unresolved site and
// keeps the record (annotated "skip", since dirty stays FIX_SKIP) unless
// --discard is set, matching the C goto no_info path.
func (p *fixrefPlugin) dbsnpNoInfo() bool {
	p.nunresolved++
	return !p.discard
}
