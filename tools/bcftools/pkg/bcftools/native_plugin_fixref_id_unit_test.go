// Binary-free unit tests for the fixref id-mode pure helpers
// (native_plugin_fixref_id.go): the dbSNP map builder's filtering/first-wins
// rules (dbsnpAddRecord), the orientation decision and GT swap (dbsnpCheck), and
// the unresolved (no-info) branch. None of these need a FASTA, a dbSNP file, or
// the upstream binary: the dbSNP map is supplied directly and every case keeps
// the dbSNP position equal to the record position so fetchRefBase is never
// consulted.
package bcftools

import (
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// TestUnitFixrefDbsnpAddRecord pins the per-record filtering and first-wins
// rules of dbsnpAddRecord, mirroring fixref.c dbsnp_init.
func TestUnitFixrefDbsnpAddRecord(t *testing.T) {
	rec := func(id, ref string, alt ...string) *vcf.Variant {
		return &vcf.Variant{Chrom: "c1", Pos: 7, ID: id, Ref: ref, Alt: alt}
	}
	tests := []struct {
		name string
		recs []*vcf.Variant
		want map[string]fixrefMarker
	}{
		{
			name: "single SNP recorded with 0-based pos and ref index",
			recs: []*vcf.Variant{rec("rs1", "C", "A")},
			want: map[string]fixrefMarker{"rs1": {pos: 6, ref: 1}}, // C=1
		},
		{
			name: "non-SNP REF (indel) skipped",
			recs: []*vcf.Variant{rec("rs1", "CA", "C")},
			want: map[string]fixrefMarker{},
		},
		{
			name: "non-SNP ALT (indel) skipped",
			recs: []*vcf.Variant{rec("rs1", "C", "CA")},
			want: map[string]fixrefMarker{},
		},
		{
			name: "no ALT skipped",
			recs: []*vcf.Variant{rec("rs1", "C")},
			want: map[string]fixrefMarker{},
		},
		{
			name: "non-ACGT REF skipped",
			recs: []*vcf.Variant{rec("rs1", "N", "A")},
			want: map[string]fixrefMarker{},
		},
		{
			name: "missing ID skipped",
			recs: []*vcf.Variant{rec(".", "C", "A"), rec("", "G", "T")},
			want: map[string]fixrefMarker{},
		},
		{
			name: "duplicate ID keeps first occurrence",
			recs: []*vcf.Variant{
				{Chrom: "c1", Pos: 5, ID: "rs9", Ref: "A", Alt: []string{"G"}},
				{Chrom: "c1", Pos: 9, ID: "rs9", Ref: "T", Alt: []string{"C"}},
			},
			want: map[string]fixrefMarker{"rs9": {pos: 4, ref: 0}}, // A=0, first wins
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			m := make(map[string]fixrefMarker)
			for _, r := range tt.recs {
				dbsnpAddRecord(m, r)
			}
			if !reflect.DeepEqual(m, tt.want) {
				t.Fatalf("map = %v, want %v", m, tt.want)
			}
		})
	}
}

// TestUnitFixrefDbsnpCheckOrientation pins the orientation decision and GT swap
// of dbsnpCheck (and applyUseID via the pre-built map), mirroring fixref.c
// dbsnp_check. The dbSNP marker shares the record position so no FASTA fetch is
// triggered.
func TestUnitFixrefDbsnpCheckOrientation(t *testing.T) {
	tests := []struct {
		name      string
		marker    fixrefMarker
		id        string
		ref, alt  string
		gt        string // S1 genotype
		ia, ib    int
		ir        int
		wantKeep  bool
		wantRef   string
		wantAlt   string
		wantGT    string
		wantDirty uint32
		wantSwap  uint32
		wantUnres uint32
		discard   bool
	}{
		{
			name:   "REF already matches -> none, unchanged",
			marker: fixrefMarker{pos: 9, ref: 0}, id: "rs1",
			ref: "A", alt: "G", gt: "0/1", ia: 0, ib: 2, ir: 0,
			wantKeep: true, wantRef: "A", wantAlt: "G", wantGT: "0/1",
			wantDirty: fixFixNone,
		},
		{
			// dbSNP REF (C=1) equals the fetched forward ref ir, and equals the
			// input ALT -> swap + GT flip.
			name:   "ALT matches REF -> swap + GT flip",
			marker: fixrefMarker{pos: 9, ref: 1}, id: "rs2",
			ref: "A", alt: "C", gt: "0/1", ia: 0, ib: 1, ir: 1,
			wantKeep: true, wantRef: "C", wantAlt: "A", wantGT: "1/0",
			wantDirty: fixFixSwap, wantSwap: 1,
		},
		{
			// dbSNP REF (A=0) equals ir and the input ALT -> swap; phase kept.
			name:   "ALT matches REF, phased GT preserved",
			marker: fixrefMarker{pos: 9, ref: 0}, id: "rs3",
			ref: "C", alt: "A", gt: "0|1", ia: 1, ib: 0, ir: 0,
			wantKeep: true, wantRef: "A", wantAlt: "C", wantGT: "1|0",
			wantDirty: fixFixSwap, wantSwap: 1,
		},
		{
			name:   "unknown ID -> unresolved (skip), kept",
			marker: fixrefMarker{pos: 9, ref: 0}, id: "rs_other",
			ref: "A", alt: "G", gt: "0/1", ia: 0, ib: 2, ir: 0,
			wantKeep: true, wantRef: "A", wantAlt: "G", wantGT: "0/1",
			wantDirty: fixFixSkip, wantUnres: 1,
		},
		{
			name:   "ID=. -> unresolved (skip), kept",
			marker: fixrefMarker{pos: 9, ref: 0}, id: ".",
			ref: "A", alt: "G", gt: "0/1", ia: 0, ib: 2, ir: 0,
			wantKeep: true, wantRef: "A", wantAlt: "G", wantGT: "0/1",
			wantDirty: fixFixSkip, wantUnres: 1,
		},
		{
			name:   "unknown ID + discard -> dropped",
			marker: fixrefMarker{pos: 9, ref: 0}, id: "rs_other",
			ref: "A", alt: "G", gt: "0/1", ia: 0, ib: 2, ir: 0,
			wantKeep: false, wantUnres: 1, wantDirty: fixFixSkip, discard: true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			p := &fixrefPlugin{discard: tt.discard, dirty: fixFixSkip}
			// The map only ever contains the marker under "rs1/rs2/rs3"; the
			// unknown-ID and ID=. cases purposely miss it.
			p.dbsnpMap = map[string]fixrefMarker{}
			if tt.id == "rs1" || tt.id == "rs2" || tt.id == "rs3" {
				p.dbsnpMap[tt.id] = tt.marker
			}
			v := &vcf.Variant{
				Chrom: "c1", Pos: 10, ID: tt.id, Ref: tt.ref, Alt: []string{tt.alt},
				Format:  []string{"GT"},
				Samples: []vcf.Sample{{Name: "S1", Data: map[string]string{"GT": tt.gt}}},
			}
			keep, err := p.dbsnpCheck(v, tt.ir, tt.ia, tt.ib)
			if err != nil {
				t.Fatalf("dbsnpCheck error: %v", err)
			}
			if keep != tt.wantKeep {
				t.Fatalf("keep = %v, want %v", keep, tt.wantKeep)
			}
			if tt.wantKeep {
				if v.Ref != tt.wantRef {
					t.Errorf("REF = %q, want %q", v.Ref, tt.wantRef)
				}
				if v.Alt[0] != tt.wantAlt {
					t.Errorf("ALT = %q, want %q", v.Alt[0], tt.wantAlt)
				}
				if got := v.Samples[0].Data["GT"]; got != tt.wantGT {
					t.Errorf("GT = %q, want %q", got, tt.wantGT)
				}
			}
			if p.dirty != tt.wantDirty {
				t.Errorf("dirty = %d, want %d", p.dirty, tt.wantDirty)
			}
			if p.nswap != tt.wantSwap {
				t.Errorf("nswap = %d, want %d", p.nswap, tt.wantSwap)
			}
			if p.nunresolved != tt.wantUnres {
				t.Errorf("nunresolved = %d, want %d", p.nunresolved, tt.wantUnres)
			}
		})
	}
}

// TestUnitFixrefDbsnpRefMismatch pins the fatal "Reference base mismatch" path:
// after a position correction (marker.pos != record pos) the dbSNP REF must
// equal the fetched forward REF, else dbsnpCheck returns an error. Here the
// position already matches, so we instead verify the in-position mismatch is
// impossible to reach without a fetch; the position-correction + mismatch path
// is covered by the live oracle. This test asserts that with no position
// correction and a consistent ir, no error is produced.
func TestUnitFixrefDbsnpNoFetchWhenAligned(t *testing.T) {
	p := &fixrefPlugin{dirty: fixFixSkip}
	p.dbsnpMap = map[string]fixrefMarker{"rs1": {pos: 9, ref: 0}}
	v := &vcf.Variant{Chrom: "c1", Pos: 10, ID: "rs1", Ref: "A", Alt: []string{"G"}}
	// ir must be passed as the aligned ref (A=0); fetchRefBase (nil fai) is not
	// called because marker.pos (9) == v.Pos-1 (9).
	keep, err := p.dbsnpCheck(v, 0, 0, 2)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if !keep || p.dirty != fixFixNone {
		t.Fatalf("keep=%v dirty=%d, want keep=true dirty=none", keep, p.dirty)
	}
}
