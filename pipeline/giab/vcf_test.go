package giab

import (
	"strings"
	"testing"
)

func TestParseVCF_Fields(t *testing.T) {
	recs := mustParse(t, vcf(rec("chr1", 12345, "A", "GT", "67.5", "PASS", "GT:PL:DP", "0/1:67,0,255:30")))
	if len(recs) != 1 {
		t.Fatalf("expected 1 record, got %d", len(recs))
	}
	r := recs[0]
	if r.Chrom != "chr1" || r.Pos != 12345 || r.Ref != "A" || r.Alt != "GT" {
		t.Fatalf("bad locus: %+v", r)
	}
	if r.Qual != "67.5" || r.Filter != "PASS" {
		t.Fatalf("bad qual/filter: %+v", r)
	}
	if r.GT() != "0/1" {
		t.Fatalf("GT(): %q", r.GT())
	}
	if r.PL() != "67,0,255" {
		t.Fatalf("PL(): %q", r.PL())
	}
	if r.gtField("DP") != "30" {
		t.Fatalf("DP subfield: %q", r.gtField("DP"))
	}
	if r.IsSNV() {
		t.Fatal("A->GT is an indel, not SNV")
	}
	if r.Key() != "chr1:12345:A:GT" {
		t.Fatalf("Key(): %q", r.Key())
	}
}

func TestVCFRecord_IsSNV(t *testing.T) {
	snv := mustParse(t, vcf(rec("chr1", 1, "A", "G", ".", ".", "GT", "0/1")))[0]
	if !snv.IsSNV() {
		t.Fatal("A->G should be SNV")
	}
	multi := mustParse(t, vcf(rec("chr1", 1, "A", "G,T", ".", ".", "GT", "1/2")))[0]
	if !multi.IsSNV() {
		t.Fatal("A->G,T (both single base) should be SNV")
	}
	star := mustParse(t, vcf(rec("chr1", 1, "A", "*", ".", ".", "GT", "0/1")))[0]
	if star.IsSNV() {
		t.Fatal("spanning-deletion * is not an SNV")
	}
}

func TestParseVCF_SkipsHeaderAndBlank(t *testing.T) {
	in := "##fileformat=VCFv4.2\n#CHROM\tPOS\n\nchr1\t1\t.\tA\tG\t.\t.\t.\n"
	recs, err := ParseVCF(strings.NewReader(in))
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) != 1 {
		t.Fatalf("expected 1 data record, got %d", len(recs))
	}
}

func TestParseVCF_BadPos(t *testing.T) {
	in := "chr1\tNOTANUMBER\t.\tA\tG\t.\t.\t.\n"
	if _, err := ParseVCF(strings.NewReader(in)); err == nil {
		t.Fatal("expected error on non-numeric POS")
	}
}

func TestResolveEngine_NoneAvailable(t *testing.T) {
	// With an empty config and (assumed) no hap.py/rtg on the test PATH, the
	// engine resolves to empty. We cannot control PATH portably, so only assert
	// that an explicit config wins.
	r := &runner{cfg: &Config{HappyBin: "/usr/bin/hap.py"}}
	e := r.resolveEngine()
	if e.engine != EngineHappy || e.bin != "/usr/bin/hap.py" {
		t.Fatalf("explicit happy_bin should win: %+v", e)
	}
	r2 := &runner{cfg: &Config{VcfevalBin: "/opt/rtg/vcfeval"}}
	e2 := r2.resolveEngine()
	if e2.engine != EngineVcfeval || e2.bin != "/opt/rtg/vcfeval" {
		t.Fatalf("explicit vcfeval_bin should win: %+v", e2)
	}
}

func TestSanitize(t *testing.T) {
	if sanitize("*") != "all" {
		t.Fatalf("star should sanitize to all")
	}
	if sanitize("low/map region") != "low_map_region" {
		t.Fatalf("sanitize: %q", sanitize("low/map region"))
	}
}
