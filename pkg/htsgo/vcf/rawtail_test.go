package vcf

import (
	"bytes"
	"strings"
	"testing"
)

// TestKeepRawSamplesByteIdentical pins the shallow-sample property: reading with
// KeepRawSamples on (FORMAT + sample columns kept verbatim as RawTail) and
// writing back is byte-identical to a normal parse-into-maps + re-serialise, for
// well-formed records. isec relies on this to skip the per-sample map round-trip.
func TestKeepRawSamplesByteIdentical(t *testing.T) {
	const vcf = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S0	S1	S2
chr1	100	rs1	A	T	60	PASS	DP=20	GT:AD:PL	0/1:12,9:30,0,40	1/1:0,22:99,40,0	0/0:18,0:0,30,99
chr1	200	.	C	G,T	.	q10	DP=5	GT:AD	1/2:0,3,4	0/1:5,2,0	./.:.,.,.
chr1	300	.	G	A	40	PASS	.	GT	0|1	1|1	0/0
`
	// Header-only with no samples (sites-only-ish) to confirm RawTail stays "".
	const sitesOnly = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	60	PASS	DP=20
`
	roundtrip := func(src string, shallow bool) string {
		r := NewReader(strings.NewReader(src))
		hdr, err := r.ReadHeader()
		if err != nil {
			t.Fatalf("ReadHeader: %v", err)
		}
		r.KeepRawSamples(shallow)
		var out bytes.Buffer
		w := NewWriter(&out, hdr)
		if err := w.WriteHeader(); err != nil {
			t.Fatalf("WriteHeader: %v", err)
		}
		for {
			v, err := r.Read()
			if err != nil {
				break
			}
			if shallow && len(hdr.Samples) > 0 && v.RawTail == "" {
				t.Errorf("KeepRawSamples on: record %s:%d has empty RawTail", v.Chrom, v.Pos)
			}
			if err := w.Write(v); err != nil {
				t.Fatalf("Write: %v", err)
			}
		}
		if err := w.Flush(); err != nil {
			t.Fatalf("Flush: %v", err)
		}
		return out.String()
	}

	if got, want := roundtrip(vcf, true), roundtrip(vcf, false); got != want {
		t.Errorf("shallow != normal round-trip:\n--- shallow ---\n%s\n--- normal ---\n%s", got, want)
	}
	// Sites-only records must round-trip identically and never set RawTail.
	if got, want := roundtrip(sitesOnly, true), roundtrip(sitesOnly, false); got != want {
		t.Errorf("sites-only shallow != normal:\n%s\nvs\n%s", got, want)
	}
}
