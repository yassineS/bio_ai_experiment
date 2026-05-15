package samtools

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
)

// TestCalmd_BasicMDNM exercises the four common code paths (match, mismatch,
// deletion, insertion+softclip) against a tiny hand-built reference. The
// expected MD/NM values are computed manually in the test table and verified
// against upstream `samtools calmd` semantics in bam_md.c.
func TestCalmd_BasicMDNM(t *testing.T) {
	in := openParity(t, "calmd/basic.sam")
	defer in.Close()
	refPath := parityPath(t, "calmd/ref.fa")

	var buf bytes.Buffer
	var warn bytes.Buffer
	if err := Calmd(in, &buf, refPath, CalmdOptions{}, &warn); err != nil {
		t.Fatalf("Calmd: %v", err)
	}

	// Parse the SAM text we just emitted and look up MD/NM by QNAME.
	got := indexCalmdSAM(t, buf.String())

	cases := []struct {
		qname string
		md    string
		nm    int64
	}{
		{"r_perfect", "10", 0},
		{"r_mismatch", "1C5T2", 2},
		{"r_with_del", "5^CG0T0A0C0G0T0", 7},
		{"r_with_ins", "10", 2},
		{"r_softclip", "7", 0},
		{"r_chr2", "10", 0},
		// End-of-contig: pos=48 on chr1 (length 52). CIGAR 3M5D3I:
		//   - 3M consumes ref pos 48..50 (T,A,C) → all match (NM=0).
		//   - 5D requests ref pos 51..55; only 51,52 exist (G,T) → MD `3^GT0`.
		//   - Upstream bam_md.c:121 breaks the CIGAR loop entirely when a D
		//     truncates; the trailing 3I is NEVER reached → NM stays at 2.
		// Locks in the rpos-accounting + break fix.
		{"r_end_d", "3^GT0", 2},
	}
	for _, c := range cases {
		t.Run(c.qname, func(t *testing.T) {
			rec, ok := got[c.qname]
			if !ok {
				t.Fatalf("record %s not in calmd output", c.qname)
			}
			mdA, ok := rec.GetAux("MD")
			if !ok {
				t.Fatalf("%s missing MD tag", c.qname)
			}
			md, _ := mdA.String()
			if md != c.md {
				t.Errorf("%s MD = %q, want %q", c.qname, md, c.md)
			}
			nmA, ok := rec.GetAux("NM")
			if !ok {
				t.Fatalf("%s missing NM tag", c.qname)
			}
			nm, _ := nmA.Int()
			if nm != c.nm {
				t.Errorf("%s NM = %d, want %d", c.qname, nm, c.nm)
			}
		})
	}

	// Unmapped records pass through and must NOT get an NM/MD tag.
	if rec, ok := got["r_unmapped"]; !ok {
		t.Errorf("unmapped record dropped from output")
	} else {
		if _, ok := rec.GetAux("MD"); ok {
			t.Errorf("unmapped record got an MD tag (should be unchanged)")
		}
		if _, ok := rec.GetAux("NM"); ok {
			t.Errorf("unmapped record got an NM tag (should be unchanged)")
		}
	}
}

// TestCalmd_UseEqualRewritesMatches verifies the -e flag rewrites the SEQ
// match positions to '='. NM/MD remain computed against the original bases.
func TestCalmd_UseEqualRewritesMatches(t *testing.T) {
	in := openParity(t, "calmd/basic.sam")
	defer in.Close()
	refPath := parityPath(t, "calmd/ref.fa")

	var buf bytes.Buffer
	if err := Calmd(in, &buf, refPath, CalmdOptions{UseEqual: true}, nil); err != nil {
		t.Fatalf("Calmd: %v", err)
	}
	got := indexCalmdSAM(t, buf.String())
	// r_mismatch ATGTACGAAC has mismatches at pos 1 (T vs ref C) and pos 7
	// (A vs ref T). With -e all match positions become '=', so we expect
	// "=T=====A==".
	rec := got["r_mismatch"]
	if rec == nil {
		t.Fatalf("r_mismatch missing")
	}
	if rec.Seq != "=T=====A==" {
		t.Errorf("r_mismatch -e SEQ = %q, want %q", rec.Seq, "=T=====A==")
	}
}

// TestCalmd_BAMRoundtrip writes BAM output and re-reads it, asserting the
// MD/NM tags survive the binary round-trip.
func TestCalmd_BAMRoundtrip(t *testing.T) {
	in := openParity(t, "calmd/basic.sam")
	defer in.Close()
	refPath := parityPath(t, "calmd/ref.fa")

	var buf bytes.Buffer
	if err := Calmd(in, &buf, refPath, CalmdOptions{OutputBAM: true}, nil); err != nil {
		t.Fatalf("Calmd BAM: %v", err)
	}
	// Read it back via the BAM reader.
	br, err := sam.NewReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("re-open BAM: %v", err)
	}
	hits := 0
	for {
		rec, err := br.Read()
		if err != nil {
			break
		}
		if rec.QName == "r_mismatch" {
			md, _ := rec.GetAux("MD")
			s, _ := md.String()
			if s != "1C5T2" {
				t.Errorf("BAM round-trip MD = %q, want %q", s, "1C5T2")
			}
			nm, _ := rec.GetAux("NM")
			v, _ := nm.Int()
			if v != 2 {
				t.Errorf("BAM round-trip NM = %d, want 2", v)
			}
			hits++
		}
	}
	if hits == 0 {
		t.Errorf("r_mismatch not found in BAM output")
	}
}

// TestCalmd_OverwritesDifferingTag verifies an existing MD/NM that differs
// from the recomputed value is replaced and a warning is emitted.
func TestCalmd_OverwritesDifferingTag(t *testing.T) {
	samText := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:52
r_stale	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	MD:Z:5A4	NM:i:99
`
	refPath := parityPath(t, "calmd/ref.fa")
	var buf, warn bytes.Buffer
	if err := Calmd(strings.NewReader(samText), &buf, refPath, CalmdOptions{}, &warn); err != nil {
		t.Fatalf("Calmd: %v", err)
	}
	got := indexCalmdSAM(t, buf.String())
	rec := got["r_stale"]
	if rec == nil {
		t.Fatalf("r_stale missing")
	}
	md, _ := rec.GetAux("MD")
	s, _ := md.String()
	if s != "10" {
		t.Errorf("MD = %q, want 10", s)
	}
	nm, _ := rec.GetAux("NM")
	v, _ := nm.Int()
	if v != 0 {
		t.Errorf("NM = %d, want 0", v)
	}
	if !strings.Contains(warn.String(), "different MD") {
		t.Errorf("missing 'different MD' warning, got: %q", warn.String())
	}
	if !strings.Contains(warn.String(), "different NM") {
		t.Errorf("missing 'different NM' warning, got: %q", warn.String())
	}
}

// TestCalmd_QuietSuppressesWarnings confirms Quiet=true silences the
// "different MD/NM" stderr line.
func TestCalmd_QuietSuppressesWarnings(t *testing.T) {
	samText := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:52
r_stale	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	MD:Z:5A4	NM:i:99
`
	refPath := parityPath(t, "calmd/ref.fa")
	var buf, warn bytes.Buffer
	if err := Calmd(strings.NewReader(samText), &buf, refPath, CalmdOptions{Quiet: true}, &warn); err != nil {
		t.Fatalf("Calmd: %v", err)
	}
	if warn.Len() != 0 {
		t.Errorf("Quiet emitted warnings: %q", warn.String())
	}
}

// TestParity_Calmd_SkipUpstream marks the upstream `samtools calmd -uAr
// mpileup.1.sam mpileup.ref.fa` regression case as a deferred parity test:
// upstream emits BGZF-compressed BAM. We can't byte-diff because the BGZF
// EOF block / deflate level differs by libdeflate version. Logical parity
// (MD + NM correctness on the same input) is exercised in the table tests
// above with hand-computed expectations.
func TestParity_Calmd_UpstreamCorpus(t *testing.T) {
	t.Skip("BGZF byte-identical output requires upstream's libdeflate; logical MD/NM parity covered by TestCalmd_BasicMDNM; tracked in docs/PARITY_ROADMAP.md#samtools")
}

// indexCalmdSAM parses SAM body lines into a QNAME → *sam.Record map.
// Helper used across calmd tests; header lines are silently skipped.
func indexCalmdSAM(t *testing.T, text string) map[string]*sam.Record {
	t.Helper()
	out := make(map[string]*sam.Record)
	r, err := sam.NewSAMReader(strings.NewReader(text))
	if err != nil {
		t.Fatalf("re-open SAM: %v", err)
	}
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		out[rec.QName] = rec
	}
	return out
}
