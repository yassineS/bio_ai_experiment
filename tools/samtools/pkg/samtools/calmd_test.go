package samtools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestCalmd_ThreadsDeterministic verifies that -@/--threads only parallelises
// the BGZF I/O (input inflate + compressed-BAM deflate): the decoded records
// and their recomputed MD/NM are identical regardless of the worker count, as
// in upstream bam_md.c where the fill-md compute itself is single-threaded.
// Guards the threaded input-decode and threaded BAM-output paths wired for -@.
func TestCalmd_ThreadsDeterministic(t *testing.T) {
	refPath := parityPath(t, "calmd/ref.fa")
	decode := func(threads int) string {
		in := openParity(t, "calmd/basic.sam")
		defer in.Close()
		var buf bytes.Buffer
		if err := Calmd(in, &buf, refPath, CalmdOptions{OutputBAM: true, Threads: threads}, nil); err != nil {
			t.Fatalf("Calmd -@%d: %v", threads, err)
		}
		br, err := sam.NewReader(bytes.NewReader(buf.Bytes()))
		if err != nil {
			t.Fatalf("re-open BAM -@%d: %v", threads, err)
		}
		var sb strings.Builder
		n := 0
		for {
			rec, err := br.Read()
			if err != nil {
				break
			}
			md, _ := rec.GetAux("MD")
			mds, _ := md.String()
			nm, _ := rec.GetAux("NM")
			nmv, _ := nm.Int()
			sb.WriteString(rec.QName)
			sb.WriteString("|MD=")
			sb.WriteString(mds)
			sb.WriteString("|NM=")
			sb.WriteString(strconv.Itoa(int(nmv)))
			sb.WriteByte('\n')
			n++
		}
		if n == 0 {
			t.Fatalf("no records decoded at -@%d", threads)
		}
		return sb.String()
	}
	if one, many := decode(1), decode(8); one != many {
		t.Fatalf("calmd -@1 vs -@8 decoded output differs:\n-@1:\n%s\n-@8:\n%s", one, many)
	}
}

// TestCalmdFile_ThreadsDeterministic drives the path-based CalmdFile entry
// with a real BGZF-framed BAM input, so it exercises the parallel BGZF input
// inflate that -@ >= 2 engages through CalmdFile's raw opener (OpenRaw). A
// prior wiring bug opened the input through the decompressing opener even with
// -@, so the parallel input decode never engaged; this pins the fix. The
// decoded MD/NM must be identical for -@1 and -@8.
func TestCalmdFile_ThreadsDeterministic(t *testing.T) {
	refPath := parityPath(t, "calmd/ref.fa")

	// Materialise basic.sam as a BGZF-framed BAM on disk so CalmdFile's raw
	// opener has BGZF blocks to inflate in parallel.
	src := openParity(t, "calmd/basic.sam")
	defer src.Close()
	sr, err := sam.NewSAMReader(src)
	if err != nil {
		t.Fatalf("open basic.sam: %v", err)
	}
	bamPath := filepath.Join(t.TempDir(), "basic.bam")
	bf, err := os.Create(bamPath)
	if err != nil {
		t.Fatalf("create bam: %v", err)
	}
	bw := sam.NewBAMWriter(bf)
	if err := bw.WriteHeader(sr.Header()); err != nil {
		t.Fatalf("write bam header: %v", err)
	}
	for {
		rec, rerr := sr.Read()
		if rerr != nil {
			break
		}
		if err := bw.Write(rec); err != nil {
			t.Fatalf("write bam record: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("close bam writer: %v", err)
	}
	if err := bf.Close(); err != nil {
		t.Fatalf("close bam file: %v", err)
	}

	decode := func(threads int) string {
		var buf bytes.Buffer
		if err := CalmdFile(bamPath, &buf, refPath, CalmdOptions{OutputBAM: true, Threads: threads}, nil); err != nil {
			t.Fatalf("CalmdFile -@%d: %v", threads, err)
		}
		br, err := sam.NewReader(bytes.NewReader(buf.Bytes()))
		if err != nil {
			t.Fatalf("re-open BAM -@%d: %v", threads, err)
		}
		var sb strings.Builder
		n := 0
		for {
			rec, err := br.Read()
			if err != nil {
				break
			}
			md, _ := rec.GetAux("MD")
			mds, _ := md.String()
			nm, _ := rec.GetAux("NM")
			nmv, _ := nm.Int()
			sb.WriteString(rec.QName)
			sb.WriteString("|MD=")
			sb.WriteString(mds)
			sb.WriteString("|NM=")
			sb.WriteString(strconv.Itoa(int(nmv)))
			sb.WriteByte('\n')
			n++
		}
		if n == 0 {
			t.Fatalf("no records decoded at -@%d", threads)
		}
		return sb.String()
	}
	if one, many := decode(1), decode(8); one != many {
		t.Fatalf("CalmdFile -@1 vs -@8 decoded output differs:\n-@1:\n%s\n-@8:\n%s", one, many)
	}
}

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

// TestParity_Calmd_UpstreamCorpus is the LIVE parity gate for `samtools
// calmd`. Upstream's default output is BGZF-compressed BAM, which cannot be
// byte-identical to ours without upstream's libdeflate (a forbidden dep) — so
// instead BOTH sides emit plain SAM text (upstream stays in SAM mode when
// given no -b/-u), and the streams are compared BYTE-FOR-BYTE modulo the @PG
// line. This is strictly stronger than a decoded-tag comparison: it pins the
// recomputed NM:i and MD:Z values AND their exact aux-list ordering (upstream
// bam_md.c removes a differing MD/NM and re-appends it at the end, leaving an
// unchanged tag in place — placement our port now reproduces). Per the
// project rules the upstream binary is built on demand; a build failure is
// fatal, never a skip.
func TestParity_Calmd_UpstreamCorpus(t *testing.T) {
	bin := upstreamSamtools(t)

	cases := []struct {
		name      string
		flags     []string
		sam, ref  string
		goOptions CalmdOptions
	}{
		{
			name: "basic_default",
			sam:  parityPath(t, "calmd/basic.sam"), ref: parityPath(t, "calmd/ref.fa"),
		},
		{
			name: "basic_useequal_-e", flags: []string{"-e"},
			sam: parityPath(t, "calmd/basic.sam"), ref: parityPath(t, "calmd/ref.fa"),
			goOptions: CalmdOptions{UseEqual: true},
		},
		{
			name: "realn01_default",
			sam:  parityPath(t, "calmd/realn01.sam"), ref: parityPath(t, "calmd/realn01.fa"),
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// Upstream → plain SAM on stdout.
			args := append([]string{"calmd"}, c.flags...)
			args = append(args, c.sam, c.ref)
			cmd := exec.Command(bin, args...)
			var upOut, upErr bytes.Buffer
			cmd.Stdout = &upOut
			cmd.Stderr = &upErr
			if err := cmd.Run(); err != nil {
				t.Fatalf("upstream samtools %v: %v\n%s", args, err, upErr.String())
			}

			// Our port → plain SAM into a buffer.
			f, err := os.Open(c.sam)
			if err != nil {
				t.Fatalf("open %s: %v", c.sam, err)
			}
			defer f.Close()
			var got, warn bytes.Buffer
			if err := Calmd(f, &got, c.ref, c.goOptions, &warn); err != nil {
				t.Fatalf("Calmd: %v", err)
			}

			want := dropPGLines(upOut.Bytes())
			have := dropPGLines(got.Bytes())
			if have != want {
				t.Fatalf("calmd %s differs from upstream (@PG excluded):\n--- ours ---\n%s\n--- upstream ---\n%s",
					c.name, have, want)
			}
		})
	}
}

// TestCalmd_BinQual verifies the -q flag bins base qualities: every value
// >= 3 maps to qual/10*10+7 (integer division); values below 3 are kept.
func TestCalmd_BinQual(t *testing.T) {
	// Qual string with Phred values 0,2,3,9,10,13,25,40,7,17 (ASCII '!'+v).
	// '!'=0, '#'=2, '$'=3, '*'=9, '+'=10, '.'=13, ':'=25, 'I'=40, '('=7, '2'=17.
	samText := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:52
r_q	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	!#$*+.:I(2
`
	refPath := parityPath(t, "calmd/ref.fa")
	var buf bytes.Buffer
	if err := Calmd(strings.NewReader(samText), &buf, refPath, CalmdOptions{BinQual: true}, nil); err != nil {
		t.Fatalf("Calmd: %v", err)
	}
	rec := indexCalmdSAM(t, buf.String())["r_q"]
	if rec == nil {
		t.Fatalf("r_q missing")
	}
	// Expected per upstream bam_md.c:204-208:
	//   0->0, 2->2 (both < 3, unchanged); 3->7, 9->7, 10->17, 13->17,
	//   25->27, 40->47, 7->7, 17->17.
	want := []byte{0, 2, 7, 7, 17, 17, 27, 47, 7, 17}
	if !bytes.Equal(rec.Qual, want) {
		t.Errorf("binned qual = %v, want %v", rec.Qual, want)
	}
}

// TestCalmd_DropTags verifies the -d flag keeps only the RG aux tag and, for
// a record without RG, drops every aux tag (including the freshly computed
// NM/MD, mirroring upstream's post-fill DROP_TAG ordering).
func TestCalmd_DropTags(t *testing.T) {
	samText := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:52
r_rg	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	RG:Z:grp1	XX:i:5	MD:Z:1A8	NM:i:9
r_norg	0	chr1	1	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XX:i:5	NM:i:9
`
	refPath := parityPath(t, "calmd/ref.fa")
	var buf bytes.Buffer
	if err := Calmd(strings.NewReader(samText), &buf, refPath, CalmdOptions{DropTags: true}, nil); err != nil {
		t.Fatalf("Calmd: %v", err)
	}
	got := indexCalmdSAM(t, buf.String())

	rg := got["r_rg"]
	if rg == nil {
		t.Fatalf("r_rg missing")
	}
	if len(rg.Aux) != 1 || rg.Aux[0].Tag != "RG" {
		t.Errorf("r_rg aux = %v, want only RG", rg.Aux)
	}
	if v, _ := rg.Aux[0].String(); v != "grp1" {
		t.Errorf("r_rg RG = %q, want grp1", v)
	}
	if _, ok := rg.GetAux("MD"); ok {
		t.Errorf("r_rg kept MD (should be dropped after fill)")
	}
	if _, ok := rg.GetAux("NM"); ok {
		t.Errorf("r_rg kept NM (should be dropped after fill)")
	}

	norg := got["r_norg"]
	if norg == nil {
		t.Fatalf("r_norg missing")
	}
	if len(norg.Aux) != 0 {
		t.Errorf("r_norg aux = %v, want empty (no RG to keep)", norg.Aux)
	}
}

// TestCalmd_MaxNM verifies the -n flag masks the matching bases of a
// high-NM read (SEQ -> '=', qual -> 0) while leaving a low-NM read and the
// emitted NM/MD untouched.
func TestCalmd_MaxNM(t *testing.T) {
	// chr1 ref = ACGTACGTAC...
	//   r_hi  ATGTACGAAC: mismatches at idx 1 (T/C) and idx 7 (A/T) -> NM=2.
	//   r_lo  ACGTTCGTAC: mismatch at idx 4 (T/A) -> NM=1.
	samText := `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:52
r_hi	0	chr1	1	60	10M	*	0	0	ATGTACGAAC	IIIIIIIIII
r_lo	0	chr1	1	60	10M	*	0	0	ACGTTCGTAC	IIIIIIIIII
`
	refPath := parityPath(t, "calmd/ref.fa")
	var buf bytes.Buffer
	if err := Calmd(strings.NewReader(samText), &buf, refPath, CalmdOptions{MaxNM: 2}, nil); err != nil {
		t.Fatalf("Calmd: %v", err)
	}
	got := indexCalmdSAM(t, buf.String())

	// r_hi has NM=2 >= 2: matching bases masked, mismatches (idx 1,7) kept.
	hi := got["r_hi"]
	if hi == nil {
		t.Fatalf("r_hi missing")
	}
	if hi.Seq != "=T=====A==" {
		t.Errorf("r_hi masked SEQ = %q, want %q", hi.Seq, "=T=====A==")
	}
	wantQual := []byte{0, 40, 0, 0, 0, 0, 0, 40, 0, 0}
	if !bytes.Equal(hi.Qual, wantQual) {
		t.Errorf("r_hi masked qual = %v, want %v", hi.Qual, wantQual)
	}
	if md, _ := hi.GetAux("MD"); true {
		s, _ := md.String()
		if s != "1C5T2" {
			t.Errorf("r_hi MD = %q, want 1C5T2 (unchanged by masking)", s)
		}
	}
	if nm, _ := hi.GetAux("NM"); true {
		v, _ := nm.Int()
		if v != 2 {
			t.Errorf("r_hi NM = %d, want 2 (unchanged by masking)", v)
		}
	}

	// r_lo has NM=1 < 2: untouched.
	lo := got["r_lo"]
	if lo == nil {
		t.Fatalf("r_lo missing")
	}
	if lo.Seq != "ACGTTCGTAC" {
		t.Errorf("r_lo SEQ = %q, want unchanged ACGTTCGTAC", lo.Seq)
	}
	loQual := []byte{40, 40, 40, 40, 40, 40, 40, 40, 40, 40}
	if !bytes.Equal(lo.Qual, loQual) {
		t.Errorf("r_lo qual = %v, want all 40 (unchanged)", lo.Qual)
	}
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

// extractTag re-parses SAM text and returns a qname->aux-string map for the
// named Z-type tag. Records lacking the tag are simply absent from the map.
func extractTag(t *testing.T, text, tag string) map[string]string {
	t.Helper()
	out := make(map[string]string)
	for q, rec := range indexCalmdSAM(t, text) {
		if a, ok := rec.GetAux(tag); ok {
			if s, ok := a.String(); ok {
				out[q] = s
			}
		}
	}
	return out
}

// upstreamCalmd runs the live `samtools calmd` binary on the realn01
// fixture with the given args and returns its SAM-text output. The upstream
// binary is built on demand; a build failure is fatal. This replaces
// reading committed realn01_exp*.sam golden files.
func upstreamCalmd(t *testing.T, args ...string) string {
	t.Helper()
	bin := upstreamSamtools(t)
	sam := parityPath(t, "calmd/realn01.sam")
	ref := parityPath(t, "calmd/realn01.fa")
	cmdArgs := append([]string{"calmd"}, args...)
	cmdArgs = append(cmdArgs, sam, ref)
	cmd := exec.Command(bin, cmdArgs...)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream samtools calmd %v: %v\n%s", cmdArgs, err, errBuf.String())
	}
	return out.String()
}

// TestCalmd_RealignBAQ verifies the -r flag drives BAQ realignment: each
// record in the output carries a BQ:Z tag byte-identical to the live
// upstream `samtools calmd -r` output on the same input.
func TestCalmd_RealignBAQ(t *testing.T) {
	in := openParity(t, "calmd/realn01.sam")
	defer in.Close()
	refPath := parityPath(t, "calmd/realn01.fa")

	var buf, warn bytes.Buffer
	if err := Calmd(in, &buf, refPath, CalmdOptions{RealignBAQ: true}, &warn); err != nil {
		t.Fatalf("Calmd -r: %v", err)
	}
	got := extractTag(t, buf.String(), "BQ")

	want := extractTag(t, upstreamCalmd(t, "-r"), "BQ")
	if len(want) == 0 {
		t.Fatal("upstream output has no BQ tags")
	}
	for q, w := range want {
		if got[q] != w {
			t.Errorf("read %s: BQ\n  got  %q\n  want %q", q, got[q], w)
		}
	}
}

// TestCalmd_ApplyBAQ verifies the -rA combination applies BAQ to the base
// qualities and writes a ZQ:Z tag. The adjusted qualities and ZQ tag are
// checked byte-for-byte against the live upstream `samtools calmd -rA`
// output.
func TestCalmd_ApplyBAQ(t *testing.T) {
	in := openParity(t, "calmd/realn01.sam")
	defer in.Close()
	refPath := parityPath(t, "calmd/realn01.fa")

	var buf, warn bytes.Buffer
	if err := Calmd(in, &buf, refPath, CalmdOptions{RealignBAQ: true, AdjustCapQ: true}, &warn); err != nil {
		t.Fatalf("Calmd -rA: %v", err)
	}
	gotRecs := indexCalmdSAM(t, buf.String())

	goldRecs := indexCalmdSAM(t, upstreamCalmd(t, "-r", "-A"))
	checked := 0
	for q, gold := range goldRecs {
		goldZQ, hasZQ := gold.GetAux("ZQ")
		if !hasZQ {
			continue
		}
		got, ok := gotRecs[q]
		if !ok {
			t.Errorf("read %s missing from output", q)
			continue
		}
		gotZQ, ok := got.GetAux("ZQ")
		if !ok {
			t.Errorf("read %s missing ZQ tag", q)
			continue
		}
		gz, _ := goldZQ.String()
		gotz, _ := gotZQ.String()
		if gotz != gz {
			t.Errorf("read %s: ZQ\n  got  %q\n  want %q", q, gotz, gz)
		}
		if !bytes.Equal(got.Qual, gold.Qual) {
			t.Errorf("read %s: adjusted QUAL mismatch\n  got  %v\n  want %v", q, got.Qual, gold.Qual)
		}
		checked++
	}
	if checked == 0 {
		t.Fatal("no ZQ-bearing records checked")
	}
}

// TestCalmd_CapMapQ verifies the -C flag caps MAPQ via baq.SamCapMapq. A
// 10M read with three quality-13 mismatches caps to 36 at threshold 40;
// thresholds <= 10 leave MAPQ untouched (matching bam_md.c's `capQ > 10`).
func TestCalmd_CapMapQ(t *testing.T) {
	// '.' is Phred 13 (ASCII 46-33). Reference CCCAAAAAAA gives the read
	// three quality-13 mismatches.
	sam10 := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:10\n" +
		"r\t0\tc\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t..........\n"
	refFA := ">c\nCCCAAAAAAA\n"

	dir := t.TempDir()
	refPath := filepath.Join(dir, "c.fa")
	if err := os.WriteFile(refPath, []byte(refFA), 0o644); err != nil {
		t.Fatalf("write ref: %v", err)
	}

	var buf, warn bytes.Buffer
	if err := Calmd(strings.NewReader(sam10), &buf, refPath, CalmdOptions{CapMapQ: 40}, &warn); err != nil {
		t.Fatalf("Calmd -C40: %v", err)
	}
	recs := indexCalmdSAM(t, buf.String())
	if recs["r"].MapQ != 36 {
		t.Errorf("MAPQ with -C40 = %d, want 36", recs["r"].MapQ)
	}

	// Threshold <= 10 must be ignored: MAPQ stays at the original 60.
	var buf2 bytes.Buffer
	if err := Calmd(strings.NewReader(sam10), &buf2, refPath, CalmdOptions{CapMapQ: 5}, &warn); err != nil {
		t.Fatalf("Calmd -C5: %v", err)
	}
	if r := indexCalmdSAM(t, buf2.String())["r"]; r.MapQ != 60 {
		t.Errorf("MAPQ with -C5 = %d, want 60 (unchanged)", r.MapQ)
	}
}
