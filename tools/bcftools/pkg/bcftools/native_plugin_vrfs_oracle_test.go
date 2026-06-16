package bcftools

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// Live-oracle parity tests for the `+vrfs` native plugin. They build small
// deterministic BAM fixtures with the upstream samtools, then drive BOTH the
// genuine upstream bcftools (BCFTOOLS_PLUGINS -> vendored vrfs.so) and OUR port
// over the same inputs, diffing stdout byte-for-byte.
//
// The vrfs pileup runs in mpileup2 LEGACY_MODE, whose realignment step is a
// stub (no BAQ), so pure-SNV pileups — and even the indel cases here, since the
// counts come straight from the unrealigned CIGAR — are byte-exact. We
// therefore assert exact stdout equality (assertStdoutParity) rather than the
// proximity comparison; see docs/UPSTREAM_BUGS.md for the parity boundary.

// upstreamSamtoolsForVrfs resolves the vendored samtools binary used to build
// BAM fixtures. It t.Skip's (not Fatal) when absent so the test degrades
// gracefully on a checkout without the samtools submodule binary, while the
// binary-free TestUnitVrfs* tests still cover the pure logic.
func upstreamSamtoolsForVrfs(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code", "samtools", "samtools"))
	if err != nil {
		t.Skipf("samtools abs: %v", err)
	}
	if fi, err := os.Stat(abs); err != nil || fi.IsDir() {
		t.Skipf("upstream samtools binary not available at %s; build it to enable the vrfs oracle", abs)
	}
	return abs
}

// vrfsRunSamtools runs the upstream samtools with argv, optionally piping
// stdinData, returning stdout. It fails the test on a non-zero exit.
func vrfsRunSamtools(t *testing.T, sam string, stdin []byte, argv ...string) []byte {
	t.Helper()
	cmd := exec.Command(sam, argv...)
	if stdin != nil {
		cmd.Stdin = bytes.NewReader(stdin)
	}
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("samtools %v: %v\n%s", argv, err, errBuf.String())
	}
	return out.Bytes()
}

// vrfsBuildFixtures writes a tiny reference, builds the BAMs from SAM text, and
// writes the sites / aln-list files into dir. It returns the directory paths
// the vrfs invocations need.
type vrfsFixtures struct {
	dir       string
	ref       string
	alnList   string
	sites     string
	indelAln  string
	indelSite string
}

func vrfsBuildFixtures(t *testing.T, sam string) vrfsFixtures {
	t.Helper()
	dir := t.TempDir()
	write := func(name, content string) string {
		p := filepath.Join(dir, name)
		if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
			t.Fatalf("write %s: %v", name, err)
		}
		return p
	}

	// Reference: chr1 = ACGT repeated to 60 bp. 1-based pos11=G, pos21=A.
	ref := write("ref.fa", ">chr1\n"+strings.Repeat("ACGT", 15)+"\n")
	vrfsRunSamtools(t, sam, nil, "faidx", ref)

	refSeq := strings.Repeat("ACGT", 15)
	rd := func(qname string, start int, cigar, seq, rg string) string {
		qual := strings.Repeat("I", len(seq))
		return fmt.Sprintf("%s\t0\tchr1\t%d\t60\t%s\t*\t0\t0\t%s\t%s\tRG:Z:%s", qname, start, cigar, seq, qual, rg)
	}
	sub := func(start int, subs map[int]byte) string {
		b := []byte(refSeq[start-1 : start-1+20])
		for p, c := range subs {
			b[p-start] = c
		}
		return string(b)
	}

	// s1.bam (sample S1): 5 reads, 1 with alt A at pos11.
	var s1 strings.Builder
	s1.WriteString("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:60\n@RG\tID:rg1\tSM:S1\n")
	for i := 0; i < 5; i++ {
		subs := map[int]byte{}
		if i == 0 {
			subs[11] = 'A'
		}
		s1.WriteString(rd(fmt.Sprintf("r%d", i), 5, "20M", sub(5, subs), "rg1") + "\n")
	}
	s1bam := vrfsSamToBam(t, sam, dir, "s1", s1.String())

	// s2.bam (sample S2): 8 reads, 2 with alt T at pos11, 1 with alt C at pos21.
	var s2 strings.Builder
	s2.WriteString("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:60\n@RG\tID:rg2\tSM:S2\n")
	for i := 0; i < 8; i++ {
		subs := map[int]byte{}
		if i < 2 {
			subs[11] = 'T'
		}
		if i == 3 {
			subs[21] = 'C'
		}
		s2.WriteString(rd(fmt.Sprintf("q%d", i), 5, "20M", sub(5, subs), "rg2") + "\n")
	}
	s2bam := vrfsSamToBam(t, sam, dir, "s2", s2.String())

	alnList := write("alns.txt", s1bam+"\n"+s2bam+"\n")
	sites := write("sites.txt", "chr1\t11\tG\tA\nchr1\t11\tG\tT\nchr1\t21\tA\tC\n")

	// indel.bam (sample S1): deletion / insertion / plain reads for the
	// indel-classification path.
	var ix strings.Builder
	ix.WriteString("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:60\n@RG\tID:rg1\tSM:S1\n")
	delSeq := refSeq[4:10] + refSeq[12:24] // 6M2D12M body
	insSeq := refSeq[4:10] + "TT" + refSeq[10:22]
	for i := 0; i < 3; i++ {
		ix.WriteString(rd(fmt.Sprintf("d%d", i), 5, "6M2D12M", delSeq, "rg1") + "\n")
	}
	for i := 0; i < 3; i++ {
		ix.WriteString(rd(fmt.Sprintf("i%d", i), 5, "6M2I12M", insSeq, "rg1") + "\n")
	}
	for i := 0; i < 4; i++ {
		ix.WriteString(rd(fmt.Sprintf("p%d", i), 5, "20M", refSeq[4:24], "rg1") + "\n")
	}
	indelBam := vrfsSamToBam(t, sam, dir, "indel", ix.String())
	indelAln := write("indelalns.txt", indelBam+"\n")
	// Sites exercising the indel-after (pos10) and is_del (pos11) classification.
	indelSite := write("indelsites.txt", "chr1\t10\tC\tCGT\nchr1\t10\tC\tA\nchr1\t11\tG\tA\nchr1\t11\tG\tGTT\n")

	return vrfsFixtures{
		dir: dir, ref: ref, alnList: alnList, sites: sites,
		indelAln: indelAln, indelSite: indelSite,
	}
}

// vrfsSamToBam converts SAM text to a sorted, indexed BAM in dir and returns
// its path.
func vrfsSamToBam(t *testing.T, sam, dir, name, samText string) string {
	t.Helper()
	bamPath := filepath.Join(dir, name+".bam")
	// samtools view -b | sort -o name.bam
	viewCmd := exec.Command(sam, "view", "-b", "-")
	viewCmd.Stdin = strings.NewReader(samText)
	var rawBam, verr bytes.Buffer
	viewCmd.Stdout = &rawBam
	viewCmd.Stderr = &verr
	if err := viewCmd.Run(); err != nil {
		t.Fatalf("samtools view (%s): %v\n%s", name, err, verr.String())
	}
	vrfsRunSamtools(t, sam, rawBam.Bytes(), "sort", "-o", bamPath, "-")
	vrfsRunSamtools(t, sam, nil, "index", bamPath)
	return bamPath
}

// runVrfsBoth runs upstream and our port with the same argv (and
// BCFTOOLS_PLUGINS set) and returns their stdout.
func runVrfsBoth(t *testing.T, upstream, ours string, argv []string) (up, our []byte) {
	t.Helper()
	pluginDir := pluginDirAbs(t)
	run := func(bin string) []byte {
		cmd := exec.Command(bin, argv...)
		cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDir)
		var out, errBuf bytes.Buffer
		cmd.Stdout = &out
		cmd.Stderr = &errBuf
		if err := cmd.Run(); err != nil {
			t.Fatalf("%s %v: %v\n%s", bin, argv, err, errBuf.String())
		}
		return out.Bytes()
	}
	return run(upstream), run(ours)
}

func TestVrfsOracleParity(t *testing.T) {
	upstream, ours := requireLive(t)
	sam := upstreamSamtoolsForVrfs(t)
	fx := vrfsBuildFixtures(t, sam)

	cases := []struct {
		name string
		argv []string
	}{
		{"snv-streaming", []string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1"}},
		{"snv-index", []string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1", "-i"}},
		{"snv-default-empty", []string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites}},
		{"recalc-data", []string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1", "-r", "data"}},
		{"nbins-10", []string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1", "-n", "10"}},
		{"nbins-25", []string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1", "-n", "25"}},
		{"batch-k", []string{"+vrfs", "-a", fx.alnList, "--batch", "k=3"}},
		{"batch-1of2", []string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1", "-b", "1/2"}},
		{"batch-2of2", []string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1", "-b", "2/2"}},
		{"indel-classification", []string{"+vrfs", "-f", fx.ref, "-a", fx.indelAln, "-s", fx.indelSite, "-d", "1"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			up, our := runVrfsBoth(t, upstream, ours, tc.argv)
			if !bytes.Equal(up, our) {
				t.Fatalf("vrfs stdout diverges (argv=%v)\n--- upstream ---\n%s\n--- ours ---\n%s",
					tc.argv, up, our)
			}
		})
	}
}

// TestVrfsOracleMerge validates the --merge-batches / --merge-files paths by
// producing two batch files with upstream and merging them with both binaries.
func TestVrfsOracleMerge(t *testing.T) {
	upstream, ours := requireLive(t)
	sam := upstreamSamtoolsForVrfs(t)
	fx := vrfsBuildFixtures(t, sam)
	pluginDir := pluginDirAbs(t)

	runUp := func(argv []string) []byte {
		cmd := exec.Command(upstream, argv...)
		cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+pluginDir)
		var out, errBuf bytes.Buffer
		cmd.Stdout = &out
		cmd.Stderr = &errBuf
		if err := cmd.Run(); err != nil {
			t.Fatalf("upstream %v: %v\n%s", argv, err, errBuf.String())
		}
		return out.Bytes()
	}

	b1 := filepath.Join(fx.dir, "b1.txt")
	b2 := filepath.Join(fx.dir, "b2.txt")
	runUp([]string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1", "-b", "1/2", "-o", b1})
	runUp([]string{"+vrfs", "-f", fx.ref, "-a", fx.alnList, "-s", fx.sites, "-d", "1", "-b", "2/2", "-o", b2})
	list := filepath.Join(fx.dir, "batchlist.txt")
	if err := os.WriteFile(list, []byte(b1+"\n"+b2+"\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	mergeCases := [][]string{
		{"+vrfs", "-m", list, "-d", "1"},
		{"+vrfs", "-M", b1, "-M", b2},
	}
	for _, argv := range mergeCases {
		t.Run(strings.Join(argv[1:3], " "), func(t *testing.T) {
			up, our := runVrfsBoth(t, upstream, ours, argv)
			if !bytes.Equal(up, our) {
				t.Fatalf("vrfs merge stdout diverges (argv=%v)\n--- upstream ---\n%s\n--- ours ---\n%s",
					argv, up, our)
			}
		})
	}
}
