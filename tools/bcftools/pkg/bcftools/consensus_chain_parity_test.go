package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"sync"
	"testing"
)

// upstreamConsensusTools locates the upstream bcftools / bgzip / tabix
// binaries built from the reference_code submodules. It is resolved once
// per test binary via sync.Once. The fields are absolute paths.
type upstreamConsensusTools struct {
	bcftools string
	bgzip    string
	tabix    string
}

var (
	upstreamConsensusOnce  sync.Once
	upstreamConsensusPaths upstreamConsensusTools
	upstreamConsensusErr   error
)

// upstreamBcftoolsConsensus returns the absolute paths to the upstream
// bcftools, bgzip and tabix binaries, building nothing itself: the binaries
// must already have been compiled from the reference_code/bcftools and
// reference_code/htslib submodules (see docs/PARITY_ROADMAP.md). It is
// memoised with sync.Once so the filesystem probing happens at most once.
func upstreamBcftoolsConsensus(t *testing.T) upstreamConsensusTools {
	t.Helper()
	upstreamConsensusOnce.Do(func() {
		root, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code"))
		if err != nil {
			upstreamConsensusErr = err
			return
		}
		paths := upstreamConsensusTools{
			bcftools: filepath.Join(root, "bcftools", "bcftools"),
			bgzip:    filepath.Join(root, "htslib", "bgzip"),
			tabix:    filepath.Join(root, "htslib", "tabix"),
		}
		for _, p := range []string{paths.bcftools, paths.bgzip, paths.tabix} {
			if _, err := os.Stat(p); err != nil {
				upstreamConsensusErr = err
				return
			}
		}
		upstreamConsensusPaths = paths
	})
	if upstreamConsensusErr != nil {
		t.Fatalf("upstream bcftools/bgzip/tabix not built: %v\n"+
			"build them with: (cd reference_code/htslib && ./configure && make) "+
			"&& (cd reference_code/bcftools && make)", upstreamConsensusErr)
	}
	return upstreamConsensusPaths
}

// writeConsensusFixtures writes a FASTA reference and a bgzipped+tabixed
// VCF into dir, returning the FASTA path and the VCF.gz path. The VCF is
// bgzipped and indexed with the upstream tools so it can be fed to upstream
// bcftools consensus, which requires an indexed input.
func writeConsensusFixtures(t *testing.T, dir, refName string, ref []byte, vcfBody string) (faPath, vcfGzPath string) {
	t.Helper()
	tools := upstreamBcftoolsConsensus(t)

	faPath = filepath.Join(dir, "ref.fa")
	faContent := ">" + refName + "\n" + string(ref) + "\n"
	if err := os.WriteFile(faPath, []byte(faContent), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	// Index the FASTA (.fai); upstream bcftools consensus requires it.
	if out, err := exec.Command(tools.bgzip, "--version").CombinedOutput(); err != nil {
		t.Fatalf("bgzip not runnable: %v: %s", err, out)
	}
	faiCmd := exec.Command(tools.bcftools, "faidx", faPath)
	if out, err := faiCmd.CombinedOutput(); err != nil {
		// faidx lives in samtools; bcftools has no faidx. Fall back to a
		// hand-written .fai for a single-line sequence.
		if werr := os.WriteFile(faPath+".fai", []byte(
			refName+"\t"+strconv.Itoa(len(ref))+"\t"+strconv.Itoa(len(refName)+2)+"\t"+strconv.Itoa(len(ref))+"\t"+strconv.Itoa(len(ref)+1)+"\n",
		), 0o644); werr != nil {
			t.Fatalf("write .fai (bcftools faidx failed: %v: %s): %v", err, out, werr)
		}
	}

	vcfPath := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(vcfPath, []byte(vcfBody), 0o644); err != nil {
		t.Fatalf("write vcf: %v", err)
	}
	vcfGzPath = vcfPath + ".gz"
	bgz := exec.Command(tools.bgzip, "-c", vcfPath)
	gzOut, err := os.Create(vcfGzPath)
	if err != nil {
		t.Fatalf("create vcf.gz: %v", err)
	}
	var bgzErr bytes.Buffer
	bgz.Stdout = gzOut
	bgz.Stderr = &bgzErr
	if err := bgz.Run(); err != nil {
		gzOut.Close()
		t.Fatalf("bgzip vcf: %v: %s", err, bgzErr.String())
	}
	gzOut.Close()
	if out, err := exec.Command(tools.tabix, "-p", "vcf", vcfGzPath).CombinedOutput(); err != nil {
		t.Fatalf("tabix index: %v: %s", err, out)
	}
	return faPath, vcfGzPath
}

// runUpstreamConsensus runs upstream bcftools consensus and returns the
// consensus FASTA bytes and, when chainPath is non-empty, the chain bytes.
func runUpstreamConsensus(t *testing.T, dir, faPath, vcfGzPath, chainPath string, extra ...string) (faOut, chainOut []byte) {
	t.Helper()
	tools := upstreamBcftoolsConsensus(t)
	args := []string{"consensus", "-f", faPath}
	if chainPath != "" {
		args = append(args, "-c", chainPath)
	}
	args = append(args, extra...)
	args = append(args, vcfGzPath)
	cmd := exec.Command(tools.bcftools, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools %v failed: %v\nstderr: %s", args, err, stderr.String())
	}
	faOut = stdout.Bytes()
	if chainPath != "" {
		var err error
		chainOut, err = os.ReadFile(chainPath)
		if err != nil {
			t.Fatalf("read upstream chain: %v", err)
		}
	}
	return faOut, chainOut
}

// runPortConsensus runs the Go port (via the file-aware entry point, which
// reads the FASTA and transparently decompresses the bgzipped VCF) and
// returns the consensus FASTA bytes and, when chainPath is non-empty, the
// chain bytes.
func runPortConsensus(t *testing.T, dir, faPath, vcfGzPath, chainPath string, opts ConsensusOptions) (faOut, chainOut []byte) {
	t.Helper()
	opts.ChainFile = chainPath
	var buf bytes.Buffer
	if _, err := ConsensusFile(vcfGzPath, faPath, &buf, opts); err != nil {
		t.Fatalf("port ConsensusFile: %v", err)
	}
	faOut = buf.Bytes()
	if chainPath != "" {
		var err error
		chainOut, err = os.ReadFile(chainPath)
		if err != nil {
			t.Fatalf("read port chain: %v", err)
		}
	}
	return faOut, chainOut
}

// consensusFixtureRef is the shared reference used by the chain parity
// tests. It is long enough that the line-wrapped output and chain offsets
// exercise multiple ungapped blocks.
const consensusFixtureRef = "AAAACCCCGGGGTTTTACGTACGTACGTACGTAAAACCCCGGGGTTTT"

// TestConsensusChainParitySNPsAndIndels compares the port against upstream
// for a mix of SNPs, an insertion and a deletion, byte-for-byte on both the
// consensus FASTA and the liftover chain file.
func TestConsensusChainParitySNPsAndIndels(t *testing.T) {
	dir := t.TempDir()
	body := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=48>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t2\t.\tA\tG\t.\tPASS\t.\n" + // SNP
		"chr1\t6\t.\tC\tCTT\t.\tPASS\t.\n" + // insertion
		"chr1\t10\t.\tGGG\tG\t.\tPASS\t.\n" + // deletion
		"chr1\t20\t.\tT\tA\t.\tPASS\t.\n" // SNP
	faPath, vcfGz := writeConsensusFixtures(t, dir, "chr1", []byte(consensusFixtureRef), body)

	upFa, upChain := runUpstreamConsensus(t, dir, faPath, vcfGz, filepath.Join(dir, "up.chain"))
	goFa, goChain := runPortConsensus(t, dir, faPath, vcfGz, filepath.Join(dir, "go.chain"), ConsensusOptions{})

	if !bytes.Equal(upFa, goFa) {
		t.Fatalf("consensus FASTA mismatch:\n upstream:\n%s\n port:\n%s", upFa, goFa)
	}
	if !bytes.Equal(upChain, goChain) {
		t.Fatalf("chain mismatch:\n upstream:\n%s\n port:\n%s", upChain, goChain)
	}
}

// TestConsensusChainParityBackToBack exercises the back-to-back gap merge
// in push_chain_gap with two adjacent indels.
func TestConsensusChainParityBackToBack(t *testing.T) {
	dir := t.TempDir()
	body := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=48>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t3\t.\tA\tAGG\t.\tPASS\t.\n" + // insertion (REF A at pos 3)
		"chr1\t9\t.\tGGGG\tG\t.\tPASS\t.\n" + // deletion (REF GGGG at pos 9-12)
		"chr1\t33\t.\tA\tACCCC\t.\tPASS\t.\n" // insertion (REF A at pos 33)
	faPath, vcfGz := writeConsensusFixtures(t, dir, "chr1", []byte(consensusFixtureRef), body)

	upFa, upChain := runUpstreamConsensus(t, dir, faPath, vcfGz, filepath.Join(dir, "up.chain"))
	goFa, goChain := runPortConsensus(t, dir, faPath, vcfGz, filepath.Join(dir, "go.chain"), ConsensusOptions{})

	if !bytes.Equal(upFa, goFa) {
		t.Fatalf("consensus FASTA mismatch:\n upstream:\n%s\n port:\n%s", upFa, goFa)
	}
	if !bytes.Equal(upChain, goChain) {
		t.Fatalf("chain mismatch:\n upstream:\n%s\n port:\n%s", upChain, goChain)
	}
}

// TestConsensusHaplotypeParity compares -H selectors (phased index,
// NpIu phased-index / unphased-IUPAC, and plain IUPAC) against upstream.
func TestConsensusHaplotypeParity(t *testing.T) {
	dir := t.TempDir()
	// A diploid sample with phased and unphased SNP genotypes.
	body := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=48>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t2\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n" + // phased het
		"chr1\t6\t.\tC\tT\t.\tPASS\t.\tGT\t1|0\n" + // phased het, other order
		"chr1\t10\t.\tG\tC\t.\tPASS\t.\tGT\t0/1\n" + // unphased het
		"chr1\t20\t.\tT\tA\t.\tPASS\t.\tGT\t1|1\n" // phased hom-alt
	faPath, vcfGz := writeConsensusFixtures(t, dir, "chr1", []byte(consensusFixtureRef), body)

	cases := []struct {
		name string
		arg  string
		sel  HaplotypeSelector
		idx  int
	}{
		{"H1", "1", HapIndex, 1},
		{"H2", "2", HapIndex, 2},
		{"H1pIu", "1pIu", HapPhasedIUPAC, 1},
		{"H2pIu", "2pIu", HapPhasedIUPAC, 2},
		{"HI", "I", HapIUPAC, 0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			upFa, _ := runUpstreamConsensus(t, dir, faPath, vcfGz, "", "-s", "S1", "-H", tc.arg)
			goFa, _ := runPortConsensus(t, dir, faPath, vcfGz, "", ConsensusOptions{
				Sample:         "S1",
				Haplotype:      tc.sel,
				HaplotypeIndex: tc.idx,
			})
			if !bytes.Equal(upFa, goFa) {
				t.Fatalf("-H %s consensus mismatch:\n upstream:\n%s\n port:\n%s", tc.arg, upFa, goFa)
			}
		})
	}
}

// TestConsensusHaplotypeChainParity confirms that the chain file also
// matches upstream when a -H haplotype indel is applied for one sample.
func TestConsensusHaplotypeChainParity(t *testing.T) {
	dir := t.TempDir()
	body := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=48>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t3\t.\tA\tAGG\t.\tPASS\t.\tGT\t0|1\n" + // phased het insertion (REF A at pos 3)
		"chr1\t9\t.\tGGGG\tG\t.\tPASS\t.\tGT\t1|1\n" // hom-alt deletion (REF GGGG at pos 9-12)
	faPath, vcfGz := writeConsensusFixtures(t, dir, "chr1", []byte(consensusFixtureRef), body)

	upFa, upChain := runUpstreamConsensus(t, dir, faPath, vcfGz, filepath.Join(dir, "up.chain"), "-s", "S1", "-H", "2")
	goFa, goChain := runPortConsensus(t, dir, faPath, vcfGz, filepath.Join(dir, "go.chain"), ConsensusOptions{
		Sample:         "S1",
		Haplotype:      HapIndex,
		HaplotypeIndex: 2,
	})
	if !bytes.Equal(upFa, goFa) {
		t.Fatalf("consensus FASTA mismatch:\n upstream:\n%s\n port:\n%s", upFa, goFa)
	}
	if !bytes.Equal(upChain, goChain) {
		t.Fatalf("chain mismatch:\n upstream:\n%s\n port:\n%s", upChain, goChain)
	}
}
