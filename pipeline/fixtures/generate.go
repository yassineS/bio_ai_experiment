// Package fixtures deterministically generates real-sized, valid,
// cross-consistent bioinformatics inputs (reference FASTA, sorted+indexed BAM,
// the same alignments as CRAM, a bgzipped+tabixed VCF plus a plain VCF, and
// plain/BED12 interval files) for the parity-and-performance pipeline.
//
// Generation strategy: a seeded math/rand source produces the raw text (FASTA
// sequence, SAM records, VCF lines, BED intervals); the vendored UPSTREAM tools
// then turn that text into the binary/indexed formats so the fixtures are
// guaranteed valid and reproducible:
//
//	reference FASTA + .fai   : generated text, indexed by `samtools faidx`
//	sorted BAM + .bai/.csi   : SAM text -> `samtools sort` -> `samtools index`
//	CRAM (+ reference)       : `samtools view -C -T ref`
//	VCF.gz + .tbi            : VCF text -> `bgzip` -> `tabix -p vcf`
//	plain VCF                : the same VCF text, uncompressed
//	multi-sample VCF(.gz)    : N-sample VCF text -> `bgzip` -> `tabix -p vcf`
//	BED / BED12              : generated text over the same coordinate space
//	GFF3                     : gene/mRNA/exon/CDS rows over the same contigs
//	FASTQ (SE plain + .gz)   : seeded reads w/ adapters + low-quality tails
//	FASTQ (paired R1/R2)     : matched-name mate pairs for PE QC tools
//
// All randomness flows from a single seed (default 1) so a given scale tier is
// byte-reproducible. Fixtures are cached under pipeline/.fixtures/<scale>/ with
// a manifest; Generate regenerates only when the cache is missing or stale.
package fixtures

import (
	"bytes"
	"fmt"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// DefaultSeed is the RNG seed used unless overridden.
const DefaultSeed int64 = 1

// Options configure a Generate call.
type Options struct {
	Scale Scale
	Seed  int64  // 0 -> DefaultSeed
	Dir   string // base fixtures dir; "" -> <repo>/pipeline/.fixtures
	Force bool   // regenerate even if a valid cache exists
	Logf  func(format string, args ...any)
}

func (o Options) log(format string, args ...any) {
	if o.Logf != nil {
		o.Logf(format, args...)
	}
}

// Generate builds (or reuses) the fixture set for the requested scale and
// returns its manifest.
func Generate(opt Options) (*Manifest, error) {
	if opt.Seed == 0 {
		opt.Seed = DefaultSeed
	}
	if opt.Scale == "" {
		opt.Scale = Small
	}
	base := opt.Dir
	if base == "" {
		root, err := upstream.RepoRoot()
		if err != nil {
			return nil, err
		}
		base = filepath.Join(root, "pipeline", ".fixtures")
	}
	dir := filepath.Join(base, string(opt.Scale))

	if !opt.Force {
		if m, err := loadManifest(dir); err == nil && m.valid(opt.Scale, opt.Seed) {
			opt.log("fixtures: reusing cached %s set at %s", opt.Scale, dir)
			return m, nil
		}
	}
	opt.log("fixtures: generating %s set at %s (seed=%d)", opt.Scale, dir, opt.Seed)

	// Resolve upstream tools up front so a missing one fails before any work.
	samtools, err := upstream.Binary("samtools")
	if err != nil {
		return nil, err
	}
	bgzip, err := upstream.Binary("bgzip")
	if err != nil {
		return nil, err
	}
	tabix, err := upstream.Binary("tabix")
	if err != nil {
		return nil, err
	}

	if err := os.RemoveAll(dir); err != nil {
		return nil, err
	}
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return nil, err
	}

	p := ParamsFor(opt.Scale)
	rng := rand.New(rand.NewSource(opt.Seed))
	contigs := makeContigs(p)

	m := &Manifest{
		Version: manifestVersion,
		Scale:   string(opt.Scale),
		Seed:    opt.Seed,
		Params:  p,
		Files:   map[string]string{},
		Sizes:   map[string]int64{},
		Digests: map[string]string{},
	}

	// --- Reference FASTA + .fai ---
	fastaPath := filepath.Join(dir, "ref.fa")
	if err := writeFasta(fastaPath, contigs, rng); err != nil {
		return nil, err
	}
	if err := run(samtools, "faidx", fastaPath); err != nil {
		return nil, fmt.Errorf("samtools faidx: %w", err)
	}
	_ = m.recordFile("fasta", fastaPath, true)
	_ = m.recordFile("fai", fastaPath+".fai", false)

	// --- SAM -> sorted BAM + index ---
	samPath := filepath.Join(dir, "reads.sam")
	if err := writeSAM(samPath, contigs, p, rng); err != nil {
		return nil, err
	}
	bamPath := filepath.Join(dir, "reads.bam")
	if err := run(samtools, "sort", "-o", bamPath, samPath); err != nil {
		return nil, fmt.Errorf("samtools sort: %w", err)
	}
	if err := run(samtools, "index", bamPath); err != nil {
		return nil, fmt.Errorf("samtools index: %w", err)
	}
	// Also a CSI index (for the large-contig path / parity coverage).
	if err := run(samtools, "index", "-c", bamPath); err != nil {
		return nil, fmt.Errorf("samtools index -c: %w", err)
	}
	_ = m.recordFile("bam", bamPath, false)
	_ = m.recordFile("bai", bamPath+".bai", false)
	_ = m.recordFile("csi", bamPath+".csi", false)

	// --- CRAM (+ reference) ---
	cramPath := filepath.Join(dir, "reads.cram")
	if err := run(samtools, "view", "-C", "-T", fastaPath, "-o", cramPath, bamPath); err != nil {
		return nil, fmt.Errorf("samtools view -C: %w", err)
	}
	if err := run(samtools, "index", cramPath); err != nil {
		return nil, fmt.Errorf("samtools index cram: %w", err)
	}
	_ = m.recordFile("cram", cramPath, false)
	_ = m.recordFile("cram_crai", cramPath+".crai", false)

	// --- VCF: plain + bgzipped + tabixed ---
	vcfPath := filepath.Join(dir, "variants.vcf")
	if err := writeVCF(vcfPath, contigs, p, rng); err != nil {
		return nil, err
	}
	// bgzip writes variants.vcf.gz and (with -k) keeps the plain file.
	if err := run(bgzip, "-kf", vcfPath); err != nil {
		return nil, fmt.Errorf("bgzip vcf: %w", err)
	}
	vcfGz := vcfPath + ".gz"
	if err := run(tabix, "-p", "vcf", "-f", vcfGz); err != nil {
		return nil, fmt.Errorf("tabix vcf: %w", err)
	}
	_ = m.recordFile("vcf_plain", vcfPath, true)
	_ = m.recordFile("vcf", vcfGz, false)
	_ = m.recordFile("vcf_tbi", vcfGz+".tbi", false)

	// --- Multi-sample VCF: plain + bgzipped + tabixed ---
	// Several vcftools modes (relatedness, het, LD) need more than one sample;
	// the single-sample VCF above keeps the simpler per-site modes simple.
	vcfMultiPath := filepath.Join(dir, "variants.multi.vcf")
	if err := writeMultiSampleVCF(vcfMultiPath, contigs, p, rng); err != nil {
		return nil, err
	}
	if err := run(bgzip, "-kf", vcfMultiPath); err != nil {
		return nil, fmt.Errorf("bgzip multi vcf: %w", err)
	}
	vcfMultiGz := vcfMultiPath + ".gz"
	if err := run(tabix, "-p", "vcf", "-f", vcfMultiGz); err != nil {
		return nil, fmt.Errorf("tabix multi vcf: %w", err)
	}
	_ = m.recordFile("vcf_multi_plain", vcfMultiPath, true)
	_ = m.recordFile("vcf_multi", vcfMultiGz, false)
	_ = m.recordFile("vcf_multi_tbi", vcfMultiGz+".tbi", false)

	// --- BED (plain + BED12) ---
	bedPath := filepath.Join(dir, "intervals.bed")
	bed12Path := filepath.Join(dir, "intervals12.bed")
	genomePath := filepath.Join(dir, "genome.txt")
	if err := writeBED(bedPath, bed12Path, genomePath, contigs, p, rng); err != nil {
		return nil, err
	}
	_ = m.recordFile("bed", bedPath, true)
	_ = m.recordFile("bed12", bed12Path, true)
	_ = m.recordFile("genome", genomePath, true)

	// --- BEDPE (paired-end intervals for bedpairtobed / bedpairtopair) ---
	bedpePath := filepath.Join(dir, "pairs.bedpe")
	if err := writeBEDPE(bedpePath, contigs, p, rng); err != nil {
		return nil, err
	}
	_ = m.recordFile("bedpe", bedpePath, true)

	// --- GFF3 (genes/mRNA/exon/CDS over the same contigs) ---
	gffPath := filepath.Join(dir, "annotations.gff3")
	if err := writeGFF(gffPath, contigs, p, rng); err != nil {
		return nil, err
	}
	_ = m.recordFile("gff", gffPath, true)

	// --- FASTQ: single-end (plain + .gz) and paired-end ---
	fastqPath := filepath.Join(dir, "reads.fastq")
	fastqGzPath := fastqPath + ".gz"
	if err := writeFastqSE(fastqPath, fastqGzPath, p, rng); err != nil {
		return nil, err
	}
	_ = m.recordFile("fastq", fastqPath, true)
	_ = m.recordFile("fastq_gz", fastqGzPath, false)

	fastq1Path := filepath.Join(dir, "reads_R1.fastq")
	fastq2Path := filepath.Join(dir, "reads_R2.fastq")
	if err := writeFastqPE(fastq1Path, fastq2Path, p, rng); err != nil {
		return nil, err
	}
	_ = m.recordFile("fastq1", fastq1Path, true)
	_ = m.recordFile("fastq2", fastq2Path, true)

	if err := m.save(dir); err != nil {
		return nil, err
	}
	opt.log("fixtures: done (%d files)", len(m.Files))
	return m, nil
}

// contig is one reference sequence.
type contig struct {
	Name string
	Len  int
	Seq  []byte
}

// makeContigs returns the contig skeleton (names + lengths) for a scale.
func makeContigs(p Params) []contig {
	cs := make([]contig, p.NumContigs)
	for i := range cs {
		cs[i] = contig{Name: fmt.Sprintf("chr%d", i+1), Len: p.ContigLen}
	}
	return cs
}

// run executes a vendored binary, capturing combined output for error context.
func run(bin string, args ...string) error {
	var buf bytes.Buffer
	cmd := exec.Command(bin, args...)
	cmd.Stdout = &buf
	cmd.Stderr = &buf
	if err := cmd.Run(); err != nil {
		return fmt.Errorf("%s %s: %v\n%s", filepath.Base(bin), strings.Join(args, " "), err, buf.String())
	}
	return nil
}
