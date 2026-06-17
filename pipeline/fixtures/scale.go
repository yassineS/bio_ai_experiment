package fixtures

import (
	"fmt"
	"os"
	"strings"
)

// Scale is a fixture-size tier. The pipeline generates the same set of files at
// every tier; only the dimensions (reference length, read count, variant count,
// interval count) grow. Smoke is sized for CI (sub-MB, runs in seconds); large
// approaches half a gigabyte for realistic performance benchmarking.
type Scale string

// The supported scale tiers.
const (
	Smoke  Scale = "smoke"
	Small  Scale = "small"
	Medium Scale = "medium"
	Large  Scale = "large"
)

// Params are the per-tier generation dimensions.
//
// Approximate on-disk footprint of the generated set (reference + BAM + CRAM +
// VCFs + BED files and indexes):
//
//	smoke  : a few hundred KB   (CI-friendly; default for `-scale=smoke`)
//	small  : ~5 MB              (default tier)
//	medium : ~50 MB             (benchmarks)
//	large  : ~500 MB            (heavy benchmarks)
//
// The numbers are chosen so that read bases (~ReadLen*Reads) and reference
// bases (NumContigs*ContigLen) land near those targets; they are documented in
// pipeline/README.md and asserted loosely by the manifest size, not exactly.
type Params struct {
	NumContigs int // reference contigs
	ContigLen  int // bases per contig
	Reads      int // total read pairs' worth of single reads
	ReadLen    int // bases per read
	Variants   int // VCF records
	Intervals  int // BED records

	// FastqReads is the number of records in the single-end FASTQ fixture and
	// the number of pairs in the paired-end fixtures. The QC/adapter-trimming
	// tools (seqtk, prinseq, sickle, skewer, fastp) consume these.
	FastqReads int
	// FastqReadLen is the (mean) read length of the FASTQ fixtures; the
	// generator varies the actual length per read around this so length-filter
	// flags have something to do.
	FastqReadLen int
	// Genes is the number of gene loci in the GFF3 fixture (each gene expands
	// into an mRNA plus a few exon/CDS rows over the shared coordinate space).
	Genes int
	// MultiSamples is the number of samples in the multi-sample VCF fixture
	// (used by vcftools relatedness / LD / per-sample modes that need >1
	// sample). The single-sample VCF keeps one sample for the simpler modes.
	MultiSamples int
}

// paramsByScale maps each tier to its dimensions.
var paramsByScale = map[Scale]Params{
	Smoke:  {NumContigs: 2, ContigLen: 20_000, Reads: 2_000, ReadLen: 100, Variants: 400, Intervals: 500, FastqReads: 2_000, FastqReadLen: 100, Genes: 60, MultiSamples: 8},
	Small:  {NumContigs: 4, ContigLen: 250_000, Reads: 40_000, ReadLen: 150, Variants: 8_000, Intervals: 6_000, FastqReads: 40_000, FastqReadLen: 150, Genes: 800, MultiSamples: 12},
	Medium: {NumContigs: 8, ContigLen: 2_000_000, Reads: 300_000, ReadLen: 150, Variants: 60_000, Intervals: 40_000, FastqReads: 300_000, FastqReadLen: 150, Genes: 6_000, MultiSamples: 16},
	Large:  {NumContigs: 16, ContigLen: 12_000_000, Reads: 2_500_000, ReadLen: 150, Variants: 400_000, Intervals: 250_000, FastqReads: 2_000_000, FastqReadLen: 150, Genes: 40_000, MultiSamples: 24},
}

// AllScales lists the tiers in increasing size order.
func AllScales() []Scale { return []Scale{Smoke, Small, Medium, Large} }

// ParseScale resolves a scale name (case-insensitive), falling back to the
// PIPELINE_SCALE environment variable and finally to Small.
func ParseScale(s string) (Scale, error) {
	if s == "" {
		s = os.Getenv("PIPELINE_SCALE")
	}
	if s == "" {
		return Small, nil
	}
	switch Scale(strings.ToLower(s)) {
	case Smoke:
		return Smoke, nil
	case Small:
		return Small, nil
	case Medium:
		return Medium, nil
	case Large:
		return Large, nil
	default:
		return "", fmt.Errorf("unknown scale %q (want one of smoke|small|medium|large)", s)
	}
}

// ParamsFor returns the generation dimensions for a tier.
func ParamsFor(s Scale) Params { return paramsByScale[s] }
