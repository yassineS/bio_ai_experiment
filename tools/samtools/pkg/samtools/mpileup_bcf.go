package samtools

// mpileup_bcf.go implements `samtools mpileup -g/-u` — BCF (and
// uncompressed-BCF) genotype-likelihood output.
//
// Upstream history: `samtools mpileup` used to carry the genotype-likelihood
// caller itself (htslib's bam2bcf). As of samtools 1.10 that path was
// REMOVED from `samtools mpileup`; the modern tool (the vendored
// reference_code/samtools is 1.23.1) only emits the text pileup and prints
//
//	"Note that using \"samtools mpileup\" to generate BCF or VCF files has
//	 been removed. To output these formats, please use \"bcftools mpileup\"
//	 instead."
//
// when BCF/VCF output is requested (the `-g`/`-u` short options are no
// longer in the getopt string at all, so upstream rejects them outright).
//
// `bcftools mpileup` is upstream's sanctioned replacement and is itself a
// thin driver over the very same htslib bam2bcf genotype-likelihood
// pipeline. This repository already contains a complete, tested port of
// that pipeline under tools/bcftools (errmod MAQ model, bcf_call_glfgen /
// bcf_call_combine / bcf_call2bcf, BCF emit through pkg/htsgo/bcf). Rather
// than re-port bam2bcf a second time, `samtools mpileup -g/-u` here
// DELEGATES to the bcftools mpileup engine, translating the samtools
// mpileup options into bcftools MpileupOptions. The emitted BCF carries the
// per-site genotype likelihoods (FORMAT/PL), the `<*>` unseen allele, and
// the standard INFO/FORMAT tags exactly as `bcftools mpileup` produces them
// — which is the live upstream behaviour for this output format.
//
// The one indel-row caveat of the bcftools mpileup port is inherited here
// and documented in docs/PARITY_ROADMAP.md#samtools (mpileup): the SNP
// genotype-likelihood path is complete; full bam2bcf_indel calling is the
// shared deferred remainder tracked against bcftools mpileup.

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

// MpileupBCFOptions configures the `samtools mpileup -g/-u` BCF
// genotype-likelihood output path. It is a small, BCF-specific superset of
// the text-pileup MpileupOptions: the fields that drive the
// genotype-likelihood model and the BCF container.
//
// Defaults that differ from the text path follow the genotype-likelihood
// caller (the htslib bam2bcf defaults that `bcftools mpileup` applies):
// when a numeric field is left at its zero value the bcftools engine fills
// the upstream default (min-BQ 1, max-depth 250, max-BQ 60, delta-BQ 30).
type MpileupBCFOptions struct {
	// Inputs is the list of BAM/SAM paths to pile up (one sample column
	// per input). Mirrors the text path's MpileupOptions.Inputs.
	Inputs []string
	// FastaRef is the reference FASTA (-f). Required: the genotype-
	// likelihood caller needs the REF base at every site.
	FastaRef string
	// BamList is the -b file listing one BAM path per line; resolved paths
	// are appended to Inputs.
	BamList string
	// Regions are samtools-style "chr[:beg-end]" specifiers (-r).
	Regions []string
	// PositionsFile is the -l BED/positions file restricting emitted sites.
	PositionsFile string

	// MinMAPQ skips reads with MAPQ below this value (-q).
	MinMAPQ uint8
	// MinBaseQ skips bases with quality below this value (-Q). Zero means
	// "use the bcftools genotype-likelihood default" (1).
	MinBaseQ uint8
	// MaxDepth caps reads per position (-d). Zero means "use the bcftools
	// genotype-likelihood default" (250).
	MaxDepth int

	// CountOrphans includes anomalous-pair reads (-A).
	CountOrphans bool
	// IgnoreOverlaps disables the smart-overlap mate de-weighting (-x).
	IgnoreOverlaps bool
	// NoBAQ disables BAQ realignment (-B).
	NoBAQ bool
	// RedoBAQ recomputes BAQ from scratch (-E).
	RedoBAQ bool

	// Uncompressed selects uncompressed BCF (-u). When false the output is
	// BGZF-compressed BCF (-g).
	Uncompressed bool

	// Output is the destination path (-o); empty means the supplied writer
	// (stdout) is used. Wired by the caller, not by this struct.
	Output string
	// Threads is accepted for flag compatibility; the engine is
	// single-threaded.
	Threads int
}

// MpileupBCF runs the `samtools mpileup -g/-u` genotype-likelihood caller,
// writing BCF (or, when opts.Uncompressed is set, uncompressed BCF) to out.
//
// It delegates to the shared bcftools mpileup engine (the htslib bam2bcf
// port), so the emitted records carry FORMAT/PL genotype likelihoods, the
// `<*>` unseen allele, and the standard INFO/FORMAT tags identical to
// `bcftools mpileup`. See the file-level comment for the upstream-removal
// history and the indel-row remainder.
func MpileupBCF(opts MpileupBCFOptions, out io.Writer) error {
	if len(opts.Inputs) == 0 && opts.BamList == "" {
		return fmt.Errorf("samtools mpileup: no input files")
	}
	if opts.FastaRef == "" {
		return fmt.Errorf("samtools mpileup: -f/--fasta-ref is required for -g/-u BCF output")
	}

	format := bcftools.OutputBCF
	if opts.Uncompressed {
		format = bcftools.OutputBCFUncompressed
	}

	bopts := bcftools.MpileupOptions{
		Inputs:         opts.Inputs,
		FastaRef:       opts.FastaRef,
		BamList:        opts.BamList,
		Regions:        opts.Regions,
		MinMQ:          opts.MinMAPQ,
		MinBQ:          opts.MinBaseQ,
		MaxDepth:       opts.MaxDepth,
		CountOrphans:   opts.CountOrphans,
		IgnoreOverlaps: opts.IgnoreOverlaps,
		NoBAQ:          opts.NoBAQ,
		RedoBAQ:        opts.RedoBAQ,
		OutputFormat:   format,
		Threads:        opts.Threads,
	}
	// samtools mpileup's -l positions file maps to bcftools mpileup's
	// -T/--targets-file (a post-filter restricting emitted sites). Both
	// accept BED-style and "chrom\tpos" listings.
	if opts.PositionsFile != "" {
		bopts.TargetsFile = opts.PositionsFile
	}

	if err := bcftools.MpileupFile(bopts, out); err != nil {
		return fmt.Errorf("samtools mpileup: %w", err)
	}
	return nil
}
