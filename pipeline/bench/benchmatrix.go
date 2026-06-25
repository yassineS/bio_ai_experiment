package bench

import (
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
)

// benchPlan is one fully-resolved invocation: the binaries to time on each side
// and their argv, plus optional stdin and per-side stdout redirection (the
// streaming operations whose output is stdout write it to a temp file so the
// write cost is counted symmetrically).
type benchPlan struct {
	ourTool, upKey      string
	ourArgs, upArgs     []string
	stdin               string
	ourStdout, upStdout string
}

// BenchCell is one benchmark in the matrix: a named operation whose cost is
// measured for OUR binary against the vendored UPSTREAM binary. Group buckets
// the cell by the input file type it stresses (BAM/CRAM/VCF/BED/FASTQ) so the
// report can table the file-format coverage the manuscript needs.
type BenchCell struct {
	Name  string
	Group string
	Build func(m *fixtures.Manifest, tmp string) benchPlan
}

// sameArgs builds a plan where both sides share argv (samtools/bcftools/
// vcftools, whose CLI shape we match exactly). stdoutFile, when true, captures
// stdout to a per-side temp file.
func sameArgs(ourTool, upKey string, stdoutFile bool, args ...string) func(*fixtures.Manifest, string) benchPlan {
	return func(_ *fixtures.Manifest, _ string) benchPlan {
		p := benchPlan{ourTool: ourTool, upKey: upKey, ourArgs: args, upArgs: args}
		if stdoutFile {
			// Stream to the null device, not a temp file: the write() syscalls
			// are still issued (so output-generation cost is counted) but nothing
			// hits disk. This keeps cells whose stdout is huge — e.g. `bcftools
			// call` over a 15M-record mpileup VCF (~1.7 GB at medium, ~20 GB at
			// large) — from filling /tmp and aborting the run.
			p.ourStdout = os.DevNull
			p.upStdout = os.DevNull
		}
		return p
	}
}

// bedArgs builds a plan for a bed* tool: OUR side is the standalone binary
// (e.g. bedintersect) with no subcommand token; UPSTREAM is `bedtools <sub>`.
// bed* operations all stream to stdout.
func bedArgs(ourTool, sub string, args ...string) func(*fixtures.Manifest, string) benchPlan {
	return func(_ *fixtures.Manifest, _ string) benchPlan {
		return benchPlan{
			ourTool: ourTool, upKey: "bedtools",
			ourArgs:   args,
			upArgs:    append([]string{sub}, args...),
			ourStdout: os.DevNull,
			upStdout:  os.DevNull,
		}
	}
}

// BenchMatrix returns the heavy, streaming operations benchmarked across the
// scale tiers. The set deliberately spans every full-scale input format — BAM,
// CRAM, VCF, BED, FASTQ — and includes the same-file vs different-file overlap
// cases and the two-input set operations (bedtools intersect, bcftools isec)
// called out as load-bearing for the scalability study.
func BenchMatrix() []BenchCell {
	cells := []BenchCell{
		// ---- BAM (samtools) ----
		{"sam_view_bam2bam", "BAM", func(m *fixtures.Manifest, tmp string) benchPlan {
			out := filepath.Join(tmp, "out.bam")
			a := []string{"view", "-b", "-o", out, m.Path("bam")}
			return benchPlan{ourTool: "samtools", upKey: "samtools", ourArgs: a, upArgs: a}
		}},
		{"sam_sort", "BAM", func(m *fixtures.Manifest, tmp string) benchPlan {
			out := filepath.Join(tmp, "sorted.bam")
			a := []string{"sort", "-o", out, m.Path("bam")}
			return benchPlan{ourTool: "samtools", upKey: "samtools", ourArgs: a, upArgs: a}
		}},
		{"sam_flagstat", "BAM", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("samtools", "samtools", true, "flagstat", m.Path("bam"))(m, tmp)
		}},
		{"sam_stats", "BAM", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("samtools", "samtools", true, "stats", m.Path("bam"))(m, tmp)
		}},
		{"sam_depth", "BAM", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("samtools", "samtools", true, "depth", m.Path("bam"))(m, tmp)
		}},
		{"sam_mpileup", "BAM", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("samtools", "samtools", true, "mpileup", "-f", m.Path("fasta"), m.Path("bam"))(m, tmp)
		}},

		// ---- CRAM (samtools): the reference-compressed encode and decode ----
		{"sam_view_bam2cram", "CRAM", func(m *fixtures.Manifest, tmp string) benchPlan {
			out := filepath.Join(tmp, "out.cram")
			a := []string{"view", "-C", "-T", m.Path("fasta"), "-o", out, m.Path("bam")}
			return benchPlan{ourTool: "samtools", upKey: "samtools", ourArgs: a, upArgs: a}
		}},
		{"sam_view_cram2bam", "CRAM", func(m *fixtures.Manifest, tmp string) benchPlan {
			out := filepath.Join(tmp, "out.bam")
			a := []string{"view", "-b", "-T", m.Path("fasta"), "-o", out, m.Path("cram")}
			return benchPlan{ourTool: "samtools", upKey: "samtools", ourArgs: a, upArgs: a}
		}},

		// ---- VCF (bcftools) ----
		{"bcf_view", "VCF", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("bcftools", "bcftools", true, "view", m.Path("vcf_plain"))(m, tmp)
		}},
		{"bcf_norm", "VCF", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("bcftools", "bcftools", true, "norm", "-f", m.Path("fasta"), m.Path("vcf_plain"))(m, tmp)
		}},
		{"bcf_stats", "VCF", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("bcftools", "bcftools", true, "stats", m.Path("vcf_plain"))(m, tmp)
		}},
		{"bcf_query", "VCF", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("bcftools", "bcftools", true, "query", "-f", `%CHROM\t%POS\t%REF\t%ALT\n`, m.Path("vcf_plain"))(m, tmp)
		}},
		{"bcf_call", "VCF", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("bcftools", "bcftools", true, "call", "-m", m.Path("vcf_pl"))(m, tmp)
		}},
		// Two-input set operation: bcftools isec over the two tabix-indexed VCFs.
		{"bcf_isec", "VCF", func(m *fixtures.Manifest, tmp string) benchPlan {
			od := filepath.Join(tmp, "ourisec")
			ud := filepath.Join(tmp, "upisec")
			return benchPlan{
				ourTool: "bcftools", upKey: "bcftools",
				ourArgs: []string{"isec", "-p", od, m.Path("vcf"), m.Path("vcf_multi")},
				upArgs:  []string{"isec", "-p", ud, m.Path("vcf"), m.Path("vcf_multi")},
			}
		}},

		// ---- BED (bedtools) — incl. same-file and different-file intersect ----
		{"bed_intersect_self", "BED", func(m *fixtures.Manifest, tmp string) benchPlan {
			return bedArgs("bedintersect", "intersect", "-a", m.Path("bed"), "-b", m.Path("bed"))(m, tmp)
		}},
		{"bed_intersect_pair", "BED", func(m *fixtures.Manifest, tmp string) benchPlan {
			return bedArgs("bedintersect", "intersect", "-a", m.Path("bed"), "-b", m.Path("bed12"))(m, tmp)
		}},
		{"bed_merge", "BED", func(m *fixtures.Manifest, tmp string) benchPlan {
			return bedArgs("bedmerge", "merge", "-i", m.Path("bed"))(m, tmp)
		}},
		{"bed_coverage", "BED", func(m *fixtures.Manifest, tmp string) benchPlan {
			return bedArgs("bedcoverage", "coverage", "-a", m.Path("bed"), "-b", m.Path("bed12"))(m, tmp)
		}},
		{"bed_genomecov", "BED", func(m *fixtures.Manifest, tmp string) benchPlan {
			return bedArgs("bedgenomecov", "genomecov", "-i", m.Path("bed"), "-g", m.Path("genome"))(m, tmp)
		}},
		{"bed_sort", "BED", func(m *fixtures.Manifest, tmp string) benchPlan {
			return bedArgs("bedsort", "sort", "-i", m.Path("bed"))(m, tmp)
		}},

		// ---- FASTQ (QC tools) ----
		{"seqtk_seq", "FASTQ", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("seqtk", "seqtk", true, "seq", "-A", m.Path("fastq"))(m, tmp)
		}},
		{"sickle_se", "FASTQ", func(m *fixtures.Manifest, tmp string) benchPlan {
			out := filepath.Join(tmp, "trimmed.fq")
			a := []string{"se", "-f", m.Path("fastq"), "-t", "sanger", "-o", out}
			return benchPlan{ourTool: "sickle", upKey: "sickle", ourArgs: a, upArgs: a}
		}},
		{"sickle_pe", "FASTQ", func(m *fixtures.Manifest, tmp string) benchPlan {
			o1 := filepath.Join(tmp, "trim_R1.fq")
			o2 := filepath.Join(tmp, "trim_R2.fq")
			s := filepath.Join(tmp, "trim_singles.fq")
			a := []string{"pe", "-f", m.Path("fastq1"), "-r", m.Path("fastq2"),
				"-t", "sanger", "-o", o1, "-p", o2, "-s", s}
			return benchPlan{ourTool: "sickle", upKey: "sickle", ourArgs: a, upArgs: a}
		}},
		{"seqtk_comp", "FASTQ", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("seqtk", "seqtk", true, "comp", m.Path("fasta"))(m, tmp)
		}},
		{"seqtk_trimfq", "FASTQ", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("seqtk", "seqtk", true, "trimfq", m.Path("fastq"))(m, tmp)
		}},
		{"seqtk_fqchk", "FASTQ", func(m *fixtures.Manifest, tmp string) benchPlan {
			return sameArgs("seqtk", "seqtk", true, "fqchk", m.Path("fastq"))(m, tmp)
		}},
	}
	return cells
}
