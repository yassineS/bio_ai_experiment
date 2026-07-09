package realbench

import (
	"os"
	"path/filepath"
	"strings"
)

// InputKind is a bitmask of the inputs a cell requires. A cell whose required
// inputs are not all present SKIPs (it never crashes). The bits are OR-ed so a
// cell can require, e.g., both a BAM and the reference (CRAM paths).
type InputKind uint32

const (
	NeedBAM          InputKind = 1 << iota // the tier BAM (.bam + .bai)
	NeedCRAM                               // the tier CRAM (.cram + .crai)
	NeedRef                                // the indexed reference FASTA (.fa + .fai)
	NeedVCF                                // the bgzipped + indexed VCF (.vcf.gz + .tbi)
	NeedFastq1                             // read-1 FASTQ (R1.fq.gz)
	NeedFastq2                             // read-2 FASTQ (R2.fq.gz)
	NeedBED                                // the intervals BED
	NeedGFF                                // the genes GFF (.gff.gz)
	NeedBED4                               // a BED4 (named) derived from the intervals BED
	NeedBEDPE                              // a BEDPE (>=10 field) derived from the intervals BED
	NeedWindow                             // an 8-field paired/windowed BED derived from the BED
	NeedFastqPlain                         // a decompressed plain FASTQ derived from Fastq1
	NeedNameSortBAM                        // a name-collated BAM derived from BAM (fixmate input)
	NeedFixmateBAM                         // a name-sort|fixmate -m|coord-sort BAM (markdup input)
	NeedBedGraph                           // a 4-col BedGraph derived from the intervals BED (unionbedg)
	NeedSampleRename                       // a one-line sample-rename file (bcftools reheader -s)
	NeedNormGFF                            // a csq-normalised GFF (biotype= injected) derived from GFF
)

// PostKind names how a cell's primary output is turned into a comparable text
// stream for the byte-exact parity check.
//
//   - PostStdout: the command's stdout IS the comparison stream; the digester
//     consumes it directly (provenance-stripped).
//   - PostViewSAM: the command WRITES an alignment file (BAM/CRAM); it is
//     re-decoded through `<bin> view -h` so two outputs with different
//     BGZF/CRAM framing are compared by decoded records (the repo-documented
//     caveat), exactly like realparity's postViewSAM.
//   - PostBgzipD: the command writes a BGZF/.gz file at {out}; it is re-decoded
//     through a multistream gzip reader so the decoded payload — not the
//     framing — is compared.
//   - PostStdoutGzip: the command's stdout is a BGZF/gzip stream (e.g. `bgzip
//     -c`); it is gzip-decoded on the fly so the decompressed payload — not the
//     framing — is compared.
//   - PostFile: the command writes a plain-text file at {out}; its bytes are
//     read back and compared directly.
//   - PostNone: the command's success (exit status) is the signal (quickcheck);
//     parity is "same exit verdict", no output is compared.
//   - PostOursOnly: there is no comparable upstream pair; the cell runs ours
//     only as a perf cell and records parity SKIP.
type PostKind int

const (
	PostStdout PostKind = iota
	PostStdoutGzip
	PostViewSAM
	PostBgzipD
	PostFile
	PostNone
	PostOursOnly
)

// CellSpec is the declarative description of one matrix cell. The same spec is
// instantiated against our binary and the upstream binary; argv is built by
// substituting the resolved input paths into Args via the placeholders below.
type CellSpec struct {
	// Tool is the logical tool name (e.g. "samtools", "bcftools", "seqtk").
	// For bed* cells it is the bed* tool name (e.g. "bedintersect"); the runner
	// maps it to our per-bed binary and to `bedtools <Sub>` upstream.
	Tool string
	// Name is the unique cell identifier (tool_subcommand_variant).
	Name string
	// Subcommand is the human label for the subcommand exercised (for the report).
	Subcommand string
	// Need is the OR-ed set of required inputs.
	Need InputKind
	// Post selects how the comparable stream is produced.
	Post PostKind
	// OurArgs is the argv for our binary (after the program name). Placeholders
	// are substituted at run time.
	OurArgs []string
	// UpArgs is the argv for the upstream binary when it differs from OurArgs
	// (e.g. our per-bed binary takes no subcommand, but upstream needs
	// `bedtools <sub> ...`). When nil, OurArgs is used for both sides.
	UpArgs []string
	// WriteOut, when non-empty, is the file extension the command writes its
	// primary output to (e.g. ".bam", ".cram", ".vcf.gz", ".txt"). The runner
	// allocates a per-side temp path, substitutes it for {out}, then post-
	// processes it per Post.
	WriteOut string
	// NeedsCRAMEnv forces the offline CRAM environment (REF_PATH pointing away
	// from any network cache) for both the command and any re-decode step.
	NeedsCRAMEnv bool
	// WorkDirOut, when true, runs the command in a per-side temp working
	// directory and compares the named output file (Compare) relative to it.
	// Used by tools (fastp, mosdepth) that write fixed-named output files into
	// cwd rather than to an -o path.
	WorkDirOut bool
	// StdoutView, when true, marks a PostViewSAM cell whose comparable alignment
	// is written to STDOUT rather than to a file (WriteOut) or a work-dir file
	// (WorkDirOut). The runner captures the timed run's stdout into a temp .bam
	// and re-decodes it through `samtools view -h`, so two BAMs with different
	// BGZF framing are compared by decoded records. Used by `bedtag`, which (like
	// upstream `bedtools tag`) writes its tagged BAM only to stdout. The decode
	// samtools is resolved per side (ours vs upstream) because the producing
	// binary (bedtag/bedtools) has no `view` subcommand of its own.
	StdoutView bool
	// Compare, when WorkDirOut is set, names the output file (relative to the
	// work dir) whose decoded bytes are compared; its decode follows Post.
	Compare string
}

// placeholders the Args may contain. {out}/{outdir} are per-side temp paths
// allocated by the runner.
const (
	phBAM          = "{bam}"
	phCRAM         = "{cram}"
	phRef          = "{ref}"
	phFai          = "{fai}"
	phVCF          = "{vcf}"
	phFastq1       = "{fastq1}"
	phFastq2       = "{fastq2}"
	phFastqPlain   = "{fastqplain}"
	phBED          = "{bed}"
	phGFF          = "{gff}"
	phNormGFF      = "{normgff}"
	phBED4         = "{bed4}"
	phBEDPE        = "{bedpe}"
	phWindow       = "{window}"
	phNameBAM      = "{namebam}"
	phFixmateBAM   = "{fixmatebam}"
	phBedGraph     = "{bedgraph}"
	phSampleRename = "{samplerename}"
	phOut          = "{out}"
	phOutdir       = "{outdir}"
)

// Inputs holds the resolved input file paths (empty when not provided) for the
// active tier.
type Inputs struct {
	Ref    string
	BAM    string
	CRAM   string
	VCF    string
	Fastq1 string
	Fastq2 string
	BED    string
	GFF    string
	// BED4/BEDPE/Window are deterministic synthetic inputs derived from BED at
	// run start (see deriveInputs). They give the bed* cells whose subcommands
	// need more than BED3 (a name column, or paired/windowed records) a real,
	// honest input instead of an invalid BED3 invocation.
	BED4   string
	BEDPE  string
	Window string
	// FastqPlain is a DECOMPRESSED copy of Fastq1, written at run start (see
	// deriveInputs). prinseq-lite.pl 0.20.4 cannot read gzip (it prints
	// "UNKNOWN format" and produces no output), so the prinseq cells run
	// against this plain FASTQ instead of the bgzipped R1 — that way both
	// bin/ours/prinseq and bin/upstream/prinseq read a format both support and
	// produce a real, comparable result.
	FastqPlain string
	// NameBAM/FixmateBAM/BedGraph/SampleRename are deterministic prerequisite
	// transforms derived at run start (see deriveInputs) so the cells whose
	// UPSTREAM oracle is stricter than ours get an input the oracle accepts,
	// and BOTH sides see byte-identical inputs:
	//   - NameBAM: a name-collated (queryname-grouped) BAM — the input
	//     upstream `samtools fixmate` requires (it errors on coord-sorted).
	//   - FixmateBAM: name-sort | fixmate -m | coord-sort — the input
	//     upstream `samtools markdup` requires (needs ms + MC + coord order).
	//   - BedGraph: a 4-column BedGraph (chrom start end value) — upstream
	//     `bedtools unionbedg` SIGABRTs on a bare BED3.
	//   - SampleRename: a one-line sample-rename map for `bcftools reheader
	//     -s` (a bare `reheader -o` gives no modification directive and
	//     upstream errors with usage).
	NameBAM      string
	FixmateBAM   string
	BedGraph     string
	SampleRename string
	// NormGFF is a csq-normalised copy of GFF written at run start (see
	// deriveInputs): it injects the `biotype=` attribute (from
	// transcript_type/gene_type) that upstream `bcftools csq` requires,
	// mirroring reference_code/bcftools/misc/gff2gff. Upstream exits 255 on
	// a bare GENCODE GFF3 without it; feeding the SAME normalised GFF to
	// both sides lets the csq cell compare ours vs upstream byte-for-byte.
	NormGFF string
}

// have reports whether a single InputKind bit's file is present.
func (in Inputs) have(bit InputKind) bool {
	switch bit {
	case NeedBAM:
		return in.BAM != ""
	case NeedCRAM:
		return in.CRAM != ""
	case NeedRef:
		return in.Ref != ""
	case NeedVCF:
		return in.VCF != ""
	case NeedFastq1:
		return in.Fastq1 != ""
	case NeedFastq2:
		return in.Fastq2 != ""
	case NeedBED:
		return in.BED != ""
	case NeedGFF:
		return in.GFF != ""
	case NeedBED4:
		return in.BED4 != ""
	case NeedBEDPE:
		return in.BEDPE != ""
	case NeedWindow:
		return in.Window != ""
	case NeedFastqPlain:
		return in.FastqPlain != ""
	case NeedNameSortBAM:
		return in.NameBAM != ""
	case NeedFixmateBAM:
		return in.FixmateBAM != ""
	case NeedBedGraph:
		return in.BedGraph != ""
	case NeedSampleRename:
		return in.SampleRename != ""
	case NeedNormGFF:
		return in.NormGFF != ""
	}
	return true
}

// missing returns the first missing required input's name, or "" when all
// required inputs are present.
func (in Inputs) missing(need InputKind) string {
	for bit, name := range map[InputKind]string{
		NeedBAM:          "-bam",
		NeedCRAM:         "-cram",
		NeedRef:          "-ref",
		NeedVCF:          "-vcf",
		NeedFastq1:       "-fastq1",
		NeedFastq2:       "-fastq2",
		NeedBED:          "-bed",
		NeedGFF:          "-gff",
		NeedBED4:         "-bed (bed4)",
		NeedBEDPE:        "-bed (bedpe)",
		NeedWindow:       "-bed (window)",
		NeedFastqPlain:   "-fastq1 (plain)",
		NeedNameSortBAM:  "-bam (name-collated)",
		NeedFixmateBAM:   "-bam (fixmate'd)",
		NeedBedGraph:     "-bed (bedgraph)",
		NeedSampleRename: "-vcf (sample-rename)",
		NeedNormGFF:      "-gff (csq-normalised)",
	} {
		if need&bit != 0 && !in.have(bit) {
			return name
		}
	}
	return ""
}

// substituteArgs replaces every placeholder in args with the resolved path. The
// .fai placeholder resolves to ref+".fai". out/outDir are per-side temp targets.
func substituteArgs(args []string, in Inputs, out, outDir string) []string {
	fai := ""
	if in.Ref != "" {
		fai = in.Ref + ".fai"
	}
	repl := strings.NewReplacer(
		phBAM, in.BAM,
		phCRAM, in.CRAM,
		phRef, in.Ref,
		phFai, fai,
		phVCF, in.VCF,
		phFastq1, in.Fastq1,
		phFastq2, in.Fastq2,
		phFastqPlain, in.FastqPlain,
		phBED, in.BED,
		phGFF, in.GFF,
		phNormGFF, in.NormGFF,
		phBED4, in.BED4,
		phBEDPE, in.BEDPE,
		phWindow, in.Window,
		phNameBAM, in.NameBAM,
		phFixmateBAM, in.FixmateBAM,
		phBedGraph, in.BedGraph,
		phSampleRename, in.SampleRename,
		phOut, out,
		phOutdir, outDir,
	)
	res := make([]string, len(args))
	for i, a := range args {
		res[i] = repl.Replace(a)
	}
	return res
}

// abs resolves p to an absolute path, leaving "" untouched and falling back to
// the original on error.
func abs(p string) string {
	if p == "" {
		return ""
	}
	if a, err := filepath.Abs(p); err == nil {
		return a
	}
	return p
}

// fileExists reports whether p names an existing regular file.
func fileExists(p string) bool {
	if p == "" {
		return false
	}
	st, err := os.Stat(p)
	return err == nil && !st.IsDir()
}
