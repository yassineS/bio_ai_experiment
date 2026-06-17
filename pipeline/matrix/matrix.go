// Package matrix defines the command-combination model that drives the
// parity-and-performance pipeline. A single Entry pairs one of our ported
// tools' invocations with the equivalent upstream invocation, plus how their
// outputs should be compared and whether the run is heavy enough to time.
//
// The model is deliberately small and declarative: later agents populate the
// full per-tool matrices by appending Entry values to the package Registry via
// Register, and use the Expand helper to turn a flag-choice description into a
// curated set of Entry values without hand-writing every combination.
package matrix

import "fmt"

// InputKind names the primary fixture an Entry consumes. The runner uses it to
// pick the correct generated input path from the fixture manifest.
type InputKind string

// The set of fixture kinds the generator produces. They share one coordinate
// space so that, e.g., a BED file and a BAM file refer to the same contigs.
const (
	InputBAM         InputKind = "bam"             // sorted+indexed BAM over the reference
	InputCRAM        InputKind = "cram"            // same alignments as CRAM (+ reference)
	InputVCF         InputKind = "vcf"             // bgzipped+tabixed VCF
	InputVCFPlain    InputKind = "vcf_plain"       // uncompressed VCF
	InputVCFMulti    InputKind = "vcf_multi_plain" // uncompressed multi-sample VCF
	InputBED         InputKind = "bed"             // plain BED3/BED6
	InputBED12       InputKind = "bed12"           // BED12
	InputFASTA       InputKind = "fasta"           // reference FASTA (+ .fai)
	InputFASTQ       InputKind = "fastq"           // single-end FASTQ
	InputFASTQPaired InputKind = "fastq_paired"    // paired-end FASTQ ({fastq1}/{fastq2})
	InputGFF         InputKind = "gff"             // GFF3 annotation
	InputNone        InputKind = ""                // entry supplies its own inputs via Args
)

// CompareMode selects how the runner decides whether our output matches
// upstream's for a given Entry.
type CompareMode string

const (
	// ByteExact requires provenance-stripped stdout/output-file equality.
	// This is the default and strongest mode; use it wherever the upstream
	// output is deterministic and reproducible.
	ByteExact CompareMode = "ByteExact"

	// Similarity is for documented heuristic or non-reproducible paths
	// (adapter auto-detection, unseeded RNG, libm float scores). It compares
	// structure and numeric fields within a tolerance instead of byte-for-byte
	// and records the maximum deviation observed.
	Similarity CompareMode = "Similarity"

	// DirContents compares the set and (stripped) contents of files written to
	// an output directory, rather than a single stdout stream. Used by tools
	// that emit multiple files (e.g. mosdepth, split outputs).
	DirContents CompareMode = "DirContents"

	// BGZFDecoded gunzips each side's stdout before a byte-exact compare. Use it
	// for tools that emit a BGZF/gzip stream (e.g. `bgzip -c`): the compressed
	// bytes differ (our klauspost deflate backend vs htslib framing) but the
	// decompressed content must be identical.
	BGZFDecoded CompareMode = "BGZFDecoded"

	// BAMDecoded decodes each side's stdout BAM through the upstream `samtools
	// view -h` before a provenance-stripped byte-exact compare. Use it for tools
	// that emit BGZF BAM on stdout (e.g. `samtools cat`, `bedtobam`): the binary
	// framing is not comparable, but the decoded SAM records must match.
	BAMDecoded CompareMode = "BAMDecoded"
)

// Entry is one runnable parity case: a single invocation of our tool and the
// equivalent upstream invocation, compared per Compare.
type Entry struct {
	// Tool is our tool name; it must match a tools/<Tool>/cmd/<Tool> binary
	// and is also the lookup key for the upstream binary (see runner.Upstream).
	Tool string

	// Subcommand is the upstream subcommand (e.g. "view", "intersect"), or ""
	// for tools without subcommands. For our binary it is prepended to Args
	// only when UsesSubcommand is true; some of our CLIs are the subcommand
	// itself (e.g. bedintersect == "bedtools intersect").
	Subcommand string

	// UpstreamTool overrides the upstream binary key when our tool name differs
	// from the upstream binary (e.g. our "bedintersect" maps to upstream
	// "bedtools"). Empty means use Tool.
	UpstreamTool string

	// UsesSubcommand reports whether OUR binary also takes Subcommand as its
	// first argument. When false (the common bed* case) only the upstream
	// invocation receives the subcommand token.
	UsesSubcommand bool

	// Name is a short, unique-within-tool label for reports and -run filters.
	Name string

	// Args are the flags/positionals shared by both invocations, with fixture
	// placeholders (see runner.resolvePlaceholders) substituted by the runner.
	// When OurArgs / UpstreamArgs are set they override Args for that side.
	Args []string

	// OurArgs, when non-nil, replaces Args for OUR invocation only. UpstreamArgs
	// does the same for the upstream invocation. These exist for tools whose CLI
	// shape genuinely differs from upstream's (e.g. our skewer is subcommand-
	// based with -i/-o flags while upstream skewer is flat with positionals), so
	// a single shared Args cannot express both sides. The same {placeholder}
	// substitution applies to each. When set, the per-side subcommand prepend
	// still applies (UsesSubcommand for ours, always for upstream).
	OurArgs      []string
	UpstreamArgs []string

	// Input is the primary fixture kind this entry consumes. Placeholders in
	// Args referencing it (e.g. "{bam}") are resolved from the manifest.
	Input InputKind

	// Compare selects the comparison strategy. Defaults to ByteExact when empty.
	Compare CompareMode

	// Heavy marks entries whose timing ratio is worth reporting; the runner
	// always records wall-clock but only surfaces ratios for heavy entries.
	Heavy bool

	// Skip, when non-empty, documents why the entry is intentionally not run
	// (e.g. "upstream segfaults on empty CRAM"); it is reported as SKIP.
	Skip string

	// OutputFiles names the output files an entry writes through an output
	// PREFIX rather than to stdout, given relative to the prefix the {out}
	// placeholder resolves to (e.g. ".frq", ".mosdepth.summary.txt"). When
	// non-empty the runner gives each side its own temp directory, resolves
	// {out} to "<tmpdir>/out", runs the tool, and compares each named output
	// file between the two directories instead of comparing stdout. Files
	// ending in ".gz" are decompressed before comparison so BGZF block framing
	// differences do not cause spurious divergence. This is the mechanism the
	// vcftools and mosdepth matrices use; it pairs with Compare: DirContents.
	OutputFiles []string
}

// CompareModeOrDefault returns the entry's comparison mode, defaulting to
// ByteExact when unset.
func (e Entry) CompareModeOrDefault() CompareMode {
	if e.Compare == "" {
		return ByteExact
	}
	return e.Compare
}

// UpstreamKey returns the binary lookup key for the upstream side.
func (e Entry) UpstreamKey() string {
	if e.UpstreamTool != "" {
		return e.UpstreamTool
	}
	return e.Tool
}

// Registry is an ordered collection of entries that tools register into. The
// driver consumes the global registry (see Default) after all init functions
// have run.
type Registry struct {
	entries []Entry
}

// Add appends entries to the registry, validating each.
func (r *Registry) Add(entries ...Entry) {
	for _, e := range entries {
		if e.Tool == "" {
			panic(fmt.Sprintf("matrix: entry %q has no Tool", e.Name))
		}
		if e.Name == "" {
			panic(fmt.Sprintf("matrix: entry for tool %q has no Name", e.Tool))
		}
		r.entries = append(r.entries, e)
	}
}

// Entries returns a copy of the registered entries.
func (r *Registry) Entries() []Entry {
	out := make([]Entry, len(r.entries))
	copy(out, r.entries)
	return out
}

// FilterTools returns the entries whose Tool is in the allowed set. An empty
// allowed set returns all entries.
func (r *Registry) FilterTools(allowed map[string]bool) []Entry {
	if len(allowed) == 0 {
		return r.Entries()
	}
	var out []Entry
	for _, e := range r.entries {
		if allowed[e.Tool] {
			out = append(out, e)
		}
	}
	return out
}

// defaultRegistry is the global registry that smoke/full matrices register into
// from their package init functions.
var defaultRegistry = &Registry{}

// Default returns the process-wide registry.
func Default() *Registry { return defaultRegistry }

// Register adds entries to the global default registry. Matrix-definition
// packages call this from init().
func Register(entries ...Entry) { defaultRegistry.Add(entries...) }
