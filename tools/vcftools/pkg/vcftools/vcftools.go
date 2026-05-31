// Package vcftools provides utilities for working with VCF files
package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Params holds all parameters for vcftools operations
type Params struct {
	// Output options
	OutPrefix string
	UseStdout bool
	Recode    bool
	// RecodeBCF emits the recoded output as BGZF-compressed BCF v2.2 to
	// `<prefix>.recode.bcf` (or stdout under `--stdout` / `-c`), mirroring
	// upstream's `--recode-bcf` (parameters.cpp:317 → vcf_file.cpp:119).
	// May be combined with `Recode`; both files are written independently.
	RecodeBCF     bool
	RecodeInfoAll bool

	// BCFInputFile, when non-empty, makes `Run` ignore its `input
	// io.Reader` argument and read the BCF at the given path instead.
	// The CLI passes this through from `--bcf FILE` (parameters.cpp:173).
	// The BCF stream is BGZF-decompressed and decoded via
	// `pkg/htsgo/bcf`; each record is converted to `vcf.Variant`
	// so the entire downstream filter pipeline works unchanged.
	BCFInputFile string

	// ContigsFile, when non-empty, supplies supplemental ##contig=
	// lines for BCF header construction in --recode-bcf when the
	// source VCF lacks contig declarations of its own. Mirrors
	// upstream's `--contigs FILE` (parameters.cpp:197 →
	// variant_file.cpp:45-69). The CLI parses the file once and
	// prepends the lines to the header MetaInfo before the BCF
	// writer is constructed.
	ContigsFile string

	// Position filtering
	Chr                  string
	NotChr               string
	FromBp               int
	ToBp                 int
	PositionsFile        string
	ExcludePositionsFile string

	// --positions-overlap FILE / --exclude-positions-overlap FILE keep or
	// drop a site when ANY base of its reference allele falls on a position
	// listed in the file. The point-positions variants (--positions /
	// --exclude-positions) above test only the leading POS column. Multi-
	// base REF records (indels, MNPs) therefore behave differently between
	// the two flag pairs: a record at POS=200 with REF="ACGT" matches a
	// positions-overlap entry at any of 200..203, but matches a plain
	// positions entry only at 200. Ported from upstream
	// parameters.cpp:221,315 + entry_filters.cpp:408-548
	// (filter_sites_by_overlap_positions). File format is the same
	// two-column whitespace-separated (CHROM, POS) text used by --positions,
	// with `#` comments and blank lines tolerated.
	PositionsOverlapFile        string
	ExcludePositionsOverlapFile string

	// SNP ID filtering
	SNP         string
	SNPs        string
	ExcludeSNP  string
	ExcludeSNPs string
	Thin        int

	// Variant type filtering
	KeepOnlyIndels bool
	RemoveIndels   bool
	MinAlleles     int
	MaxAlleles     int

	// Quality filtering
	MinQ              float64
	RemoveFilteredAll bool

	// Allele frequency filtering
	Maf    float64
	MaxMaf float64
	Mac    int
	MaxMac int

	// Non-reference (per-ALT) allele frequency / count filtering. Upstream
	// vcftools registers eight flags here (parameters.cpp:287-290, 302-305):
	//
	//   --non-ref-af / --max-non-ref-af           — plain min/max AF
	//   --non-ref-af-any / --max-non-ref-af-any   — "any ALT passes" AF
	//   --non-ref-ac / --max-non-ref-ac           — plain min/max AC
	//   --non-ref-ac-any / --max-non-ref-ac-any   — "any ALT passes" AC
	//
	// Plain semantics (entry_filters.cpp:807, 905): site dropped if ANY ALT
	// falls outside [min, max]. "-any" semantics (lines 810, 908): N_failed
	// counts ALTs outside the threshold; site dropped if N_failed equals
	// (N_alleles - 1), i.e. EVERY ALT failed.
	//
	// Two important upstream asymmetries that we preserve verbatim:
	//
	//  1. AF fallback gate (line 814) is keyed on the PLAIN thresholds, not
	//     the -any thresholds. Consequence: `--non-ref-af-any` and
	//     `--max-non-ref-af-any` are effectively NO-OPS unless the plain
	//     flag is also set. We still wire and validate the flags for
	//     command-line parity, but mirror upstream's degenerate behaviour.
	//  2. AC fallback gate (line 912) is keyed on the -any thresholds.
	//     Consequence: plain `--non-ref-ac` / `--max-non-ref-ac` do NOT drop
	//     monomorphic (N_alleles == 1) sites, whereas the -any variants do
	//     (N_failed == 0 == N_alleles-1 fires the fallback). The plain AF
	//     `--non-ref-af` / `--max-non-ref-af` do drop monomorphic sites via
	//     the same N_failed == N_alleles-1 mechanism.
	//
	// All fields default to 0; the Max* defaults of 0 mean "no upper bound".
	// We adopt 0 (instead of upstream's sentinel +inf / INT_MAX) because
	// none of upstream's flags actually require setting max <= 0 to filter
	// — a 0-max ALT count or AF is degenerate. Internally we treat
	// MaxNonRefAF/MaxNonRefAC == 0 as "unset". A nonzero value applies
	// `freq > Max` / `count > Max`.
	MinNonRefAF    float64
	MinNonRefAC    int
	MaxNonRefAF    float64
	MaxNonRefAC    int
	MinNonRefAFAny float64
	MinNonRefACAny int
	MaxNonRefAFAny float64
	MaxNonRefACAny int

	// Genotype filtering
	MaxMissing float64
	// MaxMissingCount is the maximum number of missing chromosomes (i.e.
	// missing haploid alleles, NOT missing samples) tolerated per site.
	// Default 0 means "feature disabled". Mirrors upstream's
	// `--max-missing-count` (parameters.cpp:286 + entry_filters.cpp:918);
	// the upstream check is `(N_chr - N_non_missing_chr) > max`, i.e.
	// values are absolute allele counts. A diploid "./." sample
	// contributes 2; "0/." contributes 1. Set via --max-missing-count.
	// Note: upstream's default sentinel is INT_MAX. We use 0 (== disabled)
	// because a zero threshold is the only useful "drop any site with any
	// missing allele" semantics — for that use `--max-missing-count 0`
	// explicitly which sets MaxMissingCount = -1 (we re-interpret -1 as
	// "active threshold of zero" below to disambiguate from default).
	MaxMissingCount int
	// MaxMissingCountSet records whether --max-missing-count was supplied
	// on the command line. Without this flag we cannot distinguish
	// "user passed 0" (meaning: drop every site with any missing call)
	// from "user did not pass the flag at all" (meaning: no filter).
	MaxMissingCountSet bool
	MinMeanDP          float64
	MaxMeanDP          float64
	MinDP              int
	MaxDP              int
	MinGQ              int

	// MinHWEPvalue is the minimum exact-test p-value (Wigginton et al.
	// 2005) for Hardy-Weinberg equilibrium that a biallelic site must
	// pass. Default 0 means "feature disabled". Set via --hwe FLOAT.
	// Upstream (parameters.cpp:254) also forces `max_alleles = 2` when
	// `--hwe` is supplied; we mirror that in CheckUnsupported / the CLI
	// wiring so `--hwe 0.05` on a multi-allelic site drops the site
	// before the HWE test even runs.
	MinHWEPvalue float64

	// Statistics output
	Freq          bool
	Counts        bool
	Freq2         bool
	Counts2       bool
	Depth         bool
	SiteDepth     bool
	SiteMeanDepth bool
	SiteQuality   bool
	MissingIndv   bool
	MissingSite   bool
	Hardy         bool
	TsTvSummary   bool
	TsTvBinSize   int
	TsTvByCount   bool
	TsTvByQual    bool
	SitePi        bool
	Het           bool
	Singletons    bool
	HistIndelLen  bool
	GenoDepth     bool

	// Phase 2: Population genetics statistics
	WindowPi      int
	WindowPiStep  int
	TajimaD       int
	SNPDensity    int
	WeirFstPop    []string
	FstWindowSize int
	FstWindowStep int
	FilterSummary bool

	// Phase 4: Format conversions
	Output012       bool
	OutputPlink     bool
	OutputPlinkTped bool
	ChromMap        string

	// Sample filtering
	IndvList       []string
	RemoveIndvList []string
	KeepFile       string
	RemoveFile     string

	// --max-indv N caps the number of kept individuals at N. Ported from
	// upstream parameters.cpp:292 + variant_file_filters.cpp:105-147
	// (filter_individuals_randomly). MaxIndvSet records whether the flag
	// was supplied (since "0" is meaningful — drop every sample). Upstream
	// uses srand(time(NULL)) + random_shuffle so the kept-sample identity
	// is non-deterministic across runs; this port instead deterministically
	// keeps the first MaxIndv samples in input order. Tests therefore pin
	// the COUNT, not the identity (the byte-parity claim is on
	// |kept samples| = min(MaxIndv, |kept|) rather than which N samples
	// were chosen). Documented in docs/PARITY_ROADMAP.md (wave 10).
	MaxIndv    int
	MaxIndvSet bool

	// --remove-filtered-geno-all sets a per-genotype filter that drops
	// (sets GT to ./.) any genotype whose FORMAT FT field is not "PASS"
	// or ".". Ported from upstream parameters.cpp:323 +
	// vcf_entry.cpp:580-608 (filter_genotypes_by_filter_status).
	// --remove-filtered-geno NAME (repeatable) drops a genotype whose
	// FORMAT FT lists any of the named flags. Both compose: --all wins
	// (when set, the named-list path is irrelevant).
	RemoveFilteredGenoAll  bool
	RemoveFilteredGenoList []string

	// Linkage disequilibrium analysis (--geno-r2 / --hap-r2 family).
	// GenoR2 enables --geno-r2 output (<prefix>.geno.ld). HapR2 enables
	// --hap-r2 output (<prefix>.hap.ld). GenoR2Positions / HapR2Positions
	// supply chrom/pos files that restrict pairs to those where at least one
	// endpoint is listed (analogous to upstream --geno-r2-positions /
	// --hap-r2-positions). The four window fields bound the pairwise SNP /
	// bp distance: zero means "no bound" for the maxima, zero for the minima
	// means "no minimum required". MinR2 thresholds the emitted r² (default
	// 0 = emit all pairs).
	GenoR2          bool
	HapR2           bool
	GenoR2Positions string
	HapR2Positions  string
	LDWindow        int
	LDWindowBp      int
	LDWindowMin     int
	LDWindowBpMin   int
	MinR2           float64

	// BED-based site filtering. Bed keeps only sites whose 1-based POS lies
	// inside any interval (0-based half-open) in the file. ExcludeBed is the
	// inverse. Both compose with other position/quality filters.
	Bed        string
	ExcludeBed string

	// FASTA-like positional mask. Mask names a file with `>CHROM` headers
	// followed by digit lines (one per reference base). A site at (CHROM,
	// POS) is kept when its mask digit is <= MaskMin (default 0).
	// InvertMask flips keep/drop. Ported from upstream
	// parameters.cpp:262/279/280 + entry_filters.cpp:674-752 (wave 11).
	// The mask reader is forward-only — VCF records reordered relative to
	// the mask's chromosome order may be dropped (this matches upstream
	// because of its streaming ifstream state). See mask_filter.go for the
	// full algorithm and parity notes.
	Mask       string
	InvertMask bool
	MaskMin    int

	// VCF comparison (--diff family). Diff names the second VCF to compare
	// against; the boolean flags request individual output files. See
	// diff.go for the column layout of each output. DiffIndvMap is a path to
	// a two-column whitespace-separated file that renames file-2 sample IDs
	// to their file-1 equivalents (upstream variant_file_diff.cpp:11-34).
	// DiffDiscordanceMatrix emits the 4x4 genotype-by-genotype counts in
	// <prefix>.diff.discordance_matrix (upstream variant_file_diff.cpp:944).
	Diff string
	// DiffBCF is the path to a BCF-formatted second file for the
	// --diff-* family, mirroring upstream's `--diff-bcf FILE`
	// (parameters.cpp:210). Mutually exclusive with Diff at the CLI
	// layer; when set, the diff loader reads via the shared
	// `pkg/htsgo/bcf` reader stack.
	DiffBCF               string
	DiffSite              bool
	DiffIndv              bool
	DiffSiteDiscordance   bool
	DiffIndvDiscordance   bool
	DiffIndvMap           string
	DiffDiscordanceMatrix bool

	// DiffSwitchError enables --diff-switch-error: per-event log
	// <prefix>.diff.switch + per-individual summary
	// <prefix>.diff.indv.switch. Ported from upstream
	// variant_file_diff.cpp:1207 (output_switch_error). See
	// diff_switch.go for the file layouts.
	DiffSwitchError bool

	// MendelPedFile is the path supplied to --mendel <PED>. When non-empty,
	// vcftools writes <prefix>.mendel listing Mendelian inconsistencies
	// across trios defined in the PED file. Ported from upstream
	// variant_file_output.cpp:5332 (output_mendel_inconsistencies); the
	// PED columns are `family child father mother` (with a header line
	// that is always skipped). See mendel.go for the column layout.
	MendelPedFile string

	// BEAGLE genotype-likelihood output. BEAGLEGL writes log10-scale GL
	// triplets derived from FORMAT/PL; BEAGLEPL writes the raw PL triplets.
	// Both are biallelic-SNP only.
	BEAGLEGL bool
	BEAGLEPL bool

	// Inter-chromosomal LD outputs. InterchromGenoR2 and InterchromHapR2
	// emit `<prefix>.interchrom.geno.ld` / `<prefix>.interchrom.hap.ld`
	// (only cross-chromosome pairs). GenoChiSq emits `<prefix>.geno.chisq`
	// (per-pair Pearson chi-square test across all pairs, same- and
	// cross-chromosome).
	InterchromGenoR2 bool
	InterchromHapR2  bool
	GenoChiSq        bool

	// Relatedness statistics. Relatedness enables <prefix>.relatedness
	// (Yang 2010 unadjusted A_jk). Relatedness2 enables
	// <prefix>.relatedness2 (KING-robust kinship; Manichaikul 2010).
	Relatedness  bool
	Relatedness2 bool

	// PhasedBlocks enables <prefix>.blocks reporting per-individual
	// contiguous runs of phased ("a|b") diploid genotypes.
	PhasedBlocks bool

	// Runs of homozygosity. LROH enables <prefix>.LROH. LROHMinVariants is
	// the minimum number of consecutive homozygous variants for a run to
	// be emitted (default 10 when LROH is true and the value is zero).
	LROH            bool
	LROHMinVariants int

	// Filter-name include/exclude (operate on the FILTER column). Both are
	// comma-separated lists. RemoveFiltered drops sites whose FILTER lists
	// any of the named tags; KeepFiltered keeps only sites that list at
	// least one of the named tags.
	RemoveFiltered string
	KeepFiltered   string

	// KeepINFO is a comma-separated list of INFO Flag-type tag names. It
	// is a SITE FILTER (NOT a recode-column selector): a site passes only
	// if at least one of the listed flag tags is present in the variant's
	// INFO field. Mirrors upstream parameters.cpp:266
	// (`site_INFO_flags_to_keep`) and the filter algorithm at
	// entry_filters.cpp:1033-1063. Upstream errors out if a named tag is
	// not declared as Type=Flag in the header; this port preserves that
	// behaviour. Multiple tags compose via OR ("any present").
	//
	// To restrict the INFO column in `.recode.vcf` output to a subset of
	// tags, use `RecodeINFO` (upstream `--recode-INFO`,
	// parameters.cpp:319 → `recode_INFO_to_keep`) instead. Pre-wave-17
	// versions of this port misaligned `KeepINFO` with that recode-column
	// semantic; see docs/UPSTREAM_BUGS.md for the migration note.
	KeepINFO string

	// RemoveINFO is a comma-separated list of INFO Flag-type tag names. It
	// is a SITE FILTER: a site is DROPPED if any of the listed flag tags
	// is present in the variant's INFO field. Mirrors upstream
	// parameters.cpp:328 (`site_INFO_flags_to_remove`) and the filter
	// algorithm at entry_filters.cpp:1068-1086. Upstream errors out if a
	// named tag is not declared as Type=Flag in the header; this port
	// preserves that behaviour. Multiple tags compose via OR ("any
	// present drops").
	//
	// Composes with `KeepINFO` per upstream's filter_sites_by_INFO
	// ordering (entry_filters.cpp:1033-1086): keep narrows first, remove
	// vetoes the survivors.
	//
	// Pre-wave-18 versions of this port wired `RemoveINFO` as a
	// recode-column stripper (a port-only invention with no upstream
	// equivalent); wave 18 (see docs/UPSTREAM_BUGS.md) repointed it at
	// the upstream site-filter semantic.
	RemoveINFO string

	// RecodeINFO is a comma-separated list of INFO tag names to retain in
	// the recoded `.recode.vcf` output. Mirrors upstream's
	// `--recode-INFO TAG` (parameters.cpp:319 → `recode_INFO_to_keep`).
	// Composes with `RecodeInfoAll` (which preserves all INFO).
	RecodeINFO string

	// --get-INFO TAG[,TAG]... extracts the named INFO tags as a TSV file
	// <prefix>.INFO with columns CHROM POS REF ALT <tags...>. The flag is
	// comma-separated; the upstream CLI accepts the same value-style.
	GetINFO string

	// --phased filters sites to those where every kept-individual GT is
	// phased (separator '|', or upstream's haploid case which is treated
	// as phased). Mirrors upstream's parameters.cpp:311 / entry_filters.cpp
	// filter_sites_by_phase (lines 989-1010).
	Phased bool

	// --ldhat / --ldhat-geno emit a paired (<prefix>.ldhat.sites,
	// <prefix>.ldhat.locs) bundle in LDhat's input format. Both require
	// exactly one chromosome to be selected via --chr (upstream errors
	// otherwise: parameters.cpp:717). --ldhat additionally implies
	// --phased (the phased-only site filter); --ldhat-geno does not.
	// See ldhat.go for the file layout.
	LDhat     bool
	LDhatGeno bool

	// --ldhelmet emits a paired (<prefix>.ldhelmet.snps,
	// <prefix>.ldhelmet.pos) bundle in LDhelmet's input format. Upstream
	// (parameters.cpp:275) sets phased_only = true and remove_indels =
	// true; --ldhelmet also shares the LDhat --chr requirement
	// (parameters.cpp:717). See ldhelmet.go for the file layout.
	LDhelmet bool

	// --IMPUTE (case-sensitive) emits a three-file bundle in IMPUTE
	// reference-panel format:
	//   <prefix>.impute.legend, <prefix>.impute.hap, <prefix>.impute.hap.indv
	// Upstream (parameters.cpp:255) sets phased_only = true,
	// min_site_call_rate = 1.0, min_alleles = 2, max_alleles = 2. See
	// impute.go for the file layout.
	IMPUTE bool

	// PCA / PCANoNorm / PCASNPLoadings hold the upstream `--pca`,
	// `--pca-no-norm`, and `--pca-snp-loadings INT` flags
	// (parameters.cpp:308-310). Wave-19 ports the eigendecomposition via
	// gonum (`gonum.org/v1/gonum/mat`'s symmetric eigensolver), the
	// project's second sanctioned third-party-dep zone after the CRAM
	// codec carveout (see CLAUDE.md). PCA writes `<prefix>.pca`;
	// PCASNPLoadings > 0 additionally writes `<prefix>.pca.loadings`
	// with the first K SNP loadings. Setting PCANoNorm disables
	// per-SNP variance normalisation (still mean-centres). See pca.go
	// for the algorithm and docs/UPSTREAM_BUGS.md for the missing-data
	// fix-on-port.
	PCA            bool
	PCANoNorm      bool
	PCASNPLoadings int

	// KeptSites enables --kept-sites: writes `<prefix>.kept.sites`, a
	// two-column (CHROM, POS) TSV listing every site that survived all
	// filters. Mirrors upstream parameters.cpp:268 +
	// variant_file_output.cpp:4285-4326 (output_kept_sites). The file
	// has a `CHROM\tPOS` header followed by one row per kept site, in
	// input order. Composes with every other filter — the contents
	// match the sites that --recode would emit.
	KeptSites bool

	// RemovedSites enables --removed-sites: writes `<prefix>.removed.sites`,
	// a two-column (CHROM, POS) TSV listing every site that was dropped
	// by any filter. Mirrors upstream parameters.cpp:330 +
	// variant_file_output.cpp:4328-4373 (output_removed_sites). Header is
	// `CHROM\tPOS`; rows appear in input order. With no filters set,
	// the file contains only the header (every site is kept).
	RemovedSites bool

	// Derived enables --derived: when --freq / --counts is active, the
	// allele columns in <prefix>.frq / <prefix>.frq.count are reordered
	// so the ancestral allele (INFO/AA, case-insensitive) appears first.
	// Sites where AA is missing, ".", "?", or does not match REF/ALT are
	// dropped (upstream prints a one-off warning and continues). Mirrors
	// upstream parameters.cpp:201 + variant_file_output.cpp:67-159
	// (output_frequency). Multi-allelic sites are already dropped by our
	// existing --freq biallelic restriction, so --derived only affects
	// the biallelic subset (matches the subset the port emits at all).
	Derived bool

	// ExtractFormatInfo names the FORMAT field to extract per-genotype
	// across all kept samples, writing a tab-separated file
	// `<prefix>.<NAME>.FORMAT`. Sites whose FORMAT column does not list
	// NAME are skipped. For samples whose colon-separated value vector
	// is too short to contain NAME, a literal "." is emitted (matches
	// upstream `read_indv_generic_entry` in vcf_entry.cpp:610-639).
	// Ported from upstream parameters.cpp:222 +
	// variant_file_format_convert.cpp:1204-1263
	// (output_FORMAT_information).
	ExtractFormatInfo string

	// IndvBurden enables --indv-burden: writes <prefix>.iburden with one
	// row per kept individual giving the per-individual count of
	// hom-ref / het / hom-alt / missing diploid genotypes across the
	// kept sites. With --derived the "ref" / "alt" columns become "anc"
	// / "der" and sites without a resolvable AA are skipped. Mirrors
	// upstream parameters.cpp:257 + variant_file_output.cpp:378-498
	// (output_indv_burden).
	IndvBurden bool

	// IndvFreqBurden enables --indv-freq-burden: writes
	// <prefix>.ifreqburden, a per-individual × per-allele-count
	// matrix where cell [i, j] is the number of non-reference (or
	// non-ancestral, with --derived) alleles that individual i carries
	// at sites where the chosen allele's per-site count across the
	// kept individuals is j. Mirrors upstream parameters.cpp:258 +
	// variant_file_output.cpp:501-627 (output_indv_freq_burden with
	// double_count_hom_alt=0).
	IndvFreqBurden bool

	// IndvFreqBurden2 enables --indv-freq-burden2: same as
	// --indv-freq-burden but with double_count_hom_alt=1 — a hom-alt
	// genotype contributes 1 (not 2) to the per-allele-count bin.
	// Upstream and the port both write to <prefix>.ifreqburden, so
	// supplying both is meaningless (only one wins). Mirrors upstream
	// parameters.cpp:259 + variant_file_output.cpp:501-627.
	IndvFreqBurden2 bool

	// HapcountBED names the BED file passed to `--hapcount`. When
	// non-empty, the runner writes `<prefix>.hapcount` with per-BED-bin
	// haplotype-multiplicity tables. Implies `--phased` (upstream
	// parameters.cpp:248 sets `phased_only = true`). Three upstream
	// defects in `output_haplotype_count` are fixed on port; see
	// hapcount.go and docs/UPSTREAM_BUGS.md for the writeup.
	HapcountBED string

	// TempDir mirrors upstream's `--temp DIR` flag
	// (parameters.cpp:341). Upstream uses it as the base directory for
	// `mkstemp` spill files in the LD / format-convert paths
	// (variant_file_output.cpp:1441, variant_file_format_convert.cpp).
	// This port does not spill to disk for any of those paths, so the
	// flag is accepted-and-ignored for command-line parity. When
	// non-empty, a stderr line is emitted noting the flag was parsed
	// but unused. Documented in docs/PARITY_ROADMAP.md#vcftools.
	TempDir string
}

// positionSet represents a set of positions to include/exclude
type positionSet map[string]map[int]bool

// variantSource abstracts over `*vcf.Reader` and a BCF-backed adapter so
// `Run` can iterate variants from either a text VCF or a binary BCF via
// the same loop. The interface deliberately mirrors the small subset of
// `*vcf.Reader` that Run uses.
type variantSource interface {
	Header() *vcf.Header
	Read() (*vcf.Variant, error)
	Close() error
}

// vcfVariantSource is the trivial adapter over a `*vcf.Reader`. The
// reader itself doesn't own an underlying file (Run's caller does),
// so Close is a no-op.
type vcfVariantSource struct {
	r   *vcf.Reader
	hdr *vcf.Header
}

func (s *vcfVariantSource) Header() *vcf.Header         { return s.hdr }
func (s *vcfVariantSource) Read() (*vcf.Variant, error) { return s.r.Read() }
func (s *vcfVariantSource) Close() error                { return nil }

// bcfVariantSource adapts a `*bcf.Reader` to the variantSource
// interface. It carries `closers` for the file handle (and any
// intermediate BGZF reader) so the caller can clean up after the
// last Read.
type bcfVariantSource struct {
	r       *bcf.Reader
	hdr     *vcf.Header
	closers []io.Closer
}

func (s *bcfVariantSource) Header() *vcf.Header { return s.hdr }
func (s *bcfVariantSource) Read() (*vcf.Variant, error) {
	rec, err := s.r.Read()
	if err != nil {
		return nil, err
	}
	return rec.ToVariant(s.r.Header()), nil
}
func (s *bcfVariantSource) Close() error {
	// closers is stored in inner-to-outer order ([bz, f]): bz reads
	// through f, so bz must close first to flush/release its own
	// state before the underlying file descriptor goes away.
	var firstErr error
	for i := 0; i < len(s.closers); i++ {
		if err := s.closers[i].Close(); err != nil && firstErr == nil {
			firstErr = err
		}
	}
	return firstErr
}

// augmentHeaderContigs reads a `--contigs FILE` and prepends `##contig=`
// MetaInfo lines to header when it has none of its own. Mirrors upstream
// `get_contigs` (variant_file.cpp:45-69) which is only consulted from
// `vcf_file::print_bcf` when `meta_data.has_contigs == false`.
//
// Two accepted line formats per upstream:
//
//   - A bare contig identifier on its own line, e.g. `chr1` — turned into
//     `##contig=<ID=chr1>`.
//   - A full meta-info form starting with `##contig=<` — used verbatim.
//
// Blank lines and lines starting with `#` (i.e. neither bare names nor
// the meta-info form) are skipped silently. Upstream is similarly
// forgiving (line 60).
func augmentHeaderContigs(hdr *vcf.Header, path string) error {
	hasContig := false
	for _, line := range hdr.MetaInfo {
		if strings.HasPrefix(line, "##contig=") {
			hasContig = true
			break
		}
	}
	if hasContig {
		// Upstream gates the file read on `has_contigs == false`. We do
		// the same so a stale `--contigs FILE` paired with a VCF that
		// already declares contigs is silently ignored, matching
		// upstream behaviour.
		return nil
	}

	f, err := iohelper.OpenReader(path)
	if err != nil {
		return fmt.Errorf("opening --contigs %s: %w", path, err)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 64*1024), 1024*1024)
	var added []string
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		switch {
		case strings.HasPrefix(line, "##contig=<"):
			added = append(added, line)
		case strings.HasPrefix(line, "#"):
			// Comment or unexpected meta line — skip.
			continue
		default:
			// Bare contig name. Wrap it in the canonical structured form.
			added = append(added, "##contig=<ID="+line+">")
		}
	}
	if err := scanner.Err(); err != nil {
		return fmt.Errorf("reading --contigs %s: %w", path, err)
	}
	if len(added) == 0 {
		return nil
	}
	// Insert the contig lines immediately AFTER the ##fileformat line.
	// VCFv4.2 mandates ##fileformat be the first line of the header;
	// putting contigs ahead of it produces a spec-violating header
	// that stricter parsers (e.g. recent htslib) reject. If no
	// ##fileformat line is present we append, mirroring upstream's
	// "tack on after meta_data.lines" pattern (variant_file.cpp:154).
	insertAt := -1
	for i, line := range hdr.MetaInfo {
		if strings.HasPrefix(line, "##fileformat=") {
			insertAt = i + 1
			break
		}
	}
	if insertAt < 0 {
		hdr.MetaInfo = append(hdr.MetaInfo, added...)
		return nil
	}
	merged := make([]string, 0, len(hdr.MetaInfo)+len(added))
	merged = append(merged, hdr.MetaInfo[:insertAt]...)
	merged = append(merged, added...)
	merged = append(merged, hdr.MetaInfo[insertAt:]...)
	hdr.MetaInfo = merged
	return nil
}

// newBCFVariantSource opens the BCF at path, BGZF-decompresses, and
// returns a variantSource over the resulting record stream. The
// returned source owns the file handle and BGZF reader; Close
// releases bz before f to match the documented inner-first ordering.
//
// Error messages use the bare path (no flag prefix) so callers that
// surface them under different flag names (e.g. --bcf vs --diff-bcf)
// can prepend their own context without producing
// "opening --diff-bcf X: opening --bcf X: …" chains.
func newBCFVariantSource(path string) (variantSource, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("opening BCF %s: %w", path, err)
	}
	bz, err := bgzip.NewReader(f)
	if err != nil {
		_ = f.Close()
		return nil, fmt.Errorf("BGZF reader on BCF %s: %w", path, err)
	}
	r, err := bcf.NewReader(bz)
	if err != nil {
		_ = bz.Close()
		_ = f.Close()
		return nil, fmt.Errorf("BCF reader on %s: %w", path, err)
	}
	return &bcfVariantSource{r: r, hdr: r.Header().VCF, closers: []io.Closer{bz, f}}, nil
}

// Run executes vcftools with the given parameters
func Run(input io.Reader, params *Params) error {
	// Reject requested features that this port does not implement yet, instead
	// of silently producing no output.
	if err := checkUnsupported(params); err != nil {
		return err
	}

	// Open the input as either a VCF stream (the caller-supplied reader)
	// or a BCF file (when --bcf is set; the io.Reader argument is then
	// ignored to match upstream's mutually-exclusive --vcf / --bcf
	// flags).
	var src variantSource
	if params.BCFInputFile != "" {
		s, err := newBCFVariantSource(params.BCFInputFile)
		if err != nil {
			return err
		}
		src = s
	} else {
		vr := vcf.NewReader(input)
		hdr, err := vr.ReadHeader()
		if err != nil {
			return fmt.Errorf("reading VCF header: %w", err)
		}
		src = &vcfVariantSource{r: vr, hdr: hdr}
	}
	defer src.Close()
	reader := src
	header := src.Header()
	// `--contigs FILE` augments the header with `##contig=<...>` lines
	// when the input lacks any of its own (upstream variant_file.cpp:45-69
	// + vcf_file.cpp:154). Only consulted before downstream code that
	// might want contig declarations (BCF writer, ##contig-aware
	// filters).
	if params.ContigsFile != "" {
		if err := augmentHeaderContigs(header, params.ContigsFile); err != nil {
			return err
		}
	}
	// err is declared here so sub-blocks below can reuse the variable
	// via `var x, err = …` without redeclaring at every call site.
	var err error

	// Load position filters if needed
	var includePositions, excludePositions positionSet
	if params.PositionsFile != "" {
		includePositions, err = loadPositions(params.PositionsFile)
		if err != nil {
			return fmt.Errorf("loading positions file: %w", err)
		}
	}
	if params.ExcludePositionsFile != "" {
		excludePositions, err = loadPositions(params.ExcludePositionsFile)
		if err != nil {
			return fmt.Errorf("loading exclude positions file: %w", err)
		}
	}

	// Load overlap position filters. Same file format as --positions but
	// the match check sweeps every base in [POS, POS+len(REF)-1] against the
	// set — see entry_filters.cpp:513-547. Stored separately from the plain
	// include/exclude sets so the two flag pairs can coexist (the upstream
	// implementation reuses one map and consequently has a documented
	// quirk where combining them silently degrades to overlap semantics for
	// whichever was supplied second; we keep the two filters independent).
	var includePositionsOverlap, excludePositionsOverlap positionSet
	if params.PositionsOverlapFile != "" {
		includePositionsOverlap, err = loadPositions(params.PositionsOverlapFile)
		if err != nil {
			return fmt.Errorf("loading positions-overlap file: %w", err)
		}
	}
	if params.ExcludePositionsOverlapFile != "" {
		excludePositionsOverlap, err = loadPositions(params.ExcludePositionsOverlapFile)
		if err != nil {
			return fmt.Errorf("loading exclude-positions-overlap file: %w", err)
		}
	}

	// Load BED-based filters (--bed / --exclude-bed). Both are optional;
	// supplying both composes them (a site must pass include AND not be
	// excluded).
	var includeBed, excludeBed *bedRegions
	if params.Bed != "" {
		includeBed, err = loadBedRegions(params.Bed)
		if err != nil {
			return fmt.Errorf("loading --bed file: %w", err)
		}
	}
	if params.ExcludeBed != "" {
		excludeBed, err = loadBedRegions(params.ExcludeBed)
		if err != nil {
			return fmt.Errorf("loading --exclude-bed file: %w", err)
		}
	}

	// Load --mask / --invert-mask FASTA-style positional mask. We surface
	// load-time errors (missing file, --mask-min out of range) before we
	// stream the VCF. The mask filter is stateful and forward-only; see
	// mask_filter.go for the parity notes.
	var mask *maskFilter
	if params.Mask != "" {
		mask, err = loadMaskFilter(params.Mask, params.InvertMask, params.MaskMin)
		if err != nil {
			return fmt.Errorf("loading --mask file: %w", err)
		}
	}

	// Load --weir-fst-pop population files if requested. We validate here so
	// that errors (missing file, sample appearing in multiple populations,
	// fewer than 2 populations) are surfaced before we start streaming the
	// VCF.
	var weirFstPops [][]string
	if len(params.WeirFstPop) > 0 {
		weirFstPops, err = loadPopulationFiles(params.WeirFstPop)
		if err != nil {
			return fmt.Errorf("loading --weir-fst-pop files: %w", err)
		}
	}

	// Load SNP ID filters
	var includeSNPs, excludeSNPs map[string]bool
	if params.SNP != "" {
		includeSNPs = make(map[string]bool)
		includeSNPs[params.SNP] = true
	}
	if params.SNPs != "" {
		includeSNPs, err = loadSNPIDs(params.SNPs)
		if err != nil {
			return fmt.Errorf("loading SNPs file: %w", err)
		}
	}
	if params.ExcludeSNP != "" {
		excludeSNPs = make(map[string]bool)
		excludeSNPs[params.ExcludeSNP] = true
	}
	if params.ExcludeSNPs != "" {
		excludeSNPs, err = loadSNPIDs(params.ExcludeSNPs)
		if err != nil {
			return fmt.Errorf("loading exclude SNPs file: %w", err)
		}
	}

	// Build sample filter set
	keepSamples, err := buildSampleFilter(header, params)
	if err != nil {
		return fmt.Errorf("building sample filter: %w", err)
	}

	// Filter header samples if needed
	filteredHeader := filterHeaderSamples(header, keepSamples)

	// Initialize statistics
	stats := newStatistics(filteredHeader)
	if len(weirFstPops) >= 2 {
		stats.weirFst = newWeirFstAccumulator(weirFstPops)
	}

	// Initialise LD runner (no-op when no LD flag is set).
	var ldRun *ldRunner
	if params.GenoR2 || params.HapR2 || params.GenoR2Positions != "" || params.HapR2Positions != "" {
		ldRun, err = newLDRunner(params)
		if err != nil {
			return fmt.Errorf("initialising LD analysis: %w", err)
		}
	}

	// Initialise --diff runner (no-op when --diff isn't set or no diff
	// sub-output is requested).
	diffRun, err := newDiffRunner(params, filteredHeader.Samples)
	if err != nil {
		return fmt.Errorf("initialising --diff analysis: %w", err)
	}

	// Initialise BEAGLE output writers. They're created lazily here so a
	// failure to open the file surfaces before we stream any variants.
	var beagleGL, beaglePL *beagleWriter
	if params.BEAGLEGL {
		beagleGL, err = newBEAGLEWriter(params.OutPrefix, beagleGLMode())
		if err != nil {
			return fmt.Errorf("initialising --BEAGLE-GL: %w", err)
		}
	}
	if params.BEAGLEPL {
		beaglePL, err = newBEAGLEWriter(params.OutPrefix, beaglePLMode())
		if err != nil {
			return fmt.Errorf("initialising --BEAGLE-PL: %w", err)
		}
	}

	// Inter-chromosomal LD / chi-square buffer (all-pairs after streaming).
	var interLD *interchromLDRunner
	if params.InterchromGenoR2 || params.InterchromHapR2 || params.GenoChiSq {
		interLD, err = newInterchromLDRunner(params)
		if err != nil {
			return fmt.Errorf("initialising interchrom LD: %w", err)
		}
	}

	// --relatedness accumulator.
	var rel *relatednessRunner
	if params.Relatedness {
		rel = newRelatednessRunner(filteredHeader.Samples)
	}

	// --relatedness2 accumulator.
	var rel2 *relatedness2Runner
	if params.Relatedness2 {
		rel2 = newRelatedness2Runner(filteredHeader.Samples)
	}

	// --LROH runner.
	var lroh *lrohRunner
	if params.LROH {
		lroh = newLROHRunner(filteredHeader.Samples, params.LROHMinVariants)
	}

	// --phased-blocks runner.
	var phasedBlocks *phasedBlockRunner
	if params.PhasedBlocks {
		phasedBlocks = newPhasedBlockRunner(filteredHeader.Samples)
	}

	// --get-INFO writer.
	var getInfo *getInfoRunner
	if params.GetINFO != "" {
		tags := splitCSV(params.GetINFO)
		getInfo, err = newGetInfoRunner(params.OutPrefix, tags)
		if err != nil {
			return fmt.Errorf("initialising --get-INFO: %w", err)
		}
	}

	// --extract-FORMAT-info NAME writer (one file per requested tag).
	// Upstream's `output_FORMAT_information` is a per-FORMAT-name TSV;
	// the flag is single-valued, so we only need one runner. See
	// extract_format.go for the row-emission rules.
	var extractFmt *extractFormatRunner
	if params.ExtractFormatInfo != "" {
		extractFmt, err = newExtractFormatRunner(params.OutPrefix, params.ExtractFormatInfo, filteredHeader.Samples)
		if err != nil {
			return fmt.Errorf("initialising --extract-FORMAT-info: %w", err)
		}
	}

	// --indv-burden / --indv-freq-burden / --indv-freq-burden2: per-
	// individual diploid-burden accumulators. Both runners stream all
	// kept sites (the filter chain ahead handles --chr / --maf / etc.)
	// and the writers flush at end-of-stream. Upstream dispatches the
	// three flags through two routines (vcftools.cpp:61-64) — we mirror
	// that exactly: --indv-freq-burden2 just sets doubleCountHomAlt.
	var indvBurden *indvBurdenRunner
	if params.IndvBurden {
		indvBurden = newIndvBurdenRunner(filteredHeader.Samples, params.Derived)
	}
	var indvFreqBurden *indvFreqBurdenRunner
	if params.IndvFreqBurden || params.IndvFreqBurden2 {
		// header.Samples carries the original (pre-sample-filter) VCF
		// kept (post-sample-filter) sample list for BOTH the row order
		// AND the leading INDV column label — fixing the upstream
		// `meta_data.indv[indv_count]` bug at line 621. See burden.go
		// and docs/UPSTREAM_BUGS.md for the writeup.
		indvFreqBurden = newIndvFreqBurdenRunner(
			filteredHeader.Samples,
			params.IndvFreqBurden2,
			params.Derived,
		)
	}

	// --hapcount BED: parse BED bins up front so missing-file /
	// overlap errors surface before we stream the VCF. Upstream errors
	// with "Could not open BED file" / "BED file must be
	// non-overlapping." at variant_file_output.cpp:1178 / :1214. The
	// runner is fed by the main per-variant loop below. Three upstream
	// bugs are fixed-on-port — see hapcount.go and
	// docs/UPSTREAM_BUGS.md.
	var hapcount *hapcountRunner
	if params.HapcountBED != "" {
		hapcount, err = newHapcountRunner(params.OutPrefix, params.HapcountBED, len(filteredHeader.Samples))
		if err != nil {
			return fmt.Errorf("initialising --hapcount: %w", err)
		}
	}

	// --pca / --pca-no-norm / --pca-snp-loadings: build the GRM from
	// centred (and optionally variance-normalised) genotypes, then
	// eigendecompose it via gonum. Mirrors upstream
	// variant_file_output.cpp:4871-5246. The runner accumulates per-site
	// rows during the variant loop and emits `<prefix>.pca` (and
	// optionally `<prefix>.pca.loadings`) at end-of-stream. Wave 19 — see
	// pca.go header and docs/UPSTREAM_BUGS.md (missing-data fix-on-port).
	var pca *pcaRunner
	// Upstream `--pca-no-norm` (parameters.cpp:298) implicitly enables
	// PCA output; `--pca-snp-loadings INT` (parameters.cpp:299) by
	// itself only writes `.pca.loadings`, NOT `.pca` (vcftools.cpp:110-111
	// dispatches the two outputs independently).
	if params.PCA || params.PCANoNorm || params.PCASNPLoadings != 0 {
		writePCA := params.PCA || params.PCANoNorm
		pca = newPCARunner(filteredHeader.Samples, writePCA, !params.PCANoNorm, params.PCASNPLoadings, params.UseStdout)
	}

	// --temp DIR: upstream parameters.cpp:341 stores DIR as the base
	// path for mkstemp spill files in the LD / format-convert paths.
	// This port does not spill, so the flag is accepted-and-noop. Emit
	// a stderr line so a confused user doesn't think it had effect.
	if params.TempDir != "" {
		fmt.Fprintf(os.Stderr, "--temp %q: accepted for CLI parity; this port does not spill to disk so the flag has no effect.\n", params.TempDir)
	}

	// --ldhat / --ldhat-geno: buffer per-site genotype rows and emit the
	// paired (.ldhat.sites, .ldhat.locs) bundle at end-of-stream. Only one
	// flag can be active at a time on the upstream CLI (each increments
	// num_outputs in parameters.cpp); if both are requested here we honour
	// the phased layout for parity with the upstream order-of-evaluation
	// in vcftools.cpp:90-91.
	var ldhat *ldhatRunner
	switch {
	case params.LDhat:
		ldhat = newLDhatRunner(ldhatPhased, params.OutPrefix, filteredHeader.Samples)
	case params.LDhatGeno:
		ldhat = newLDhatRunner(ldhatUnphased, params.OutPrefix, filteredHeader.Samples)
	}

	// --ldhelmet: buffer per-site haplotype contributions and emit the
	// paired (.ldhelmet.snps, .ldhelmet.pos) bundle at end-of-stream.
	var ldhelmet *ldhelmetRunner
	if params.LDhelmet {
		ldhelmet = newLDhelmetRunner(params.OutPrefix, filteredHeader.Samples)
	}

	// --IMPUTE: buffer per-site rows for the IMPUTE (.legend, .hap,
	// .hap.indv) bundle.
	var impute *imputeRunner
	if params.IMPUTE {
		impute = newImputeRunner(params.OutPrefix, filteredHeader.Samples)
	}

	// --mendel: parse the PED file, intersect with kept samples, open the
	// output file. Errors here (missing PED, no trios) fail fast.
	var mendel *mendelRunner
	if params.MendelPedFile != "" {
		mendel, err = newMendelRunner(params.OutPrefix, params.MendelPedFile, filteredHeader.Samples)
		if err != nil {
			return fmt.Errorf("initialising --mendel: %w", err)
		}
	}

	// Pre-parse filter-name sets and INFO-tag sets so we don't re-tokenise
	// every line in the hot path.
	removeFilteredSet := parseFilterList(params.RemoveFiltered)
	keepFilteredSet := parseFilterList(params.KeepFiltered)
	// recodeInfoSet drives the `--recode-INFO TAG` recode-column selector
	// (upstream parameters.cpp:319 → recode_INFO_to_keep). It restricts
	// the INFO column in `.recode.vcf` output to the listed tags. It is
	// NOT a site filter.
	recodeInfoSet := parseInfoTagList(params.RecodeINFO)

	// keepInfoSiteSet / removeInfoSiteSet drive the upstream
	// filter_sites_by_INFO routine (entry_filters.cpp:1033-1086): keep
	// requires at least one named Flag present (OR), remove drops the
	// site if any named Flag is present (OR-veto). Both sets require
	// every tag to be declared `Type=Flag` in the header — upstream
	// errors otherwise (entry_filters.cpp:1053/1072); we mirror that
	// error once at Run start since the check is header-invariant.
	keepInfoSiteSet := parseInfoTagList(params.KeepINFO)
	removeInfoSiteSet := parseInfoTagList(params.RemoveINFO)
	if err := validateFlagTypeINFO("--keep-INFO", keepInfoSiteSet, filteredHeader); err != nil {
		return err
	}
	if err := validateFlagTypeINFO("--remove-INFO", removeInfoSiteSet, filteredHeader); err != nil {
		return err
	}

	// Set up output writer for recode
	var recodeWriter *vcf.Writer
	if params.Recode {
		var w io.Writer
		if params.UseStdout {
			w = os.Stdout
		} else {
			outFile := params.OutPrefix + ".recode.vcf"
			f, err := iohelper.OpenWriter(outFile)
			if err != nil {
				return fmt.Errorf("opening output file: %w", err)
			}
			defer f.Close()
			w = f
		}
		recodeWriter = vcf.NewWriter(w, filteredHeader)
		if err := recodeWriter.WriteHeader(); err != nil {
			return fmt.Errorf("writing VCF header: %w", err)
		}
	}

	// Set up output writer for --recode-bcf. The BCF spec requires the
	// output stream to be BGZF-compressed (upstream `vcf_file.cpp:127`
	// calls `bgzf_open`). We layer:
	//
	//	file (or os.Stdout)  →  *bgzip.Writer  →  *bcf.Writer
	//
	// The BCF writer takes care of magic, length-prefixed text header,
	// and the per-record (l_shared,l_indiv,body) framing. Header is
	// synthesised from filteredHeader, so any prior --keep / --remove
	// sample filtering or contig dropping is reflected in the BCF
	// dictionary indices.
	var (
		recodeBCFWriter *bcf.Writer
		recodeBCFGzip   *bgzip.Writer
	)
	if params.RecodeBCF {
		var w io.Writer
		if params.UseStdout {
			w = os.Stdout
		} else {
			outFile := params.OutPrefix + ".recode.bcf"
			f, err := iohelper.OpenWriter(outFile)
			if err != nil {
				return fmt.Errorf("opening BCF output file: %w", err)
			}
			defer f.Close()
			w = f
		}
		recodeBCFGzip = bgzip.NewWriter(w)
		recodeBCFWriter, err = bcf.NewWriterFromVCFHeader(recodeBCFGzip, filteredHeader)
		if err != nil {
			return fmt.Errorf("initialising BCF writer: %w", err)
		}
		if err := recodeBCFWriter.WriteHeader(); err != nil {
			return fmt.Errorf("writing BCF header: %w", err)
		}
	}

	// --kept-sites / --removed-sites trace writers. The tracker is
	// non-nil even when both flags are off so the recordKept /
	// recordRemoved calls in the hot loop stay branchless aside from
	// the writer nil-check inside the method. See sitetrace.go for the
	// upstream byte-for-byte format derivation.
	siteTrace, err := newSiteTracker(params.OutPrefix, params.KeptSites, params.RemovedSites)
	if err != nil {
		return fmt.Errorf("initialising site-trace output: %w", err)
	}
	// Reviewer-flagged PR #133 nit: Run has ~20 early-return paths; without
	// a deferred close the buffered writers leak FDs and leave truncated
	// .kept.sites / .removed.sites files on disk. The explicit close at
	// the bottom of Run still runs first on the success path; the defer
	// is a safety net for all error paths.
	defer func() { _ = siteTrace.close() }()

	// Process variants
	keptSites := 0
	totalSites := 0
	thinChrom := ""
	thinPos := 0
	var allVariants []*vcf.Variant // For format conversions that need all data

	for {
		variant, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("reading variant: %w", err)
		}

		totalSites++

		// Apply --thin: drop sites that fall within `Thin` bp of the
		// most recently-kept site on the same chromosome. Mirrors
		// upstream entry::filter_sites_by_thinning (entry_filters.cpp:
		// 1012-1031): the first site on a chromosome is always kept;
		// later sites need POS - thin_pos >= Thin.
		if params.Thin > 0 {
			if variant.Chrom == thinChrom && variant.Pos-thinPos < params.Thin {
				if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
					return err
				}
				continue
			}
			thinChrom = variant.Chrom
			thinPos = variant.Pos
		}

		// Apply filters
		if !passFilters(variant, params, includePositions, excludePositions, includePositionsOverlap, excludePositionsOverlap, includeSNPs, excludeSNPs, includeBed, excludeBed) {
			if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
				return err
			}
			continue
		}

		// --mask / --invert-mask: forward-only FASTA-positional filter. The
		// cursor mutates with every call so this must run for every input
		// site, after the cheaper position/chr filters that don't move the
		// cursor. Upstream's filter ordering (entry_filters.cpp:47) puts
		// mask between mean-depth and phase; our equivalent slot is here
		// (the depth filters live inside passFilters already).
		if mask != nil {
			if !mask.passes(variant.Chrom, variant.Pos) {
				if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
					return err
				}
				continue
			}
		}

		// --remove-filtered <names>: drop sites listing any of the named
		// FILTERs. --keep-filtered <names>: keep only sites listing at
		// least one of the named FILTERs. Both compose with
		// --remove-filtered-all (which is the union of all non-PASS sites).
		if !passRemoveFilteredNames(variant, removeFilteredSet) {
			if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
				return err
			}
			continue
		}
		if !passKeepFilteredNames(variant, keepFilteredSet) {
			if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
				return err
			}
			continue
		}

		// --keep-INFO TAG / --remove-INFO TAG: site filters.
		// keep narrows first (drop site if NONE of the named Flags is
		// present, OR semantics), then remove vetoes the survivors
		// (drop site if ANY of the named Flags is present, OR-veto).
		// Mirrors upstream's `filter_sites_by_INFO`
		// (entry_filters.cpp:1033-1086) in upstream order.
		if len(keepInfoSiteSet) > 0 {
			if !passKeepINFOSite(variant, keepInfoSiteSet) {
				if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
					return err
				}
				continue
			}
		}
		if len(removeInfoSiteSet) > 0 {
			if !passRemoveINFOSite(variant, removeInfoSiteSet) {
				if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
					return err
				}
				continue
			}
		}

		// Filter samples
		filteredVariant := filterVariantSamples(variant, keepSamples)

		// Apply genotype-level filters
		filteredVariant = filterGenotypes(filteredVariant, params)

		// --phased (also implied by --ldhat, --ldhelmet, --IMPUTE,
		// --hapcount): drop sites with any unphased kept-individual
		// genotype. Mirrors upstream's filter_sites_by_phase
		// (entry_filters.cpp:989-1010), which iterates over kept
		// individuals only — so we run after the sample filter. The
		// --hapcount implication comes from parameters.cpp:248
		// (`phased_only = true`).
		if params.Phased || params.LDhat || params.LDhelmet || params.IMPUTE || params.HapcountBED != "" {
			if !isPhasedSite(filteredVariant) {
				if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
					return err
				}
				continue
			}
		}

		// --ldhelmet implies --remove-indels (parameters.cpp:275). We've
		// already passed the generic remove-indels check above via
		// passFilters when params.RemoveIndels is set; here we apply it
		// implicitly for --ldhelmet even if the user didn't set
		// --remove-indels explicitly.
		if params.LDhelmet && isIndelVariant(filteredVariant) {
			if err := siteTrace.recordRemoved(variant.Chrom, variant.Pos); err != nil {
				return err
			}
			continue
		}

		// --IMPUTE implies min_alleles == max_alleles == 2 (biallelic
		// only; parameters.cpp:255). The generic --min-alleles /
		// --max-alleles defaults in Params already cover the biallelic
		// case (MinAlleles=2 default), but the multi-allelic guard is
		// enforced inside imputeRunner.addVariant for the warning
		// behaviour to match upstream.

		// Update statistics
		stats.addVariant(filteredVariant, params)

		// Feed LD runner (writes pairwise output incrementally).
		if ldRun != nil {
			ldRun.addVariant(filteredVariant)
		}

		// Feed --diff runner.
		if diffRun != nil {
			if err := diffRun.addVariant(filteredVariant); err != nil {
				return fmt.Errorf("writing diff output: %w", err)
			}
		}

		// Emit BEAGLE rows (header is emitted lazily on the first call).
		if beagleGL != nil {
			if err := beagleGL.write(filteredVariant, filteredHeader.Samples); err != nil {
				return fmt.Errorf("writing BEAGLE-GL output: %w", err)
			}
		}
		if beaglePL != nil {
			if err := beaglePL.write(filteredVariant, filteredHeader.Samples); err != nil {
				return fmt.Errorf("writing BEAGLE-PL output: %w", err)
			}
		}

		// Inter-chromosomal LD: buffer for end-of-stream pair emission.
		if interLD != nil {
			interLD.addVariant(filteredVariant)
		}

		// Per-variant relatedness contribution.
		if rel != nil {
			rel.addVariant(filteredVariant)
		}
		if rel2 != nil {
			rel2.addVariant(filteredVariant)
		}

		// Per-variant LROH state update.
		if lroh != nil {
			lroh.addVariant(filteredVariant)
		}

		// Per-variant phased-block state update.
		if phasedBlocks != nil {
			phasedBlocks.addVariant(filteredVariant)
		}

		// Per-variant INFO extraction (--get-INFO).
		if getInfo != nil {
			if err := getInfo.addVariant(filteredVariant); err != nil {
				return fmt.Errorf("writing --get-INFO output: %w", err)
			}
		}

		// Per-variant FORMAT extraction (--extract-FORMAT-info NAME).
		if extractFmt != nil {
			if err := extractFmt.addVariant(filteredVariant); err != nil {
				return fmt.Errorf("writing --extract-FORMAT-info output: %w", err)
			}
		}

		// Per-individual burden accumulators (--indv-burden,
		// --indv-freq-burden, --indv-freq-burden2). Both runners drop
		// non-diploid sites internally to match upstream.
		if indvBurden != nil {
			indvBurden.addVariant(filteredVariant)
		}
		if indvFreqBurden != nil {
			indvFreqBurden.addVariant(filteredVariant)
		}

		// --hapcount: per-BED-bin haplotype-count update. Output is
		// emitted incrementally on each chromosome transition and at
		// close().
		if hapcount != nil {
			if err := hapcount.addVariant(filteredVariant); err != nil {
				return fmt.Errorf("writing --hapcount output: %w", err)
			}
		}

		// --pca: accumulate centred / normalised genotypes for the GRM.
		// The eigendecomposition runs at end-of-stream.
		if pca != nil {
			if err := pca.addVariant(filteredVariant); err != nil {
				return fmt.Errorf("--pca: %w", err)
			}
		}

		// LDhat buffer (biallelic-only filtering is applied inside).
		if ldhat != nil {
			ldhat.addVariant(filteredVariant)
		}

		// LDhelmet buffer (no biallelic filter — upstream just looks up
		// alleles[geno]).
		if ldhelmet != nil {
			ldhelmet.addVariant(filteredVariant)
		}

		// IMPUTE buffer (biallelic-only + per-site missing/unphased
		// guard applied inside).
		if impute != nil {
			impute.addVariant(filteredVariant)
		}

		// --mendel: per-trio Mendelian inconsistency check.
		if mendel != nil {
			if err := mendel.addVariant(filteredVariant); err != nil {
				return fmt.Errorf("writing --mendel output: %w", err)
			}
		}

		// Collect variants for format conversions
		if params.Output012 || params.OutputPlink || params.OutputPlinkTped {
			allVariants = append(allVariants, filteredVariant)
		}

		// Write to output if recoding. Both --recode (text VCF) and
		// --recode-bcf (binary BCF) may be active at the same time;
		// the outVariant body is shared.
		if params.Recode || params.RecodeBCF {
			var outInfo map[string]string
			switch {
			case len(recodeInfoSet) > 0:
				// --recode-INFO TAG (recode-column selector): restrict
				// the INFO column to the listed tags. Composes with
				// --recode-INFO-all by overriding it.
				outInfo = filterRecodeInfo(filteredVariant.Info, recodeInfoSet)
			case params.RecodeInfoAll:
				outInfo = filteredVariant.Info
			default:
				outInfo = make(map[string]string)
			}
			outVariant := &vcf.Variant{
				Chrom:   filteredVariant.Chrom,
				Pos:     filteredVariant.Pos,
				ID:      filteredVariant.ID,
				Ref:     filteredVariant.Ref,
				Alt:     filteredVariant.Alt,
				Qual:    filteredVariant.Qual,
				Filter:  filteredVariant.Filter,
				Info:    outInfo,
				Format:  filteredVariant.Format,
				Samples: filteredVariant.Samples,
			}
			if params.Recode {
				if err := recodeWriter.Write(outVariant); err != nil {
					return fmt.Errorf("writing variant: %w", err)
				}
			}
			if params.RecodeBCF {
				if err := recodeBCFWriter.Write(outVariant); err != nil {
					return fmt.Errorf("writing BCF variant: %w", err)
				}
			}
		}

		if err := siteTrace.recordKept(variant.Chrom, variant.Pos); err != nil {
			return err
		}
		keptSites++
	}

	// Flush recode writer if needed
	if recodeWriter != nil {
		if err := recodeWriter.Flush(); err != nil {
			return fmt.Errorf("flushing output: %w", err)
		}
	}

	// Flush + close the BCF writer chain. Order matters: BCF buffered
	// frames first, then the BGZF stream (which appends an EOF block
	// on Close per the BGZF spec). The underlying file handle is closed
	// by the deferred f.Close() from the writer setup above.
	if recodeBCFWriter != nil {
		if err := recodeBCFWriter.Flush(); err != nil {
			return fmt.Errorf("flushing BCF output: %w", err)
		}
		if err := recodeBCFGzip.Close(); err != nil {
			return fmt.Errorf("closing BGZF stream: %w", err)
		}
	}

	// Print summary to stderr
	fmt.Fprintf(os.Stderr, "\nAfter filtering, kept %d out of a possible %d Sites\n", keptSites, totalSites)

	// Output statistics
	if err := outputStatistics(stats, params); err != nil {
		return fmt.Errorf("outputting statistics: %w", err)
	}

	// Output format conversions
	if err := outputFormatConversions(allVariants, filteredHeader, params); err != nil {
		return fmt.Errorf("outputting format conversions: %w", err)
	}

	// Flush LD outputs.
	if err := ldRun.close(); err != nil {
		return fmt.Errorf("closing LD output: %w", err)
	}

	// Flush --diff outputs (also emits file-2-only sites and per-individual
	// reports).
	if err := diffRun.close(); err != nil {
		return fmt.Errorf("closing --diff output: %w", err)
	}

	// Flush BEAGLE outputs.
	if err := beagleGL.close(); err != nil {
		return fmt.Errorf("closing --BEAGLE-GL output: %w", err)
	}
	if err := beaglePL.close(); err != nil {
		return fmt.Errorf("closing --BEAGLE-PL output: %w", err)
	}

	// Inter-chromosomal LD requires all sites in memory; emit pairs now.
	if interLD != nil {
		if err := interLD.flush(); err != nil {
			return fmt.Errorf("flushing interchrom LD: %w", err)
		}
		if err := interLD.close(); err != nil {
			return fmt.Errorf("closing interchrom LD output: %w", err)
		}
	}

	// Relatedness output.
	if rel != nil {
		if err := rel.writeOutput(params.OutPrefix); err != nil {
			return fmt.Errorf("writing --relatedness output: %w", err)
		}
	}
	if rel2 != nil {
		if err := rel2.writeOutput(params.OutPrefix); err != nil {
			return fmt.Errorf("writing --relatedness2 output: %w", err)
		}
	}

	// LROH output.
	if lroh != nil {
		if err := lroh.writeOutput(params.OutPrefix); err != nil {
			return fmt.Errorf("writing --LROH output: %w", err)
		}
	}

	// Phased blocks output.
	if phasedBlocks != nil {
		if err := phasedBlocks.writeOutput(params.OutPrefix); err != nil {
			return fmt.Errorf("writing --phased-blocks output: %w", err)
		}
	}

	// --get-INFO output.
	if err := getInfo.close(); err != nil {
		return fmt.Errorf("closing --get-INFO output: %w", err)
	}

	// --extract-FORMAT-info output.
	if err := extractFmt.close(); err != nil {
		return fmt.Errorf("closing --extract-FORMAT-info output: %w", err)
	}

	// --indv-burden output.
	if err := indvBurden.writeOutput(params.OutPrefix); err != nil {
		return fmt.Errorf("writing --indv-burden output: %w", err)
	}

	// --indv-freq-burden / --indv-freq-burden2 output. Both flags
	// share the .ifreqburden suffix upstream, so we only call once.
	if err := indvFreqBurden.writeOutput(params.OutPrefix); err != nil {
		return fmt.Errorf("writing --indv-freq-burden output: %w", err)
	}

	// --ldhat / --ldhat-geno output.
	if err := ldhat.close(); err != nil {
		return fmt.Errorf("closing --ldhat output: %w", err)
	}

	// --ldhelmet output.
	if err := ldhelmet.close(); err != nil {
		return fmt.Errorf("closing --ldhelmet output: %w", err)
	}

	// --IMPUTE output.
	if err := impute.close(); err != nil {
		return fmt.Errorf("closing --IMPUTE output: %w", err)
	}

	// --mendel output.
	if err := mendel.close(); err != nil {
		return fmt.Errorf("closing --mendel output: %w", err)
	}

	// --hapcount output. The close() handler emits the last seen
	// chromosome unconditionally (FIX for upstream bug #2 — see
	// hapcount.go).
	if err := hapcount.close(); err != nil {
		return fmt.Errorf("closing --hapcount output: %w", err)
	}

	// --pca / --pca-snp-loadings: end-of-stream eigendecomposition + write.
	if pca != nil {
		if err := pca.computeAndWrite(params.OutPrefix); err != nil {
			return fmt.Errorf("--pca: %w", err)
		}
	}

	// --kept-sites / --removed-sites output.
	if err := siteTrace.close(); err != nil {
		return err
	}

	return nil
}

// splitCSV splits a comma-separated string into trimmed non-empty tokens, in
// order. Used for --get-INFO and similar list-valued flags.
func splitCSV(s string) []string {
	if s == "" {
		return nil
	}
	parts := strings.Split(s, ",")
	out := make([]string, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" {
			out = append(out, p)
		}
	}
	return out
}

// beagleGLMode and beaglePLMode are tiny helpers so we don't need to expose
// the unexported beagleMode constants outside the package.
func beagleGLMode() beagleMode { return beagleGL }
func beaglePLMode() beagleMode { return beaglePL }

// passFilters checks if a variant passes all filters
func passFilters(v *vcf.Variant, params *Params, includePos, excludePos, includePosOverlap, excludePosOverlap positionSet, includeSNPs, excludeSNPs map[string]bool, includeBed, excludeBed *bedRegions) bool {
	// SNP ID filters
	if includeSNPs != nil && len(includeSNPs) > 0 {
		if !includeSNPs[v.ID] {
			return false
		}
	}
	if excludeSNPs != nil && len(excludeSNPs) > 0 {
		if excludeSNPs[v.ID] {
			return false
		}
	}

	// Position filters
	if params.Chr != "" && v.Chrom != params.Chr {
		return false
	}
	if params.NotChr != "" && v.Chrom == params.NotChr {
		return false
	}
	if params.FromBp > 0 && v.Pos < params.FromBp {
		return false
	}
	if params.ToBp > 0 && v.Pos > params.ToBp {
		return false
	}

	// Position include/exclude
	if includePos != nil {
		if chromPos, ok := includePos[v.Chrom]; ok {
			if !chromPos[v.Pos] {
				return false
			}
		} else {
			return false
		}
	}
	if excludePos != nil {
		if chromPos, ok := excludePos[v.Chrom]; ok {
			if chromPos[v.Pos] {
				return false
			}
		}
	}

	// --positions-overlap / --exclude-positions-overlap: sweep every base
	// in [POS, POS+len(REF)-1] against the set. Ported from upstream
	// entry_filters.cpp:513-547. Length is taken from REF (the loop is
	// `for ui=POS; ui<POS+REF.size(); ui++`); we guard against an empty
	// REF (defensive — valid VCFs always have len(REF) >= 1) so the loop
	// always considers at least the POS itself, matching what upstream
	// does when REF.size() >= 1.
	if includePosOverlap != nil {
		chromPos, ok := includePosOverlap[v.Chrom]
		if !ok {
			return false
		}
		refLen := len(v.Ref)
		if refLen < 1 {
			refLen = 1
		}
		found := false
		for p := v.Pos; p < v.Pos+refLen; p++ {
			if chromPos[p] {
				found = true
				break
			}
		}
		if !found {
			return false
		}
	}
	if excludePosOverlap != nil {
		if chromPos, ok := excludePosOverlap[v.Chrom]; ok {
			refLen := len(v.Ref)
			if refLen < 1 {
				refLen = 1
			}
			for p := v.Pos; p < v.Pos+refLen; p++ {
				if chromPos[p] {
					return false
				}
			}
		}
	}

	// BED-based include/exclude. Sites must be inside any --bed interval and
	// must not be inside any --exclude-bed interval. Each is independent of
	// the other (no implicit subtraction).
	if includeBed != nil {
		if !includeBed.containsVCFPos(v.Chrom, v.Pos) {
			return false
		}
	}
	if excludeBed != nil {
		if excludeBed.containsVCFPos(v.Chrom, v.Pos) {
			return false
		}
	}

	// Variant type filters
	isIndel := isIndelVariant(v)
	if params.KeepOnlyIndels && !isIndel {
		return false
	}
	if params.RemoveIndels && isIndel {
		return false
	}

	// Allele count filter
	numAlleles := len(v.Alt) + 1 // +1 for reference
	if params.MinAlleles > 0 && numAlleles < params.MinAlleles {
		return false
	}
	if params.MaxAlleles > 0 && numAlleles > params.MaxAlleles {
		return false
	}

	// Quality filter
	if params.MinQ > 0 && (v.Qual < 0 || v.Qual < params.MinQ) {
		return false
	}

	// Filter flag
	if params.RemoveFilteredAll {
		if len(v.Filter) == 0 || (len(v.Filter) == 1 && v.Filter[0] != "PASS") {
			if len(v.Filter) == 0 || v.Filter[0] != "PASS" {
				return false
			}
		}
	}

	// Allele frequency filters
	if params.Maf > 0 || params.MaxMaf > 0 || params.Mac > 0 || params.MaxMac > 0 {
		maf, mac := calculateMAF(v)
		if params.Maf > 0 && maf < params.Maf {
			return false
		}
		if params.MaxMaf > 0 && maf > params.MaxMaf {
			return false
		}
		if params.Mac > 0 && mac < params.Mac {
			return false
		}
		if params.MaxMac > 0 && mac > params.MaxMac {
			return false
		}
	}

	// Non-reference allele filters (--non-ref-af, --non-ref-ac, and the
	// matching --max-* / -any variants). Mirrors upstream
	// entry_filters.cpp:770-823 (frequency) and 869-919 (count).
	//
	// For every ALT allele we evaluate FOUR independent checks per branch:
	//   plain min  — drop site immediately if `value < MinPlain` (line 807/905)
	//   plain max  — drop site immediately if `value > MaxPlain` (line 807/905)
	//   any  min   — increment N_failed if `value < MinAny`     (line 810/908)
	//   any  max   — increment N_failed if `value > MaxAny`     (line 810/908)
	//
	// After the loop, upstream applies a fallback that drops the site when
	// N_failed equals (N_alleles - 1). The gate that activates this fallback
	// differs between AF and AC:
	//   AF (line 814): gate uses PLAIN thresholds (min_non_ref_af > 0 OR
	//                  max_non_ref_af < +inf). The -any thresholds never
	//                  trigger the fallback on their own — so the AF -any
	//                  flags are no-ops alone. We mirror this.
	//   AC (line 912): gate uses -any thresholds (min_non_ref_ac_any > 0 OR
	//                  max_non_ref_ac_any < INT_MAX). The plain AC flags
	//                  never trigger the fallback — so plain --non-ref-ac
	//                  keeps monomorphic sites whereas --non-ref-af drops
	//                  them.
	//
	// We collapse this into a single pass per branch that accumulates
	// N_failed across all ALTs (lifting the per-ALT early-return so the
	// _any-flavoured semantics can decide post-loop), exactly matching
	// upstream's structure.
	afPlainOn := params.MinNonRefAF > 0 || params.MaxNonRefAF > 0
	afAnyOn := params.MinNonRefAFAny > 0 || params.MaxNonRefAFAny > 0
	acPlainOn := params.MinNonRefAC > 0 || params.MaxNonRefAC > 0
	acAnyOn := params.MinNonRefACAny > 0 || params.MaxNonRefACAny > 0
	if afPlainOn || afAnyOn || acPlainOn || acAnyOn {
		altCounts, totalCalled := calculateAlleleCounts(v)
		// nAlt is the count of real ALTs (upstream's N_alleles - 1).
		// ALT == "." or empty is upstream's "no ALT" sentinel.
		nAlt := 0
		for _, alt := range v.Alt {
			if alt != "" && alt != "." {
				nAlt++
			}
		}
		// Per-branch N_failed accumulators for the -any fallback.
		afNFailed := 0
		acNFailed := 0
		for i, alt := range v.Alt {
			if alt == "" || alt == "." {
				continue
			}
			c := 0
			if i < len(altCounts) {
				c = altCounts[i]
			}
			// AC branch — plain immediate-drop, then -any tally.
			if acPlainOn {
				if params.MinNonRefAC > 0 && c < params.MinNonRefAC {
					return false
				}
				if params.MaxNonRefAC > 0 && c > params.MaxNonRefAC {
					return false
				}
			}
			if acAnyOn {
				failed := false
				if params.MinNonRefACAny > 0 && c < params.MinNonRefACAny {
					failed = true
				}
				if params.MaxNonRefACAny > 0 && c > params.MaxNonRefACAny {
					failed = true
				}
				if failed {
					acNFailed++
				}
			}
			// AF branch — same shape, freq instead of count. Upstream
			// uses freq = count / N_non_missing_chr; if there are no
			// called chromosomes the divide is degenerate (NaN), and
			// upstream's `freq < threshold` becomes false (NaN compares
			// false). We skip the AF checks entirely when totalCalled
			// == 0 to mirror that behaviour without producing NaN.
			if (afPlainOn || afAnyOn) && totalCalled > 0 {
				freq := float64(c) / float64(totalCalled)
				if afPlainOn {
					if params.MinNonRefAF > 0 && freq < params.MinNonRefAF {
						return false
					}
					if params.MaxNonRefAF > 0 && freq > params.MaxNonRefAF {
						return false
					}
				}
				if afAnyOn {
					failed := false
					if params.MinNonRefAFAny > 0 && freq < params.MinNonRefAFAny {
						failed = true
					}
					if params.MaxNonRefAFAny > 0 && freq > params.MaxNonRefAFAny {
						failed = true
					}
					if failed {
						afNFailed++
					}
				}
			}
		}
		// AF fallback (entry_filters.cpp:814). Gate keyed on PLAIN
		// thresholds. Upstream's check is `N_failed == (N_alleles -
		// 1)`; N_alleles includes REF, so (N_alleles - 1) == nAlt
		// here. A monomorphic site (nAlt == 0) satisfies 0 == 0 and
		// is dropped via this fallback — matching the wave-6 behaviour
		// that pinned the AF flag's mono-drop quirk.
		if afPlainOn && afNFailed == nAlt {
			return false
		}
		// AC fallback (entry_filters.cpp:912). Gate keyed on `-any`
		// thresholds — this is the upstream asymmetry vs the AF
		// branch's plain-keyed gate above.
		if acAnyOn && acNFailed == nAlt {
			return false
		}
	}

	// Genotype filters. Upstream's --max-missing is the MIN fraction of
	// non-missing genotypes (0.0 = allow all, 1.0 = require all
	// non-missing). The Params field is the same semantics: 0 means
	// "feature disabled" (no filter), >0 means apply. We guard against
	// the zero default explicitly so `--max-missing 1.0` (require all
	// non-missing) still applies — the previous `< 1` guard mistakenly
	// dropped that exact case.
	if params.MaxMissing > 0 {
		missingRate := calculateMissingRate(v)
		if missingRate > (1 - params.MaxMissing) {
			return false
		}
	}

	if params.MinMeanDP > 0 || params.MaxMeanDP > 0 {
		meanDP := calculateMeanDepth(v)
		if params.MinMeanDP > 0 && meanDP < params.MinMeanDP {
			return false
		}
		if params.MaxMeanDP > 0 && meanDP > params.MaxMeanDP {
			return false
		}
	}

	// --max-missing-count: drop the site if the total number of missing
	// chromosomes (haploid alleles) exceeds the threshold. Mirrors
	// upstream `(N_chr - N_non_missing_chr) > max_missing_call_count`
	// (entry_filters.cpp:918). Only active when the flag was explicitly
	// set on the command line.
	if params.MaxMissingCountSet {
		missing := countMissingChromosomes(v)
		if missing > params.MaxMissingCount {
			return false
		}
	}

	// --hwe FLOAT: drop biallelic sites whose exact-test HWE p-value is
	// below the threshold. Upstream pairs this with `max_alleles = 2`
	// (parameters.cpp:254). The `min_alleles`/`max_alleles` check runs
	// earlier above; here we only need to skip non-biallelic sites
	// quietly (they've already been rejected by max_alleles=2) and skip
	// sites with zero called diploid genotypes (upstream's SNPHWE
	// returns p_hwe=1.0 for an empty site, never failing the filter).
	if params.MinHWEPvalue > 0 {
		hom1, het, hom2, biallelic := countDiploidGenotypes(v)
		if !biallelic {
			// Multi-allelic or REF-only: upstream's max_alleles=2 will
			// have already dropped most of these; for REF-only "no
			// ALT" sites SNPHWE returns 1.0 (every observation is a
			// hom1), so we treat them as passing.
			return true
		}
		p := snpHWE(het, hom1, hom2)
		if p < params.MinHWEPvalue {
			return false
		}
	}

	return true
}

// isIndelVariant checks if a variant is an indel
func isIndelVariant(v *vcf.Variant) bool {
	refLen := len(v.Ref)
	for _, alt := range v.Alt {
		if len(alt) != refLen {
			return true
		}
	}
	return false
}

// isSNPVariant mirrors upstream `entry::is_SNP` (entry_getters.cpp): a
// site is a SNP when REF is a single base and at least one ALT is a
// single base distinct from REF, both restricted to A/C/G/T.
func isSNPVariant(v *vcf.Variant) bool {
	if len(v.Ref) != 1 || !allACGT(v.Ref) {
		return false
	}
	for _, alt := range v.Alt {
		if len(alt) == 1 && allACGT(alt) && alt != v.Ref {
			return true
		}
	}
	return false
}

// calculateAlleleCounts returns the per-ALT allele counts and the total
// number of non-missing called chromosomes across all samples. altCounts[i]
// is the count of ALT allele i (using v.Alt's 0-based indexing; upstream
// indexes ALT alleles at 1..N where 0 is REF, so the caller must add the
// REF count separately if needed). Missing alleles ('.') are skipped.
// Mirrors entry_getters.cpp:389-422 (entry::get_allele_counts) but
// returns just the ALT slice for the --non-ref-af / --non-ref-ac filters.
func calculateAlleleCounts(v *vcf.Variant) (altCounts []int, totalCalled int) {
	altCounts = make([]int, len(v.Alt))
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "" || gt == "." {
			continue
		}
		// Split on '/' or '|'. The vcf package leaves GT as the raw
		// string ("0|0", "1/2", "./.").
		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})
		for _, a := range alleles {
			if a == "" || a == "." {
				continue
			}
			// Each called chromosome contributes 1 to the total.
			totalCalled++
			// Allele 0 is REF; >0 is ALT index 1..N which maps to
			// v.Alt[idx-1].
			idx, err := strconv.Atoi(a)
			if err != nil || idx <= 0 {
				continue
			}
			if idx-1 < len(altCounts) {
				altCounts[idx-1]++
			}
		}
	}
	return altCounts, totalCalled
}

// calculateMAF calculates minor allele frequency and count
func calculateMAF(v *vcf.Variant) (maf float64, mac int) {
	if len(v.Samples) == 0 {
		return 0, 0
	}

	// Count alleles
	alleleCounts := make(map[string]int)
	totalAlleles := 0

	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "." || gt == "./." || gt == ".|." {
			continue
		}

		// Parse genotype
		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})

		for _, allele := range alleles {
			if allele != "." {
				alleleCounts[allele]++
				totalAlleles++
			}
		}
	}

	if totalAlleles == 0 {
		return 0, 0
	}

	// Find minor allele (least frequent)
	minCount := totalAlleles
	for _, count := range alleleCounts {
		if count < minCount {
			minCount = count
		}
	}

	// If all alleles are the same, minor allele count is 0
	if len(alleleCounts) == 1 {
		return 0, 0
	}

	mac = minCount
	maf = float64(minCount) / float64(totalAlleles)
	return maf, mac
}

// calculateMissingRate calculates the proportion of missing genotypes
func calculateMissingRate(v *vcf.Variant) float64 {
	if len(v.Samples) == 0 {
		return 0
	}

	missing := 0
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "." || gt == "./." || gt == ".|." {
			missing++
		}
	}

	return float64(missing) / float64(len(v.Samples))
}

// calculateMeanDepth calculates mean depth across samples
func calculateMeanDepth(v *vcf.Variant) float64 {
	if len(v.Samples) == 0 {
		return 0
	}

	totalDepth := 0
	count := 0

	for _, sample := range v.Samples {
		dpStr, ok := sample.Data["DP"]
		if !ok {
			continue
		}
		dp, err := strconv.Atoi(dpStr)
		if err != nil {
			continue
		}
		totalDepth += dp
		count++
	}

	if count == 0 {
		return 0
	}

	return float64(totalDepth) / float64(count)
}

// loadPositions loads a positions file
func loadPositions(filename string) (positionSet, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	positions := make(positionSet)
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}

		chrom := fields[0]
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			continue
		}

		if positions[chrom] == nil {
			positions[chrom] = make(map[int]bool)
		}
		positions[chrom][pos] = true
	}

	return positions, scanner.Err()
}

// loadSNPIDs loads SNP IDs from a file
func loadSNPIDs(filename string) (map[string]bool, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	snpIDs := make(map[string]bool)
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		// SNP ID is the first field (or whole line if no whitespace)
		fields := strings.Fields(line)
		if len(fields) > 0 {
			snpIDs[fields[0]] = true
		}
	}

	return snpIDs, scanner.Err()
}

// buildSampleFilter builds a set of samples to keep
func buildSampleFilter(header *vcf.Header, params *Params) (map[string]bool, error) {
	// Decide whether the user touched any sample-affecting flag. --max-indv
	// is its own gate because it's not an identity-based filter — it caps
	// the count, so it can be active even when --indv/--keep aren't.
	hasIdentityFilter := len(params.IndvList) > 0 || len(params.RemoveIndvList) > 0 ||
		params.KeepFile != "" || params.RemoveFile != ""
	hasMaxIndv := params.MaxIndvSet && params.MaxIndv >= 0
	if !hasIdentityFilter && !hasMaxIndv {
		return nil, nil
	}

	// Start with all samples
	keep := make(map[string]bool)
	for _, sample := range header.Samples {
		keep[sample] = true
	}

	// Apply keep file
	if params.KeepFile != "" {
		samples, err := loadSampleFile(params.KeepFile)
		if err != nil {
			return nil, err
		}
		newKeep := make(map[string]bool)
		for _, sample := range samples {
			if keep[sample] {
				newKeep[sample] = true
			}
		}
		keep = newKeep
	}

	// Apply keep list
	if len(params.IndvList) > 0 {
		newKeep := make(map[string]bool)
		for _, sample := range params.IndvList {
			if keep[sample] {
				newKeep[sample] = true
			}
		}
		keep = newKeep
	}

	// Apply remove file
	if params.RemoveFile != "" {
		samples, err := loadSampleFile(params.RemoveFile)
		if err != nil {
			return nil, err
		}
		for _, sample := range samples {
			delete(keep, sample)
		}
	}

	// Apply remove list
	for _, sample := range params.RemoveIndvList {
		delete(keep, sample)
	}

	// --max-indv N: cap the number of kept individuals at N. Ported from
	// upstream parameters.cpp:292 + variant_file_filters.cpp:105-147
	// (filter_individuals_randomly). Upstream uses srand(time(NULL)) +
	// random_shuffle, so the kept-sample identity varies across runs and
	// cannot be byte-matched. This port instead deterministically keeps
	// the first N kept samples in input (header) order — see the comment
	// on Params.MaxIndv for the parity claim. N < 0 is "no cap"; N == 0
	// is a valid drop-everything request (matches upstream's behaviour
	// when max_N_indv == 0: filter_individuals_randomly returns early
	// only when max_N_indv < 0).
	if hasMaxIndv {
		ordered := make([]string, 0, len(keep))
		for _, sample := range header.Samples {
			if keep[sample] {
				ordered = append(ordered, sample)
			}
		}
		if len(ordered) > params.MaxIndv {
			capped := make(map[string]bool, params.MaxIndv)
			for i := 0; i < params.MaxIndv; i++ {
				capped[ordered[i]] = true
			}
			keep = capped
		}
	}

	return keep, nil
}

// loadSampleFile loads a file with one sample name per line
func loadSampleFile(filename string) ([]string, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var samples []string
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" && !strings.HasPrefix(line, "#") {
			samples = append(samples, line)
		}
	}

	return samples, scanner.Err()
}

// filterHeaderSamples filters header samples based on keep set
func filterHeaderSamples(header *vcf.Header, keepSamples map[string]bool) *vcf.Header {
	if keepSamples == nil {
		return header
	}

	filtered := &vcf.Header{
		MetaInfo: header.MetaInfo,
		Samples:  []string{},
	}

	for _, sample := range header.Samples {
		if keepSamples[sample] {
			filtered.Samples = append(filtered.Samples, sample)
		}
	}

	return filtered
}

// filterVariantSamples filters variant samples based on keep set
func filterVariantSamples(v *vcf.Variant, keepSamples map[string]bool) *vcf.Variant {
	if keepSamples == nil {
		return v
	}

	filtered := &vcf.Variant{
		Chrom:  v.Chrom,
		Pos:    v.Pos,
		ID:     v.ID,
		Ref:    v.Ref,
		Alt:    v.Alt,
		Qual:   v.Qual,
		Filter: v.Filter,
		Info:   v.Info,
		Format: v.Format,
	}

	for _, sample := range v.Samples {
		if keepSamples[sample.Name] {
			filtered.Samples = append(filtered.Samples, sample)
		}
	}

	return filtered
}

// filterGenotypes applies genotype-level filters (sets genotypes to missing if they fail filters)
func filterGenotypes(v *vcf.Variant, params *Params) *vcf.Variant {
	// Build the FT-flag set once per call so the hot loop just probes a map.
	var ftDropSet map[string]struct{}
	if len(params.RemoveFilteredGenoList) > 0 {
		ftDropSet = make(map[string]struct{}, len(params.RemoveFilteredGenoList))
		for _, name := range params.RemoveFilteredGenoList {
			name = strings.TrimSpace(name)
			if name == "" {
				continue
			}
			ftDropSet[name] = struct{}{}
		}
	}
	ftFilterActive := params.RemoveFilteredGenoAll || len(ftDropSet) > 0

	// If no genotype filters specified, return as-is.
	if params.MinDP == 0 && params.MaxDP == 0 && params.MinGQ == 0 && !ftFilterActive {
		return v
	}

	// Create a copy to avoid modifying original
	filtered := &vcf.Variant{
		Chrom:  v.Chrom,
		Pos:    v.Pos,
		ID:     v.ID,
		Ref:    v.Ref,
		Alt:    v.Alt,
		Qual:   v.Qual,
		Filter: v.Filter,
		Info:   v.Info,
		Format: v.Format,
	}

	// Upstream's filter_genotypes_by_depth (vcf_entry.cpp:540-557)
	// activates only when the site's FORMAT carries DP at all
	// (DP_idx != -1). When it does, every sample is checked: a missing
	// or unparseable per-sample DP is set to -1 by set_indv_DEPTH
	// (entry_setters.cpp:31-35) and will fail any positive MinDP.
	// Same shape applies to GQ.
	hasDPFormat := false
	hasGQFormat := false
	for _, fld := range v.Format {
		if fld == "DP" {
			hasDPFormat = true
		}
		if fld == "GQ" {
			hasGQFormat = true
		}
	}
	for _, sample := range v.Samples {
		// Check DP (depth) filter
		if (params.MinDP > 0 || params.MaxDP > 0) && hasDPFormat {
			dp := -1
			if dpStr, ok := sample.Data["DP"]; ok && dpStr != "" && dpStr != "." {
				if parsed, err := strconv.Atoi(dpStr); err == nil {
					dp = parsed
				}
			}
			if (params.MinDP > 0 && dp < params.MinDP) || (params.MaxDP > 0 && dp > params.MaxDP) {
				filtered.Samples = append(filtered.Samples, sampleWithMissingGT(sample))
				continue
			}
		}

		// Check GQ (genotype quality) filter
		if params.MinGQ > 0 && hasGQFormat {
			gq := -1
			if gqStr, ok := sample.Data["GQ"]; ok && gqStr != "" && gqStr != "." {
				if parsed, err := strconv.Atoi(gqStr); err == nil {
					gq = parsed
				}
			}
			if gq < params.MinGQ {
				filtered.Samples = append(filtered.Samples, sampleWithMissingGT(sample))
				continue
			}
		}

		// --remove-filtered-geno-all / --remove-filtered-geno NAME:
		// inspect the FORMAT FT field. Upstream parses FT as a
		// semicolon-separated list (vcf_entry_setters.cpp:188-212);
		// entries equal to "" or "." mean "unfiltered" and are
		// dropped from the list. --all keeps a genotype only if its
		// first FT entry is "PASS" or the list is effectively empty
		// (mirrors upstream vcf_entry.cpp:594-599: it inspects
		// GFILTERs[0]); --remove-filtered-geno NAME drops the
		// genotype if any FT entry matches a named flag
		// (vcf_entry.cpp:601-605). Sites without an FT FORMAT column
		// are left untouched (upstream's filter_genotypes_by_filter_status
		// in entry_filters.cpp:94-108 returns early when FT_idx == -1).
		if ftFilterActive {
			if ftEntries, hasFT := parseSampleFT(sample); hasFT {
				if shouldDropByFT(ftEntries, params.RemoveFilteredGenoAll, ftDropSet) {
					filtered.Samples = append(filtered.Samples, sampleWithMissingGT(sample))
					continue
				}
			}
		}

		// Genotype passes filters, keep as-is
		filtered.Samples = append(filtered.Samples, sample)
	}

	return filtered
}

// sampleWithMissingGT clones a sample and replaces its GT field with "./."
// while preserving every other FORMAT field. Mirrors upstream's recode
// emission for include_genotype[ui]==false (vcf_entry.cpp:341 — only the
// GT slot is rewritten; FT/DP/GQ/... are passed through unchanged).
func sampleWithMissingGT(sample vcf.Sample) vcf.Sample {
	out := vcf.Sample{
		Name: sample.Name,
		Data: make(map[string]string, len(sample.Data)),
	}
	for k, val := range sample.Data {
		if k == "GT" {
			out.Data[k] = "./."
		} else {
			out.Data[k] = val
		}
	}
	return out
}

// parseSampleFT splits the sample's FT FORMAT field on ';' and returns the
// non-trivial entries (empty or "." are dropped, matching upstream's
// vcf_entry_setters.cpp:188-212). The second return value is false when
// the sample has no FT field at all.
func parseSampleFT(sample vcf.Sample) ([]string, bool) {
	raw, ok := sample.Data["FT"]
	if !ok {
		return nil, false
	}
	if raw == "" || raw == "." {
		return nil, true
	}
	parts := strings.Split(raw, ";")
	out := parts[:0]
	for _, p := range parts {
		if p == "" || p == "." {
			continue
		}
		out = append(out, p)
	}
	return out, true
}

// shouldDropByFT returns true if a genotype should be set to missing given
// its FT entries. Mirrors the two upstream branches at
// vcf_entry.cpp:594-605:
//   - removeAll == true: drop unless the first FT entry is "PASS" (an
//     empty list — i.e. the unfiltered case — is implicitly kept).
//   - removeAll == false: drop if ANY FT entry matches a named flag.
func shouldDropByFT(ftEntries []string, removeAll bool, namedDrops map[string]struct{}) bool {
	if removeAll {
		if len(ftEntries) == 0 {
			return false
		}
		return ftEntries[0] != "PASS"
	}
	if len(namedDrops) == 0 {
		return false
	}
	for _, e := range ftEntries {
		if _, hit := namedDrops[e]; hit {
			return true
		}
	}
	return false
}

// outputStatistics outputs all requested statistics
func outputStatistics(stats *statistics, params *Params) error {
	if params.Freq {
		if err := stats.outputFrequency(params.OutPrefix, false, false); err != nil {
			return err
		}
	}

	if params.Counts {
		if err := stats.outputFrequency(params.OutPrefix, true, false); err != nil {
			return err
		}
	}

	if params.Freq2 {
		if err := stats.outputFrequency(params.OutPrefix, false, true); err != nil {
			return err
		}
	}

	if params.Counts2 {
		if err := stats.outputFrequency(params.OutPrefix, true, true); err != nil {
			return err
		}
	}

	if params.SiteMeanDepth {
		if err := stats.outputSiteMeanDepth(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.SiteDepth {
		if err := stats.outputSiteDepth(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.SiteQuality {
		if err := stats.outputSiteQuality(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.MissingSite {
		if err := stats.outputMissingSite(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.MissingIndv {
		if err := stats.outputMissingIndv(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.Hardy {
		if err := stats.outputHWE(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.TsTvSummary {
		if err := stats.outputTsTvSummary(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.TsTvBinSize > 0 {
		if err := stats.outputTsTvByBin(params.OutPrefix, params.TsTvBinSize); err != nil {
			return err
		}
	}

	if params.TsTvByCount {
		if err := stats.outputTsTvByCount(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.TsTvByQual {
		if err := stats.outputTsTvByQual(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.Depth {
		if err := stats.outputDepth(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.SitePi {
		if err := stats.outputSitePi(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.WindowPi > 0 {
		if err := stats.outputWindowedPi(params.OutPrefix, params.WindowPi, params.WindowPiStep); err != nil {
			return err
		}
	}

	// Phase 2: Population genetics statistics

	if params.Het {
		if err := stats.outputHet(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.Singletons {
		if err := stats.outputSingletons(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.FilterSummary {
		if err := stats.outputFilterSummary(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.SNPDensity > 0 {
		if err := stats.outputSNPDensity(params.OutPrefix, params.SNPDensity); err != nil {
			return err
		}
	}

	if params.GenoDepth {
		if err := stats.outputGenoDepth(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.HistIndelLen {
		if err := stats.outputIndelHist(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.TajimaD > 0 {
		if err := stats.outputTajimaD(params.OutPrefix, params.TajimaD); err != nil {
			return err
		}
	}

	// Weir & Cockerham 1984 Fst (per-site, plus optional windowed output).
	if stats.weirFst != nil {
		if err := stats.weirFst.outputWeirFst(params.OutPrefix); err != nil {
			return err
		}
		if params.FstWindowSize > 0 {
			if err := stats.weirFst.outputWindowedWeirFst(params.OutPrefix, params.FstWindowSize, params.FstWindowStep); err != nil {
				return err
			}
		}
	}

	return nil
}

// outputFormatConversions outputs requested format conversions
func outputFormatConversions(variants []*vcf.Variant, header *vcf.Header, params *Params) error {
	if params.Output012 {
		if err := output012Matrix(variants, header, params.OutPrefix); err != nil {
			return err
		}
	}

	if params.OutputPlink || params.OutputPlinkTped {
		// Load chromosome map if provided
		chromMap, err := loadChromMap(params.ChromMap)
		if err != nil {
			return fmt.Errorf("loading chromosome map: %w", err)
		}

		if params.OutputPlink {
			if err := outputPlink(variants, header, params.OutPrefix, chromMap); err != nil {
				return err
			}
		}

		if params.OutputPlinkTped {
			if err := outputPlinkTped(variants, header, params.OutPrefix, chromMap); err != nil {
				return err
			}
		}
	}

	return nil
}

// Helper function to get output file path
func getOutputPath(prefix, suffix string) string {
	return filepath.Join(".", prefix+suffix)
}

// checkUnsupported returns an error if the parameters request a feature that
// this Go port does not implement yet. Previously these options were accepted
// and silently ignored, which produced no output and looked like success.
func checkUnsupported(params *Params) error {
	var missing []string

	// --ldhat / --ldhat-geno / --ldhelmet: upstream errors at
	// parameters.cpp:717 unless exactly one chromosome is selected via
	// --chr. Mirror that check here so misuse fails fast rather than
	// producing a malformed file.
	if (params.LDhat || params.LDhatGeno || params.LDhelmet) && params.Chr == "" {
		return fmt.Errorf("Require a chromosome (--chr) when outputting LDhat format.")
	}

	if len(missing) > 0 {
		return fmt.Errorf("not implemented in this Go port yet: %s (see tools/vcftools/ROADMAP.md)", strings.Join(missing, ", "))
	}
	return nil
}
