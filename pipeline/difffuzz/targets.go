package difffuzz

// DefaultTargets returns a curated set of fuzz targets spanning the formats and
// tool families the manuscript cares about (VCF, SAM/BAM-text, BED, gzip). Each
// target's two sides take the SAME input via the "{in}" token (a temp file)
// unless noted; both binaries are vendored locally so these run offline.
//
// The targets are deliberately ones whose CLI shape matches between our port
// and upstream so a single Args template drives both, keeping the harness's
// comparison apples-to-apples.
func DefaultTargets() []Target {
	return []Target{
		{
			// Our bcftools is a subcommand dispatcher like upstream, so BOTH sides
			// take the "view" subcommand (UsesSubcommand: true).
			Name:           "bcftools-view",
			Tool:           "bcftools",
			Subcommand:     "view",
			UsesSubcommand: true,
			Args:           []string{"{in}"},
			Format:         FormatVCF,
			SeedFixture:    "vcf_plain",
		},
		{
			Name:           "bcftools-query",
			Tool:           "bcftools",
			Subcommand:     "query",
			UsesSubcommand: true,
			Args:           []string{"-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", "{in}"},
			Format:         FormatVCF,
			SeedFixture:    "vcf_plain",
		},
		{
			Name:           "samtools-view",
			Tool:           "samtools",
			Subcommand:     "view",
			UsesSubcommand: true,
			Args:           []string{"{in}"},
			Format:         FormatSAM,
			// No SAM-text fixture in the manifest; seed via the structured SAM
			// generator + raw/mutation strategies.
			SeedFixture: "",
		},
		{
			Name:           "samtools-flagstat",
			Tool:           "samtools",
			Subcommand:     "flagstat",
			UsesSubcommand: true,
			Args:           []string{"{in}"},
			Format:         FormatSAM,
			SeedFixture:    "",
		},
		{
			// our bedmerge == upstream `bedtools merge`. Our binary is the
			// subcommand itself, so UsesSubcommand stays false and only the
			// upstream side gets the "merge" token.
			Name:        "bedtools-merge",
			Tool:        "bedmerge",
			UpstreamKey: "bedtools",
			Subcommand:  "merge",
			Args:        []string{"-i", "{in}"},
			Format:      FormatBED,
			SeedFixture: "bed",
		},
		{
			Name:        "bedtools-intersect",
			Tool:        "bedintersect",
			UpstreamKey: "bedtools",
			Subcommand:  "intersect",
			// Self-intersection: a single fuzzed file used for both -a and -b.
			Args:        []string{"-a", "{in}", "-b", "{in}"},
			Format:      FormatBED,
			SeedFixture: "bed",
		},
		{
			// bgzip -d -c reads a gzip stream on stdin and writes plain bytes.
			// Mutated/raw inputs exercise the decompressor's error parity.
			Name:        "bgzip-decompress",
			Tool:        "bgzip",
			Args:        []string{"-d", "-c"},
			Format:      FormatGzip,
			SeedFixture: "vcf", // the bgzipped VCF is a valid gzip member stream
		},
	}
}

// QuickTargets returns a small subset suitable for a seconds-long -quick run.
func QuickTargets() []Target {
	all := DefaultTargets()
	pick := map[string]bool{"bcftools-view": true, "samtools-flagstat": true, "bedtools-merge": true}
	var out []Target
	for _, t := range all {
		if pick[t.Name] {
			out = append(out, t)
		}
	}
	return out
}
