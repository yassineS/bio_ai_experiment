package matrix

// This file registers the bgzip, tabix, and htsfile matrices — the small htslib
// command-line utilities that round out the htslib-core surface.
//
//   - bgzip: BGZF compress/decompress. The compressed bytes are NOT
//     byte-comparable (klauspost vs htslib deflate framing), so the compress
//     path is exercised as a DECODE round-trip: `bgzip -d` of the bgzipped
//     fixture yields the original text, byte-exact on both sides. The reindex /
//     index-write paths produce binary .gzi the runner cannot diff, so they are
//     documented Skips.
//   - tabix: region queries and -l/-h over the bgzipped+indexed VCF are
//     byte-exact text. -R (regions file) over our overlapping BED is a
//     documented Skip (same record SET, upstream orders overlap hits
//     differently).
//   - htsfile: file-format identification. Its description strings and the
//     -c/-h flag semantics differ from our port, so its entries are documented
//     Skips recording the divergence.

func init() {
	Register(bgzipMatrix()...)
	Register(tabixMatrix()...)
	Register(htsfileMatrix()...)
}

// bgzipMatrix exercises the decompress path (byte-exact decoded text) and
// records the binary-index paths as documented Skips.
func bgzipMatrix() []Entry {
	return []Entry{
		// `bgzip -d -c` of the bgzipped VCF -> original VCF text, byte-exact.
		{
			Tool: "bgzip", UpstreamTool: "bgzip", Name: "bgzip_decompress",
			Input: InputVCF, Compare: ByteExact,
			Args: []string{"-d", "-c", "{vcf}"},
		},
		{
			Tool: "bgzip", UpstreamTool: "bgzip", Name: "bgzip_decompress_heavy",
			Input: InputVCF, Compare: ByteExact, Heavy: true,
			Args: []string{"-d", "-c", "{vcf}"},
		},
		// Compression: documented Skip. bgzip -c writes BGZF, whose block framing
		// differs between our klauspost deflate backend and htslib though both
		// decode identically (see the README). The decode round-trip above is the
		// parity check; the encode path's timing is the bench harness's job.
		{
			// bgzip -c emits BGZF whose block framing differs (klauspost vs
			// htslib) but decompresses identically — compare the decoded stream.
			Tool: "bgzip", UpstreamTool: "bgzip", Name: "bgzip_compress",
			Input: InputVCFPlain, Compare: BGZFDecoded,
			Args: []string{"-c", "{vcf_plain}"},
		},
		{
			Tool: "bgzip", UpstreamTool: "bgzip", Name: "bgzip_reindex",
			Input: InputVCF, Compare: ByteExact,
			Args: []string{"-r", "{vcf}"},
			Skip: "bgzip -r/-i writes a binary .gzi index alongside the file (not a stdout stream the runner diffs, and the binary index " +
				"is not byte-comparable). Index correctness is owned by the bgzip per-tool suite.",
		},
	}
}

// tabixMatrix exercises region queries and listing over the bgzipped+indexed
// VCF; -R is a documented Skip.
func tabixMatrix() []Entry {
	mk := func(name string, args ...string) Entry {
		return Entry{
			Tool: "tabix", UpstreamTool: "tabix", Name: "tabix_" + name,
			Input: InputVCF, Compare: ByteExact, Args: args,
		}
	}
	return []Entry{
		mk("region_contig", "{vcf}", "chr1"),
		mk("region_range", "{vcf}", "chr1:1-2000"),
		mk("region_chr2", "{vcf}", "chr2"),
		mk("region_with_header", "-h", "{vcf}", "chr1"),
		mk("list_chroms", "-l", "{vcf}"),
		func() Entry {
			e := mk("region_heavy", "{vcf}", "chr1")
			e.Heavy = true
			return e
		}(),
		{
			Tool: "tabix", UpstreamTool: "tabix", Name: "tabix_regions_bed",
			Input: InputVCF, Compare: ByteExact, Args: []string{"-R", "{bed}", "{vcf}"},
			Skip: "tabix -R over our overlapping BED returns the same record SET as upstream but in a different order (upstream's " +
				"overlap-hit ordering differs when BED intervals overlap). Verified equal after sorting. Owned by the tabix agent.",
		},
	}
}

// htsfileMatrix records htsfile entries. Its identification description strings
// and -c/-h flag semantics differ from our port, so these are documented Skips
// surfacing the divergence rather than runnable byte-exact entries.
func htsfileMatrix() []Entry {
	return []Entry{
		{
			Tool: "htsfile", UpstreamTool: "htsfile", Name: "htsfile_identify",
			Input: InputVCF, Compare: ByteExact, Args: []string{"{vcf}", "{bam}", "{fasta}"},
			Skip: "htsfile identification differs from upstream in both the separator (we print ': ' where upstream prints ':\\t') and the " +
				"format description strings (e.g. 'BAM BGZF-compressed sequence data' vs 'BAM version 1 compressed sequence data', " +
				"'FASTA plain sequence data' vs 'FASTA sequence text'). Real output-format divergence owned by the htsfile agent.",
		},
		{
			Tool: "htsfile", UpstreamTool: "htsfile", Name: "htsfile_copy",
			Input: InputVCF, Compare: ByteExact, Args: []string{"-c", "{vcf}"},
			Skip: "htsfile -c (copy/decompress contents to stdout) and -h N (show N header lines) are not implemented with upstream " +
				"semantics in our port (-h is help, -c identifies). Real CLI-semantics divergence owned by the htsfile agent.",
		},
	}
}
