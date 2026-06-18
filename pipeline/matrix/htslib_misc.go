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
//   - htsfile: file-format identification (byte-exact) and the -c/--view
//     format-aware re-serialisation, which for a VCF/BCF round-trips through the
//     VCF reader/writer to canonical VCF text (byte-exact, provenance-stripped).

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
			// bgzip -r writes a binary .gzi index alongside the file. Copy the
			// bgzipped VCF into each side's out dir, reindex the copy, and
			// byte-compare the sibling .gzi — our index is byte-identical to
			// htslib's.
			Tool: "bgzip", UpstreamTool: "bgzip", Name: "bgzip_reindex",
			Input: InputVCF, Compare: ByteExact,
			Args:        []string{"-r", "{out}.vcf.gz"},
			CopyToOut:   map[string]string{"vcf": ".vcf.gz"},
			OutputFiles: []string{".vcf.gz.gzi"},
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
			// tabix -R now orders the regions exactly as htslib's regidx does —
			// chromosomes in first-appearance order, then by (start asc, end
			// DESC), so two regions sharing a start emit the longer one first.
			// Byte-exact against upstream over the overlapping-interval fixture.
			Tool: "tabix", UpstreamTool: "tabix", Name: "tabix_regions_bed",
			Input: InputVCF, Compare: ByteExact, Args: []string{"-R", "{bed}", "{vcf}"},
		},
	}
}

// htsfileMatrix records htsfile entries. Its identification description strings
// and -c/-h flag semantics differ from our port, so these are documented Skips
// surfacing the divergence rather than runnable byte-exact entries.
func htsfileMatrix() []Entry {
	return []Entry{
		{
			// Identification now matches hts_format_description byte-for-byte:
			// the path/description separator is a TAB, and the description
			// strings (version word, BGZF vs "compressed" for BAM/BCF, the
			// " text"/" data" suffix) follow upstream exactly.
			Tool: "htsfile", UpstreamTool: "htsfile", Name: "htsfile_identify",
			Input: InputVCF, Compare: ByteExact, Args: []string{"{vcf}", "{bam}", "{fasta}"},
		},
		{
			// htsfile -c/--view re-serialises a VCF/BCF through the format reader
			// (the htslib round-trip); our port does the same via the VCF
			// reader/writer. The implicit ##FILTER=<ID=PASS> line htslib injects
			// is provenance boilerplate stripped by the comparison, so the
			// decoded VCF text matches byte-for-byte over the bgzipped fixture.
			Tool: "htsfile", UpstreamTool: "htsfile", Name: "htsfile_copy",
			Input: InputVCF, Compare: ByteExact, Args: []string{"-c", "{vcf}"},
		},
	}
}
