package main

import (
	"bufio"
	"bytes"
	"os"
	"strings"
)

// describeInputs builds the report's input table: for each provided file it
// records the role, absolute path, size, and (BAM/VCF) contig count. The contig
// count is read by inspecting the file header — for the VCF via the bgzip-gzip
// stream's ##contig lines, for the BAM by asking a resolved binary for its @SQ
// header (cheap, and avoids re-implementing BAM header parsing here). A failure
// to count contigs leaves Contigs at 0 (rendered "—"); it never aborts the run.
func describeInputs(cfg config) []inputInfo {
	var out []inputInfo
	add := func(role, path string, contigs int) {
		if path == "" {
			return
		}
		var size int64
		if st, err := os.Stat(path); err == nil {
			size = st.Size()
		}
		out = append(out, inputInfo{Role: role, Path: path, SizeB: size, Contigs: contigs})
	}
	add("ref", cfg.in.ref, faiContigs(cfg.in.ref))
	add("bam", cfg.in.bam, bamContigs(cfg))
	add("vcf", cfg.in.vcf, vcfContigs(cfg.in.vcf))
	return out
}

// faiContigs counts contigs from a FASTA index (.fai) sibling when present
// (one line per contig). It is a cheap, binary-free count for the reference.
func faiContigs(ref string) int {
	if ref == "" {
		return 0
	}
	b, err := os.ReadFile(ref + ".fai")
	if err != nil {
		return 0
	}
	return countNonEmptyLines(b)
}

// bamContigs counts the BAM's reference sequences (@SQ lines) by asking a
// resolved samtools (ours, else upstream) for its header. Returns 0 if neither
// binary is available or the call fails.
func bamContigs(cfg config) int {
	if cfg.in.bam == "" {
		return 0
	}
	bin := cfg.bins.oursSamtools
	if bin == "" {
		bin = cfg.bins.upSamtools
	}
	if bin == "" {
		return 0
	}
	// The header is small, so buffering it here is fine (unlike the body-producing
	// battery cells, which stream).
	var hdr bytes.Buffer
	_, err := runOnce(bin, []string{"view", "-H", cfg.in.bam}, "", nil, &hdr)
	if err != nil {
		return 0
	}
	return countPrefixLines(hdr.Bytes(), "@SQ\t")
}

// vcfContigs counts ##contig header lines in a VCF, transparently decompressing
// a bgzipped (.gz) file. Returns 0 on any error.
func vcfContigs(vcf string) int {
	if vcf == "" {
		return 0
	}
	b, err := readVCFHeader(vcf)
	if err != nil {
		return 0
	}
	return countPrefixLines(b, "##contig")
}

// readVCFHeader returns the header bytes of a VCF, decompressing .gz transparently
// and stopping at the #CHROM column line so a huge body is not read into memory.
func readVCFHeader(path string) ([]byte, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var rd *bufio.Reader
	if strings.HasSuffix(path, ".gz") {
		gz, gerr := newGzipReader(f)
		if gerr != nil {
			return nil, gerr
		}
		defer gz.Close()
		rd = bufio.NewReader(gz)
	} else {
		rd = bufio.NewReader(f)
	}
	var buf bytes.Buffer
	for {
		line, err := rd.ReadBytes('\n')
		buf.Write(line)
		if len(line) > 0 && line[0] != '#' {
			break
		}
		if bytes.HasPrefix(line, []byte("#CHROM")) {
			break
		}
		if err != nil {
			break
		}
	}
	return buf.Bytes(), nil
}

func countNonEmptyLines(b []byte) int {
	n := 0
	for _, ln := range bytes.Split(b, []byte("\n")) {
		if len(bytes.TrimSpace(ln)) > 0 {
			n++
		}
	}
	return n
}

func countPrefixLines(b []byte, prefix string) int {
	p := []byte(prefix)
	n := 0
	for _, ln := range bytes.Split(b, []byte("\n")) {
		if bytes.HasPrefix(ln, p) {
			n++
		}
	}
	return n
}
