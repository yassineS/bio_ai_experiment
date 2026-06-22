package roundtrip

// This file adds the explicitly *bidirectional* container-format interop
// checks. The base round-trip checks in roundtrip.go prove a single tool can
// encode then decode a format (and, where an upstream binary is present,
// cross-check the result). The checks here go further: for every container
// format they prove BOTH directions of producer/consumer interop on
// multi-contig fixtures —
//
//  1. ours-writes / upstream-reads : our binary produces the file, the UPSTREAM
//     binary reads it back, and the decoded payload (provenance-stripped)
//     equals the source; and
//  2. upstream-writes / ours-reads : the upstream binary produces the file, OUR
//     binary reads it back, and the decoded payload equals the source.
//
// They additionally prove *index* interop: our .bai/.csi/.tbi must let the
// upstream tool answer a region query, and upstream's index must let ours
// answer the same query (the queried records must match).
//
// Every check needs the upstream binary; when it is unavailable the check
// returns Skip (never Fail). The fixtures are the shared multi-contig set
// (≥2 contigs at every scale), so cross-contig bins, RNEXT, and
// coordinate-sort ordering are all exercised.

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
	"github.com/yassineS/bio_ai_experiment/pipeline/runner"
)

// interopChecks runs every bidirectional interop check and returns their
// results. Each returns Skip (not Fail) when its upstream binary is missing.
func (e *env) interopChecks() []Result {
	return []Result{
		e.bgzfInterop(),
		e.bamInterop(),
		e.cramInterop(),
		e.vcfGzInterop(),
		e.bcfInterop(),
		e.fastqInterop(),
		e.baiInterop(),
		e.csiInterop(),
		e.tbiInterop(),
	}
}

// samRecords decodes an alignment file (BAM/CRAM/SAM) through bin to normalized
// SAM text (records sorted, @PG/provenance stripped) so two producers compare
// equal despite cosmetic header/provenance differences. ref may be "" for BAM.
// trailing is the path optionally followed by a region (e.g. "chr2").
func samRecords(bin, ref string, trailing ...string) (string, error) {
	args := []string{"view", "-h"}
	if ref != "" {
		args = append(args, "-T", ref)
	}
	args = append(args, trailing...)
	out, err := runCmd(bin, args...)
	if err != nil {
		return "", err
	}
	return upstream.NormalizeSAM(string(out), true), nil
}

// sortAuxTags canonicalises a normalized-SAM block so two producers that emit
// the same optional (aux) tags in a different order compare equal. The order of
// SAM optional fields (columns 12+) is not semantically significant, and CRAM
// encoders legitimately differ in it (e.g. our encoder orders RG/MD/NM
// differently from upstream's). It sorts the aux fields of every record line
// while leaving header lines and the 11 mandatory columns untouched.
func sortAuxTags(sam string) string {
	lines := strings.Split(sam, "\n")
	for i, ln := range lines {
		if ln == "" || strings.HasPrefix(ln, "@") {
			continue
		}
		f := strings.Split(ln, "\t")
		if len(f) <= 12 {
			continue // no, or a single, aux tag: nothing to reorder
		}
		aux := f[11:]
		sort.Strings(aux)
		lines[i] = strings.Join(append(f[:11:11], aux...), "\t")
	}
	return strings.Join(lines, "\n")
}

// vcfRecords decodes a VCF/BCF/VCF.gz through the bcftools-style bin to
// provenance-stripped text so two producers compare equal.
func vcfRecords(bin, path string) ([]byte, error) {
	out, err := runCmd(bin, "view", path)
	if err != nil {
		return nil, err
	}
	return runner.StripProvenance(out), nil
}

// bgzfInterop: BGZF block-gzip interop, both directions. Our bgzip-compressed
// payload must decompress byte-identically under upstream bgzip, and upstream's
// bgzip output must decompress byte-identically under ours.
func (e *env) bgzfInterop() Result {
	const name, format = "bgzf-interop-bidirectional", "BGZF"
	src := e.path("vcf_plain")
	if src == "" {
		return skip(name, format, "no vcf_plain fixture")
	}
	our, err := e.our("bgzip")
	if err != nil {
		return skip(name, format, "our bgzip unavailable: "+err.Error())
	}
	up, err := e.up("bgzip")
	if err != nil {
		return skip(name, format, "upstream bgzip unavailable: "+err.Error())
	}
	raw, err := os.ReadFile(src)
	if err != nil {
		return fail(name, format, err)
	}
	// ours-writes / upstream-reads.
	ourGz := e.out("bgzf.our.gz")
	if err := runToFile(ourGz, our, "-c", src); err != nil {
		return fail(name, format, err)
	}
	if back, err := runCmd(up, "-dc", ourGz); err != nil {
		return fail(name, format, fmt.Errorf("upstream decode of our bgzip: %w", err))
	} else if !bytes.Equal(raw, back) {
		return fail(name, format, fmt.Errorf("upstream read of our BGZF differs from source"))
	}
	// upstream-writes / ours-reads.
	upGz := e.out("bgzf.up.gz")
	if err := runToFile(upGz, up, "-c", src); err != nil {
		return fail(name, format, err)
	}
	if back, err := runCmd(our, "-dc", upGz); err != nil {
		return fail(name, format, fmt.Errorf("our decode of upstream bgzip: %w", err))
	} else if !bytes.Equal(raw, back) {
		return fail(name, format, fmt.Errorf("our read of upstream BGZF differs from source"))
	}
	return ok(name, format, "BGZF interop both ways (ours↔upstream) byte-identical")
}

// bamInterop: BAM write/read interop, both directions, on the multi-contig
// fixture. The source records are the decoded fixture; each side's BAM is
// re-decoded and compared (normalized SAM) against that source.
func (e *env) bamInterop() Result {
	const name, format = "bam-interop-bidirectional", "BAM"
	src := e.path("bam")
	if src == "" {
		return skip(name, format, "no bam fixture")
	}
	our, err := e.our("samtools")
	if err != nil {
		return skip(name, format, "our samtools unavailable: "+err.Error())
	}
	up, err := e.up("samtools")
	if err != nil {
		return skip(name, format, "upstream samtools unavailable: "+err.Error())
	}
	want, err := samRecords(our, "", src)
	if err != nil {
		return fail(name, format, err)
	}
	// ours-writes / upstream-reads.
	ourBam := e.out("bam.our.bam")
	if err := runToFile(ourBam, our, "view", "-b", src); err != nil {
		return fail(name, format, err)
	}
	if got, err := samRecords(up, "", ourBam); err != nil {
		return fail(name, format, fmt.Errorf("upstream read of our BAM: %w", err))
	} else if got != want {
		return fail(name, format, fmt.Errorf("upstream read of our BAM differs from source"))
	}
	// upstream-writes / ours-reads.
	upBam := e.out("bam.up.bam")
	if err := runToFile(upBam, up, "view", "-b", src); err != nil {
		return fail(name, format, err)
	}
	if got, err := samRecords(our, "", upBam); err != nil {
		return fail(name, format, fmt.Errorf("our read of upstream BAM: %w", err))
	} else if got != want {
		return fail(name, format, fmt.Errorf("our read of upstream BAM differs from source"))
	}
	return ok(name, format, "BAM interop both ways (ours↔upstream) records identical")
}

// cramInterop: CRAM write/read interop, both directions, against a shared
// reference. CRAM is reference-based and decoders legitimately differ in
// whether they auto-regenerate the MD/NM tags from the reference on read — that
// is a *reader* policy, not a property of the CRAM bytes. To isolate genuine
// container interop (can reader R consume a CRAM written by writer W?) the
// reader is held fixed while the writer varies:
//
//   - ours-writes / upstream-reads: the UPSTREAM reader must decode OUR CRAM to
//     exactly what it decodes from its OWN CRAM of the same source; and
//   - upstream-writes / ours-reads: OUR reader must decode the UPSTREAM CRAM to
//     exactly what it decodes from OUR OWN CRAM.
//
// Both directions thus prove the file is interchangeable for that reader, on
// the multi-contig fixture (cross-contig slices exercised).
func (e *env) cramInterop() Result {
	const name, format = "cram-interop-bidirectional", "CRAM"
	src, ref := e.path("bam"), e.path("fasta")
	if src == "" || ref == "" {
		return skip(name, format, "no bam/fasta fixture")
	}
	our, err := e.our("samtools")
	if err != nil {
		return skip(name, format, "our samtools unavailable: "+err.Error())
	}
	up, err := e.up("samtools")
	if err != nil {
		return skip(name, format, "upstream samtools unavailable: "+err.Error())
	}
	ourCram := e.out("cram.our.cram")
	upCram := e.out("cram.up.cram")
	if err := runToFile(ourCram, our, "view", "-C", "-T", ref, src); err != nil {
		return fail(name, format, fmt.Errorf("our CRAM encode: %w", err))
	}
	if err := runToFile(upCram, up, "view", "-C", "-T", ref, src); err != nil {
		return fail(name, format, fmt.Errorf("upstream CRAM encode: %w", err))
	}
	// ours-writes / upstream-reads: the upstream reader sees our CRAM == its own.
	upOnOurs, err := samRecords(up, ref, ourCram)
	if err != nil {
		return fail(name, format, fmt.Errorf("upstream read of our CRAM: %w", err))
	}
	upOnUps, err := samRecords(up, ref, upCram)
	if err != nil {
		return fail(name, format, fmt.Errorf("upstream read of upstream CRAM: %w", err))
	}
	if sortAuxTags(upOnOurs) != sortAuxTags(upOnUps) {
		return fail(name, format, fmt.Errorf("upstream reader: our CRAM decodes differently from upstream's CRAM"))
	}
	// upstream-writes / ours-reads: our reader sees the upstream CRAM == our own.
	ourOnUps, err := samRecords(our, ref, upCram)
	if err != nil {
		return fail(name, format, fmt.Errorf("our read of upstream CRAM: %w", err))
	}
	ourOnOurs, err := samRecords(our, ref, ourCram)
	if err != nil {
		return fail(name, format, fmt.Errorf("our read of our CRAM: %w", err))
	}
	if sortAuxTags(ourOnUps) != sortAuxTags(ourOnOurs) {
		return fail(name, format, fmt.Errorf("our reader: upstream CRAM decodes differently from our CRAM"))
	}
	return ok(name, format, "CRAM interop both ways (ours↔upstream) records identical per fixed reader")
}

// vcfGzInterop: bgzipped VCF (bcftools -O z) write/read interop, both
// directions, on the multi-contig multi-sample VCF.
func (e *env) vcfGzInterop() Result {
	const name, format = "vcfgz-interop-bidirectional", "VCF.gz"
	src := e.path("vcf_plain")
	if src == "" {
		return skip(name, format, "no vcf_plain fixture")
	}
	return e.bcftoolsInterop(name, format, src, "z", "vcfgz")
}

// bcfInterop: binary BCF (bcftools -O b) write/read interop, both directions.
func (e *env) bcfInterop() Result {
	const name, format = "bcf-interop-bidirectional", "BCF"
	src := e.path("vcf_plain")
	if src == "" {
		return skip(name, format, "no vcf_plain fixture")
	}
	return e.bcftoolsInterop(name, format, src, "b", "bcf")
}

// bcftoolsInterop runs the shared ours↔upstream bcftools write/read interop for
// output type ot ("z" for VCF.gz, "b" for BCF). The decoded records (provenance
// stripped) on each side must equal the decoded source.
func (e *env) bcftoolsInterop(name, format, src, ot, tag string) Result {
	our, err := e.our("bcftools")
	if err != nil {
		return skip(name, format, "our bcftools unavailable: "+err.Error())
	}
	up, err := e.up("bcftools")
	if err != nil {
		return skip(name, format, "upstream bcftools unavailable: "+err.Error())
	}
	want, err := vcfRecords(our, src)
	if err != nil {
		return fail(name, format, err)
	}
	// ours-writes / upstream-reads.
	ourFile := e.out(tag + ".our")
	if err := runToFile(ourFile, our, "view", "-O"+ot, src); err != nil {
		return fail(name, format, err)
	}
	if got, err := vcfRecords(up, ourFile); err != nil {
		return fail(name, format, fmt.Errorf("upstream read of our %s: %w", format, err))
	} else if !bytes.Equal(got, want) {
		return fail(name, format, fmt.Errorf("upstream read of our %s differs from source", format))
	}
	// upstream-writes / ours-reads.
	upFile := e.out(tag + ".up")
	if err := runToFile(upFile, up, "view", "-O"+ot, src); err != nil {
		return fail(name, format, err)
	}
	if got, err := vcfRecords(our, upFile); err != nil {
		return fail(name, format, fmt.Errorf("our read of upstream %s: %w", format, err))
	} else if !bytes.Equal(got, want) {
		return fail(name, format, fmt.Errorf("our read of upstream %s differs from source", format))
	}
	return ok(name, format, fmt.Sprintf("%s interop both ways (ours↔upstream) records identical", format))
}

// fastqInterop: FASTQ via BGZF, both directions. Our bgzip-compressed FASTQ
// must decompress byte-identically under upstream bgzip, and vice versa.
func (e *env) fastqInterop() Result {
	const name, format = "fastq-interop-bidirectional", "FASTQ"
	src := e.path("fastq")
	if src == "" {
		return skip(name, format, "no fastq fixture")
	}
	our, err := e.our("bgzip")
	if err != nil {
		return skip(name, format, "our bgzip unavailable: "+err.Error())
	}
	up, err := e.up("bgzip")
	if err != nil {
		return skip(name, format, "upstream bgzip unavailable: "+err.Error())
	}
	raw, err := os.ReadFile(src)
	if err != nil {
		return fail(name, format, err)
	}
	// ours-writes / upstream-reads.
	ourGz := e.out("fastq.our.gz")
	if err := runToFile(ourGz, our, "-c", src); err != nil {
		return fail(name, format, err)
	}
	if back, err := runCmd(up, "-dc", ourGz); err != nil {
		return fail(name, format, fmt.Errorf("upstream decode of our FASTQ.gz: %w", err))
	} else if !bytes.Equal(raw, back) {
		return fail(name, format, fmt.Errorf("upstream read of our FASTQ.gz differs from source"))
	}
	// upstream-writes / ours-reads.
	upGz := e.out("fastq.up.gz")
	if err := runToFile(upGz, up, "-c", src); err != nil {
		return fail(name, format, err)
	}
	if back, err := runCmd(our, "-dc", upGz); err != nil {
		return fail(name, format, fmt.Errorf("our decode of upstream FASTQ.gz: %w", err))
	} else if !bytes.Equal(raw, back) {
		return fail(name, format, fmt.Errorf("our read of upstream FASTQ.gz differs from source"))
	}
	return ok(name, format, "FASTQ-over-BGZF interop both ways (ours↔upstream) byte-identical")
}

// baiInterop: .bai index interop. The region query is answered across a contig
// boundary (the second contig) so cross-contig bins are exercised. Our index
// must let upstream answer the query, and upstream's index must let ours.
func (e *env) baiInterop() Result {
	return e.indexInterop("bai-index-interop", "BAI", "bai", []string{"index", "-b"})
}

// csiInterop: .csi index interop (coordinate-sorted BAM, CSI binning).
func (e *env) csiInterop() Result {
	return e.indexInterop("csi-index-interop", "CSI", "csi", []string{"index", "-c"})
}

// indexInterop runs the shared samtools .bai/.csi index interop. It copies the
// fixture BAM into the scratch dir (so the index lands next to a writable copy),
// then for each index producer (ours, upstream) it builds the index with the
// given indexArgs and has the *other* tool answer a region query on the second
// contig; the queried record sets must match the reference query.
func (e *env) indexInterop(name, format, suffix string, indexArgs []string) Result {
	src := e.path("bam")
	if src == "" {
		return skip(name, format, "no bam fixture")
	}
	our, err := e.our("samtools")
	if err != nil {
		return skip(name, format, "our samtools unavailable: "+err.Error())
	}
	up, err := e.up("samtools")
	if err != nil {
		return skip(name, format, "upstream samtools unavailable: "+err.Error())
	}
	region, err := secondContig(up, src)
	if err != nil {
		return fail(name, format, err)
	}
	// Reference answer: query the original (already-indexed) fixture.
	wantRecords, err := samRecords(up, "", src, region)
	if err != nil {
		// Fall back to building a fresh reference query below if the fixture's own
		// index is absent; but the fixture ships with one, so treat this as fatal.
		return fail(name, format, fmt.Errorf("reference region query: %w", err))
	}
	check := func(indexer, querier string, dir string) error {
		bam := e.out(fmt.Sprintf("%s.%s.bam", suffix, dir))
		if err := copyFile(bam, src); err != nil {
			return err
		}
		// Remove any co-copied index so only the freshly built one is used.
		_ = os.Remove(bam + ".bai")
		_ = os.Remove(bam + ".csi")
		args := append(append([]string{}, indexArgs...), bam)
		if _, err := runCmd(indexer, args...); err != nil {
			return fmt.Errorf("%s build: %w", dir, err)
		}
		if _, err := os.Stat(bam + "." + suffix); err != nil {
			return fmt.Errorf("%s: expected %s index not written: %w", dir, suffix, err)
		}
		got, err := samRecords(querier, "", bam, region)
		if err != nil {
			return fmt.Errorf("%s query: %w", dir, err)
		}
		if got != wantRecords {
			return fmt.Errorf("%s: region %s records differ from reference", dir, region)
		}
		return nil
	}
	// ours-indexes / upstream-queries, then upstream-indexes / ours-queries.
	if err := check(our, up, "ours-writes-upstream-reads"); err != nil {
		return fail(name, format, err)
	}
	if err := check(up, our, "upstream-writes-ours-reads"); err != nil {
		return fail(name, format, err)
	}
	return ok(name, format, fmt.Sprintf("%s index interop both ways: region %s queried identically", format, region))
}

// tbiInterop: .tbi (tabix) index interop on the bgzipped multi-contig VCF. Our
// tabix index must let upstream tabix answer a region query and vice versa.
func (e *env) tbiInterop() Result {
	const name, format = "tbi-index-interop", "TBI"
	src := e.path("vcf_plain")
	if src == "" {
		return skip(name, format, "no vcf_plain fixture")
	}
	ourBgzip, err := e.our("bgzip")
	if err != nil {
		return skip(name, format, "our bgzip unavailable: "+err.Error())
	}
	ourTabix, err := e.our("tabix")
	if err != nil {
		return skip(name, format, "our tabix unavailable: "+err.Error())
	}
	upTabix, err := e.up("tabix")
	if err != nil {
		return skip(name, format, "upstream tabix unavailable: "+err.Error())
	}
	// bgzip the plain VCF once (the index is what we exercise, not the framing —
	// BGZF interop is covered by bgzfInterop).
	vcfgz := e.out("tbi.vcf.gz")
	if err := runToFile(vcfgz, ourBgzip, "-c", src); err != nil {
		return fail(name, format, err)
	}
	// Pick the second contig from the VCF body so the query crosses a contig.
	region, err := secondContigVCF(src)
	if err != nil {
		return fail(name, format, err)
	}
	check := func(indexer, querier, dir string) error {
		_ = os.Remove(vcfgz + ".tbi")
		if _, err := runCmd(indexer, "-p", "vcf", "-f", vcfgz); err != nil {
			return fmt.Errorf("%s tabix index: %w", dir, err)
		}
		got, err := runCmd(querier, vcfgz, region)
		if err != nil {
			return fmt.Errorf("%s tabix query: %w", dir, err)
		}
		want, err := runCmd(indexer, vcfgz, region)
		if err != nil {
			return fmt.Errorf("%s tabix self-query: %w", dir, err)
		}
		if !bytes.Equal(runner.StripProvenance(got), runner.StripProvenance(want)) {
			return fmt.Errorf("%s: region %s records differ", dir, region)
		}
		return nil
	}
	if err := check(ourTabix, upTabix, "ours-writes-upstream-reads"); err != nil {
		return fail(name, format, err)
	}
	if err := check(upTabix, ourTabix, "upstream-writes-ours-reads"); err != nil {
		return fail(name, format, err)
	}
	return ok(name, format, fmt.Sprintf("TBI index interop both ways: region %s queried identically", region))
}

// secondContig returns the name of the second @SQ contig in an alignment file
// (read through bin), so region queries cross a contig boundary on multi-contig
// fixtures. It falls back to the first contig if only one exists.
func secondContig(bin, path string) (string, error) {
	out, err := runCmd(bin, "view", "-H", path)
	if err != nil {
		return "", err
	}
	var names []string
	for _, ln := range bytes.Split(out, []byte("\n")) {
		if !bytes.HasPrefix(ln, []byte("@SQ")) {
			continue
		}
		for _, f := range bytes.Split(ln, []byte("\t")) {
			if bytes.HasPrefix(f, []byte("SN:")) {
				names = append(names, string(f[3:]))
			}
		}
	}
	if len(names) == 0 {
		return "", fmt.Errorf("no @SQ contigs in %s", filepath.Base(path))
	}
	if len(names) >= 2 {
		return names[1], nil
	}
	return names[0], nil
}

// secondContigVCF returns the chromosome of the second distinct contig seen in
// a plain VCF's data rows (so a tabix region query crosses a contig). It falls
// back to the first if only one exists.
func secondContigVCF(path string) (string, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return "", err
	}
	var first, second string
	for _, ln := range bytes.Split(b, []byte("\n")) {
		if len(ln) == 0 || ln[0] == '#' {
			continue
		}
		tab := bytes.IndexByte(ln, '\t')
		if tab <= 0 {
			continue
		}
		chrom := string(ln[:tab])
		if first == "" {
			first = chrom
			continue
		}
		if chrom != first {
			second = chrom
			break
		}
	}
	if second != "" {
		return second, nil
	}
	if first != "" {
		return first, nil
	}
	return "", fmt.Errorf("no data rows in %s", filepath.Base(path))
}

// copyFile copies src to dst (overwriting), preserving content only.
func copyFile(dst, src string) error {
	b, err := os.ReadFile(src)
	if err != nil {
		return err
	}
	return os.WriteFile(dst, b, 0o644)
}
