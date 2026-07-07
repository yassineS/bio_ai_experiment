package main

import (
	"bytes"
	"os"
	"strconv"
	"strings"
	"testing"
)

// bigSampleVCF returns a sorted VCF large enough to span several BGZF blocks,
// so the parallel decode path is genuinely exercised.
func bigSampleVCF() string {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.2\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	for pos := 100; pos < 100+5000; pos++ {
		sb.WriteString("chr1\t")
		sb.WriteString(strconv.Itoa(pos))
		sb.WriteString("\t.\tA\tT\t.\t.\tDP=42;AF=0.5\n")
	}
	return sb.String()
}

// TestTabix_ThreadsByteIdentical verifies -@/--threads only parallelises the
// BGZF (de)compression: the built .tbi index and the reheadered stdout stream
// are byte-for-byte identical whether tabix runs single-threaded or with a
// worker pool. Query/list-chroms are deliberately left single-threaded.
func TestTabix_ThreadsByteIdentical(t *testing.T) {
	vcf := bigSampleVCF()

	// --- build: .tbi identical for -@1 and -@4 ---
	buildTBI := func(threads int) []byte {
		dir := t.TempDir()
		gz := writeBGZF(t, dir, "in.vcf.gz", vcf)
		var out, errb bytes.Buffer
		code := run([]string{"-p", "vcf", "-@", strconv.Itoa(threads), gz}, nil, &out, &errb)
		if code != 0 {
			t.Fatalf("build -@%d exit=%d stderr=%s", threads, code, errb.String())
		}
		b, err := os.ReadFile(gz + ".tbi")
		if err != nil {
			t.Fatalf("read tbi -@%d: %v", threads, err)
		}
		return b
	}
	if one, many := buildTBI(1), buildTBI(4); !bytes.Equal(one, many) {
		t.Errorf(".tbi differs between -@1 (%d bytes) and -@4 (%d bytes)", len(one), len(many))
	}

	// --- reheader: emitted bgzipped stdout identical for -@1 and -@4 ---
	dir := t.TempDir()
	gz := writeBGZF(t, dir, "in.vcf.gz", vcf)
	newHdr := writePlain(t, dir, "hdr.txt",
		"##fileformat=VCFv4.2\n##NEW=yes\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	reheader := func(threads int) []byte {
		var out, errb bytes.Buffer
		code := run([]string{"-p", "vcf", "-@", strconv.Itoa(threads), "-r", newHdr, gz}, nil, &out, &errb)
		if code != 0 {
			t.Fatalf("reheader -@%d exit=%d stderr=%s", threads, code, errb.String())
		}
		return out.Bytes()
	}
	if one, many := reheader(1), reheader(4); !bytes.Equal(one, many) {
		t.Errorf("reheader stdout differs between -@1 (%d bytes) and -@4 (%d bytes)", len(one), len(many))
	}
}
