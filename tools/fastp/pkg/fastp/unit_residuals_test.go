package fastp

// Binary-free unit tests (TestUnit*) for the pure helpers behind the two
// fastp residuals: the split read->file-index assignment / zero-padded
// numbering, and the JSON sub-structure builders (per-read curves, k-mer
// histogram, insert-size block, sequencing string). These pass with the
// upstream submodule UNPOPULATED — they exercise the deterministic Go logic
// directly, with hand-computed expectations derived from upstream's algorithm.

import (
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// TestUnitSplitFileName checks the zero-padded "<NUM>.<base>" naming
// (upstream ThreadConfig::initWriterForSplit), 1-based and zero-padded to the
// configured digit width.
func TestUnitSplitFileName(t *testing.T) {
	cases := []struct {
		base   string
		index  int
		digits int
		want   string
	}{
		{"out.fq", 0, 4, "0001.out.fq"},
		{"out.fq", 9, 4, "0010.out.fq"},
		{"out.fq", 123, 4, "0124.out.fq"},
		{"out.fq", 0, 0, "1.out.fq"},
		{"out.fq", 0, 1, "1.out.fq"},
		{"out.fq", 0, 3, "001.out.fq"},
		{"sample.R1.fastq.gz", 4, 4, "0005.sample.R1.fastq.gz"},
	}
	for _, c := range cases {
		if got := splitFileName(c.base, c.index, c.digits); got != c.want {
			t.Errorf("splitFileName(%q,%d,%d) = %q, want %q", c.base, c.index, c.digits, got, c.want)
		}
	}
}

// TestUnitResolveSplitConfig checks the per-file size, file-count cap, and
// worker-thread clamping for --split N and --split_by_lines L.
func TestUnitResolveSplitConfig(t *testing.T) {
	// --split N: size = total/N, capped file count N, threads clamped to N.
	opts := DefaultProcessOptions()
	opts.SplitNumber = 4
	opts.Threads = 8
	cfg := resolveSplitConfig(opts, 5000)
	if cfg.Size != 1250 || cfg.Number != 4 || cfg.Digits != 4 || cfg.ByFileLines {
		t.Fatalf("byNumber cfg = %+v", cfg)
	}
	if cfg.Threads != 4 {
		t.Fatalf("byNumber threads clamp = %d, want 4", cfg.Threads)
	}
	// Fewer reads than files -> size clamps to 1.
	if c := resolveSplitConfig(opts, 2); c.Size != 1 {
		t.Fatalf("byNumber clamp size = %d, want 1", c.Size)
	}
	// --split_by_lines L: size = L/4, no cap, no thread clamp.
	opts2 := DefaultProcessOptions()
	opts2.SplitByLines = 4000
	opts2.Threads = 8
	cfg2 := resolveSplitConfig(opts2, 9999)
	if cfg2.Size != 1000 || cfg2.Number != 0 || !cfg2.ByFileLines || cfg2.Threads != 8 {
		t.Fatalf("byLines cfg = %+v", cfg2)
	}
}

// assignFor builds a splitWriter, announces n input positions all surviving,
// and returns the per-entry file index and the opened-file set.
func assignFor(cfg SplitConfig, n int) ([]int, []int) {
	sw := newSplitWriter("out.fq", cfg, fastq.Phred33)
	for i := 0; i < n; i++ {
		sw.SetInputPos(i)
		sw.Write(&fastq.Record{ID: "r", Sequence: []byte("ACGT"), Quality: []byte("IIII")})
	}
	return sw.assignFiles()
}

// TestUnitSplitAssignment checks the upstream pack/thread round-robin file
// assignment is deterministic for a fixed thread count and matches the
// hand-computed pack-to-file boundaries. (There is NO cross-thread invariant:
// upstream assigns pack i to thread i%threads, and each thread owns a strided
// set of files, so the read->file mapping genuinely depends on the thread
// count — see the package comment in split.go and the note on
// TestUnitSplitThreadAssignmentDiffers.)
func TestUnitSplitAssignment(t *testing.T) {
	// 1000 reads, byFileNumber size 250 (N=4). Packs of 256: [256,256,256,232].
	// Single thread: thread 0 owns all packs, rolling after each 256>=250:
	//   pack0->file0, pack1->file1, pack2->file2, pack3->file3.
	cfg1 := SplitConfig{Size: 250, Digits: 4, Number: 4, Threads: 1}
	files1, opened1 := assignFor(cfg1, 1000)
	wantOpened := []int{0, 1, 2, 3}
	if !reflect.DeepEqual(opened1, wantOpened) {
		t.Fatalf("w1 opened = %v, want %v", opened1, wantOpened)
	}
	// Single-thread: read p -> file p/256 (pack index), capped at N-1.
	for p, f := range files1 {
		want := p / 256
		if want > 3 {
			want = 3
		}
		if f != want {
			t.Fatalf("w1 read %d -> file %d, want %d", p, f, want)
		}
	}

	// 2 threads: pack i -> thread i%2; thread0 owns packs 0,2 (files 0 then 2),
	// thread1 owns packs 1,3 (files 1 then 3). Same set of files opened.
	cfg2 := cfg1
	cfg2.Threads = 2
	files2, opened2 := assignFor(cfg2, 1000)
	if !reflect.DeepEqual(opened2, wantOpened) {
		t.Fatalf("w2 opened = %v, want %v", opened2, wantOpened)
	}
	// w2 read->file: pack0(0..255)->file0, pack1(256..511)->file1,
	// pack2(512..767)->file2, pack3(768..999)->file3.
	for p, f := range files2 {
		pack := p / 256
		want := []int{0, 1, 2, 3}[pack]
		if f != want {
			t.Fatalf("w2 read %d -> file %d, want %d", p, f, want)
		}
	}
}

// TestUnitSplitThreadAssignmentDiffers documents (and pins) that the read->file
// assignment is thread-count-dependent, exactly as upstream fastp's: for
// --split N the per-thread strided rollover sends different reads to file 0
// depending on -w. This is intentional — our contract is byte-parity with
// upstream PER thread count, not a (non-existent) cross-thread invariant.
func TestUnitSplitThreadAssignmentDiffers(t *testing.T) {
	// 6000 reads, --split 4 (size = 6000/4 = 1500). Pack size 256.
	cfg := func(threads int) SplitConfig {
		return SplitConfig{Size: 1500, Digits: 4, Number: 4, Threads: threads}
	}
	f1, _ := assignFor(cfg(1), 6000)
	f4, _ := assignFor(cfg(4), 6000)
	// w1: thread0 fills file0 with packs until cur>=1500 -> 6 packs (1536
	// reads), so reads 0..1535 land in file0.
	if f1[1535] != 0 || f1[1536] != 1 {
		t.Fatalf("w1 file0 boundary wrong: f1[1535]=%d f1[1536]=%d", f1[1535], f1[1536])
	}
	// w4: thread0 owns packs 0,4,8,...; file0 gets pack0 (256) then more strided
	// packs, so read 1535 (in pack5, owned by thread1) is NOT in file0.
	if f4[1535] == 0 {
		t.Fatalf("w4 unexpectedly assigned read 1535 to file0 (should be thread-dependent)")
	}
	// The two assignments must therefore differ somewhere (matching upstream).
	differ := false
	for p := range f1 {
		if f1[p] != f4[p] {
			differ = true
			break
		}
	}
	if !differ {
		t.Fatal("expected thread-count-dependent assignment for --split N (matches upstream)")
	}
}

// TestUnitSplitEmptyTrailingFiles verifies byFileNumber materializes all N
// files (trailing empties) when the input is shorter than N files of records.
func TestUnitSplitEmptyTrailingFiles(t *testing.T) {
	cfg := SplitConfig{Size: 1, Digits: 4, Number: 4, Threads: 2}
	_, opened := assignFor(cfg, 3)
	want := []int{0, 1, 2, 3}
	if !reflect.DeepEqual(opened, want) {
		t.Fatalf("opened = %v, want %v (all N files materialized)", opened, want)
	}
}

// makeRecord builds a Record with Phred33 quality string from a base sequence
// and a per-base phred-score slice.
func makeRecord(seq string, q []int) *fastq.Record {
	qb := make([]byte, len(q))
	for i, v := range q {
		qb[i] = byte(v + 33)
	}
	return &fastq.Record{ID: "r", Sequence: []byte(seq), Quality: qb}
}

// TestUnitReadCurvesContentAndQuality checks the per-cycle content fractions,
// per-base quality means (including the mean-substitution for absent bases),
// the GC curve and the q40 tally against hand-computed values.
func TestUnitReadCurvesContentAndQuality(t *testing.T) {
	rc := &readCurves{}
	// Two reads, length 4. Qualities chosen to cross the q40 boundary (>=40).
	rc.stat(makeRecord("ACGT", []int{40, 41, 10, 20}), 33)
	rc.stat(makeRecord("AGGT", []int{38, 42, 39, 5}), 33)

	if rc.cycles() != 4 {
		t.Fatalf("cycles = %d, want 4", rc.cycles())
	}
	// q40: bases with phred>=40: read1 cycle0(40),cycle1(41); read2 cycle1(42) = 3.
	if rc.q40 != 3 {
		t.Fatalf("q40 = %d, want 3", rc.q40)
	}

	content := rc.contentCurves()
	// content_curves order/keys: A,T,C,G,N,GC.
	for _, k := range []string{"A", "T", "C", "G", "N", "GC"} {
		if _, ok := content[k]; !ok {
			t.Fatalf("content_curves missing key %q", k)
		}
	}
	// Cycle 0: both reads 'A' -> A=1.0, others 0.
	if content["A"][0] != 1.0 {
		t.Fatalf("content A[0] = %v, want 1.0", content["A"][0])
	}
	// Cycle 1: bases C,G -> C=0.5, G=0.5, GC=1.0.
	if content["C"][1] != 0.5 || content["G"][1] != 0.5 || content["GC"][1] != 1.0 {
		t.Fatalf("content cycle1 C=%v G=%v GC=%v", content["C"][1], content["G"][1], content["GC"][1])
	}

	qual := rc.qualityCurves()
	for _, k := range []string{"A", "T", "C", "G", "mean"} {
		if _, ok := qual[k]; !ok {
			t.Fatalf("quality_curves missing key %q", k)
		}
	}
	// mean at cycle0 = (40+38)/2 = 39.
	if qual["mean"][0] != 39.0 {
		t.Fatalf("mean[0] = %v, want 39", qual["mean"][0])
	}
	// 'C' never occurs at cycle 0 -> qual curve uses the mean there.
	if qual["C"][0] != qual["mean"][0] {
		t.Fatalf("absent-base qual substitution failed: C[0]=%v mean[0]=%v", qual["C"][0], qual["mean"][0])
	}
	// 'A' qual at cycle0 = (40+38)/2 = 39 (both reads have A there).
	if qual["A"][0] != 39.0 {
		t.Fatalf("A qual[0] = %v, want 39", qual["A"][0])
	}
}

// TestUnitKmer3Kmer2 checks the kmer index->string mappings (bases A,T,C,G).
func TestUnitKmer3Kmer2(t *testing.T) {
	if got := kmer3(0); got != "AAA" {
		t.Fatalf("kmer3(0) = %q, want AAA", got)
	}
	// val 0x3F = 0b111111 -> (3,3,3) -> GGG.
	if got := kmer3(0x3F); got != "GGG" {
		t.Fatalf("kmer3(0x3F) = %q, want GGG", got)
	}
	if got := kmer2(0); got != "AA" {
		t.Fatalf("kmer2(0) = %q, want AA", got)
	}
	if got := kmer2(0x0F); got != "GG" {
		t.Fatalf("kmer2(0x0F) = %q, want GG", got)
	}
	// Decode order: bits map most-significant first.
	if got := kmer3(0x01); got != "AAT" {
		t.Fatalf("kmer3(1) = %q, want AAT", got)
	}
}

// TestUnitKmerCounts checks the 5-mer histogram has 1024 entries, that N
// resets the window, and that a known 5-mer is counted.
func TestUnitKmerCounts(t *testing.T) {
	rc := &readCurves{}
	// "AAAAA" -> one 5-mer AAAAA. "AAAAT" -> AAAAT. With "AAAAAT" we get AAAAA
	// (positions 0-4) then AAAAT (positions 1-5).
	rc.stat(makeRecord("AAAAAT", []int{30, 30, 30, 30, 30, 30}), 33)
	counts := rc.kmerCounts()
	if len(counts) != 1024 {
		t.Fatalf("kmer count entries = %d, want 1024", len(counts))
	}
	if counts["AAAAA"] != 1 {
		t.Fatalf("AAAAA count = %d, want 1", counts["AAAAA"])
	}
	if counts["AAAAT"] != 1 {
		t.Fatalf("AAAAT count = %d, want 1", counts["AAAAT"])
	}
	// N resets the window: "AANAAAA" yields only AAAAA from the tail run.
	rc2 := &readCurves{}
	rc2.stat(makeRecord("AANAAAAA", []int{30, 30, 30, 30, 30, 30, 30, 30}), 33)
	c2 := rc2.kmerCounts()
	// Positions after the N: A A A A A at indices 3..7 -> one AAAAA at i=7.
	if c2["AAAAA"] != 1 {
		t.Fatalf("post-N AAAAA count = %d, want 1", c2["AAAAA"])
	}
}

// TestUnitBuildInsertSize checks the insert_size block construction: peak is
// the most-populated bin below max, unknown is the last bucket, histogram
// excludes the unknown bucket.
func TestUnitBuildInsertSize(t *testing.T) {
	// max=4 -> hist length 5; buckets 0..3 are sizes, bucket 4 is unknown.
	s := &ProcessStats{InsertHist: []int64{0, 2, 5, 1, 3}}
	is := buildInsertSize(s)
	if is == nil {
		t.Fatal("buildInsertSize returned nil")
	}
	if is.Peak != 2 {
		t.Fatalf("peak = %d, want 2 (bin with 5)", is.Peak)
	}
	if is.Unknown != 3 {
		t.Fatalf("unknown = %d, want 3", is.Unknown)
	}
	if !reflect.DeepEqual(is.Histogram, []int64{0, 2, 5, 1}) {
		t.Fatalf("histogram = %v, want [0 2 5 1]", is.Histogram)
	}
	// No histogram -> nil (SE).
	if buildInsertSize(&ProcessStats{}) != nil {
		t.Fatal("buildInsertSize should be nil when InsertHist is empty")
	}
}

// TestUnitSequencingInfo checks the deterministic summary.sequencing string
// for SE and PE from the before-stream cycle counts.
func TestUnitSequencingInfo(t *testing.T) {
	s := &ProcessStats{}
	s.curvesBefore[0] = &readCurves{}
	s.curvesBefore[0].grow(100)
	if got := sequencingInfo(s, false); got != "single end (100 cycles)" {
		t.Fatalf("SE sequencing = %q", got)
	}
	s.curvesBefore[1] = &readCurves{}
	s.curvesBefore[1].grow(75)
	if got := sequencingInfo(s, true); got != "paired end (100 cycles + 75 cycles)" {
		t.Fatalf("PE sequencing = %q", got)
	}
}
