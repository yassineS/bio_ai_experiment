// Output splitting (upstream fastp -s / -S / -d).
//
// Splitting routes the post-filter output stream across multiple numbered
// files. The naming scheme mirrors upstream's ThreadConfig::initWriterForSplit
// (threadconfig.cpp:106-125): each split file is named
// "<dir>/<NUM>.<base>" where NUM is the 1-based split index zero-padded to
// SplitPrefixDigits (e.g. 0001.out.fq, 0002.out.fq).
//
// The file-boundary distribution reproduces upstream's *multi-threaded* pack
// assignment, not just the single-thread case. Upstream reads the input in
// fixed-size packs of PACK_SIZE (256) reads and hands pack i to worker thread
// i % thread (seprocessor.cpp / peprocessor.cpp readerTask). Each thread t owns
// a disjoint set of split files: it starts at split index t and, on rollover,
// advances by +thread (ThreadConfig: mWorkingSplit starts at threadId and
// markProcessed does mWorkingSplit += thread). A thread rolls its current file
// at pack boundaries once its accumulated count reaches split.size:
//
//   - byFileNumber (-s N): the accumulated count is the *input* read count of
//     the packs the thread consumed (markProcessed(pack->count)), and rollover
//     is capped so a thread never advances past split.number-1. Threads that
//     finish early emit empty trailing files (writeEmptyFilesForSplitting).
//   - byFileLines (-S L): the accumulated count is the number of reads that
//     *passed* filtering (markProcessed(readPassed)), with no file-count cap.
//
// Because the reader assigns packs to threads by a deterministic counter (not
// by which worker happens to be free) and each thread consumes its packs
// strictly in FIFO order, the per-file contents are fully deterministic for a
// fixed thread count -- there is no thread-race dependence. We therefore match
// upstream byte-for-byte for any -w, not just -w 1.

package fastp

import (
	"fmt"
	"path/filepath"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// SplitConfig holds the resolved output-splitting parameters. Size is the
// number of records per split file (computed from --split or
// --split_by_lines); Digits is the zero-pad width (--split_prefix_digits);
// Number, when > 0, caps the number of split files (--split mode); Threads is
// the worker-thread count (-w/--thread) that governs the pack-to-thread
// distribution. ByFileLines selects the --split_by_lines rollover rule.
type SplitConfig struct {
	Size        int
	Digits      int
	Number      int // 0 means split-by-lines (no file-count cap).
	Threads     int
	ByFileLines bool
}

// resolveSplitConfig computes the SplitConfig from the user options and the
// known total record count. For --split N (SplitNumber > 0) the per-file
// size is totalRecords/N (at least 1), matching upstream main.cpp:493-497.
// For --split_by_lines L (SplitByLines > 0) the per-file size is L/4
// records (upstream main.cpp:373), with no file-count cap. Threads mirrors the
// upstream worker count, clamped to split.number in byFileNumber mode exactly
// as options.cpp does (thread = min(thread, split.number)).
func resolveSplitConfig(opts ProcessOptions, totalRecords int) SplitConfig {
	digits := opts.SplitPrefixDigits
	threads := opts.Threads
	if threads < 1 {
		threads = 1
	}
	if opts.SplitByLines > 0 {
		size := opts.SplitByLines / 4
		if size < 1 {
			size = 1
		}
		return SplitConfig{Size: size, Digits: digits, Number: 0, Threads: threads, ByFileLines: true}
	}
	// --split N by file number.
	n := opts.SplitNumber
	if n < 1 {
		n = 1
	}
	// Upstream caps the worker count at the number of files (options.cpp:
	// "thread number cannot be more than the number of file to split").
	if threads > n {
		threads = n
	}
	size := totalRecords / n
	if size <= 0 {
		size = 1
	}
	return SplitConfig{Size: size, Digits: digits, Number: n, Threads: threads, ByFileLines: false}
}

// splitPackSize is upstream fastp's PACK_SIZE (common.h:34). The reader emits
// packs of this many input reads and assigns pack i to thread i % thread, and
// each thread's split rollover is only checked once per pack. The per-file
// boundaries are therefore quantized to (and interleaved by) this value.
const splitPackSize = 256

// splitWriter routes the surviving FASTQ records to a sequence of numbered
// output files, reproducing upstream's multi-threaded pack/thread split
// distribution. Callers must announce the current input read index with
// SetInputPos before processing each input read (or pair); the writer uses it
// to determine which worker thread owns the read and thus which split file the
// record belongs to. Records are buffered (tagged with their input position)
// and flushed on Close so the rollover can be resolved exactly as upstream
// does, independent of how many reads pass filtering.
type splitWriter struct {
	basePath string
	cfg      SplitConfig
	encoding fastq.QualityEncoding

	// pos is the input read index announced by the most recent SetInputPos.
	pos int
	// entries collects, in call order, the records that survived filtering
	// (one per Write call), each tagged with the input position it came from.
	entries []splitEntry
	// maxPos tracks the highest input position announced, so Close knows the
	// total input read count (needed for the byFileNumber rollover cap and the
	// trailing empty files).
	maxPos int
}

// splitEntry is one surviving record tagged with the input read index it came
// from, buffered until Close resolves its destination file.
type splitEntry struct {
	pos    int
	record *fastq.Record
}

// newSplitWriter creates a splitWriter for the given base output path and
// split configuration.
func newSplitWriter(basePath string, cfg SplitConfig, encoding fastq.QualityEncoding) *splitWriter {
	if cfg.Threads < 1 {
		cfg.Threads = 1
	}
	return &splitWriter{
		basePath: basePath,
		cfg:      cfg,
		encoding: encoding,
		pos:      -1,
		maxPos:   -1,
	}
}

// SetInputPos records the input read index of the read about to be processed.
// It must be called before the corresponding processOneSE/processPairOnce so a
// subsequent Write is attributed to the right thread and split file.
func (sw *splitWriter) SetInputPos(pos int) {
	sw.pos = pos
	if pos > sw.maxPos {
		sw.maxPos = pos
	}
}

// splitFileName returns the split file name for the given 0-based index,
// matching upstream's "<dir>/<NUM>.<base>" with 1-based, zero-padded NUM.
func splitFileName(basePath string, index, digits int) string {
	num := fmt.Sprintf("%d", index+1)
	if digits > 0 {
		for len(num) < digits {
			num = "0" + num
		}
	}
	dir := filepath.Dir(basePath)
	base := filepath.Base(basePath)
	name := num + "." + base
	if dir == "." && !strings.HasPrefix(basePath, "."+string(filepath.Separator)) {
		return name
	}
	return filepath.Join(dir, name)
}

// Write buffers a surviving record, tagging it with the current input
// position. The actual file assignment is deferred to Close. It satisfies the
// recordWriter interface so a splitWriter is a drop-in for *fastq.Writer.
func (sw *splitWriter) Write(record *fastq.Record) error {
	sw.entries = append(sw.entries, splitEntry{pos: sw.pos, record: record})
	return nil
}

// Flush is a no-op: records are buffered until Close. It satisfies the
// recordWriter interface.
func (sw *splitWriter) Flush() error {
	return nil
}

// recordWriter is the minimal sink used by the processing loops: a
// *fastq.Writer (single output file) or a *splitWriter (numbered files).
type recordWriter interface {
	Write(record *fastq.Record) error
	Flush() error
}

// Close resolves each buffered record's destination split file using the
// upstream multi-threaded rollover simulation, writes every split file that
// upstream would open (including the empty ones), and returns the first error.
func (sw *splitWriter) Close() error {
	fileOf, opened := sw.assignFiles()

	// Group records by destination file, preserving call (input) order. A file
	// is owned by a single thread that consumes its packs in increasing index
	// order, so the records of a file come out in increasing input-position
	// order -- exactly the order entries were buffered.
	perFile := make(map[int][]*fastq.Record)
	for i, e := range sw.entries {
		f := fileOf[i]
		perFile[f] = append(perFile[f], e.record)
	}

	for _, f := range opened {
		if err := sw.writeFile(f, perFile[f]); err != nil {
			return err
		}
	}
	return nil
}

// assignFiles runs the upstream pack/thread rollover simulation. It returns,
// for each buffered entry (by buffer index), the 0-based split file index it
// belongs to, and the sorted list of every file index upstream would open
// (initWriterForSplit) -- including the empty files created by a rollover on
// the last pack and by writeEmptyFilesForSplitting in byFileNumber mode.
func (sw *splitWriter) assignFiles() ([]int, []int) {
	total := sw.maxPos + 1
	if total < 0 {
		total = 0
	}
	threads := sw.cfg.Threads
	if threads < 1 {
		threads = 1
	}

	// Per-thread state: the split file currently open and the accumulated
	// count used to decide rollover. Each thread t starts on file index t,
	// which upstream opens via initConfig before any pack is processed.
	workingSplit := make([]int, threads)
	cur := make([]int, threads)
	openedSet := make(map[int]bool)
	for t := 0; t < threads; t++ {
		workingSplit[t] = t
		// initConfig -> initWriterForSplit opens each thread's base file, but
		// only when that base index is a valid file. In byFileNumber mode a
		// thread whose base index is >= split.number cannot exist (threads are
		// clamped to split.number), so t is always a real file here.
		openedSet[t] = true
	}

	// fileForPos[p] is the destination file for any record at input position p.
	fileForPos := make([]int, total)

	// passedAt[p] counts the surviving records buffered at input position p
	// (used only for the byFileLines passed-count rollover). For SE each
	// surviving read contributes 1; for PE each surviving mate contributes 1
	// to its own writer, mirroring upstream's per-pack readPassed count (which
	// is identical for the R1 and R2 writers since a pair passes together).
	var passedAt []int
	if sw.cfg.ByFileLines {
		passedAt = make([]int, total)
		for _, e := range sw.entries {
			if e.pos >= 0 && e.pos < total {
				passedAt[e.pos]++
			}
		}
	}

	packCount := 0
	if total > 0 {
		packCount = (total + splitPackSize - 1) / splitPackSize
	}

	for pi := 0; pi < packCount; pi++ {
		t := pi % threads
		start := pi * splitPackSize
		end := start + splitPackSize
		if end > total {
			end = total
		}
		packPassed := 0
		for p := start; p < end; p++ {
			fileForPos[p] = workingSplit[t]
			if sw.cfg.ByFileLines {
				packPassed += passedAt[p]
			}
		}
		// markProcessed: advance the thread's accumulated count, then roll the
		// file at this pack boundary if it has reached split.size. A rollover
		// opens (initWriterForSplit) the next file, even on the final pack --
		// which is why upstream can emit a trailing empty file.
		if sw.cfg.ByFileLines {
			cur[t] += packPassed
		} else {
			cur[t] += end - start // input read count of this pack
		}
		if cur[t] >= sw.cfg.Size {
			if sw.cfg.ByFileLines || workingSplit[t]+threads < sw.cfg.Number {
				workingSplit[t] += threads
				cur[t] = 0
				openedSet[workingSplit[t]] = true
			}
		}
	}

	// byFileNumber: writeEmptyFilesForSplitting opens each remaining file in a
	// thread's residue class up to split.number-1, so every index 0..N-1 is
	// materialized even when the input is short.
	if !sw.cfg.ByFileLines && sw.cfg.Number > 0 {
		for t := 0; t < threads; t++ {
			for ws := t; ws < sw.cfg.Number; ws += threads {
				openedSet[ws] = true
			}
		}
	}

	out := make([]int, len(sw.entries))
	for i, e := range sw.entries {
		if e.pos >= 0 && e.pos < total {
			out[i] = fileForPos[e.pos]
		}
	}

	opened := make([]int, 0, len(openedSet))
	for f := range openedSet {
		opened = append(opened, f)
	}
	sort.Ints(opened)
	return out, opened
}

// writeFile writes the records destined for split file index f.
func (sw *splitWriter) writeFile(f int, records []*fastq.Record) error {
	name := splitFileName(sw.basePath, f, sw.cfg.Digits)
	wc, err := iohelper.OpenWriter(name)
	if err != nil {
		return fmt.Errorf("open split file %q: %w", name, err)
	}
	w := fastq.NewWriter(wc, sw.encoding)
	for _, rec := range records {
		if err := w.Write(rec); err != nil {
			wc.Close()
			return err
		}
	}
	if err := w.Flush(); err != nil {
		wc.Close()
		return err
	}
	return wc.Close()
}
