// Output splitting (upstream fastp -s / -S / -d).
//
// Splitting routes the post-filter output stream across multiple numbered
// files. The naming scheme mirrors upstream's ThreadConfig::initWriterForSplit
// (threadconfig.cpp:106-125): each split file is named
// "<dir>/<NUM>.<base>" where NUM is the 1-based split index zero-padded to
// SplitPrefixDigits (e.g. 0001.out.fq, 0002.out.fq). A new file is opened
// whenever the current one reaches SplitSize records, matching
// ThreadConfig::markProcessed (threadconfig.cpp:127-149) in the
// single-thread case.

package fastp

import (
	"fmt"
	"path/filepath"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// SplitConfig holds the resolved output-splitting parameters. Size is the
// number of records per split file (computed from --split or
// --split_by_lines); Digits is the zero-pad width (--split_prefix_digits);
// Number, when > 0, caps the number of split files (--split mode).
type SplitConfig struct {
	Size   int
	Digits int
	Number int // 0 means split-by-lines (no file-count cap).
}

// resolveSplitConfig computes the SplitConfig from the user options and the
// known total record count. For --split N (SplitNumber > 0) the per-file
// size is totalRecords/N (at least 1), matching upstream main.cpp:493-497.
// For --split_by_lines L (SplitByLines > 0) the per-file size is L/4
// records (upstream main.cpp:373), with no file-count cap.
func resolveSplitConfig(opts ProcessOptions, totalRecords int) SplitConfig {
	digits := opts.SplitPrefixDigits
	if opts.SplitByLines > 0 {
		size := opts.SplitByLines / 4
		if size < 1 {
			size = 1
		}
		return SplitConfig{Size: size, Digits: digits, Number: 0}
	}
	// --split N by file number.
	n := opts.SplitNumber
	if n < 1 {
		n = 1
	}
	size := totalRecords / n
	if size <= 0 {
		size = 1
	}
	return SplitConfig{Size: size, Digits: digits, Number: n}
}

// splitPackSize is upstream fastp's PACK_SIZE (common.h:34). Splitting only
// rolls over to a new file at pack boundaries (ThreadConfig::markProcessed is
// invoked once per pack with the pack's read count), so output file sizes are
// quantized to multiples of this value. We replicate that granularity for
// byte-for-byte parity of the per-file boundaries.
const splitPackSize = 256

// splitWriter routes FASTQ records to a sequence of numbered output files,
// rolling over to a new file once the current file has reached Size records,
// checked at pack (splitPackSize) boundaries to match upstream. It owns the
// file handles it opens and closes them as it rotates and on Close.
type splitWriter struct {
	basePath string
	cfg      SplitConfig
	encoding fastq.QualityEncoding

	current      *fastq.Writer
	closer       func() error
	currentReads int
	workingSplit int // 0-based index of the current split file.
}

// newSplitWriter creates a splitWriter for the given base output path and
// split configuration. The first split file is opened lazily on the first
// Write so that a zero-record run still produces 0001.<base> via Close.
func newSplitWriter(basePath string, cfg SplitConfig, encoding fastq.QualityEncoding) *splitWriter {
	return &splitWriter{
		basePath:     basePath,
		cfg:          cfg,
		encoding:     encoding,
		workingSplit: 0,
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

// openCurrent opens the split file for the current workingSplit index.
func (sw *splitWriter) openCurrent() error {
	name := splitFileName(sw.basePath, sw.workingSplit, sw.cfg.Digits)
	wc, err := iohelper.OpenWriter(name)
	if err != nil {
		return fmt.Errorf("open split file %q: %w", name, err)
	}
	sw.current = fastq.NewWriter(wc, sw.encoding)
	sw.closer = wc.Close
	sw.currentReads = 0
	return nil
}

// rotate flushes and closes the current split file and advances to the next.
func (sw *splitWriter) rotate() error {
	if sw.current != nil {
		if err := sw.current.Flush(); err != nil {
			return err
		}
		if sw.closer != nil {
			if err := sw.closer(); err != nil {
				return err
			}
		}
		sw.current = nil
		sw.closer = nil
	}
	sw.workingSplit++
	return nil
}

// Write appends a record to the current split file, rotating first if the
// current file is full. The Size==0 guard treats every record as a single
// file (defensive; callers always set Size >= 1).
func (sw *splitWriter) Write(record *fastq.Record) error {
	if sw.current == nil {
		if err := sw.openCurrent(); err != nil {
			return err
		}
	}
	// Roll over before writing when the current file is full. Upstream checks
	// this only at pack boundaries (markProcessed is per-pack), so we likewise
	// only consider rotating when currentReads is a positive multiple of the
	// pack size. In file-number mode we never roll past the last allowed file
	// (upstream keeps appending to it).
	if sw.cfg.Size > 0 && sw.currentReads > 0 &&
		sw.currentReads%splitPackSize == 0 && sw.currentReads >= sw.cfg.Size {
		atLastFile := sw.cfg.Number > 0 && sw.workingSplit+1 >= sw.cfg.Number
		if !atLastFile {
			if err := sw.rotate(); err != nil {
				return err
			}
			if err := sw.openCurrent(); err != nil {
				return err
			}
		}
	}
	if err := sw.current.Write(record); err != nil {
		return err
	}
	sw.currentReads++
	return nil
}

// Flush flushes the current split file's buffer, if any. It satisfies the
// recordWriter interface so a splitWriter is a drop-in for *fastq.Writer.
func (sw *splitWriter) Flush() error {
	if sw.current != nil {
		return sw.current.Flush()
	}
	return nil
}

// recordWriter is the minimal sink used by the processing loops: a
// *fastq.Writer (single output file) or a *splitWriter (numbered files).
type recordWriter interface {
	Write(record *fastq.Record) error
	Flush() error
}

// Close flushes and closes the final split file.
func (sw *splitWriter) Close() error {
	if sw.current == nil {
		// No records written: still create an empty first split file so the
		// output set is non-empty, matching upstream's behaviour of always
		// opening 0001.<base>.
		if err := sw.openCurrent(); err != nil {
			return err
		}
	}
	if err := sw.current.Flush(); err != nil {
		return err
	}
	if sw.closer != nil {
		return sw.closer()
	}
	return nil
}
