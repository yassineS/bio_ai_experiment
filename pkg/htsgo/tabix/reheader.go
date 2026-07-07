package tabix

import (
	"bytes"
	"io"
	"os"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// Reheader rewrites a bgzipped tab-delimited file, replacing its leading
// header — the contiguous run of lines at the top whose first byte equals the
// configured meta character — with the bytes from headerPath, and emits the
// result as a fresh bgzipped stream to out.
//
// The data body (every line after the original header) is preserved verbatim.
// A trailing newline is appended to the replacement header if it does not end
// in one, so the first data line is never merged onto the last header line.
// This mirrors htslib tabix's `reheader` behavior for text formats.
func Reheader(dataPath, headerPath string, meta byte, out io.Writer) error {
	return ReheaderThreaded(dataPath, headerPath, meta, out, 0)
}

// ReheaderThreaded is Reheader with the BGZF (de)compression wired to a worker
// count (upstream tabix -@/--threads, which shares one pool across the input
// and output bgzf streams in reheader_file). When threads > 1 the input body is
// inflated across the pool and the output is deflated across it; at
// DefaultCompression the framed output is byte-identical to the serial writer,
// so the reheadered stream is stable for any thread count. threads < 2 keeps
// the single-threaded reader/writer.
func ReheaderThreaded(dataPath, headerPath string, meta byte, out io.Writer, threads int) error {
	body, err := readBodyAfterHeader(dataPath, meta, threads)
	if err != nil {
		return err
	}
	header, err := os.ReadFile(headerPath)
	if err != nil {
		return err
	}
	if len(header) > 0 && header[len(header)-1] != '\n' {
		header = append(header, '\n')
	}

	var bw io.WriteCloser
	if threads > 1 {
		mw, werr := bgzip.NewMultiWriter(out, bgzip.DefaultCompression, threads)
		if werr != nil {
			return werr
		}
		bw = mw
	} else {
		bw = bgzip.NewWriter(out)
	}
	if _, err := bw.Write(header); err != nil {
		bw.Close()
		return err
	}
	if _, err := bw.Write(body); err != nil {
		bw.Close()
		return err
	}
	return bw.Close()
}

// readBodyAfterHeader decodes the bgzipped file at path and returns the bytes
// of every line after the leading meta-character header. If the file does not
// begin with a header line, the entire decoded content is returned unchanged.
// threads >= 2 inflates the BGZF blocks across a worker pool (byte-identical
// output for any count); < 2 stays single-threaded.
func readBodyAfterHeader(path string, meta byte, threads int) ([]byte, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	br, err := bgzip.NewMultiReader(f, threads)
	if err != nil {
		return nil, err
	}
	defer br.Close()
	data, err := io.ReadAll(br)
	if err != nil {
		return nil, err
	}
	return data[headerLen(data, meta):], nil
}

// headerLen returns the byte length of the leading header in data: the prefix
// consisting of consecutive lines whose first byte equals meta. The returned
// length includes the terminating newline of the last header line. When data
// does not start with a header line the result is 0.
func headerLen(data []byte, meta byte) int {
	pos := 0
	for pos < len(data) && data[pos] == meta {
		nl := bytes.IndexByte(data[pos:], '\n')
		if nl < 0 {
			// Header line runs to EOF with no newline: the whole file is
			// header.
			return len(data)
		}
		pos += nl + 1
	}
	return pos
}
