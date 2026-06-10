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
	body, err := readBodyAfterHeader(dataPath, meta)
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

	bw := bgzip.NewWriter(out)
	if _, err := bw.Write(header); err != nil {
		return err
	}
	if _, err := bw.Write(body); err != nil {
		return err
	}
	return bw.Close()
}

// readBodyAfterHeader decodes the bgzipped file at path and returns the bytes
// of every line after the leading meta-character header. If the file does not
// begin with a header line, the entire decoded content is returned unchanged.
func readBodyAfterHeader(path string, meta byte) ([]byte, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	br, err := bgzip.NewReader(f)
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
