package bedunionbedg

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// ReadChromSizes parses a chrom-sizes file (one `chrom\tsize` per line), used
// with -empty to know each chromosome's full extent. It also accepts samtools
// .fai files (first two whitespace-separated columns). Blank lines and comments
// (lines starting with '#') are skipped.
func ReadChromSizes(r io.Reader) (map[string]int64, error) {
	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 64*1024), 16*1024*1024)
	out := make(map[string]int64)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("chrom-sizes line %q must have at least 2 fields", line)
		}
		size, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("invalid size %q for chromosome %q: %v", fields[1], fields[0], err)
		}
		out[fields[0]] = size
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}
