// Shared helpers for the FASTA-backed native bcftools plugins
// (fill-from-fasta, fixref). They cover reading appended header-lines files and
// resolving a declared INFO tag's Type from the header meta lines.
package bcftools

import (
	"bufio"
	"os"
	"strings"
)

// readHeaderLinesFile reads a file of additional VCF header lines (the
// fill-from-fasta -h/--header-lines option). Blank lines are skipped; every
// other non-empty line is returned verbatim (with any trailing CR stripped),
// to be appended to the output header.
func readHeaderLinesFile(path string) ([]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var lines []string
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 8*1024*1024)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r")
		if line == "" {
			continue
		}
		lines = append(lines, line)
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return lines, nil
}

// infoTagType returns the declared Type (e.g. "Integer", "String", "Float") of
// the INFO tag id as found in the header meta lines, plus whether such an INFO
// header line was found. The lookup is the equivalent of
// bcf_hdr_id2type(hdr,BCF_HL_INFO,id).
func infoTagType(meta []string, id string) (string, bool) {
	for _, line := range meta {
		if !strings.HasPrefix(line, "##INFO=<") {
			continue
		}
		if headerID(line) != id {
			continue
		}
		return headerField(line, "Type"), true
	}
	return "", false
}

// headerField extracts the value of key=... from a structured header line
// (e.g. Type=Integer). It returns "" when the key is absent.
func headerField(line, key string) string {
	needle := key + "="
	i := strings.Index(line, needle)
	if i < 0 {
		return ""
	}
	rest := line[i+len(needle):]
	if len(rest) > 0 && rest[0] == '"' {
		end := strings.IndexByte(rest[1:], '"')
		if end < 0 {
			return rest[1:]
		}
		return rest[1 : 1+end]
	}
	end := strings.IndexAny(rest, ",>")
	if end < 0 {
		return rest
	}
	return rest[:end]
}
