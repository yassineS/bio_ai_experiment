package main

import "bytes"

// detectTextual reports whether the peeked input bytes look like one of the
// uncompressed textual bioinformatics formats that htslib bgzip treats with its
// text-mode block-flush heuristic: plain text, SAM, VCF, BED, FASTA, FASTQ, and
// the FAI/FQI index variants. It mirrors the relevant branches of htslib's
// hts_detect_format (hts.c): any binary magic (gzip/BGZF/BAM/BCF/CRAM/...) or
// non-text content yields false, in which case bgzip streams the data straight
// through without aligning blocks to line boundaries.
//
// The peek window is whatever Peek returned (htslib peeks up to 1024 bytes);
// detection works on that prefix exactly as upstream does.
func detectTextual(s []byte) bool {
	if len(s) == 0 {
		// empty_format: upstream leaves format unset (not textual).
		return false
	}

	// Compressed streams are never textual for bgzip's purposes: gzip/BGZF,
	// bzip2, xz, and zstd all carry distinctive leading magic. (htslib would
	// decompress-and-peek, but bgzip only sets textual=1 when compression ==
	// no_compression, so any compressed input takes the binary path.)
	if len(s) >= 2 && s[0] == 0x1f && s[1] == 0x8b {
		return false // gzip / BGZF
	}
	if len(s) >= 3 && bytes.Equal(s[:3], []byte("BZh")) {
		return false // bzip2
	}
	if len(s) >= 6 && bytes.Equal(s[:6], []byte{0xfd, '7', 'z', 'X', 'Z', 0x00}) {
		return false // xz
	}
	if len(s) >= 4 && bytes.Equal(s[:4], []byte{0x28, 0xb5, 0x2f, 0xfd}) {
		return false // zstd
	}

	// Binary magic for the htslib container formats (s[3] <= 4 guards the
	// block in upstream; the individual magics are unambiguous regardless).
	if len(s) >= 6 && bytes.Equal(s[:4], []byte("CRAM")) && s[4] >= 1 && s[4] <= 7 && s[5] <= 7 {
		return false
	}
	if len(s) >= 4 {
		switch {
		case bytes.Equal(s[:4], []byte("BAM\x01")),
			bytes.Equal(s[:4], []byte("BAI\x01")),
			bytes.Equal(s[:4], []byte("BCF\x04")),
			bytes.Equal(s[:4], []byte("BCF\x02")),
			bytes.Equal(s[:4], []byte("CSI\x01")),
			bytes.Equal(s[:4], []byte("TBI\x01")):
			return false
		}
	}

	// VCF: "##fileformat=VCF".
	if len(s) >= 16 && bytes.Equal(s[:16], []byte("##fileformat=VCF")) {
		return true
	}

	// SAM header line.
	if len(s) >= 4 && s[0] == '@' {
		switch {
		case bytes.Equal(s[:4], []byte("@HD\t")),
			bytes.Equal(s[:4], []byte("@SQ\t")),
			bytes.Equal(s[:4], []byte("@RG\t")),
			bytes.Equal(s[:4], []byte("@PG\t")),
			bytes.Equal(s[:4], []byte("@CO\t")):
			return true
		}
	}

	// Binary/non-text container formats with their own magic.
	if len(s) >= 8 && bytes.Equal(s[:4], []byte("d4\xdd\xdd")) {
		return false // d4
	}
	if cmpNonblank("{\"htsget\":", s) == 0 {
		return false // htsget JSON
	}
	if len(s) > 8 && bytes.Equal(s[:8], []byte("crypt4gh")) {
		return false
	}

	// FASTA / FASTQ.
	if s[0] == '>' && isFastaq(s) {
		return true
	}
	if s[0] == '@' && isFastaq(s) {
		return true
	}

	// Tab-delimited text: SAM body, BED, FAI/FQI. We do not have the filename
	// extension here, so FAI/FQI cannot be distinguished from BED — but they
	// are all textual, and any tab-separated text that is not SAM/BED still
	// falls through to the plain-text check below, which classifies it as
	// textual. So matching SAM or BED column shapes is sufficient.
	if cols, n := parseTabbedText(s); n > 0 {
		switch {
		case colmatchAtLeast(cols, "ZiZiiCZiiZZOOOOOOOOOOOOOOOOOOOO+", 9): // SAM
			return true
		case colmatchAtLeast(cols, "Zii+", 3): // BED
			return true
		}
	}

	// Arbitrary text file: every byte is printable/whitespace.
	return isTextOnly(s)
}

// isTextOnly reports whether every byte is a printable character (>= ' ') or
// one of tab/CR/LF, mirroring htslib's is_text_only.
func isTextOnly(s []byte) bool {
	for _, c := range s {
		if !(c >= ' ' || c == '\t' || c == '\r' || c == '\n') {
			return false
		}
	}
	return true
}

// cmpNonblank compares key against s, skipping whitespace in s, mirroring
// htslib's cmp_nonblank. It returns 0 on a full prefix match.
func cmpNonblank(key string, s []byte) int {
	ki, si := 0, 0
	for ki < len(key) {
		if si >= len(s) {
			return +1
		}
		c := s[si]
		if c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f' {
			si++
			continue
		}
		if c != key[ki] {
			if key[ki] < c {
				return -1
			}
			return +1
		}
		si++
		ki++
	}
	return 0
}

// isFastaq mirrors htslib's is_fastaq: the first line must be textual, and the
// second line (if present) must be base-encoding letters (incl. N, excl. '=').
func isFastaq(s []byte) bool {
	nl := bytes.IndexByte(s, '\n')
	var firstLim []byte
	if nl >= 0 {
		firstLim = s[:nl]
	} else {
		firstLim = s
	}
	if !isTextOnly(firstLim) {
		return false
	}
	if nl < 0 {
		// Very long first line: treat as FASTA/Q.
		return true
	}
	u := s[nl+1:]
	i := 0
	for i < len(u) {
		c := u[i]
		if seqNT16Table[c] == 15 && toUpper(c) != 'N' {
			break
		}
		if c == '=' {
			return false
		}
		i++
	}
	if i == len(u) {
		return true
	}
	return u[i] == '\r' || u[i] == '\n'
}

func toUpper(c byte) byte {
	if c >= 'a' && c <= 'z' {
		return c - ('a' - 'A')
	}
	return c
}

// bamCigarStr is the set of CIGAR operation characters (htslib BAM_CIGAR_STR).
const bamCigarStr = "MIDNSHP=XB"

// parseTabbedText mirrors htslib's parse_tabbed_text: it classifies tab- (or
// newline-) delimited columns by type, returning the type string and the column
// count, or -1 if a non-printable, non-delimiter byte is seen.
func parseTabbedText(s []byte) (string, int) {
	const (
		digit       = 1
		leadingSign = 2
		cigarOp     = 4
		other       = 8
		maxColumns  = 23 // column_len(24) - 1
	)
	var cols []byte
	start := 0
	seen := 0
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch {
		case c >= ' ':
			switch {
			case c >= '0' && c <= '9':
				seen |= digit
			case (c == '+' || c == '-') && i == start:
				seen |= leadingSign
			case isCigarOp(c) && i > start && s[i-1] >= '0' && s[i-1] <= '9':
				seen |= cigarOp
			default:
				seen |= other
			}
		case c == '\t' || c == '\r' || c == '\n':
			length := i - start
			var typ byte
			switch {
			case seen == digit || seen == (leadingSign|digit):
				typ = 'i'
			case seen == (digit | cigarOp):
				typ = 'C'
			case length == 1:
				switch s[start] {
				case '*':
					typ = 'C'
				case '+', '-', '.':
					typ = 's'
				default:
					typ = 'Z'
				}
			case length >= 5 && s[start+2] == ':' && s[start+4] == ':':
				typ = 'O'
			default:
				typ = 'Z'
			}
			cols = append(cols, typ)
			if c != '\t' || len(cols) >= maxColumns {
				return string(cols), len(cols)
			}
			start = i + 1
			seen = 0
		default:
			return "", -1
		}
	}
	return string(cols), len(cols)
}

func isCigarOp(c byte) bool {
	for i := 0; i < len(bamCigarStr); i++ {
		if bamCigarStr[i] == c {
			return true
		}
	}
	return false
}

// colmatchAtLeast reports whether columns matches pattern as a prefix (mirroring
// htslib's colmatch) and the matched length is >= min. '+' in the pattern marks
// a "match rest" position, 'Z' matches any type.
func colmatchAtLeast(columns, pattern string, min int) bool {
	i := 0
	for i < len(columns) {
		if i < len(pattern) && pattern[i] == '+' {
			return i >= min
		}
		if i >= len(pattern) {
			return false
		}
		if !(columns[i] == pattern[i] || pattern[i] == 'Z') {
			return false
		}
		i++
	}
	return i >= min
}
