package main

import (
	"bufio"
	"io"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// windowSize mirrors htslib bgzip.c's WINDOW_SIZE (== BGZF_BLOCK_SIZE): the
// chunk size in which the compressor reads input and the upper bound on a
// single BGZF block's uncompressed payload.
const windowSize = bgzip.MaxBlockSize

// detectPeekLen is the number of leading bytes hts_detect_format peeks at when
// classifying an uncompressed stream (hts.c uses a 1024-byte window).
const detectPeekLen = 1024

// seqNT16Table maps a byte to its IUPAC 4-bit nucleotide code, or 15 for
// characters that are not nucleotide letters. It mirrors htslib's
// seq_nt16_table and is used by isFastaq to decide whether the second line of
// a '>'/'@' record looks like sequence data.
var seqNT16Table = buildSeqNT16Table()

func buildSeqNT16Table() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = 15
	}
	// "=ACMGRSVTWYHKDBN" indexed 0..15 in htslib's seq_nt16_str.
	const code = "=ACMGRSVTWYHKDBN"
	for i := 0; i < len(code); i++ {
		c := code[i]
		t[c] = byte(i)
		// Lower-case bases map identically.
		if c >= 'A' && c <= 'Z' {
			t[c+('a'-'A')] = byte(i)
		}
	}
	return t
}

// compressStream reads uncompressed bytes from src and writes a BGZF stream to
// bw, replicating htslib bgzip's framing exactly. When the input is detected as
// a textual bioinformatics format (and forceBinary is false) it applies the
// text-mode block-flush heuristic from bgzip.c: the header occupies its own
// block(s) and each WINDOW_SIZE chunk of records is flushed at the last newline
// boundary. Binary or undetected input streams straight through, letting the
// writer pack natural 65280-byte blocks. compressStream does not close bw.
func compressStream(bw *bgzip.Writer, src io.Reader, forceBinary bool) error {
	br := bufio.NewReaderSize(src, windowSize)

	textual := false
	if !forceBinary {
		peek, _ := br.Peek(detectPeekLen)
		textual = detectTextual(peek)
	}

	if !textual {
		return streamBinary(bw, br)
	}
	return streamText(bw, br)
}

// streamBinary mirrors the binary branch of bgzip.c: read WINDOW_SIZE chunks
// and write them through with no forced flushing.
func streamBinary(bw *bgzip.Writer, br *bufio.Reader) error {
	buf := make([]byte, windowSize)
	for {
		c, err := io.ReadFull(br, buf)
		if c > 0 {
			if _, werr := bw.Write(buf[:c]); werr != nil {
				return werr
			}
		}
		if err == io.EOF || err == io.ErrUnexpectedEOF {
			return nil
		}
		if err != nil {
			return err
		}
	}
}

// streamText mirrors the text-mode loop of bgzip.c (lines ~435-498): it keeps a
// WINDOW_SIZE buffer, accumulates the header into its own block, then flushes a
// block boundary after the last newline of every chunk so records align to line
// boundaries. Leftover bytes after the flush point carry to the next iteration.
func streamText(bw *bgzip.Writer, br *bufio.Reader) error {
	buf := make([]byte, windowSize)
	inHeader := true
	n := 0
	longLine := false

	for {
		// Read up to WINDOW_SIZE-n bytes into buf+n. Like htslib's hread on
		// a regular file, fill the request completely unless EOF intervenes,
		// so chunk boundaries (and therefore block boundaries) match upstream.
		c, rerr := io.ReadFull(br, buf[n:])
		if rerr == io.ErrUnexpectedEOF {
			rerr = io.EOF
		}
		if c == 0 {
			if rerr == io.EOF {
				break
			}
			if rerr != nil {
				return rerr
			}
			break
		}

		c2 := c + n
		flush := false

		if inHeader && (longLine || buf[0] == '@' || buf[0] == '#') {
			// Scan forward to the last header line.
			lastStart := 0
			n = 0
			for n < c2 {
				if buf[n] != '\n' {
					n++
					continue
				}
				n++
				lastStart = n
				if n < c2 && !(buf[n] == '@' || buf[n] == '#') {
					inHeader = false
					break
				}
			}
			if lastStart == 0 {
				n = c2
				longLine = true
			} else {
				n = lastStart
				flush = true
				longLine = false
			}
		} else {
			// Scan backwards to the last newline.
			n += c
			for {
				n--
				if n < 0 || buf[n] == '\n' {
					break
				}
			}
			if n >= 0 {
				flush = true
				n++
			} else {
				n = c2
			}
		}

		if _, err := bw.Write(buf[:n]); err != nil {
			return err
		}
		if flush {
			if err := bw.Flush(); err != nil {
				return err
			}
		}

		copy(buf, buf[n:c2])
		n = c2 - n

		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return rerr
		}
	}

	// Trailing data.
	if n > 0 {
		if _, err := bw.Write(buf[:n]); err != nil {
			return err
		}
	}
	return nil
}
