// seqreg.go: a faithful port of upstream seqtk's region-list parsing and
// masking used by `seqtk seq -M` (the stk_reg_read / stk_mask functions
// in reference_code/seqtk/seqtk.c). The file may be a BED file
// (name<TAB>beg<TAB>end) or a plain name list (one sequence name per
// line); when only a name is given the whole sequence is masked.

package seqtk

import (
	"bufio"
	"io"
	"math"
)

// regInterval is one [beg, end) interval (0-based half-open) belonging to
// a sequence name.
type regInterval struct {
	beg int64
	end int64
}

// regHash maps a sequence name to its list of mask intervals, mirroring
// upstream's khash_t(reg).
type regHash map[string][]regInterval

func isDigitByte(b byte) bool { return b >= '0' && b <= '9' }

// atol parses a leading signed integer the way C's atol does: it consumes
// an optional sign and the leading run of digits and ignores the rest.
func atol(s []byte) int64 {
	i := 0
	neg := false
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		neg = s[i] == '-'
		i++
	}
	var v int64
	for i < len(s) && isDigitByte(s[i]) {
		v = v*10 + int64(s[i]-'0')
		i++
	}
	if neg {
		return -v
	}
	return v
}

// LoadMaskFile parses a BED / name-list stream from r into the option's
// internal region table, so that Seq applies -M masking. Call it once
// before Seq when MaskFile is set.
func (o *SeqOptions) LoadMaskFile(r io.Reader) error {
	h, err := readRegions(r)
	if err != nil {
		return err
	}
	o.maskRegions = h
	return nil
}

// readRegions parses a BED / name-list stream into a regHash, mirroring
// stk_reg_read (seqtk.c:115). The tokenisation matches kseq's
// ks_getuntil with KS_SEP_SPACE: fields are separated by runs of
// whitespace, and the first field on each line is the name.
func readRegions(r io.Reader) (regHash, error) {
	h := make(regHash)
	br := bufio.NewReader(r)

	// Tokeniser state: read one whitespace-delimited token at a time,
	// remembering whether the delimiter that ended it was a newline.
	readToken := func() (tok []byte, dret byte, eof bool) {
		// Skip leading separators except we must report the delimiter; in
		// kseq, ks_getuntil with KS_SEP_SPACE skips leading whitespace.
		for {
			c, err := br.ReadByte()
			if err != nil {
				return nil, 0, true
			}
			if !isSpaceByte(c) {
				tok = append(tok, c)
				break
			}
		}
		for {
			c, err := br.ReadByte()
			if err != nil {
				return tok, 0, false
			}
			if isSpaceByte(c) {
				return tok, c, false
			}
			tok = append(tok, c)
		}
	}

	for {
		name, dret, eof := readToken()
		if eof && len(name) == 0 {
			break
		}
		var beg, end int64 = -1, -1
		if dret != '\n' {
			tok2, dret2, eof2 := readToken()
			if !eof2 && len(tok2) > 0 && isDigitByte(tok2[0]) {
				beg = atol(tok2)
				if dret2 != '\n' {
					tok3, dret3, _ := readToken()
					if len(tok3) > 0 && isDigitByte(tok3[0]) {
						end = atol(tok3)
						if end < 0 {
							end = -1
						}
						dret = dret3
					} else {
						dret = dret3
					}
				} else {
					dret = dret2
				}
			} else {
				dret = dret2
			}
		}
		// Skip the rest of the line.
		if dret != '\n' {
			for {
				c, err := br.ReadByte()
				if err != nil || c == '\n' {
					break
				}
			}
		}
		if end < 0 && beg > 0 {
			end = beg
			beg = beg - 1
		}
		if beg < 0 {
			beg = 0
			end = math.MaxInt64
		}
		key := string(name)
		h[key] = append(h[key], regInterval{beg: beg, end: end})
		if eof {
			break
		}
	}
	return h, nil
}

// maskRegion masks rec.seq according to h, mirroring stk_mask
// (seqtk.c:1337). When isComplement is true the regions are preserved and
// everything else is masked; otherwise the listed regions are masked.
// maskChr == 0 means lowercase the masked bases.
func maskRegion(rec *kseqRecord, h regHash, isComplement bool, maskChr byte) {
	intervals, found := h[string(rec.name)]
	seqLen := int64(len(rec.seq))
	if !found {
		if isComplement {
			if maskChr != 0 {
				for j := int64(0); j < seqLen; j++ {
					rec.seq[j] = maskChr
				}
			} else {
				for j := int64(0); j < seqLen; j++ {
					rec.seq[j] = toLowerByte(rec.seq[j])
				}
			}
		}
		return
	}
	if !isComplement {
		for _, iv := range intervals {
			beg, end := iv.beg, iv.end
			if beg >= seqLen {
				continue
			}
			if end > seqLen {
				end = seqLen
			}
			if maskChr == 0 {
				for j := beg; j < end; j++ {
					rec.seq[j] = toLowerByte(rec.seq[j])
				}
			} else {
				for j := beg; j < end; j++ {
					rec.seq[j] = maskChr
				}
			}
		}
		return
	}
	mask := make([]bool, seqLen)
	for _, iv := range intervals {
		beg, end := iv.beg, iv.end
		if end >= seqLen {
			end = seqLen
		}
		for j := beg; j < end && j < seqLen; j++ {
			if j >= 0 {
				mask[j] = true
			}
		}
	}
	if maskChr != 0 {
		for j := int64(0); j < seqLen; j++ {
			if !mask[j] {
				rec.seq[j] = maskChr
			}
		}
	} else {
		for j := int64(0); j < seqLen; j++ {
			if !mask[j] {
				rec.seq[j] = toLowerByte(rec.seq[j])
			}
		}
	}
}
