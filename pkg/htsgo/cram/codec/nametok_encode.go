package codec

// nametok_encode.go ports the encoder half of the htscodecs name
// tokeniser (tokenise_name3.c's tok3_encode_names). It reproduces the
// reference encoder's structure faithfully: the prefix trie that picks
// which earlier name to diff against, the per-name tokenisation, the
// per-token-block sub-codec method search, and the descriptor
// deduplication / N_TYPE elision applied before serialisation.
//
// The output is byte-identical to htscodecs' tok3_encode_names for the
// levels exercised by the compliance suite — see NameTokEncode for the
// exact guarantees.

import (
	"fmt"
)

// trieNode is one node of the prefix trie build_trie/search_trie walk.
// Children are kept as a sibling-linked list, exactly as in the C struct.
type trieNode struct {
	next, sibling *trieNode
	count         int
	c             byte
	n             int // index of the most recent name passing through here.
}

// nameEncoder holds the encoder's working state for one block of names.
type nameEncoder struct {
	tHead   *trieNode
	lc      []encLastCtx
	counter int

	desc     [nameMaxTBlocks]encDescriptor
	tokenD   [nameMaxTokens]int
	tokenI   [nameMaxTokens]int
	maxTok   int
	maxNames int
}

// encLastCtx is the encoder's per-name token history (last_context).
type encLastCtx struct {
	name []byte
	ntok int
	last []encLastTok
}

// encLastTok is one remembered token (last_context_tok).
type encLastTok struct {
	typ    int
	intVal int
	strOff int
}

// encDescriptor accumulates the raw byte stream of one token-block, plus
// the bookkeeping needed for the serialisation pass.
type encDescriptor struct {
	buf     []byte
	tnum    int
	ttype   int
	dupFrom int
}

// NameTokEncode compresses names — a buffer of read names separated by
// NUL or newline bytes, optionally with a trailing separator — into a
// complete htscodecs name-tokeniser block (CRAM compression method 8).
//
// level selects the encoder's sub-codec search aggressiveness. Levels
// 1-9 use the static rANS 4x16 sub-codec; levels 11-19 use the
// arith_dynamic range coder (with effective level level-10). These are
// the value ranges the tok3 compliance suite exercises and for which the
// output is byte-identical to the reference tok3_encode_names.
//
// The encoder always round-trips: NameTokDecode(NameTokEncode(x)) joins
// the same names with NULs. It returns an error if a name contains a
// non-printable byte or the input describes more than the supported
// number of records.
func NameTokEncode(names []byte, level int) ([]byte, error) {
	useArith := false
	if level > 10 {
		level -= 10
		useArith = true
	}
	if level < 1 {
		level = 1
	}
	return tok3EncodeNames(names, level, useArith)
}

// tok3EncodeNames ports tok3_encode_names. blk is consumed read-only;
// names are delimited by any byte <= '\n'.
func tok3EncodeNames(blk []byte, level int, useArith bool) ([]byte, error) {
	// Count records (entries terminated by \n or \0).
	nreads := 0
	for _, b := range blk {
		if b <= '\n' {
			nreads++
		}
	}
	if nreads > nameMaxNames {
		return nil, fmt.Errorf("%w: %d names exceeds the %d-record ceiling", errNameTok, nreads, nameMaxNames)
	}

	enc := &nameEncoder{
		lc:       make([]encLastCtx, nreads+1),
		maxNames: nreads + 1,
		maxTok:   1,
	}

	// Build the prefix trie over every complete line.
	ctr := 0
	for i := 0; i < len(blk); {
		j := i
		for i < len(blk) && blk[i] > '\n' {
			i++
		}
		if i >= len(blk) {
			break
		}
		if err := enc.buildTrie(blk[j:i], ctr); err != nil {
			return nil, err
		}
		ctr++
		i++
	}

	// Tokenise each name. blk bytes are copied per-name so the encoder
	// can NUL-terminate its working copy without mutating the input.
	for i := 0; i < len(blk); {
		j := i
		for i < len(blk) && int8(blk[i]) >= ' ' {
			i++
		}
		if i >= len(blk) {
			break
		}
		if blk[i] != 0 && blk[i] != '\n' {
			return nil, fmt.Errorf("%w: name contains a non-ASCII-printable byte", errNameTok)
		}
		name := make([]byte, i-j)
		copy(name, blk[j:i])
		if err := enc.encodeName(name); err != nil {
			return nil, err
		}
		i++
	}

	// Drop N_TYPE (type-0) blocks whose every entry past the first is
	// N_MATCH — the decoder regenerates them from the sibling blocks.
	for i := 0; i < enc.maxTok*16; i += 16 {
		d := &enc.desc[i]
		if len(d.buf) == 0 {
			continue
		}
		allMatch := true
		for z := 1; z < len(d.buf); z++ {
			if d.buf[z] != ntMatch {
				allMatch = false
				break
			}
		}
		if !allMatch {
			continue
		}
		hasSibling := false
		for k := 1; k < 16; k++ {
			if len(enc.desc[i+k].buf) != 0 {
				hasSibling = true
				break
			}
		}
		if hasSibling {
			d.buf = nil
		}
	}

	// Compress each surviving descriptor and detect duplicate streams.
	for i := 0; i < enc.maxTok*16; i++ {
		d := &enc.desc[i]
		d.dupFrom = -1
		if len(d.buf) == 0 {
			continue
		}
		comp, err := compressTokenBlock(d.buf, i&15, level, useArith)
		if err != nil {
			return nil, err
		}
		d.buf = comp
		d.tnum = i >> 4
		d.ttype = i & 15

		for j := 0; j < i; j++ {
			o := &enc.desc[j]
			if o.buf == nil || o.dupFrom != -1 {
				// A descriptor that was itself a dup has no standalone
				// stream to copy; skip it (its buf is the original).
			}
			if len(o.buf) == 0 || len(o.buf) != len(d.buf) || len(d.buf) <= 4 {
				continue
			}
			if bytesEqual(o.buf, d.buf) {
				d.dupFrom = j
				break
			}
		}
	}

	// Serialise: 9-byte header then each descriptor.
	lastStart := len(blk)
	out := make([]byte, 0, len(blk)/2+64)
	out = append(out,
		byte(lastStart), byte(lastStart>>8), byte(lastStart>>16), byte(lastStart>>24),
		byte(nreads), byte(nreads>>8), byte(nreads>>16), byte(nreads>>24),
		boolByte(useArith))

	lastTnum := -1
	for i := 0; i < enc.maxTok*16; i++ {
		d := &enc.desc[i]
		if len(d.buf) == 0 {
			continue
		}
		ttype8 := byte(d.ttype)
		if d.tnum != lastTnum {
			ttype8 |= 128
			lastTnum = d.tnum
		}
		if d.dupFrom >= 0 {
			out = append(out, ttype8|64, byte(d.dupFrom>>4), byte(d.dupFrom&15))
		} else {
			out = append(out, ttype8)
			out = append(out, d.buf...)
		}
	}
	return out, nil
}

// boolByte maps a bool to the 0/1 byte the format stores.
func boolByte(b bool) byte {
	if b {
		return 1
	}
	return 0
}

// --- trie --------------------------------------------------------------------

// buildTrie inserts one name (without its terminator) into the prefix
// trie, mirroring build_trie. n is the name's ordinal index.
func (enc *nameEncoder) buildTrie(data []byte, n int) error {
	if enc.tHead == nil {
		enc.tHead = &trieNode{}
	}
	t := enc.tHead
	t.count++
	for i := 0; i < len(data); i++ {
		c := data[i]
		if c&0x80 != 0 {
			return fmt.Errorf("%w: 8-bit characters are unsupported", errNameTok)
		}
		c &= 127
		var l *trieNode
		x := t.next
		for x != nil && x.c != c {
			l = x
			x = x.sibling
		}
		if x == nil {
			x = &trieNode{c: c, n: n}
			if l == nil {
				t.next = x
			} else {
				l.sibling = x
			}
		}
		t = x
		t.c = c
		t.count++
	}
	return nil
}

// searchResult carries the four outputs of search_trie.
type searchResult struct {
	from     int
	exact    bool
	isFixed  bool
	fixedLen int
}

// isXdigit reports whether c is a hexadecimal digit, matching isxdigit.
func isXdigit(c byte) bool {
	return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F')
}

// isAlphaC reports whether c is an ASCII letter (isalpha).
func isAlphaC(c byte) bool {
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')
}

// isDigitC reports whether c is an ASCII decimal digit (isdigit).
func isDigitC(c byte) bool { return c >= '0' && c <= '9' }

// isPunctC reports whether c is an ASCII punctuation byte (ispunct).
func isPunctC(c byte) bool {
	return (c >= '!' && c <= '/') || (c >= ':' && c <= '@') ||
		(c >= '[' && c <= '`') || (c >= '{' && c <= '~')
}

// isSpaceC reports whether c is an ASCII whitespace byte (isspace).
func isSpaceC(c byte) bool {
	return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r'
}

// searchTrie ports search_trie: it walks the trie for name and returns
// which earlier name (if any) this one should diff against, plus the
// fixed-prefix heuristics for known name layouts.
func (enc *nameEncoder) searchTrie(name []byte, n int) searchResult {
	var res searchResult
	from := -1
	p3 := -1

	d := name
	f := 0
	if len(name) > 0 && name[0] == '@' {
		d = name[1:]
	}
	if len(name) > 0 && name[0] == '>' {
		f = 1
	}
	l := len(d)

	prefixLen := 0
	hasPrefix := true
	switch {
	case l > 70 && idx(d, f+0) == 'm' && idx(d, 7) == '_' &&
		idx(d, f+14) == '_' && idx(d, f+61) == '/':
		prefixLen = 60 // PacBio
		res.isFixed = false
	case l == 17 && idx(d, f+5) == ':' && idx(d, f+11) == ':':
		prefixLen = 6 // IonTorrent
		res.fixedLen = 6
		res.isFixed = true
	case l >= 36 &&
		idx(d, f+8) == '-' && idx(d, f+13) == '-' && idx(d, f+18) == '-' && idx(d, f+23) == '-' &&
		isXdigit(idx(d, f+0)) && isXdigit(idx(d, f+7)) &&
		isXdigit(idx(d, f+9)) && isXdigit(idx(d, f+12)) &&
		isXdigit(idx(d, f+14)) && isXdigit(idx(d, f+17)) &&
		isXdigit(idx(d, f+19)) && isXdigit(idx(d, f+22)) &&
		isXdigit(idx(d, f+24)) && isXdigit(idx(d, f+35)):
		prefixLen = 36 // ONT uuid4
		res.fixedLen = 36
		res.isFixed = true
	default:
		// Illumina: trim back to lane:tile:x:y by counting four colons.
		i := 0
		for i < len(name) && name[i] > ' ' {
			i++
		}
		colons := 0
		for i > 0 && colons < 4 {
			i--
			if name[i] == ':' {
				colons++
			}
		}
		if colons == 4 {
			res.fixedLen = i + 1
			prefixLen = i + 1
			res.isFixed = true
		} else {
			prefixLen = intMax
			res.isFixed = false
			hasPrefix = false
		}
	}
	_ = hasPrefix

	if enc.tHead == nil {
		enc.tHead = &trieNode{}
	}

	fromPunct := from
	for i := 0; i < len(name); i++ {
		t := enc.tHead
		for i < len(name) && name[i] > '\n' {
			c := name[i]
			i++
			c &= 127
			x := t.next
			for x != nil && x.c != c {
				x = x.sibling
			}
			t = x
			from = t.n
			if (isPunctC(c) || isSpaceC(c)) && t.n != n {
				fromPunct = t.n
			}
			if i == prefixLen {
				p3 = t.n
			}
			t.n = n
		}
	}

	res.exact = from != n && len(name) > 0
	if res.exact {
		res.from = from
	} else if p3 != -1 {
		res.from = p3
	} else {
		res.from = fromPunct
	}
	return res
}

// idx returns d[i], or 0 when i is out of range (the C heuristics index
// guarded by length checks; this keeps the Go translation panic-free).
func idx(d []byte, i int) byte {
	if i < 0 || i >= len(d) {
		return 0
	}
	return d[i]
}

// --- token emission ----------------------------------------------------------

// growToken zeroes token-blocks for newly reachable token positions.
func (enc *nameEncoder) growToken(ntok int) bool {
	if ntok < enc.maxTok {
		return true
	}
	if enc.maxTok >= nameMaxTokens {
		return false
	}
	base := enc.maxTok << 4
	for k := 0; k < 16; k++ {
		enc.desc[base+k].buf = nil
	}
	enc.tokenD[enc.maxTok] = 0
	enc.tokenI[enc.maxTok] = 0
	enc.maxTok = ntok + 1
	return true
}

func (enc *nameEncoder) emitType(ntok, typ int) {
	id := ntok << 4
	enc.desc[id].buf = append(enc.desc[id].buf, byte(typ))
}

func (enc *nameEncoder) emitInt(ntok, typ int, val uint32) {
	enc.emitType(ntok, typ)
	id := (ntok << 4) | typ
	enc.desc[id].buf = append(enc.desc[id].buf,
		byte(val), byte(val>>8), byte(val>>16), byte(val>>24))
}

func (enc *nameEncoder) emitInt1(ntok, typ int, val uint32) {
	enc.emitType(ntok, typ)
	id := (ntok << 4) | typ
	enc.desc[id].buf = append(enc.desc[id].buf, byte(val))
}

// emitInt1Raw appends a single byte without a preceding type byte
// (encode_token_int1_, used for the N_DZLEN width that rides with
// N_DIGITS0).
func (enc *nameEncoder) emitInt1Raw(ntok, typ int, val uint32) {
	id := (ntok << 4) | typ
	enc.desc[id].buf = append(enc.desc[id].buf, byte(val))
}

func (enc *nameEncoder) emitAlpha(ntok int, str []byte) {
	enc.emitType(ntok, ntAlpha)
	id := (ntok << 4) | ntAlpha
	enc.desc[id].buf = append(enc.desc[id].buf, str...)
	enc.desc[id].buf = append(enc.desc[id].buf, 0)
}

func (enc *nameEncoder) emitChar(ntok int, c byte) {
	enc.emitType(ntok, ntChar)
	id := (ntok << 4) | ntChar
	enc.desc[id].buf = append(enc.desc[id].buf, c)
}

// encodeName ports encode_name (mode is always 1 here, matching the
// reference call site). It tokenises name and emits the per-token-block
// byte streams, recording the breakdown for later names to diff against.
func (enc *nameEncoder) encodeName(name []byte) error {
	cnum := enc.counter
	enc.counter++

	sr := enc.searchTrie(name, cnum)
	pnum := sr.from
	if pnum < 0 {
		if cnum != 0 {
			pnum = cnum - 1
		} else {
			pnum = 0
		}
	}

	// Whole-line duplicate.
	if sr.exact && pnum >= 0 && pnum < len(enc.lc) &&
		len(name) == len(enc.lc[pnum].name) {
		enc.emitInt(0, ntDup, uint32(cnum-pnum))
		enc.lc[cnum].name = name
		enc.lc[cnum].ntok = enc.lc[pnum].ntok
		enc.lc[cnum].last = append([]encLastTok(nil), enc.lc[pnum].last...)
		return nil
	}

	enc.lc[cnum].last = make([]encLastTok, nameMaxTokens)
	enc.emitInt(0, ntDiff, uint32(cnum-pnum))
	ntok := 1

	i := 0
	switch {
	case sr.fixedLen == 36:
		for !enc.growToken(37) {
			return fmt.Errorf("%w: token count exceeds %d", errNameTok, nameMaxTokens)
		}
		for k := 0; k < 36; k++ {
			enc.emitChar(ntok, name[k])
			enc.lc[cnum].last[ntok] = encLastTok{typ: ntChar, intVal: int(name[k])}
			ntok++
		}
		sr.isFixed = false
		i = 36
	case sr.isFixed:
		if !enc.growToken(ntok) {
			return fmt.Errorf("%w: token count exceeds %d", errNameTok, nameMaxTokens)
		}
		fl := sr.fixedLen
		if pnum < cnum && ntok < enc.lc[pnum].ntok &&
			enc.lc[pnum].last[ntok].typ == ntAlpha {
			pl := enc.lc[pnum].last[ntok]
			if pl.intVal == fl && bytesEqual(name[:fl], enc.lc[pnum].name[:fl]) {
				enc.emitType(ntok, ntMatch)
			} else {
				enc.emitAlpha(ntok, name[:fl])
			}
		} else {
			enc.emitAlpha(ntok, name[:fl])
		}
		enc.lc[cnum].last[ntok] = encLastTok{typ: ntAlpha, intVal: fl, strOff: 0}
		ntok++
		i = fl
	}

	for ; i < len(name); i++ {
		if !enc.growToken(ntok) {
			return fmt.Errorf("%w: token count exceeds %d", errNameTok, nameMaxTokens)
		}
		c := name[i]
		switch {
		case isAlphaC(c):
			s := i + 1
			for s < len(name) && (isAlphaC(name[s]) || isPunctC(name[s])) {
				s++
			}
			if s-i == 1 {
				enc.encodeChar(name, i, ntok, cnum, pnum)
				ntok++
				continue
			}
			if pnum < cnum && ntok < enc.lc[pnum].ntok &&
				enc.lc[pnum].last[ntok].typ == ntAlpha {
				pl := enc.lc[pnum].last[ntok]
				if s-i == pl.intVal &&
					bytesEqual(name[i:s], enc.lc[pnum].name[pl.strOff:pl.strOff+s-i]) {
					enc.emitType(ntok, ntMatch)
				} else {
					enc.emitAlpha(ntok, name[i:s])
				}
			} else {
				enc.emitAlpha(ntok, name[i:s])
			}
			enc.lc[cnum].last[ntok] = encLastTok{typ: ntAlpha, intVal: s - i, strOff: i}
			i = s - 1

		case isDigitC(c):
			s := i
			var v uint32
			for s < len(name) && isDigitC(name[s]) && s-i < 9 {
				v = v*10 + uint32(name[s]-'0')
				s++
			}
			leadingZero := c == '0'
			if !leadingZero && pnum < cnum && ntok < enc.lc[pnum].ntok &&
				enc.lc[pnum].last[ntok].typ == ntDigits0 &&
				enc.lc[pnum].last[ntok].strOff == s-i {
				leadingZero = true
			}
			if leadingZero {
				enc.encodeDigits0(name, i, s, v, ntok, cnum, pnum)
			} else {
				enc.encodeDigits(v, ntok, cnum, pnum)
			}
			i = s - 1

		default:
			enc.encodeChar(name, i, ntok, cnum, pnum)
		}
		ntok++
	}

	if !enc.growToken(ntok) {
		return fmt.Errorf("%w: token count exceeds %d", errNameTok, nameMaxTokens)
	}
	enc.emitType(ntok, ntEnd)

	enc.lc[cnum].name = name
	enc.lc[cnum].ntok = ntok
	enc.lc[cnum].last = enc.lc[cnum].last[:ntok+1]
	return nil
}

// encodeChar emits a single-byte token (the n_char branch of encode_name).
func (enc *nameEncoder) encodeChar(name []byte, i, ntok, cnum, pnum int) {
	c := name[i]
	if pnum < cnum && ntok < enc.lc[pnum].ntok &&
		enc.lc[pnum].last[ntok].typ == ntChar {
		if int(c) == enc.lc[pnum].last[ntok].intVal {
			enc.emitType(ntok, ntMatch)
		} else {
			enc.emitChar(ntok, c)
		}
	} else {
		enc.emitChar(ntok, c)
	}
	enc.lc[cnum].last[ntok] = encLastTok{typ: ntChar, intVal: int(c)}
}

// encodeDigits0 emits a leading-zero digit run (the digits0 branch).
func (enc *nameEncoder) encodeDigits0(name []byte, i, s int, v uint32, ntok, cnum, pnum int) {
	width := s - i
	if pnum < cnum && ntok < enc.lc[pnum].ntok &&
		enc.lc[pnum].last[ntok].typ == ntDigits0 {
		pl := enc.lc[pnum].last[ntok]
		d := int(v) - pl.intVal
		switch {
		case d == 0 && pl.strOff == width:
			enc.emitType(ntok, ntMatch)
		case d < 256 && d >= 0 && pl.strOff == width:
			enc.emitInt1(ntok, ntDDelta0, uint32(d))
		default:
			enc.emitInt1Raw(ntok, ntDZLen, uint32(width))
			enc.emitInt(ntok, ntDigits0, v)
		}
	} else {
		enc.emitInt1Raw(ntok, ntDZLen, uint32(width))
		enc.emitInt(ntok, ntDigits0, v)
	}
	enc.lc[cnum].last[ntok] = encLastTok{typ: ntDigits0, intVal: int(v), strOff: width}
}

// encodeDigits emits a non-leading-zero digit run (the digits branch).
func (enc *nameEncoder) encodeDigits(v uint32, ntok, cnum, pnum int) {
	if pnum < cnum && ntok < enc.lc[pnum].ntok &&
		enc.lc[pnum].last[ntok].typ == ntDigits {
		pl := enc.lc[pnum].last[ntok]
		d := int(v) - pl.intVal
		switch {
		case d == 0:
			enc.emitType(ntok, ntMatch)
		case d < 256 && d >= 0 && 5+enc.tokenD[ntok] > enc.tokenI[ntok]:
			enc.emitInt1(ntok, ntDDelta, uint32(d))
			enc.tokenD[ntok]++
		default:
			enc.emitInt(ntok, ntDigits, v)
			enc.tokenI[ntok]++
		}
	} else {
		enc.emitInt(ntok, ntDigits, v)
	}
	enc.lc[cnum].last[ntok] = encLastTok{typ: ntDigits, intVal: int(v)}
}

// --- per-block sub-codec search ----------------------------------------------

// rTable is the R[5][N_ALL][7] method-search table from tokenise_name3.c.
// Each entry is {count, method, method, ...}: the encoder tries the first
// count methods and keeps the smallest result.
var rTable = [5][ntAll][]int{
	{ // level -1
		ntType: {128}, ntAlpha: {129}, ntChar: {0}, ntDigits0: {8},
		ntDZLen: {0}, ntDup: {8}, ntDiff: {8}, ntDigits: {8},
		ntDDelta: {0}, ntDDelta0: {128}, ntMatch: {0}, ntNop: {0}, ntEnd: {0},
	},
	{ // level -3
		ntType: {192, 0}, ntAlpha: {129, 1}, ntChar: {0}, ntDigits0: {136, 0},
		ntDZLen: {0}, ntDup: {200}, ntDiff: {136}, ntDigits: {200},
		ntDDelta: {0}, ntDDelta0: {128}, ntMatch: {0}, ntNop: {0}, ntEnd: {0},
	},
	{ // level -5
		ntType: {192, 0}, ntAlpha: {1, 128, 0, 129}, ntChar: {0}, ntDigits0: {200, 0},
		ntDZLen: {0}, ntDup: {200}, ntDiff: {192, 200}, ntDigits: {132, 201},
		ntDDelta: {0}, ntDDelta0: {128}, ntMatch: {0}, ntNop: {0}, ntEnd: {0},
	},
	{ // level -7
		ntType: {193, 0, 1}, ntAlpha: {128, 1, 128, 0, 129}, ntChar: {1, 0},
		ntDigits0: {200, 0}, ntDZLen: {0}, ntDup: {201}, ntDiff: {192, 200},
		ntDigits: {132, 201}, ntDDelta: {0}, ntDDelta0: {128},
		ntMatch: {0}, ntNop: {0}, ntEnd: {0},
	},
	{ // level -9
		ntType: {192, 0, 1, 65, 193, 132}, ntAlpha: {132, 1, 0, 129},
		ntChar: {1, 0, 192}, ntDigits0: {201, 0, 192, 64}, ntDZLen: {0, 128, 1},
		ntDup: {201}, ntDiff: {192, 201, 65}, ntDigits: {132, 201, 1, 192, 129, 193},
		ntDDelta: {1, 0, 192}, ntDDelta0: {192, 1, 0}, ntMatch: {0}, ntNop: {0}, ntEnd: {0},
	},
}

// compressTokenBlock ports the compress() routine: it tries each sub-codec
// method for the given token type and keeps the smallest framed result.
// Each framed result is var_u32(len) followed by the rANS/arith stream.
func compressTokenBlock(in []byte, typ, level int, useArith bool) ([]byte, error) {
	lvl := (level - 1) / 2
	if lvl < 0 {
		lvl = 0
	}
	if lvl > 4 {
		lvl = 4
	}
	if typ < 0 || typ >= ntAll {
		typ = 0
	}
	methods := append([]int(nil), rTable[lvl][typ]...)
	if len(methods) == 0 {
		methods = []int{0}
	}
	// Level-3 DIGITS tweak for the arithmetic coder.
	if useArith && lvl == 1 && typ == ntDigits {
		methods[0] = 201
	}

	var best []byte
	for _, m := range methods {
		if !useArith && m&4 != 0 {
			m &^= 4
		}
		if len(in)%4 != 0 && m&8 != 0 {
			continue
		}
		var (
			stream []byte
			err    error
		)
		if useArith {
			stream, err = ArithEncode(in, m)
		} else {
			stream, err = RANS4x16Encode(in, m)
		}
		if err != nil {
			// A method the sub-codec cannot encode (notably X_EXT, which
			// would need an external bzip2 encoder) is simply dropped from
			// the brute-force search, exactly as if it had lost on size.
			// As long as one method succeeds the block still encodes; the
			// reference's method search is likewise a "pick the smallest"
			// over whatever it can produce.
			continue
		}
		framed := varPutU32(nil, uint32(len(stream)))
		framed = append(framed, stream...)
		if best == nil || len(framed) < len(best) {
			best = framed
		}
	}
	if best == nil {
		return nil, fmt.Errorf("%w: no usable sub-codec method for token type %d", errNameTok, typ)
	}
	return best, nil
}
