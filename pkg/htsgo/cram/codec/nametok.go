package codec

// nametok.go ports the htscodecs name tokeniser (tokenise_name3.c) — the
// codec behind CRAM block compression method 8. It splits each read name
// into a sequence of typed tokens (alpha runs, digit runs, leading-zero
// digit runs, single characters and delimiters), models each name as a
// diff against an earlier name, and routes each per-token-position byte
// stream through a sub-codec (rANS 4x16, the arith_dynamic range coder,
// or a verbatim copy). The compressed block is fully self-describing, so
// a method-8 CRAM block decodes standalone.
//
// Only the decoder is required to be byte-exact; the encoder here
// round-trips correctly but does not attempt to reproduce the reference
// encoder's per-level sub-codec search (see NameTokEncode).

import (
	"errors"
	"fmt"
)

// Name-token types, mirroring enum name_type in tokenise_name3.c. They
// double as the low four bits of a token-block identifier.
const (
	ntType    = 0  // N_TYPE: the per-token-position stream of token types.
	ntAlpha   = 1  // N_ALPHA: a NUL-terminated run of letters/punctuation.
	ntChar    = 2  // N_CHAR: a single byte.
	ntDigits0 = 3  // N_DIGITS0: a fixed-width (leading-zero) decimal value.
	ntDZLen   = 4  // N_DZLEN: the width byte accompanying an N_DIGITS0 value.
	ntDup     = 5  // N_DUP: the whole name duplicates an earlier one.
	ntDiff    = 6  // N_DIFF: distance back to the name this one diffs against.
	ntDigits  = 7  // N_DIGITS: a variable-width decimal value.
	ntDDelta  = 8  // N_DDELTA: an 8-bit delta against the previous N_DIGITS.
	ntDDelta0 = 9  // N_DDELTA0: an 8-bit delta against the previous N_DIGITS0.
	ntMatch   = 10 // N_MATCH: this token equals the corresponding earlier token.
	ntNop     = 11 // N_NOP: an unused token slot.
	ntEnd     = 12 // N_END: terminates a name.
	ntAll     = 13 // N_ALL: count of token types.
)

// nametok hard limits, copied verbatim from tokenise_name3.c. They cap
// every attacker-controlled count so a malformed block cannot drive an
// unbounded allocation or an unterminated loop.
const (
	nameMaxTokens  = 128 // MAX_TOKENS.
	nameMaxTBlocks = nameMaxTokens * 16
	nameMaxNames   = 10_000_000 // create_context's 10M-record ceiling.
	// nameMaxNameLen bounds a single decoded name. The reference grants
	// decode_name a ulen-sized scratch buffer; we cap each name instead so
	// a corrupt stream cannot make one token-block describe a giant name.
	nameMaxNameLen = 1 << 20
)

// errNameTok wraps every failure from the name-tokeniser decoder.
var errNameTok = errors.New("nametok")

// lastTok records, per token position, what the previous name decoded to
// at that position. The current name diffs against it (mirroring
// last_context_tok in the C source).
type lastTok struct {
	typ    int // token type (ntChar, ntAlpha, ntDigits, ntDigits0, ...).
	intVal int // numeric value, or the char/alpha length depending on typ.
	strOff int // for ntAlpha: byte offset of the run within lastName.
}

// lastCtx is the decoded state of one name: its bytes and its per-token
// breakdown, kept so later names can diff against it.
type lastCtx struct {
	lastName []byte
	lastNtok int
	last     []lastTok
}

// nameDescriptor is one decoded token-block: the uncompressed bytes for a
// single (token-position, token-type) pair and a read cursor into them.
type nameDescriptor struct {
	buf []byte
	pos int // read cursor; bufL in the C source.
}

// nameContext holds the decoder's working state: the decoded token-blocks
// and the per-name diff history.
type nameContext struct {
	desc     [nameMaxTBlocks]nameDescriptor
	lc       []lastCtx
	counter  int
	maxTok   int
	maxNames int // create_context's nreads+1.
	nreads   int // declared record count; sizes synthesised type blocks.
}

// NameTokDecode decompresses a complete htscodecs name-tokeniser block
// (CRAM compression method 8) and returns the decoded read names joined
// by NUL bytes, with a trailing NUL after the final name. It is the
// inverse of NameTokEncode and of htscodecs' tok3_decode_names.
//
// The input is untrusted: every token count, stream length and decoded
// name length is bounds-checked, allocations are bounded, and a
// malformed block yields an error rather than a panic or a runaway loop.
func NameTokDecode(in []byte) ([]byte, error) {
	if len(in) < 9 {
		return nil, fmt.Errorf("%w: block shorter than the 9-byte header", errNameTok)
	}

	ulen := uint32(in[0]) | uint32(in[1])<<8 | uint32(in[2])<<16 | uint32(in[3])<<24
	if ulen >= 1<<31-1024 {
		return nil, fmt.Errorf("%w: declared uncompressed size %d implausibly large", errNameTok, ulen)
	}
	nreads := int(uint32(in[4]) | uint32(in[5])<<8 | uint32(in[6])<<16 | uint32(in[7])<<24)
	if nreads < 0 || nreads > nameMaxNames {
		return nil, fmt.Errorf("%w: name count %d out of range", errNameTok, nreads)
	}
	useArith := in[8] != 0

	ctx := &nameContext{
		// create_context allocates nreads+1 last-contexts.
		lc:       make([]lastCtx, nreads+1),
		maxNames: nreads + 1,
		nreads:   nreads,
		maxTok:   1,
	}

	if err := ctx.unpackDescriptors(in, useArith); err != nil {
		return nil, err
	}

	// Decode names until decode_name signals the end of the stream. The
	// initial capacity is hinted by the declared size but capped so a
	// corrupt ulen cannot drive a huge up-front allocation — the slice
	// still grows on demand, bounded by the ulen check inside the loop.
	hint := int(ulen) + 1
	if hint > 1<<20 {
		hint = 1 << 20
	}
	out := make([]byte, 0, hint)
	for {
		name, done, err := ctx.decodeName()
		if err != nil {
			return nil, err
		}
		if done {
			break
		}
		out = append(out, name...)
		if len(out) > int(ulen)+1024 {
			return nil, fmt.Errorf("%w: decoded output exceeds declared size %d", errNameTok, ulen)
		}
	}
	return out, nil
}

// unpackDescriptors reads the serialised token-blocks that follow the
// 9-byte header, decompressing each through the rANS/arith sub-codec (or
// resolving it as a duplicate of an earlier block). It mirrors the
// descriptor-unpacking loop of tok3_decode_names.
func (ctx *nameContext) unpackDescriptors(in []byte, useArith bool) error {
	sz := len(in)
	o := 9
	tnum := -1

	for o < sz {
		ttype := in[o]
		o++

		if ttype&64 != 0 {
			// A duplicate of an earlier token-block.
			if o+2 > sz {
				return fmt.Errorf("%w: truncated duplicate-block reference", errNameTok)
			}
			j := int(in[o])<<4 + int(in[o+1])
			o += 2

			if ttype&128 != 0 {
				tnum++
				if tnum >= nameMaxTokens {
					return fmt.Errorf("%w: token count exceeds %d", errNameTok, nameMaxTokens)
				}
				ctx.maxTok = tnum + 1
				ctx.clearToken(tnum)
			}
			if tnum < 0 {
				return fmt.Errorf("%w: token-block before its token header", errNameTok)
			}
			if err := ctx.synthesiseTypeBlock(tnum, ttype); err != nil {
				return err
			}

			i := (tnum << 4) | int(ttype&15)
			if j >= i {
				return fmt.Errorf("%w: duplicate block %d references later block %d", errNameTok, i, j)
			}
			if ctx.desc[j].buf == nil {
				return fmt.Errorf("%w: duplicate block %d references absent block %d", errNameTok, i, j)
			}
			src := ctx.desc[j].buf
			cp := make([]byte, len(src))
			copy(cp, src)
			ctx.desc[i].buf = cp
			ctx.desc[i].pos = 0
			continue
		}

		if ttype&128 != 0 {
			tnum++
			if tnum >= nameMaxTokens {
				return fmt.Errorf("%w: token count exceeds %d", errNameTok, nameMaxTokens)
			}
			ctx.maxTok = tnum + 1
			ctx.clearToken(tnum)
		}
		if tnum < 0 {
			return fmt.Errorf("%w: token-block before its token header", errNameTok)
		}
		if err := ctx.synthesiseTypeBlock(tnum, ttype); err != nil {
			return err
		}

		// A self-describing compressed sub-stream follows. The reference
		// frames it as var_u32(clen) then the rANS/arith stream itself,
		// whose own header carries the format byte and raw-size varint.
		clen, nb, ok := varGetU32(in, o)
		if !ok {
			return fmt.Errorf("%w: truncated compressed-length varint", errNameTok)
		}
		streamStart := nb
		streamEnd := streamStart + int(clen)
		if int(clen) < 0 || streamEnd < streamStart || streamEnd > sz {
			return fmt.Errorf("%w: compressed sub-stream length %d overruns the block", errNameTok, clen)
		}

		i := (tnum << 4) | int(ttype&15)
		if i < 0 || i >= nameMaxTBlocks {
			return fmt.Errorf("%w: token-block index %d out of range", errNameTok, i)
		}

		var (
			dec []byte
			err error
		)
		if useArith {
			dec, err = ArithDecode(in[streamStart:streamEnd])
		} else {
			dec, err = RANS4x16Decode(in[streamStart:streamEnd])
		}
		if err != nil {
			return fmt.Errorf("%w: token-block %d sub-codec: %v", errNameTok, i, err)
		}
		if len(dec) > nameMaxNameLen*uint64Cap(ctx.maxNames) {
			return fmt.Errorf("%w: token-block %d decoded to %d bytes", errNameTok, i, len(dec))
		}
		ctx.desc[i].buf = dec
		ctx.desc[i].pos = 0
		o = streamEnd
	}
	return nil
}

// uint64Cap returns max(n, 1) so the decoded-stream sanity bound never
// collapses to zero when a block declares no names.
func uint64Cap(n int) int {
	if n < 1 {
		return 1
	}
	return n
}

// clearToken zeroes the sixteen token-blocks belonging to token tnum,
// matching the memset that precedes each new token in the C decoder.
func (ctx *nameContext) clearToken(tnum int) {
	base := tnum << 4
	for k := 0; k < 16; k++ {
		ctx.desc[base+k] = nameDescriptor{}
	}
}

// synthesiseTypeBlock regenerates an elided N_TYPE (type-0) block. The
// encoder drops a token's type stream when every entry past the first is
// N_MATCH; the decoder rebuilds it as [type, N_MATCH, N_MATCH, ...] when
// it meets a token whose first descriptor carries a non-zero type.
func (ctx *nameContext) synthesiseTypeBlock(tnum int, ttype byte) error {
	if ttype&15 == 0 || ttype&128 == 0 {
		return nil
	}
	base := tnum << 4
	// The reference sizes the synthesised type stream at nreads bytes —
	// one type entry per name — not maxNames; an off-by-one here would
	// invent a phantom token for a non-existent (nreads+1)th name.
	buf := make([]byte, ctx.nreads)
	if len(buf) > 0 {
		buf[0] = ttype & 15
		for k := 1; k < len(buf); k++ {
			buf[k] = ntMatch
		}
	}
	ctx.desc[base].buf = buf
	ctx.desc[base].pos = 0
	return nil
}

// --- per-stream readers, mirroring the decode_token_* helpers ----------------

// readType reads one token-type byte from token ntok's type block.
func (ctx *nameContext) readType(ntok int) (int, bool) {
	d := &ctx.desc[ntok<<4]
	if d.pos >= len(d.buf) {
		return -1, false
	}
	v := int(d.buf[d.pos])
	d.pos++
	return v, true
}

// readInt reads a little-endian uint32 from token ntok's type-typ block.
func (ctx *nameContext) readInt(ntok, typ int) (uint32, bool) {
	d := &ctx.desc[(ntok<<4)|typ]
	if d.pos+4 > len(d.buf) {
		return 0, false
	}
	v := uint32(d.buf[d.pos]) | uint32(d.buf[d.pos+1])<<8 |
		uint32(d.buf[d.pos+2])<<16 | uint32(d.buf[d.pos+3])<<24
	d.pos += 4
	return v, true
}

// readInt1 reads a single byte from token ntok's type-typ block.
func (ctx *nameContext) readInt1(ntok, typ int) (uint32, bool) {
	d := &ctx.desc[(ntok<<4)|typ]
	if d.pos >= len(d.buf) {
		return 0, false
	}
	v := uint32(d.buf[d.pos])
	d.pos++
	return v, true
}

// readChar reads one byte from token ntok's N_CHAR block.
func (ctx *nameContext) readChar(ntok int) (byte, bool) {
	d := &ctx.desc[(ntok<<4)|ntChar]
	if d.pos >= len(d.buf) {
		return 0, false
	}
	c := d.buf[d.pos]
	d.pos++
	return c, true
}

// readAlpha reads one NUL-terminated alpha run from token ntok's N_ALPHA
// block and appends it (without the terminator) to dst, returning the
// extended slice and the run length. It bounds the run by nameMaxNameLen.
func (ctx *nameContext) readAlpha(ntok int, dst []byte) ([]byte, int, bool) {
	d := &ctx.desc[(ntok<<4)|ntAlpha]
	if d.pos >= len(d.buf) {
		return dst, 0, false
	}
	n := 0
	for d.pos < len(d.buf) {
		c := d.buf[d.pos]
		d.pos++
		if c == 0 {
			return dst, n, true
		}
		dst = append(dst, c)
		n++
		if n > nameMaxNameLen {
			return dst, 0, false
		}
	}
	// Stream ended without a terminator: the C reader still returns the
	// bytes consumed so far.
	return dst, n, true
}

// --- integer formatting, mirroring append_uint32_* ---------------------------

// appendUintFixed appends v as a zero-padded decimal of exactly width
// digits (append_uint32_fixed). width is clamped to 0..9.
func appendUintFixed(dst []byte, v uint32, width int) []byte {
	if width < 0 {
		width = 0
	}
	if width > 9 {
		width = 9
	}
	var buf [9]byte
	for k := width - 1; k >= 0; k-- {
		buf[k] = byte('0' + v%10)
		v /= 10
	}
	return append(dst, buf[:width]...)
}

// appendUintVar appends v as a minimal-width decimal (append_uint32_var).
func appendUintVar(dst []byte, v uint32) []byte {
	if v == 0 {
		return append(dst, '0')
	}
	var buf [10]byte
	n := 0
	for v > 0 {
		buf[n] = byte('0' + v%10)
		v /= 10
		n++
	}
	for k := n - 1; k >= 0; k-- {
		dst = append(dst, buf[k])
	}
	return dst
}

// --- name decoding -----------------------------------------------------------

// decodeName decodes one read name. It returns the name with its trailing
// NUL, or done=true once the stream is exhausted, mirroring decode_name.
func (ctx *nameContext) decodeName() (name []byte, done bool, err error) {
	cnum := ctx.counter
	ctx.counter++
	if cnum >= ctx.maxNames {
		return nil, false, fmt.Errorf("%w: name index %d exceeds declared count", errNameTok, cnum)
	}

	t0, ok := ctx.readType(0)
	if !ok {
		// The type-0 stream is exhausted: the block is fully decoded.
		return nil, true, nil
	}
	if t0 < 0 || t0 >= ctx.maxTok*16 {
		return nil, true, nil
	}

	dist, ok := ctx.readInt(0, t0)
	if !ok || dist > uint32(cnum) {
		return nil, false, fmt.Errorf("%w: bad diff/dup distance", errNameTok)
	}
	pnum := cnum - int(dist)
	if pnum < 0 {
		pnum = 0
	}

	if t0 == ntDup {
		if pnum == cnum {
			return nil, false, fmt.Errorf("%w: name duplicates itself", errNameTok)
		}
		src := ctx.lc[pnum].lastName
		buf := make([]byte, len(src))
		copy(buf, src)
		ctx.lc[cnum].lastName = buf
		ctx.lc[cnum].lastNtok = ctx.lc[pnum].lastNtok
		ctx.lc[cnum].last = append([]lastTok(nil), ctx.lc[pnum].last...)
		return buf, false, nil
	}

	out := make([]byte, 0, 64)
	ctx.lc[cnum].last = make([]lastTok, nameMaxTokens)

	for ntok := 1; ntok < nameMaxTokens && ntok < ctx.maxTok; ntok++ {
		tok, ok := ctx.readType(ntok)
		if !ok {
			return nil, false, fmt.Errorf("%w: token %d type stream exhausted", errNameTok, ntok)
		}

		switch tok {
		case ntChar:
			c, ok := ctx.readChar(ntok)
			if !ok {
				return nil, false, fmt.Errorf("%w: token %d char stream exhausted", errNameTok, ntok)
			}
			out = append(out, c)
			ctx.lc[cnum].last[ntok] = lastTok{typ: ntChar, intVal: int(c)}

		case ntAlpha:
			var n int
			out, n, ok = ctx.readAlpha(ntok, out)
			if !ok {
				return nil, false, fmt.Errorf("%w: token %d alpha stream exhausted", errNameTok, ntok)
			}
			ctx.lc[cnum].last[ntok] = lastTok{typ: ntAlpha, strOff: len(out) - n, intVal: n}

		case ntDigits0:
			vl, ok := ctx.readInt1(ntok, ntDZLen)
			if !ok {
				return nil, false, fmt.Errorf("%w: token %d dzlen stream exhausted", errNameTok, ntok)
			}
			v, ok := ctx.readInt(ntok, ntDigits0)
			if !ok {
				return nil, false, fmt.Errorf("%w: token %d digits0 stream exhausted", errNameTok, ntok)
			}
			out = appendUintFixed(out, v, int(vl))
			ctx.lc[cnum].last[ntok] = lastTok{typ: ntDigits0, intVal: int(v), strOff: int(vl)}

		case ntDDelta0:
			if ntok >= ctx.lc[pnum].lastNtok {
				return nil, false, fmt.Errorf("%w: token %d delta0 with no prior token", errNameTok, ntok)
			}
			d, ok := ctx.readInt1(ntok, ntDDelta0)
			if !ok {
				return nil, false, fmt.Errorf("%w: token %d delta0 stream exhausted", errNameTok, ntok)
			}
			prev := ctx.lc[pnum].last[ntok]
			v := d + uint32(prev.intVal)
			out = appendUintFixed(out, v, prev.strOff)
			ctx.lc[cnum].last[ntok] = lastTok{typ: ntDigits0, intVal: int(v), strOff: prev.strOff}

		case ntDigits:
			v, ok := ctx.readInt(ntok, ntDigits)
			if !ok {
				return nil, false, fmt.Errorf("%w: token %d digits stream exhausted", errNameTok, ntok)
			}
			out = appendUintVar(out, v)
			ctx.lc[cnum].last[ntok] = lastTok{typ: ntDigits, intVal: int(v)}

		case ntDDelta:
			if ntok >= ctx.lc[pnum].lastNtok {
				return nil, false, fmt.Errorf("%w: token %d delta with no prior token", errNameTok, ntok)
			}
			d, ok := ctx.readInt1(ntok, ntDDelta)
			if !ok {
				return nil, false, fmt.Errorf("%w: token %d delta stream exhausted", errNameTok, ntok)
			}
			v := d + uint32(ctx.lc[pnum].last[ntok].intVal)
			out = appendUintVar(out, v)
			ctx.lc[cnum].last[ntok] = lastTok{typ: ntDigits, intVal: int(v)}

		case ntNop:
			ctx.lc[cnum].last[ntok] = lastTok{typ: ntNop}

		case ntMatch:
			if ntok >= ctx.lc[pnum].lastNtok {
				return nil, false, fmt.Errorf("%w: token %d match with no prior token", errNameTok, ntok)
			}
			prev := ctx.lc[pnum].last[ntok]
			switch prev.typ {
			case ntChar:
				out = append(out, byte(prev.intVal))
				ctx.lc[cnum].last[ntok] = lastTok{typ: ntChar, intVal: prev.intVal}
			case ntAlpha:
				if prev.intVal < 0 || prev.strOff < 0 ||
					prev.strOff+prev.intVal > len(ctx.lc[pnum].lastName) {
					return nil, false, fmt.Errorf("%w: token %d match alpha out of range", errNameTok, ntok)
				}
				start := len(out)
				out = append(out, ctx.lc[pnum].lastName[prev.strOff:prev.strOff+prev.intVal]...)
				ctx.lc[cnum].last[ntok] = lastTok{typ: ntAlpha, strOff: start, intVal: prev.intVal}
			case ntDigits:
				out = appendUintVar(out, uint32(prev.intVal))
				ctx.lc[cnum].last[ntok] = lastTok{typ: ntDigits, intVal: prev.intVal}
			case ntDigits0:
				out = appendUintFixed(out, uint32(prev.intVal), prev.strOff)
				ctx.lc[cnum].last[ntok] = lastTok{typ: ntDigits0, intVal: prev.intVal, strOff: prev.strOff}
			default:
				return nil, false, fmt.Errorf("%w: token %d match against unmatchable type %d", errNameTok, ntok, prev.typ)
			}

		default: // an elided N_END or an explicit N_END.
			out = append(out, 0)
			ctx.lc[cnum].last[ntok] = lastTok{typ: ntEnd}
			ctx.lc[cnum].lastName = out
			ctx.lc[cnum].lastNtok = ntok
			ctx.lc[cnum].last = ctx.lc[cnum].last[:ntok+1]
			return out, false, nil
		}

		if len(out) > nameMaxNameLen {
			return nil, false, fmt.Errorf("%w: decoded name exceeds %d bytes", errNameTok, nameMaxNameLen)
		}
	}

	return nil, false, fmt.Errorf("%w: name not terminated within %d tokens", errNameTok, nameMaxTokens)
}
