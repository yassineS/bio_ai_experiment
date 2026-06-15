// Pure-Go port of the VariantKey encoder used by the add-variantkey and
// variantkey-hex bcftools plugins. It mirrors reference_code/bcftools/variantkey.h
// (the tecnickcom/variantkey library, MIT) bit-for-bit so the 64-bit codes and
// their hexadecimal renderings are byte-identical to the upstream plugins.
//
// A VariantKey packs CHROM (5 bits), POS (28 bits, 0-based) and a REF+ALT code
// (31 bits) into a single sortable uint64. The REF+ALT code is reversible when
// REF and ALT together fit in 11 ACGT bases, otherwise a MurmurHash3-derived
// 31-bit hash (with the low bit set) is used.
package bcftools

const (
	vkShiftChrom = 59
	vkShiftPos   = 31
	vkMaxUint32  = uint32(0xFFFFFFFF)
)

// encodeNumericChrom encodes a purely numeric chromosome string (e.g. "12").
// It returns 0 if any non-digit character is encountered, matching upstream's
// encode_numeric_chrom.
func encodeNumericChrom(chrom string) uint8 {
	if len(chrom) == 0 {
		return 0
	}
	v := uint8(chrom[0] - '0')
	for i := 1; i < len(chrom); i++ {
		if chrom[i] > '9' || chrom[i] < '0' {
			return 0
		}
		v = v*10 + (chrom[i] - '0')
	}
	return v
}

// hasChromChrPrefix reports whether chrom starts with a case-insensitive "chr"
// prefix and is longer than 3 characters, matching has_chrom_chr_prefix.
func hasChromChrPrefix(chrom string) bool {
	return len(chrom) > 3 &&
		(chrom[0] == 'c' || chrom[0] == 'C') &&
		(chrom[1] == 'h' || chrom[1] == 'H') &&
		(chrom[2] == 'r' || chrom[2] == 'R')
}

// encodeChrom returns the 5-bit chromosome code (1-25, with X=23, Y=24, M/MT=25)
// or 0 for unrecognised input, porting encode_chrom including the one-character
// X/Y/M(/MT) map.
func encodeChrom(chrom string) uint8 {
	if hasChromChrPrefix(chrom) {
		chrom = chrom[3:]
	}
	if len(chrom) == 0 {
		return 0
	}
	if chrom[0] >= '0' && chrom[0] <= '9' {
		return encodeNumericChrom(chrom)
	}
	if len(chrom) == 1 || (len(chrom) == 2 && (chrom[1] == 'T' || chrom[1] == 't')) {
		switch chrom[0] {
		case 'X', 'x':
			return 23
		case 'Y', 'y':
			return 24
		case 'M', 'm':
			return 25
		}
		return 0
	}
	return 0
}

// encodeBase maps an ACGT nucleotide (any case) to 0-3, or 4 for anything else,
// porting the encode_base lookup table.
func encodeBase(c byte) uint32 {
	switch c {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	}
	return 4
}

// encodeAllele packs the bases of str into h starting at *bitpos (2 bits each),
// returning false if a non-ACGT base is found. It mutates h and bitpos exactly
// as the C encode_allele does.
func encodeAllele(h *uint32, bitpos *uint8, str string) bool {
	for i := 0; i < len(str); i++ {
		v := encodeBase(str[i])
		if v > 3 {
			return false
		}
		*bitpos -= 2
		*h |= v << *bitpos
	}
	return true
}

// encodeRefaltRev produces the reversible 31-bit REF+ALT code (used when the
// combined length is <= 11 and all bases are ACGT), porting encode_refalt_rev.
// It returns vkMaxUint32 on a non-ACGT base to signal a fall-through to hashing.
func encodeRefaltRev(ref, alt string) uint32 {
	var h uint32
	h |= uint32(len(ref)) << 27
	h |= uint32(len(alt)) << 23
	bitpos := uint8(23)
	if !encodeAllele(&h, &bitpos, ref) || !encodeAllele(&h, &bitpos, alt) {
		return vkMaxUint32
	}
	return h
}

// muxhash mixes two 32-bit values with a MurmurHash3-like step, porting muxhash.
func muxhash(k, h uint32) uint32 {
	k *= 0xcc9e2d51
	k = (k >> 17) | (k << 15)
	k *= 0x1b873593
	h ^= k
	h = (h >> 19) | (h << 13)
	return (h * 5) + 0xe6546b64
}

// encodePackchar maps a character to its 5-bit packed code, porting
// encode_packchar (non-letters map to 27).
func encodePackchar(c byte) uint32 {
	if c < 'A' {
		return 27
	}
	if c >= 'a' {
		return uint32(c-'a') + 1
	}
	return uint32(c-'A') + 1
}

// packCharsTail packs a trailing block of 1-5 characters, porting
// pack_chars_tail with its fall-through switch.
func packCharsTail(str string) uint32 {
	var h uint32
	size := len(str)
	pos := size - 1
	switch size {
	case 5:
		h ^= encodePackchar(str[pos]) << (1 + 5*1)
		pos--
		fallthrough
	case 4:
		h ^= encodePackchar(str[pos]) << (1 + 5*2)
		pos--
		fallthrough
	case 3:
		h ^= encodePackchar(str[pos]) << (1 + 5*3)
		pos--
		fallthrough
	case 2:
		h ^= encodePackchar(str[pos]) << (1 + 5*4)
		pos--
		fallthrough
	case 1:
		h ^= encodePackchar(str[pos]) << (1 + 5*5)
	}
	return h
}

// packChars packs a full block of 6 characters from the start of str, porting
// pack_chars. str must be at least 6 bytes long.
func packChars(str string) uint32 {
	return (encodePackchar(str[5]) << 1) ^
		(encodePackchar(str[4]) << (1 + 5*1)) ^
		(encodePackchar(str[3]) << (1 + 5*2)) ^
		(encodePackchar(str[2]) << (1 + 5*3)) ^
		(encodePackchar(str[1]) << (1 + 5*4)) ^
		(encodePackchar(str[0]) << (1 + 5*5))
}

// hash32 returns the 32-bit hash of a nucleotide string, porting hash32.
func hash32(str string) uint32 {
	var h uint32
	for len(str) >= 6 {
		h = muxhash(packChars(str), h)
		str = str[6:]
	}
	if len(str) > 0 {
		h = muxhash(packCharsTail(str), h)
	}
	return h
}

// encodeRefaltHash produces the 31-bit hashed REF+ALT code (low bit set to mark
// HASH mode), porting encode_refalt_hash including the MurmurHash3 finaliser.
func encodeRefaltHash(ref, alt string) uint32 {
	h := muxhash(hash32(alt), muxhash(0x3, hash32(ref)))
	h ^= h >> 16
	h *= 0x85ebca6b
	h ^= h >> 13
	h *= 0xc2b2ae35
	h ^= h >> 16
	return (h >> 1) | 0x1
}

// encodeRefalt returns the REF+ALT code, preferring the reversible encoding when
// the combined length is <= 11 and all bases are ACGT, porting encode_refalt.
func encodeRefalt(ref, alt string) uint32 {
	if len(ref)+len(alt) <= 11 {
		if h := encodeRefaltRev(ref, alt); h != vkMaxUint32 {
			return h
		}
	}
	return encodeRefaltHash(ref, alt)
}

// encodeVariantkey packs the pre-encoded components into the 64-bit key,
// porting encode_variantkey.
func encodeVariantkey(chrom uint8, pos uint32, refalt uint32) uint64 {
	return (uint64(chrom) << vkShiftChrom) | (uint64(pos) << vkShiftPos) | uint64(refalt)
}

// variantKey returns the 64-bit VariantKey for the given CHROM, 0-based POS,
// REF and ALT, porting the top-level variantkey() function. POS must already be
// 0-based (upstream passes rec->pos, which htslib stores 0-based).
func variantKey(chrom string, pos uint32, ref, alt string) uint64 {
	return encodeVariantkey(encodeChrom(chrom), pos, encodeRefalt(ref, alt))
}

const vkHexDigits = "0123456789abcdef"

// variantKeyHex renders a 64-bit VariantKey as a 16-character lowercase
// hexadecimal string, matching variantkey_hex / hex_uint64_t (which uses the
// "%016" PRIx64 lowercase formatting).
func variantKeyHex(vk uint64) string {
	var buf [16]byte
	for i := 15; i >= 0; i-- {
		buf[i] = vkHexDigits[vk&0xF]
		vk >>= 4
	}
	return string(buf[:])
}
