package cram

import (
	"fmt"
	"sort"
)

// huffmanCode is one canonical Huffman code word: the symbol it stands
// for, its bit length, and the integer value of its code bits.
type huffmanCode struct {
	symbol int32
	length int32
	code   uint32
}

// huffmanTable is a canonical Huffman decoder built from an alphabet
// and a parallel array of per-symbol bit lengths, exactly as a CRAM
// HUFFMAN encoding stores them.
type huffmanTable struct {
	codes []huffmanCode // sorted by (length, code) for the decode walk.
	// degenerate holds the single symbol when the alphabet has exactly
	// one code word of length zero; such a series consumes no bits.
	degenerate    bool
	degenerateSym int32
}

// newHuffmanTable builds a canonical Huffman decoder. CRAM stores only
// the alphabet and the bit lengths; the canonical code words are
// derived here using the standard algorithm (RFC 1951 §3.2.2): sort by
// (bit length, symbol order) and assign consecutively increasing codes,
// left-shifting at each length boundary.
//
// A one-symbol alphabet with bit length zero is the common CRAM idiom
// for a constant data series; it is represented as a degenerate table
// that yields its symbol without reading any bits.
func newHuffmanTable(symbols, lengths []int32) (*huffmanTable, error) {
	if len(symbols) != len(lengths) {
		return nil, fmt.Errorf("cram: huffman alphabet/length size mismatch (%d vs %d)",
			len(symbols), len(lengths))
	}
	if len(symbols) == 0 {
		return nil, fmt.Errorf("cram: huffman encoding has an empty alphabet")
	}
	t := &huffmanTable{}
	if len(symbols) == 1 {
		// A single code word: it is the whole alphabet. CRAM writes
		// such constant series with bit length 0, meaning "no bits on
		// the wire". Treat any single-symbol table as degenerate.
		t.degenerate = true
		t.degenerateSym = symbols[0]
		return t, nil
	}
	// Canonical ordering: ascending bit length, then ascending symbol
	// value. The tie-break must be the symbol value, not the input
	// alphabet order — htslib's writer happens to emit each length
	// group already sorted, but a spec-conformant file from another
	// encoder need not, and sorting by input order would then assign
	// the wrong codes and decode silently wrong.
	type entry struct {
		sym, length int32
	}
	entries := make([]entry, len(symbols))
	for i := range symbols {
		if lengths[i] < 0 || lengths[i] > 32 {
			return nil, fmt.Errorf("cram: huffman symbol %d has invalid bit length %d", symbols[i], lengths[i])
		}
		entries[i] = entry{sym: symbols[i], length: lengths[i]}
	}
	sort.SliceStable(entries, func(a, b int) bool {
		if entries[a].length != entries[b].length {
			return entries[a].length < entries[b].length
		}
		return entries[a].sym < entries[b].sym
	})
	var code uint32
	var prevLen int32 = -1
	t.codes = make([]huffmanCode, 0, len(entries))
	for _, e := range entries {
		if e.length == 0 {
			return nil, fmt.Errorf("cram: huffman symbol %d has zero bit length in a multi-symbol alphabet", e.sym)
		}
		if prevLen < 0 {
			prevLen = e.length
		}
		for prevLen < e.length {
			code <<= 1
			prevLen++
		}
		t.codes = append(t.codes, huffmanCode{symbol: e.sym, length: e.length, code: code})
		code++
	}
	return t, nil
}

// decode reads one symbol from the bit reader by walking the canonical
// code lengths shortest-first. It returns the degenerate symbol without
// consuming any bits when the table is degenerate.
func (t *huffmanTable) decode(br *bitReader) (int32, error) {
	if t.degenerate {
		return t.degenerateSym, nil
	}
	var acc uint32
	var nbits int32
	idx := 0
	for {
		b, err := br.readBit()
		if err != nil {
			return 0, err
		}
		acc = (acc << 1) | b
		nbits++
		if nbits > 32 {
			return 0, fmt.Errorf("cram: huffman code longer than 32 bits (corrupt bitstream)")
		}
		// Codes are sorted by (length, code); all words of the current
		// length form a contiguous run. Check every word whose length
		// equals nbits.
		for idx < len(t.codes) && t.codes[idx].length < nbits {
			idx++
		}
		for i := idx; i < len(t.codes) && t.codes[i].length == nbits; i++ {
			if t.codes[i].code == acc {
				return t.codes[i].symbol, nil
			}
		}
	}
}
