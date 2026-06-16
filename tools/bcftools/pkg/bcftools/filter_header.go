package bcftools

import (
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// tagClass describes how a bare identifier in a filter expression should be
// resolved against the VCF header. It mirrors the precedence rules of upstream
// bcftools filter.c::filters_init1 (reference_code/bcftools/filter.c:3429-3440):
// a bare tag that is declared only as FORMAT resolves to FORMAT, a tag declared
// only as INFO resolves to INFO, and a tag declared as BOTH is an ambiguous
// expression that upstream rejects with a fatal error.
type tagClass int

const (
	// tagUnknown means the identifier is not declared in the header as either
	// an INFO or FORMAT tag (it may still be a builtin column or a constant).
	tagUnknown tagClass = iota
	// tagInfo means the identifier is declared only as an INFO tag.
	tagInfo
	// tagFormat means the identifier is declared only as a FORMAT tag.
	tagFormat
	// tagBoth means the identifier is declared as both INFO and FORMAT, which
	// makes a bare reference ambiguous.
	tagBoth
)

// headerTags is the parsed INFO/FORMAT tag inventory of a VCF header, used to
// resolve bare identifiers in filter expressions the way upstream does. A nil
// *headerTags resolves every identifier to tagUnknown, which reproduces the
// header-less behaviour of CompileFilter (bare names fall back to INFO lookup
// at evaluation time, then to a string constant).
type headerTags struct {
	info   map[string]bool
	format map[string]bool
}

// newHeaderTags parses the ##INFO and ##FORMAT meta lines of hdr into a
// headerTags inventory. A nil header yields a nil inventory.
func newHeaderTags(hdr *vcf.Header) *headerTags {
	if hdr == nil {
		return nil
	}
	ht := &headerTags{
		info:   make(map[string]bool),
		format: make(map[string]bool),
	}
	for _, m := range hdr.MetaInfo {
		kind, id := structuredID(m)
		if id == "" {
			continue
		}
		switch kind {
		case "INFO":
			ht.info[id] = true
		case "FORMAT":
			ht.format[id] = true
		}
	}
	return ht
}

// classify returns how the bare identifier tag should be resolved against the
// header, following the upstream precedence rules. A nil inventory always
// returns tagUnknown.
func (ht *headerTags) classify(tag string) tagClass {
	if ht == nil {
		return tagUnknown
	}
	isInfo := ht.info[tag]
	isFmt := ht.format[tag]
	switch {
	case isInfo && isFmt:
		return tagBoth
	case isFmt:
		return tagFormat
	case isInfo:
		return tagInfo
	default:
		return tagUnknown
	}
}

// ambiguousTagError builds the fatal-error message upstream emits for a bare
// identifier that is declared as both INFO and FORMAT, byte-matching
// filter.c:3437 so the CLI surfaces the same diagnostic.
func ambiguousTagError(tag string) error {
	return fmt.Errorf("Error: ambiguous filtering expression, both INFO/%s and FORMAT/%s are defined in the VCF header.", tag, tag)
}
