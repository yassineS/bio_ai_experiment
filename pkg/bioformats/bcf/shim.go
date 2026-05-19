// Package bcf is a re-export shim for the relocated
// `github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf`. See PR-B of
// the htsgo migration roadmap.
//
// Deprecated: import `github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf`.
package bcf

import (
	htsgo "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
)

type (
	DictEntry  = htsgo.DictEntry
	DictKind   = htsgo.DictKind
	Header     = htsgo.Header
	Reader     = htsgo.Reader
	Record     = htsgo.Record
	TypedValue = htsgo.TypedValue
	Writer     = htsgo.Writer
)

const (
	TypeMissing = htsgo.TypeMissing
	TypeInt8    = htsgo.TypeInt8
	TypeInt16   = htsgo.TypeInt16
	TypeInt32   = htsgo.TypeInt32
	TypeInt64   = htsgo.TypeInt64
	TypeFloat   = htsgo.TypeFloat
	TypeChar    = htsgo.TypeChar

	MissingInt8    = htsgo.MissingInt8
	MissingInt16   = htsgo.MissingInt16
	MissingInt32   = htsgo.MissingInt32
	MissingFloat   = htsgo.MissingFloat
	MissingFloat32 = htsgo.MissingFloat32

	EndOfVectorInt8  = htsgo.EndOfVectorInt8
	EndOfVectorInt16 = htsgo.EndOfVectorInt16
	EndOfVectorInt32 = htsgo.EndOfVectorInt32
	EndOfVectorFloat = htsgo.EndOfVectorFloat

	DictContig  = htsgo.DictContig
	DictTagInfo = htsgo.DictTagInfo
	DictTagFmt  = htsgo.DictTagFmt
)

var (
	Magic        = htsgo.Magic
	ErrBadMagic  = htsgo.ErrBadMagic
	ErrTruncated = htsgo.ErrTruncated

	NewReader              = htsgo.NewReader
	NewReaderWithHeader    = htsgo.NewReaderWithHeader
	NewWriter              = htsgo.NewWriter
	NewWriterFromVCFHeader = htsgo.NewWriterFromVCFHeader
	ReadHeader             = htsgo.ReadHeader

	DecodeTyped       = htsgo.DecodeTyped
	DecodeFormatTyped = htsgo.DecodeFormatTyped
	DecodeTypedInt    = htsgo.DecodeTypedInt
	DecodeTypedInts   = htsgo.DecodeTypedInts
	DecodeTypedString = htsgo.DecodeTypedString

	EncodeMissing       = htsgo.EncodeMissing
	EncodeTypedFloat    = htsgo.EncodeTypedFloat
	EncodeTypedFloatVec = htsgo.EncodeTypedFloatVec
	EncodeTypedInt8     = htsgo.EncodeTypedInt8
	EncodeTypedInt16    = htsgo.EncodeTypedInt16
	EncodeTypedInt32    = htsgo.EncodeTypedInt32
	EncodeTypedInt32Vec = htsgo.EncodeTypedInt32Vec
	EncodeTypedString   = htsgo.EncodeTypedString

	IsMissingFloat     = htsgo.IsMissingFloat
	IsEndOfVectorFloat = htsgo.IsEndOfVectorFloat
)
