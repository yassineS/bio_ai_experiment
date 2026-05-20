package cram

// CRAM v3.0 read-feature codes. A mapped record stores its differences
// from the reference (or, in a reference-free file, its whole sequence)
// as a list of read features, each tagged with one of these one-byte
// codes. The codes are the ASCII letters the CRAM specification assigns
// in its "read feature records" table.
const (
	featBases        = 'b' // a stretch of bases (BB data series).
	featScores       = 'q' // a stretch of quality scores (QQ data series).
	featBase         = 'B' // a single base with its quality (BA + QS).
	featSubst        = 'X' // a base substitution (BS substitution code).
	featInsertion    = 'I' // an insertion of several bases (IN).
	featDeletion     = 'D' // a deletion of reference bases (DL length).
	featInsertBase   = 'i' // an insertion of a single base (BA).
	featQualityScore = 'Q' // a single quality score (QS).
	featRefSkip      = 'N' // a reference skip (RS length).
	featSoftClip     = 'S' // a soft clip of several bases (SC).
	featPadding      = 'P' // padding (PD length).
	featHardClip     = 'H' // a hard clip (HC length).
)

// readFeature is one decoded CRAM read feature: its one-byte code, the
// in-read position it applies at (1-based, relative to the start of the
// read), and the code-specific payload. Only the field relevant to the
// code is populated.
type readFeature struct {
	// code is one of the feat* constants.
	code byte
	// pos is the 1-based in-read position the feature starts at.
	pos int32

	// bases holds the inserted/soft-clipped/stretch bytes for the
	// multi-byte features (b, I, S) and the QQ stretch.
	bases []byte
	// base is the single base of a B or i feature.
	base byte
	// quality is the single quality score of a B or Q feature.
	quality byte
	// substCode is the substitution code of an X feature, indexing the
	// preservation map's substitution matrix.
	substCode byte
	// length is the operation length of a D, N, P or H feature.
	length int32
}

// decodeFeatures decodes one record's list of nFeatures read features.
// It threads the slice's data-series cursors so that successive features
// (and successive records) draw their payloads from the right place in
// each external block. The feature-position series (FP) is stored as a
// delta from the previous feature's position within the same record, so
// the running position is accumulated here.
func (rd *recordDecoder) decodeFeatures(nFeatures int32) ([]readFeature, error) {
	if nFeatures < 0 {
		return nil, errFormat("record declares a negative read-feature count %d", nFeatures)
	}
	feats := make([]readFeature, 0)
	var prevPos int32
	prev := rd.src.s.consumed()
	for i := int32(0); i < nFeatures; i++ {
		code, err := rd.byteSeries("FC")
		if err != nil {
			return nil, wrapf(err, "read feature %d code", i)
		}
		delta, err := rd.intSeries("FP")
		if err != nil {
			return nil, wrapf(err, "read feature %d position", i)
		}
		// FP is a delta from the previous feature's position; the first
		// feature's delta is relative to zero.
		prevPos += delta
		f := readFeature{code: code, pos: prevPos}
		if err := rd.decodeFeaturePayload(&f); err != nil {
			return nil, wrapf(err, "read feature %d (%c)", i, code)
		}
		feats = append(feats, f)
		// Every feature must consume series input; if one did not, the
		// declared feature count has outrun the data — stop rather than
		// loop to nFeatures (a crafted FN value could be ~2^31).
		c := rd.src.s.consumed()
		if c == prev && i+1 < nFeatures {
			return nil, errFormat("record declares %d read features but the series data is exhausted after feature %d",
				nFeatures, i)
		}
		prev = c
	}
	return feats, nil
}

// decodeFeaturePayload reads the code-specific payload of a single read
// feature, dispatching on its already-decoded code.
func (rd *recordDecoder) decodeFeaturePayload(f *readFeature) error {
	switch f.code {
	case featBases:
		b, err := rd.byteArraySeries("BB")
		if err != nil {
			return err
		}
		f.bases = b
	case featScores:
		b, err := rd.byteArraySeries("QQ")
		if err != nil {
			return err
		}
		f.bases = b
	case featBase:
		base, err := rd.byteSeries("BA")
		if err != nil {
			return err
		}
		q, err := rd.byteSeries("QS")
		if err != nil {
			return err
		}
		f.base, f.quality = base, q
	case featSubst:
		code, err := rd.byteSeries("BS")
		if err != nil {
			return err
		}
		f.substCode = code
	case featInsertion:
		b, err := rd.byteArraySeries("IN")
		if err != nil {
			return err
		}
		f.bases = b
	case featSoftClip:
		b, err := rd.byteArraySeries("SC")
		if err != nil {
			return err
		}
		f.bases = b
	case featInsertBase:
		base, err := rd.byteSeries("BA")
		if err != nil {
			return err
		}
		f.base = base
	case featQualityScore:
		q, err := rd.byteSeries("QS")
		if err != nil {
			return err
		}
		f.quality = q
	case featDeletion:
		n, err := rd.intSeries("DL")
		if err != nil {
			return err
		}
		f.length = n
	case featRefSkip:
		n, err := rd.intSeries("RS")
		if err != nil {
			return err
		}
		f.length = n
	case featPadding:
		n, err := rd.intSeries("PD")
		if err != nil {
			return err
		}
		f.length = n
	case featHardClip:
		n, err := rd.intSeries("HC")
		if err != nil {
			return err
		}
		f.length = n
	default:
		return errFormat("unknown read-feature code %#02x (%q)", f.code, string(f.code))
	}
	return nil
}
