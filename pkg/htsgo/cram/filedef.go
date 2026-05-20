package cram

import (
	"fmt"
	"io"
)

// fileDefMagic is the 4-byte signature at the start of every CRAM file.
var fileDefMagic = [4]byte{'C', 'R', 'A', 'M'}

// fileDefSize is the on-disk size of a CRAM file definition: 4-byte
// magic, 1-byte major version, 1-byte minor version and a 20-byte file
// identifier.
const fileDefSize = 26

// FileDefinition is the fixed-size header that opens every CRAM file. It
// records the format version and a free-form file identifier.
type FileDefinition struct {
	// Major is the CRAM major version number (2 or 3 in practice).
	Major uint8
	// Minor is the CRAM minor version number (0 or 1 in practice).
	Minor uint8
	// FileID is the 20-byte file identifier; it is conventionally the
	// originating file name, NUL-padded. Trailing NULs are not trimmed.
	FileID [20]byte
}

// supported reports whether the file definition is a CRAM version this
// package can parse. CRAM major version 3 (v3.0 and v3.1) is the primary
// target; major version 2 is recognised structurally but its containers
// omit the CRC32 fields, which changes parsing.
func (d FileDefinition) supported() bool {
	return d.Major == 2 || d.Major == 3
}

// hasCRC reports whether containers and blocks in this CRAM version carry
// the trailing 4-byte CRC32. CRAM v3.0 and later embed it; v2 does not.
func (d FileDefinition) hasCRC() bool {
	return d.Major >= 3
}

// VersionString returns the version as a "major.minor" string.
func (d FileDefinition) VersionString() string {
	return fmt.Sprintf("%d.%d", d.Major, d.Minor)
}

// FileIDString returns the file identifier with trailing NUL padding
// removed, as a string.
func (d FileDefinition) FileIDString() string {
	id := d.FileID[:]
	for len(id) > 0 && id[len(id)-1] == 0 {
		id = id[:len(id)-1]
	}
	return string(id)
}

// readFileDefinition reads and validates the 26-byte CRAM file
// definition from r. It returns an error if the magic bytes are not
// "CRAM" or the major version is one this package cannot parse.
func readFileDefinition(r io.Reader) (FileDefinition, error) {
	var buf [fileDefSize]byte
	if _, err := io.ReadFull(r, buf[:]); err != nil {
		if err == io.EOF || err == io.ErrUnexpectedEOF {
			return FileDefinition{}, fmt.Errorf("cram: truncated file definition (need %d bytes): %w",
				fileDefSize, io.ErrUnexpectedEOF)
		}
		return FileDefinition{}, fmt.Errorf("cram: reading file definition: %w", err)
	}
	var d FileDefinition
	if [4]byte{buf[0], buf[1], buf[2], buf[3]} != fileDefMagic {
		return FileDefinition{}, fmt.Errorf("cram: not a CRAM file: magic is %q, want %q",
			buf[0:4], fileDefMagic[:])
	}
	d.Major = buf[4]
	d.Minor = buf[5]
	copy(d.FileID[:], buf[6:26])
	if !d.supported() {
		return d, fmt.Errorf("cram: unsupported CRAM major version %d (this parser supports 2 and 3)",
			d.Major)
	}
	return d, nil
}
