package samtools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// writeSAMFile writes SAM text to a temp file and returns its path.
func writeSAMFile(t *testing.T, text string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(p, []byte(text), 0o600); err != nil {
		t.Fatalf("write sam: %v", err)
	}
	return p
}

// toReaderSlice wraps a single reader in the []io.Reader slice Depth/DepthFile
// expect.
func toReaderSlice(r io.Reader) []io.Reader { return []io.Reader{r} }

// TestDepthFileIndexedMatchesStreaming pins the core guarantee of the indexed
// region scan: DepthFile (which seeks to the .bai chunks) must emit output
// byte-identical to the linear streaming Depth for every region query, since
// the index only changes which BGZF blocks are inflated, never which positions
// are reported.
func TestDepthFileIndexedMatchesStreaming(t *testing.T) {
	bam := makeIndexedBAM(t, writeSAMFile(t, depthSAM))

	cases := []struct {
		name string
		opts DepthOptions
	}{
		{"whole-ref-a", DepthOptions{AllPositions: true, Regions: []string{"chr1"}, ExcludeFlags: 0x4}},
		{"whole-ref-noa", DepthOptions{Regions: []string{"chr1"}, ExcludeFlags: 0x4}},
		{"sub-region-a", DepthOptions{AllPositions: true, Regions: []string{"chr1:10-20"}, ExcludeFlags: 0x4}},
		{"sub-region-noa", DepthOptions{Regions: []string{"chr1:10-20"}, ExcludeFlags: 0x4}},
		{"second-ref-a", DepthOptions{AllPositions: true, Regions: []string{"chr2"}, ExcludeFlags: 0x4}},
		{"empty-region-a", DepthOptions{AllPositions: true, Regions: []string{"chr1:480-490"}, ExcludeFlags: 0x4}},
		{"minmapq-region", DepthOptions{AllPositions: true, MinMAPQ: 1, Regions: []string{"chr1:195-210"}, ExcludeFlags: 0x4}},
		{"minbaseq-region", DepthOptions{AllPositions: true, MinBaseQ: 20, Regions: []string{"chr1:395-410"}, ExcludeFlags: 0x4}},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			// Reference output via the linear streaming reader.
			fStream, err := os.Open(bam)
			if err != nil {
				t.Fatal(err)
			}
			var want bytes.Buffer
			if err := Depth(toReaderSlice(fStream), &want, tc.opts); err != nil {
				t.Fatalf("Depth (stream): %v", err)
			}
			_ = fStream.Close()

			// Indexed output via DepthFile.
			fIdx, err := os.Open(bam)
			if err != nil {
				t.Fatal(err)
			}
			var got bytes.Buffer
			if err := DepthFile([]string{bam}, toReaderSlice(fIdx), &got, tc.opts); err != nil {
				t.Fatalf("DepthFile (indexed): %v", err)
			}
			_ = fIdx.Close()

			if got.String() != want.String() {
				t.Errorf("indexed output differs from streaming output\n--- indexed ---\n%s\n--- streaming ---\n%s",
					got.String(), want.String())
			}
		})
	}
}

// TestDepthFileFallsBackWithoutIndex verifies DepthFile produces identical
// output to streaming Depth when the input has no sibling index (the indexed
// fast path is not taken and the linear scan is used instead).
func TestDepthFileFallsBackWithoutIndex(t *testing.T) {
	bam := makeIndexedBAM(t, writeSAMFile(t, depthSAM))
	// Drop the sibling index so DepthFile must fall back to streaming.
	if err := os.Remove(bam + ".bai"); err != nil {
		t.Fatalf("remove bai: %v", err)
	}

	opts := DepthOptions{AllPositions: true, Regions: []string{"chr1:10-20"}, ExcludeFlags: 0x4}

	fStream, err := os.Open(bam)
	if err != nil {
		t.Fatal(err)
	}
	var want bytes.Buffer
	if err := Depth(toReaderSlice(fStream), &want, opts); err != nil {
		t.Fatalf("Depth: %v", err)
	}
	_ = fStream.Close()

	fIdx, err := os.Open(bam)
	if err != nil {
		t.Fatal(err)
	}
	var got bytes.Buffer
	if err := DepthFile([]string{bam}, toReaderSlice(fIdx), &got, opts); err != nil {
		t.Fatalf("DepthFile: %v", err)
	}
	_ = fIdx.Close()

	if got.String() != want.String() {
		t.Errorf("fallback output differs from streaming\n got: %q\nwant: %q", got.String(), want.String())
	}
}

// TestBAMReadDepthIntoMatchesFullDecode checks that the depth-tailored BAM
// decode populates exactly the fields depth reads (RName, Pos, Flag, MapQ,
// CIGAR, and QUAL when requested) identically to the full record decode, and
// clears the fields depth never reads.
func TestBAMReadDepthIntoMatchesFullDecode(t *testing.T) {
	bam := makeIndexedBAM(t, writeSAMFile(t, depthSAM))

	// Full-decode reference records.
	fFull, err := os.Open(bam)
	if err != nil {
		t.Fatal(err)
	}
	defer fFull.Close()
	full, err := sam.NewBAMReader(fFull)
	if err != nil {
		t.Fatal(err)
	}
	var fullRecs []sam.Record
	for {
		rec, err := full.Read()
		if err != nil {
			break
		}
		fullRecs = append(fullRecs, *rec)
	}

	// Depth-tailored decode (needQual=true to exercise the QUAL branch).
	fDepth, err := os.Open(bam)
	if err != nil {
		t.Fatal(err)
	}
	defer fDepth.Close()
	dr, err := sam.NewBAMReader(fDepth)
	if err != nil {
		t.Fatal(err)
	}
	var rec sam.Record
	i := 0
	for {
		if err := dr.ReadDepthInto(&rec, true); err != nil {
			break
		}
		if i >= len(fullRecs) {
			t.Fatalf("depth decode yielded more records than full decode")
		}
		ref := fullRecs[i]
		if rec.RName != ref.RName || rec.Pos != ref.Pos || rec.Flag != ref.Flag || rec.MapQ != ref.MapQ {
			t.Errorf("record %d prefix mismatch: got {%s %d %x %d} want {%s %d %x %d}",
				i, rec.RName, rec.Pos, rec.Flag, rec.MapQ, ref.RName, ref.Pos, ref.Flag, ref.MapQ)
		}
		if len(rec.Cigar) != len(ref.Cigar) {
			t.Errorf("record %d CIGAR len: got %d want %d", i, len(rec.Cigar), len(ref.Cigar))
		} else {
			for k := range rec.Cigar {
				if rec.Cigar[k] != ref.Cigar[k] {
					t.Errorf("record %d CIGAR[%d]: got %v want %v", i, k, rec.Cigar[k], ref.Cigar[k])
				}
			}
		}
		if !bytes.Equal(rec.Qual, ref.Qual) {
			t.Errorf("record %d QUAL: got %v want %v", i, rec.Qual, ref.Qual)
		}
		// Fields depth never reads must be cleared on the reused record.
		if rec.QName != "" {
			t.Errorf("record %d QName not cleared: %q", i, rec.QName)
		}
		if rec.Seq != "" {
			t.Errorf("record %d Seq not cleared: %q", i, rec.Seq)
		}
		if len(rec.Aux) != 0 {
			t.Errorf("record %d Aux not cleared: %v", i, rec.Aux)
		}
		i++
	}
	if i != len(fullRecs) {
		t.Errorf("depth decode yielded %d records, full decode %d", i, len(fullRecs))
	}
}

// TestBAMReadDepthIntoNoQualSkipsQual confirms that with needQual=false the
// depth decode leaves QUAL empty (the SEQ/QUAL block is skipped entirely).
func TestBAMReadDepthIntoNoQualSkipsQual(t *testing.T) {
	bam := makeIndexedBAM(t, writeSAMFile(t, depthSAM))
	f, err := os.Open(bam)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	dr, err := sam.NewBAMReader(f)
	if err != nil {
		t.Fatal(err)
	}
	var rec sam.Record
	for {
		if err := dr.ReadDepthInto(&rec, false); err != nil {
			break
		}
		if len(rec.Qual) != 0 {
			t.Fatalf("needQual=false should leave QUAL empty, got %v", rec.Qual)
		}
		if len(rec.Cigar) == 0 && rec.Flag&sam.FlagUnmapped == 0 {
			t.Fatalf("mapped record decoded with empty CIGAR")
		}
	}
}
