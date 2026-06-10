package samtools

// Live-binary oracle test for POSIX getopt-style short-flag bundling in
// `samtools view`.
//
// The point under test is the new cliflag.Parse wiring in runView: command
// lines that bundle short flags (`-bS` == `-b -S`) and concatenate a value
// (`-q20` == `-q 20`) must now parse and behave exactly like upstream
// samtools' getopt parser. This test proves that two ways by:
//
//  1. Running the genuine upstream `samtools view -bS in.sam` AND our port
//     `samtools view -bS in.sam` on the same fixture, decoding both BAM
//     outputs, and asserting the surviving records match field-for-field.
//  2. Asserting that, within our own binary, the bundled/value-concat forms
//     (`-bS`, `-bSq20`, `-hb`) produce records identical to the canonical
//     spelled-out forms (`-b -S`, `-b -S -q 20`, `-h -b`).
//
// Per the project's testing rules the helpers t.Fatalf rather than t.Skip
// when the binaries cannot be built.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// viewPosixFixtureSAM has a spread of MAPQ values so a `-q20` filter has a
// visible effect (read5 with MAPQ 10 must be dropped, the rest kept).
const viewPosixFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:1000\n" +
	"@SQ\tSN:chr2\tLN:500\n" +
	"read1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
	"read2\t0\tchr1\t200\t40\t5M\t*\t0\t0\tTGCAT\tIIIII\n" +
	"read3\t0\tchr2\t10\t25\t3M2I\t*\t0\t0\tACGTA\tIIIII\n" +
	"read4\t0\tchr1\t300\t30\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
	"read5\t0\tchr1\t400\t10\t5M\t*\t0\t0\tACGTA\tIIIII\n"

// decodeBAMRecords parses a BAM byte slice into a stable, comparable view of
// each record. It is used to compare upstream and our output independently of
// any header provenance differences (@PG lines etc.).
type decodedRec struct {
	QName string
	Flag  uint16
	RName string
	Pos   int32
	MapQ  uint8
	Cigar string
	Seq   string
	Qual  string
}

func decodeBAMRecords(t *testing.T, b []byte) []decodedRec {
	t.Helper()
	br, err := sam.NewBAMReader(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	var recs []decodedRec
	for {
		rec, err := br.Read()
		if err != nil {
			break
		}
		recs = append(recs, decodedRec{
			QName: rec.QName,
			Flag:  rec.Flag,
			RName: rec.RName,
			Pos:   rec.Pos,
			MapQ:  rec.MapQ,
			Cigar: rec.Cigar.String(),
			Seq:   rec.Seq,
			Qual:  string(rec.Qual),
		})
	}
	return recs
}

// TestLiveViewPosixBundling is the upstream-parity gate: `samtools view -bS`
// (bundled) must yield the same records from our port as from upstream.
func TestLiveViewPosixBundling(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)

	dir := t.TempDir()
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, []byte(viewPosixFixtureSAM), 0o644); err != nil {
		t.Fatal(err)
	}

	// The whole point: a bundled `-bS` cluster, exactly as upstream getopt
	// accepts it, must now parse in our port too.
	upstream := decodeBAMRecords(t, runSamtools(t, live, "view", "-bS", samPath))
	mine := decodeBAMRecords(t, runSamtools(t, ours, "view", "-bS", samPath))

	if len(upstream) != len(mine) {
		t.Fatalf("record count: upstream=%d ours=%d", len(upstream), len(mine))
	}
	for i := range upstream {
		if upstream[i] != mine[i] {
			t.Errorf("record %d differs:\nupstream=%+v\nours    =%+v", i, upstream[i], mine[i])
		}
	}
}

// TestLiveViewPosixValueConcat exercises a value-concatenated, mixed cluster
// `-bSq20` against upstream: bool b, bool S, then q taking "20" inline. Only
// records with MAPQ >= 20 survive in both.
func TestLiveViewPosixValueConcat(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)

	dir := t.TempDir()
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, []byte(viewPosixFixtureSAM), 0o644); err != nil {
		t.Fatal(err)
	}

	upstream := decodeBAMRecords(t, runSamtools(t, live, "view", "-bSq20", samPath))
	mine := decodeBAMRecords(t, runSamtools(t, ours, "view", "-bSq20", samPath))

	// Sanity: the q20 filter must actually have dropped read5 (MAPQ 10).
	if len(upstream) != 4 {
		t.Fatalf("expected upstream -q20 to keep 4 records, got %d", len(upstream))
	}
	if len(upstream) != len(mine) {
		t.Fatalf("record count: upstream=%d ours=%d", len(upstream), len(mine))
	}
	for i := range upstream {
		if upstream[i] != mine[i] {
			t.Errorf("record %d differs:\nupstream=%+v\nours    =%+v", i, upstream[i], mine[i])
		}
	}
}

// TestViewPosixBundlingEquivalentToCanonical proves, within our own binary,
// that the bundled / value-concatenated spellings are equivalent to the
// canonical spelled-out forms. This isolates the parser change from any
// upstream behavioural quirks.
func TestViewPosixBundlingEquivalentToCanonical(t *testing.T) {
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, []byte(viewPosixFixtureSAM), 0o644); err != nil {
		t.Fatal(err)
	}

	cases := []struct {
		name      string
		bundled   []string
		canonical []string
	}{
		{"bS == b S", []string{"view", "-bS", samPath}, []string{"view", "-b", "-S", samPath}},
		{"bSq20 == b S q 20", []string{"view", "-bSq20", samPath}, []string{"view", "-b", "-S", "-q", "20", samPath}},
		{"hb == h b", []string{"view", "-hb", samPath}, []string{"view", "-h", "-b", samPath}},
		{"value flag at end -bSq 20", []string{"view", "-bSq", "20", samPath}, []string{"view", "-b", "-S", "-q", "20", samPath}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := runSamtools(t, ours, tc.bundled...)
			want := runSamtools(t, ours, tc.canonical...)
			if !bytes.Equal(got, want) {
				t.Errorf("bundled %v != canonical %v\nbundled records=%+v\ncanonical records=%+v",
					tc.bundled, tc.canonical,
					decodeBAMRecords(t, got), decodeBAMRecords(t, want))
			}
		})
	}
}
