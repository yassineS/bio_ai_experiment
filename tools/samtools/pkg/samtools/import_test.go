package samtools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestImport_PairedR1R2 verifies the two-file paired shape:
//   - both records get FPAIRED + FUNMAP + FMUNMAP
//   - first record gets FREAD1, second gets FREAD2
//   - /1 /2 suffix is stripped from QNAME (default)
//   - records are emitted R1, R2, R1, R2... in pair order
func TestImport_PairedR1R2(t *testing.T) {
	var buf bytes.Buffer
	n, err := FastqImportFiles(nil, &buf, FastqImportOptions{
		Read1Path:       parityPath(t, "import/r1.fq"),
		Read2Path:       parityPath(t, "import/r2.fq"),
		StripPairSuffix: true,
		OutputBAM:       false,
	})
	if err != nil {
		t.Fatalf("FastqImport: %v", err)
	}
	if n != 4 {
		t.Errorf("emitted %d records, want 4", n)
	}
	recs := readSAMRecords(t, buf.String())
	if len(recs) != 4 {
		t.Fatalf("parsed %d records, want 4", len(recs))
	}
	// Pair order: pairA/1, pairA/2, pairB/1, pairB/2.
	wantFlag1 := uint16(sam.FlagPaired | sam.FlagUnmapped | sam.FlagMateUnmapped | sam.FlagRead1)
	wantFlag2 := uint16(sam.FlagPaired | sam.FlagUnmapped | sam.FlagMateUnmapped | sam.FlagRead2)
	wantQNames := []string{"pairA", "pairA", "pairB", "pairB"}
	wantFlags := []uint16{wantFlag1, wantFlag2, wantFlag1, wantFlag2}
	for i, rec := range recs {
		if rec.QName != wantQNames[i] {
			t.Errorf("rec %d QName = %q, want %q", i, rec.QName, wantQNames[i])
		}
		if rec.Flag != wantFlags[i] {
			t.Errorf("rec %d Flag = 0x%x, want 0x%x", i, rec.Flag, wantFlags[i])
		}
	}
}

// TestImport_ThreadsByteIdentical verifies -@/--threads only parallelises the
// output BAM's BGZF deflate: because BGZF blocks are independent gzip members,
// the compressed output is byte-for-byte identical for any worker count. Uses
// OutputBAM:true so the threaded writer branch is actually exercised (SAM text
// output would not touch it).
func TestImport_ThreadsByteIdentical(t *testing.T) {
	run := func(threads int) []byte {
		var buf bytes.Buffer
		if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
			Read1Path:       parityPath(t, "import/r1.fq"),
			Read2Path:       parityPath(t, "import/r2.fq"),
			StripPairSuffix: true,
			OutputBAM:       true,
			Threads:         threads,
		}); err != nil {
			t.Fatalf("FastqImport -@%d: %v", threads, err)
		}
		return buf.Bytes()
	}
	if one, many := run(1), run(4); !bytes.Equal(one, many) {
		t.Errorf("import -@1 (%d bytes) vs -@4 (%d bytes) differ", len(one), len(many))
	}
}

// TestImport_SingleUnpaired exercises the -0 unpaired-file shape: each
// record gets just FUNMAP, no pair bits.
func TestImport_SingleUnpaired(t *testing.T) {
	var buf bytes.Buffer
	if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
		UnpairedPath: parityPath(t, "import/single.fq"),
		OutputBAM:    false,
	}); err != nil {
		t.Fatalf("FastqImport: %v", err)
	}
	recs := readSAMRecords(t, buf.String())
	if len(recs) != 2 {
		t.Fatalf("got %d records, want 2", len(recs))
	}
	for i, rec := range recs {
		if rec.Flag != uint16(sam.FlagUnmapped) {
			t.Errorf("rec %d flag = 0x%x, want 0x4 (just unmapped)", i, rec.Flag)
		}
	}
}

// TestImport_InterleavedSingleFile verifies the -s shape: a single file
// with /1 /2 suffixes that's reinterpreted as interleaved paired data.
// Matching upstream (htslib's FASTQ reader, sam.c fastq_parse1), every
// record whose name ends in /1 or /2 gets FUNMAP|FMUNMAP|FPAIRED plus
// FREAD1 or FREAD2 — the FMUNMAP bit is set from the suffix alone, exactly
// as the live `samtools import -s` binary does (flags 77 / 141). A prior
// version of this test asserted FMUNMAP was absent, which diverged from
// upstream; that was a port bug, now fixed in fastq2bam.go.
func TestImport_InterleavedSingleFile(t *testing.T) {
	var buf bytes.Buffer
	if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
		SinglePath:      parityPath(t, "import/interleaved.fq"),
		StripPairSuffix: true,
		OutputBAM:       false,
	}); err != nil {
		t.Fatalf("FastqImport: %v", err)
	}
	recs := readSAMRecords(t, buf.String())
	if len(recs) != 4 {
		t.Fatalf("got %d records, want 4", len(recs))
	}
	r1Flag := uint16(sam.FlagUnmapped | sam.FlagMateUnmapped | sam.FlagPaired | sam.FlagRead1)
	r2Flag := uint16(sam.FlagUnmapped | sam.FlagMateUnmapped | sam.FlagPaired | sam.FlagRead2)
	wantFlags := []uint16{r1Flag, r2Flag, r1Flag, r2Flag}
	for i, rec := range recs {
		if rec.Flag != wantFlags[i] {
			t.Errorf("rec %d flag = 0x%x, want 0x%x", i, rec.Flag, wantFlags[i])
		}
	}
}

// TestImport_AuxTagsExtraction verifies -T "*" captures every aux field
// from the FASTQ description, and -T "XZ" captures just that one.
func TestImport_AuxTagsExtraction(t *testing.T) {
	t.Run("star", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
			SinglePath: parityPath(t, "import/aux.fq"),
			AuxTags:    "*",
			OutputBAM:  false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		recs := readSAMRecords(t, buf.String())
		if len(recs) != 2 {
			t.Fatalf("got %d records, want 2", len(recs))
		}
		// readX: XX:i:10  XY:i:-257  XZ:Z:foo
		xx, ok := recs[0].GetAux("XX")
		if !ok {
			t.Errorf("readX missing XX")
		} else if v, _ := xx.Int(); v != 10 {
			t.Errorf("readX XX = %d, want 10", v)
		}
		xz, ok := recs[0].GetAux("XZ")
		if !ok {
			t.Errorf("readX missing XZ")
		} else if v, _ := xz.String(); v != "foo" {
			t.Errorf("readX XZ = %q, want foo", v)
		}
		// readY: AA:Z:  (empty)  XZ:Z:bar
		xz2, ok := recs[1].GetAux("XZ")
		if !ok {
			t.Errorf("readY missing XZ")
		} else if v, _ := xz2.String(); v != "bar" {
			t.Errorf("readY XZ = %q, want bar", v)
		}
	})
	t.Run("subset", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
			SinglePath: parityPath(t, "import/aux.fq"),
			AuxTags:    "XZ",
			OutputBAM:  false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		recs := readSAMRecords(t, buf.String())
		if _, ok := recs[0].GetAux("XX"); ok {
			t.Errorf("XX leaked through despite -T XZ")
		}
		if _, ok := recs[0].GetAux("XZ"); !ok {
			t.Errorf("XZ missing despite -T XZ")
		}
	})
	t.Run("empty_disables", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
			SinglePath: parityPath(t, "import/aux.fq"),
			AuxTags:    "",
			OutputBAM:  false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		recs := readSAMRecords(t, buf.String())
		for _, rec := range recs {
			if _, ok := rec.GetAux("XX"); ok {
				t.Errorf("XX present despite -T \"\"")
			}
			if _, ok := rec.GetAux("XZ"); ok {
				t.Errorf("XZ present despite -T \"\"")
			}
		}
	})
}

// TestImport_ReadGroup verifies both forms: -R "id" and -r "ID:id\tSM:foo".
// Both produce an @RG header line and an RG:Z aux on every record.
func TestImport_ReadGroup(t *testing.T) {
	t.Run("short_R", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
			SinglePath: parityPath(t, "import/single.fq"),
			ReadGroup:  "rgid",
			OutputBAM:  false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		if !strings.Contains(buf.String(), "@RG\tID:rgid\n") {
			t.Errorf("missing @RG header line, got: %q", headerOf(buf.String()))
		}
		recs := readSAMRecords(t, buf.String())
		for _, rec := range recs {
			rg, ok := rec.GetAux("RG")
			if !ok {
				t.Errorf("record %s missing RG", rec.QName)
				continue
			}
			s, _ := rg.String()
			if s != "rgid" {
				t.Errorf("RG = %q, want rgid", s)
			}
		}
	})
	t.Run("long_r", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
			SinglePath:    parityPath(t, "import/single.fq"),
			ReadGroupLine: "ID:rg42\tSM:sampleX",
			OutputBAM:     false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		if !strings.Contains(buf.String(), "@RG\tID:rg42\tSM:sampleX\n") {
			t.Errorf("missing fully-fleshed @RG line in header: %q", headerOf(buf.String()))
		}
		recs := readSAMRecords(t, buf.String())
		rg, ok := recs[0].GetAux("RG")
		if !ok {
			t.Fatalf("RG aux missing")
		}
		s, _ := rg.String()
		if s != "rg42" {
			t.Errorf("RG ID = %q, want rg42", s)
		}
	})
}

// TestImport_OrderTag verifies the --order TAG forms: bare TAG → i:N
// counter, TAG:WIDTH → zero-padded Z string.
func TestImport_OrderTag(t *testing.T) {
	t.Run("int", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
			SinglePath: parityPath(t, "import/single.fq"),
			OrderTag:   "oi",
			OutputBAM:  false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		recs := readSAMRecords(t, buf.String())
		for i, rec := range recs {
			oi, ok := rec.GetAux("oi")
			if !ok {
				t.Errorf("rec %d missing oi tag", i)
				continue
			}
			v, _ := oi.Int()
			if v != int64(i) {
				t.Errorf("rec %d oi = %d, want %d", i, v, i)
			}
		}
	})
	t.Run("padded", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
			SinglePath: parityPath(t, "import/single.fq"),
			OrderTag:   "oi:5",
			OutputBAM:  false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		recs := readSAMRecords(t, buf.String())
		oi, ok := recs[0].GetAux("oi")
		if !ok {
			t.Fatalf("rec 0 missing oi")
		}
		s, _ := oi.String()
		if s != "00000" {
			t.Errorf("oi[0] = %q, want 00000", s)
		}
		oi1, _ := recs[1].GetAux("oi")
		s1, _ := oi1.String()
		if s1 != "00001" {
			t.Errorf("oi[1] = %q, want 00001", s1)
		}
	})
}

// TestImport_KeepPairSuffix flips StripPairSuffix off and verifies /1 /2
// are preserved (upstream -N flag).
func TestImport_KeepPairSuffix(t *testing.T) {
	var buf bytes.Buffer
	if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
		Read1Path:       parityPath(t, "import/r1.fq"),
		Read2Path:       parityPath(t, "import/r2.fq"),
		StripPairSuffix: false,
		OutputBAM:       false,
	}); err != nil {
		t.Fatalf("FastqImport: %v", err)
	}
	recs := readSAMRecords(t, buf.String())
	wantQNames := []string{"pairA/1", "pairA/2", "pairB/1", "pairB/2"}
	for i, rec := range recs {
		if rec.QName != wantQNames[i] {
			t.Errorf("rec %d QName = %q, want %q", i, rec.QName, wantQNames[i])
		}
	}
}

// TestImport_BAMRoundtrip writes BAM and reads it back, verifying paired
// flags survive the binary trip.
func TestImport_BAMRoundtrip(t *testing.T) {
	var buf bytes.Buffer
	if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
		Read1Path:       parityPath(t, "import/r1.fq"),
		Read2Path:       parityPath(t, "import/r2.fq"),
		StripPairSuffix: true,
		OutputBAM:       true,
	}); err != nil {
		t.Fatalf("FastqImport BAM: %v", err)
	}
	r, err := sam.NewReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("re-open BAM: %v", err)
	}
	n := 0
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		if !rec.IsUnmapped() {
			t.Errorf("BAM rec %d not unmapped: flag=0x%x", n, rec.Flag)
		}
		if !rec.IsPaired() {
			t.Errorf("BAM rec %d not paired: flag=0x%x", n, rec.Flag)
		}
		n++
	}
	if n != 4 {
		t.Errorf("BAM emitted %d records, want 4", n)
	}
}

// TestImport_PositionalArgs verifies the positional-arg paths: one file →
// single, two files → R1+R2.
func TestImport_PositionalArgs(t *testing.T) {
	t.Run("one_positional", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles([]string{parityPath(t, "import/single.fq")}, &buf, FastqImportOptions{
			OutputBAM: false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		if !strings.Contains(buf.String(), "read1\t4\t") {
			t.Errorf("missing read1 in output: %q", buf.String())
		}
	})
	t.Run("two_positional", func(t *testing.T) {
		var buf bytes.Buffer
		if _, err := FastqImportFiles([]string{
			parityPath(t, "import/r1.fq"),
			parityPath(t, "import/r2.fq"),
		}, &buf, FastqImportOptions{
			StripPairSuffix: true,
			OutputBAM:       false,
		}); err != nil {
			t.Fatalf("FastqImport: %v", err)
		}
		recs := readSAMRecords(t, buf.String())
		if len(recs) != 4 {
			t.Errorf("got %d records from two-positional input, want 4", len(recs))
		}
	})
}

// TestImport_HeaderHasHDLine verifies the @HD line is well-formed:
// VN:1.6 SO:unsorted GO:query.
func TestImport_HeaderHasHDLine(t *testing.T) {
	var buf bytes.Buffer
	if _, err := FastqImportFiles(nil, &buf, FastqImportOptions{
		SinglePath: parityPath(t, "import/single.fq"),
		OutputBAM:  false,
	}); err != nil {
		t.Fatalf("FastqImport: %v", err)
	}
	if !strings.Contains(buf.String(), "@HD\tVN:1.6\tSO:unsorted\tGO:query\n") {
		t.Errorf("missing or malformed @HD line: %q", headerOf(buf.String()))
	}
}

// TestImport_MismatchedPairLengths verifies that an unequal number of
// records in R1 / R2 returns an error mentioning both filenames. The
// fixtures are written into a per-test temp dir so the geometry (R1 has
// three records, R2 has two) is exactly controlled.
func TestImport_MismatchedPairLengths(t *testing.T) {
	dir := t.TempDir()
	r1Path := filepath.Join(dir, "r1.fq")
	r2Path := filepath.Join(dir, "r2.fq")
	r1 := "@a/1\nACGT\n+\nIIII\n@b/1\nACGT\n+\nIIII\n@c/1\nACGT\n+\nIIII\n"
	r2 := "@a/2\nACGT\n+\nIIII\n@b/2\nACGT\n+\nIIII\n"
	if err := os.WriteFile(r1Path, []byte(r1), 0o644); err != nil {
		t.Fatalf("write r1: %v", err)
	}
	if err := os.WriteFile(r2Path, []byte(r2), 0o644); err != nil {
		t.Fatalf("write r2: %v", err)
	}

	var buf bytes.Buffer
	_, err := FastqImportFiles(nil, &buf, FastqImportOptions{
		Read1Path: r1Path,
		Read2Path: r2Path,
		OutputBAM: false,
	})
	if err == nil {
		t.Fatalf("expected an error for mismatched pair lengths, got nil")
	}
	msg := err.Error()
	if !strings.Contains(msg, r1Path) || !strings.Contains(msg, r2Path) {
		t.Errorf("error %q does not mention both filenames %q and %q", msg, r1Path, r2Path)
	}
	if !strings.Contains(msg, "differing number of records") {
		t.Errorf("error %q missing 'differing number of records' detail", msg)
	}
}

// TestParity_Import_UpstreamCorpus is the LIVE parity gate for `samtools
// import`. For each input shape it runs the upstream binary with `-O sam`
// (so the comparison is plain SAM text, sidestepping the BGZF byte-identity
// problem entirely) and our FastqImportFiles with OutputBAM=false, then
// asserts the two SAM streams are equal MODULO the @PG line. The @PG line
// is the only legitimately non-deterministic difference: it embeds the
// upstream version string and the verbatim argv (with absolute fixture
// paths), which our port deliberately does not inject. Every other byte —
// the @HD line, the @CO "Reverse with:" line, the @RG line (where present),
// and all record fields including FLAG, SEQ, QUAL and aux tags — must match
// the upstream binary exactly. Per the project rules the upstream binary is
// built on demand and a build failure is fatal, never a skip.
func TestParity_Import_UpstreamCorpus(t *testing.T) {
	bin := upstreamSamtools(t)

	cases := []struct {
		name string
		// upArgs are the flags passed to the upstream binary, in addition
		// to `import -O sam`. The trailing -o is appended by the runner.
		upArgs []string
		// opts mirrors upArgs for the Go FastqImportFiles entry point.
		opts FastqImportOptions
	}{
		{
			name:   "unpaired_-0",
			upArgs: []string{"-0", parityPath(t, "import/single.fq")},
			opts:   FastqImportOptions{UnpairedPath: parityPath(t, "import/single.fq"), StripPairSuffix: true},
		},
		{
			name:   "paired_-1_-2",
			upArgs: []string{"-1", parityPath(t, "import/r1.fq"), "-2", parityPath(t, "import/r2.fq")},
			opts: FastqImportOptions{
				Read1Path: parityPath(t, "import/r1.fq"), Read2Path: parityPath(t, "import/r2.fq"),
				StripPairSuffix: true,
			},
		},
		{
			name:   "interleaved_-s",
			upArgs: []string{"-s", parityPath(t, "import/interleaved.fq")},
			opts:   FastqImportOptions{SinglePath: parityPath(t, "import/interleaved.fq"), StripPairSuffix: true},
		},
		{
			name:   "aux_star_-T",
			upArgs: []string{"-T", "*", "-0", parityPath(t, "import/aux.fq")},
			opts: FastqImportOptions{
				UnpairedPath: parityPath(t, "import/aux.fq"), AuxTags: "*", StripPairSuffix: true,
			},
		},
		{
			name:   "readgroup_-R",
			upArgs: []string{"-R", "rgid", "-0", parityPath(t, "import/single.fq")},
			opts: FastqImportOptions{
				UnpairedPath: parityPath(t, "import/single.fq"), ReadGroup: "rgid", StripPairSuffix: true,
			},
		},
		{
			name:   "order_int_--order",
			upArgs: []string{"--order", "oi", "-0", parityPath(t, "import/single.fq")},
			opts: FastqImportOptions{
				UnpairedPath: parityPath(t, "import/single.fq"), OrderTag: "oi", StripPairSuffix: true,
			},
		},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// Upstream → a temp .sam file (CLI infers SAM from -O sam).
			upPath := filepath.Join(t.TempDir(), "up.sam")
			args := append([]string{"import", "-O", "sam"}, c.upArgs...)
			args = append(args, "-o", upPath)
			cmd := exec.Command(bin, args...)
			var stderr bytes.Buffer
			cmd.Stderr = &stderr
			if err := cmd.Run(); err != nil {
				t.Fatalf("upstream samtools %v: %v\n%s", args, err, stderr.String())
			}
			upSAM, err := os.ReadFile(upPath)
			if err != nil {
				t.Fatalf("read upstream SAM: %v", err)
			}

			// Our port → an in-memory SAM stream.
			var got bytes.Buffer
			if _, err := FastqImportFiles(nil, &got, c.opts); err != nil {
				t.Fatalf("FastqImportFiles: %v", err)
			}

			// dropPGLines (subfeatures_live_oracle_test.go) strips the @PG
			// header line — the only legitimately non-deterministic
			// difference (upstream version string + verbatim argv).
			want := dropPGLines(upSAM)
			have := dropPGLines(got.Bytes())
			if have != want {
				t.Fatalf("import %s differs from upstream (@PG excluded):\n--- ours ---\n%s\n--- upstream ---\n%s",
					c.name, have, want)
			}
		})
	}
}

// readSAMRecords parses a SAM-text blob into []*sam.Record. Header lines
// are silently skipped.
func readSAMRecords(t *testing.T, text string) []*sam.Record {
	t.Helper()
	r, err := sam.NewSAMReader(strings.NewReader(text))
	if err != nil {
		t.Fatalf("NewSAMReader: %v", err)
	}
	var out []*sam.Record
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		out = append(out, rec)
	}
	return out
}

// headerOf returns just the leading @-prefixed header lines.
func headerOf(text string) string {
	var sb strings.Builder
	for _, line := range strings.SplitAfter(text, "\n") {
		if strings.HasPrefix(line, "@") {
			sb.WriteString(line)
		} else {
			break
		}
	}
	return sb.String()
}
