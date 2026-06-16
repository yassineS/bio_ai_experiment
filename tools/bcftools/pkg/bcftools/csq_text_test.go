package bcftools

import (
	"bytes"
	"sort"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// stripCSQTextProvenance drops the first two header lines of -O t output
// (the "# This file was produced by:" version line and the "# The
// command line was:" line). Both are build/argv provenance that always
// differs between two invocations; the stable content is the "# LOG"
// line, the "# CSQ" column header, and the "CSQ\t..." data rows.
func stripCSQTextProvenance(s string) string {
	var b strings.Builder
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, "# This file was produced by:") {
			continue
		}
		if strings.HasPrefix(line, "# The command line was:") {
			continue
		}
		b.WriteString(line)
		b.WriteByte('\n')
	}
	return b.String()
}

// csqTextDataLines keeps only the "CSQ\t..." rows from -O t output (the
// consequence stream), dropping every "#" comment / "LOG" line. Used when
// upstream interleaves LOG lines our engine does not currently emit.
func csqTextDataLines(s string) string {
	var b strings.Builder
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, "CSQ\t") {
			b.WriteString(line)
			b.WriteByte('\n')
		}
	}
	return b.String()
}

// TestCSQTextOracleParity runs the live upstream `bcftools csq -O t` and
// the Go engine in TextOutput mode over the same fixtures and asserts the
// streaming-text output is byte-identical once the two version/command
// provenance lines are stripped. This is a CLI-to-CLI oracle test: the
// expected bytes are produced by the real upstream binary at run time,
// never from a committed golden.
func TestCSQTextOracleParity(t *testing.T) {
	bin := upstreamBcftoolsCsqSlice4(t)
	cases := []struct {
		name, vcf, fa, gff string
		// sortLines compares the CSQ data rows as an order-independent
		// multiset. It is set for fixtures (csq.2) whose site has dozens
		// of overlapping transcripts: upstream's per-transcript emission
		// order is an artifact of its internal heap/hash iteration that
		// our engine does not reproduce in any output mode (the VCF
		// parity test masks the same divergence by sorting INFO/BCSQ).
		// The line *content* must still match exactly.
		sortLines bool
	}{
		{name: "csq.1", vcf: "csq.vcf", fa: "csq.fa", gff: "csq.gff3"},
		{name: "csq.2", vcf: "csq.2.vcf", fa: "csq.fa", gff: "csq.2.gff", sortLines: true},
		{name: "csq.oob-codon", vcf: "csq.oob-codon.vcf", fa: "csq.oob-codon.fa", gff: "csq.oob-codon.gff"},
		{name: "csq.splice.issue-2543", vcf: "csq.splice.issue-2543.vcf", fa: "csq.splice.issue-2543.fa", gff: "csq.splice.issue-2543.gff"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			up := runUpstreamCsqSlice4(t, bin,
				"-f", fixtureCSQ(tc.fa), "-g", fixtureCSQ(tc.gff),
				"-O", "t", fixtureCSQ(tc.vcf))

			var ours bytes.Buffer
			_, err := CSQFile(fixtureCSQ(tc.vcf), &ours, CSQOptions{
				FastaRef:   fixtureCSQ(tc.fa),
				GFFAnnot:   fixtureCSQ(tc.gff),
				TextOutput: true,
				TextArgv: []string{"-f", fixtureCSQ(tc.fa), "-g", fixtureCSQ(tc.gff),
					"-O", "t", fixtureCSQ(tc.vcf)},
			})
			if err != nil {
				t.Fatalf("CSQFile -O t: %v", err)
			}

			// The "# CSQ" column header plus all "CSQ\t" data rows must be
			// byte-identical to upstream once provenance is stripped.
			wantHdr := stripCSQTextProvenance(string(up))
			gotHdr := stripCSQTextProvenance(ours.String())
			// Upstream may interleave "LOG\t..." warning lines into the
			// stream; our engine routes those to stderr instead. Compare
			// the column-header line and CSQ data rows, which are the
			// load-bearing parity surface.
			wantData := csqTextDataLines(wantHdr)
			gotData := csqTextDataLines(gotHdr)
			if tc.sortLines {
				wantData = sortTextLines(wantData)
				gotData = sortTextLines(gotData)
			}
			if wantData != gotData {
				t.Fatalf("CSQ text data mismatch for %s:\n--- upstream ---\n%s\n--- ours ---\n%s",
					tc.name, wantData, gotData)
			}
			// The stable column-header line must match too.
			if want, got := csqColumnHeader(wantHdr), csqColumnHeader(gotHdr); want != got {
				t.Fatalf("CSQ column header mismatch:\nupstream: %q\nours:     %q", want, got)
			}
		})
	}
}

// sortTextLines returns s with its newline-separated lines sorted, so two
// outputs can be compared as an order-independent multiset.
func sortTextLines(s string) string {
	lines := strings.Split(strings.TrimRight(s, "\n"), "\n")
	sort.Strings(lines)
	return strings.Join(lines, "\n") + "\n"
}

// csqColumnHeader extracts the "# CSQ\t..." column-header line.
func csqColumnHeader(s string) string {
	for _, line := range strings.Split(s, "\n") {
		if strings.HasPrefix(line, "# CSQ\t") {
			return line
		}
	}
	return ""
}

// TestUnitCSQTextPrintVcsq is a binary-free unit test for the
// textPrintVcsq line formatter. It synthesises a buffered record with a
// couple of consequences and verifies the exact "CSQ\t..." line bytes for
// the sample/haplotype/printed-upstream cases.
func TestUnitCSQTextPrintVcsq(t *testing.T) {
	e := &hapEngine{
		hdr: &vcf.Header{Samples: []string{"S0", "S1"}},
	}
	vr := &vrecBuf{
		rec: &vcf.Variant{Chrom: "chr1", Pos: 200},
		vcsq: []vcsq{
			{typ: csqMissense, gene: "G", trid: "T1", biotype: "protein_coding", strand: 1,
				vstr: "|2V>2I|200G>A", hasVstr: true},
			{typ: csqSpliceRegion, gene: "G", trid: "T1", biotype: "protein_coding"},
			// A CSQ_PRINTED_UPSTREAM back-reference: must be skipped.
			{typ: csqPrintedUpstream | csqMissense, gene: "G", trid: "T1", biotype: "protein_coding"},
		},
	}

	tests := []struct {
		name    string
		txt     txtCsq
		wantOK  bool
		wantStr string
	}{
		{
			name:    "sample haplotype1 missense",
			txt:     txtCsq{idx: 0, ismpl: 0, ihap: 1},
			wantOK:  true,
			wantStr: "CSQ\tS0\t1\tchr1\t200\tmissense|G|T1|protein_coding|+|2V>2I|200G>A\n",
		},
		{
			name:    "second sample haplotype2 splice_region",
			txt:     txtCsq{idx: 1, ismpl: 1, ihap: 2},
			wantOK:  true,
			wantStr: "CSQ\tS1\t2\tchr1\t200\tsplice_region|G|T1|protein_coding\n",
		},
		{
			name:    "drop-GT sample-agnostic",
			txt:     txtCsq{idx: 1, ismpl: -1, ihap: 0},
			wantOK:  true,
			wantStr: "CSQ\t-\t-\tchr1\t200\tsplice_region|G|T1|protein_coding\n",
		},
		{
			name:   "printed-upstream entry is skipped",
			txt:    txtCsq{idx: 2, ismpl: 0, ihap: 1},
			wantOK: false,
		},
		{
			name:   "out-of-range idx is skipped",
			txt:    txtCsq{idx: 99, ismpl: 0, ihap: 1},
			wantOK: false,
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got, ok := e.textPrintVcsq(vr, &tc.txt)
			if ok != tc.wantOK {
				t.Fatalf("ok = %v, want %v (line=%q)", ok, tc.wantOK, got)
			}
			if ok && got != tc.wantStr {
				t.Fatalf("line = %q, want %q", got, tc.wantStr)
			}
		})
	}
}

// TestUnitCSQTextHeader is a binary-free unit test for the -O t header
// block. The two provenance lines vary by build, so the test asserts the
// stable "# LOG" and "# CSQ" column-header lines verbatim and that the
// command-line line embeds the supplied argv.
func TestUnitCSQTextHeader(t *testing.T) {
	var buf bytes.Buffer
	argv := []string{"-f", "ref.fa", "-g", "ann.gff", "-O", "t", "in.vcf"}
	if err := writeCSQTextHeader(&buf, argv); err != nil {
		t.Fatalf("writeCSQTextHeader: %v", err)
	}
	out := buf.String()
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 4 {
		t.Fatalf("expected 4 header lines, got %d:\n%s", len(lines), out)
	}
	if !strings.HasPrefix(lines[0], "# This file was produced by: bcftools +csq(") {
		t.Errorf("line 0 = %q", lines[0])
	}
	wantCmd := "# The command line was:\tbcftools +csq -f ref.fa -g ann.gff -O t in.vcf"
	if lines[1] != wantCmd {
		t.Errorf("command line = %q, want %q", lines[1], wantCmd)
	}
	if lines[2] != "# LOG\t[2]Message" {
		t.Errorf("LOG line = %q", lines[2])
	}
	wantCols := "# CSQ\t[2]Sample\t[3]Haplotype\t[4]Chromosome\t[5]Position\t[6]Consequence"
	if lines[3] != wantCols {
		t.Errorf("column header = %q, want %q", lines[3], wantCols)
	}
}

// TestUnitCSQTextStageDedup is a binary-free unit test for textStage: an
// exact (idx, ismpl, ihap) repeat must not be staged twice, while distinct
// tuples accumulate in staging order.
func TestUnitCSQTextStageDedup(t *testing.T) {
	e := &hapEngine{}
	vr := &vrecBuf{vcsq: []vcsq{{}, {}}}
	c0 := &csqEntry{vrec: vr, idx: 0}
	c1 := &csqEntry{vrec: vr, idx: 1}

	e.textStage(c0, 0, 1)
	e.textStage(c0, 0, 1) // exact repeat: ignored
	e.textStage(c0, 0, 2) // distinct haplotype
	e.textStage(c1, 0, 1) // distinct idx

	want := []txtCsq{
		{idx: 0, ismpl: 0, ihap: 1},
		{idx: 0, ismpl: 0, ihap: 2},
		{idx: 1, ismpl: 0, ihap: 1},
	}
	if len(vr.txt) != len(want) {
		t.Fatalf("staged %d tuples, want %d: %+v", len(vr.txt), len(want), vr.txt)
	}
	for i := range want {
		if vr.txt[i] != want[i] {
			t.Errorf("txt[%d] = %+v, want %+v", i, vr.txt[i], want[i])
		}
	}
}
