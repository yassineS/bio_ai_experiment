package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// Binary-free unit tests for the stats-wave pure helpers closed in this change:
// the indel-stats PED parser (parseIndelStatsPED), the contrast -f threshold
// parser (parseContrastMaxAC) and its region-wide minor-allele folding, the
// trio-stats -a/--alt-trios deferred singleton/doubleton accounting, and the
// shared -o report writer (statsReportWriter). These run with NO upstream binary
// and with the submodules unpopulated (no exec.Command, no BCFTOOLS_PLUGINS).

// writePED writes a temporary PED file and returns its path.
func writePED(t *testing.T, body string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "trios.ped")
	if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestUnitParseIndelStatsPED(t *testing.T) {
	// Header sample order fixes the indices used by the min-index sort.
	idx := sampleIndex(&vcf.Header{Samples: []string{"CHILD1", "FATHER1", "MOTHER1", "CHILD2", "FATHER2", "MOTHER2"}})

	t.Run("two_trios_sorted_by_min_index", func(t *testing.T) {
		// Listed in reverse index order; the stable min-index sort must reorder so
		// the trio whose minimum sample index is smallest comes first.
		ped := writePED(t, "FAM2\tCHILD2\tFATHER2\tMOTHER2\t1\t0\nFAM1\tCHILD1\tFATHER1\tMOTHER1\t2\t0\n")
		got, err := parseIndelStatsPED(ped, idx)
		if err != nil {
			t.Fatal(err)
		}
		want := []trioStatsTrio{
			{child: 0, father: 1, mother: 2}, // CHILD1 family (min index 0)
			{child: 3, father: 4, mother: 5}, // CHILD2 family (min index 3)
		}
		if !reflect.DeepEqual(got, want) {
			t.Fatalf("got %+v want %+v", got, want)
		}
	})

	t.Run("skips_rows_with_missing_members", func(t *testing.T) {
		// First row references samples absent from the header and is skipped.
		ped := writePED(t, "FAMX\tGHOST\tNOBODY\tNOONE\t1\t0\nFAM1\tCHILD1\tFATHER1\tMOTHER1\t2\t0\n")
		got, err := parseIndelStatsPED(ped, idx)
		if err != nil {
			t.Fatal(err)
		}
		if len(got) != 1 || got[0] != (trioStatsTrio{child: 0, father: 1, mother: 2}) {
			t.Fatalf("got %+v", got)
		}
	})

	t.Run("keeps_duplicate_trios_unlike_trio_stats", func(t *testing.T) {
		// indel-stats.c's parse_ped does NOT deduplicate; the same trio listed
		// twice yields two entries (this is the documented fix-on-port difference
		// from trio-stats.c, which would reject/skip it).
		ped := writePED(t, "FAM1\tCHILD1\tFATHER1\tMOTHER1\t2\t0\nFAM1\tCHILD1\tFATHER1\tMOTHER1\t2\t0\n")
		got, err := parseIndelStatsPED(ped, idx)
		if err != nil {
			t.Fatal(err)
		}
		if len(got) != 2 {
			t.Fatalf("expected 2 (no dedup), got %d: %+v", len(got), got)
		}
	})

	t.Run("error_no_complete_trio", func(t *testing.T) {
		ped := writePED(t, "FAMX\tGHOST\tNOBODY\tNOONE\t1\t0\n")
		if _, err := parseIndelStatsPED(ped, idx); err == nil {
			t.Fatal("expected an error when no trio resolves")
		}
	})

	t.Run("error_too_few_columns", func(t *testing.T) {
		ped := writePED(t, "FAM1\tCHILD1\tFATHER1\n")
		if _, err := parseIndelStatsPED(ped, idx); err == nil {
			t.Fatal("expected an error for <4 columns")
		}
	})
}

func TestUnitParseContrastMaxAC(t *testing.T) {
	cases := []struct {
		s        string
		nsamples int
		want     int
		wantErr  bool
	}{
		{"2", 4, 2, false},     // clean integer, verbatim
		{"0", 4, 0, false},     // integer zero stays zero
		{"5", 4, 5, false},     // integer may exceed sample count
		{"0.5", 4, 2, false},   // 0.5*4 = 2
		{"0.25", 8, 2, false},  // 0.25*8 = 2
		{"0.001", 4, 1, false}, // 0.004 floors to 0 -> bumped to 1
		{"1", 4, 1, false},     // integer 1 (not the float path)
		{"1.0", 4, 4, false},   // float 1.0*4 = 4
		{"0.75", 4, 3, false},  // 0.75*4 = 3
		{"-0.5", 4, 0, true},   // float out of [0,1]
		{"1.5", 4, 0, true},    // float out of [0,1]
		{"abc", 4, 0, true},    // unparseable
	}
	for _, c := range cases {
		c := c
		t.Run(c.s, func(t *testing.T) {
			got, err := parseContrastMaxAC(c.s, c.nsamples)
			if c.wantErr {
				if err == nil {
					t.Fatalf("expected error for %q", c.s)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if got != c.want {
				t.Fatalf("parseContrastMaxAC(%q,%d) = %d, want %d", c.s, c.nsamples, got, c.want)
			}
		})
	}
}

// contrastTestVariant builds a minimal biallelic A>T variant with the given
// per-sample GT strings, in header order ctrl... then case...
func contrastTestVariant(gts []string) *vcf.Variant {
	v := &vcf.Variant{Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"T"}, Format: []string{"GT"}}
	for i, gt := range gts {
		v.Samples = append(v.Samples, vcf.Sample{Name: "S" + string(rune('1'+i)), Data: map[string]string{"GT": gt}})
	}
	return v
}

func TestUnitContrastEnrichmentFolding(t *testing.T) {
	// Samples: S1,S2 control; S3,S4 case. Allele counts per record drive the
	// region-wide enrichNals folding. We exercise process_record directly (no
	// binary) and inspect the accumulated counts.
	newPlugin := func(maxAC int) *contrastPlugin {
		return &contrastPlugin{
			annots:     contrastNASSOC, // keep annotation path simple
			controlIdx: []int{0, 1},
			caseIdx:    []int{2, 3},
			maxAC:      maxAC,
			maxACSet:   true,
		}
	}

	t.Run("alt_is_minor_added_verbatim", func(t *testing.T) {
		// One het carrier of ALT among the case group: nals = ctrlRef4? Let's make
		// the alt the rarer allele. ctrl: 0/0, 0/0 -> ref=4; case: 0/1, 0/0 ->
		// ref=3, alt=1. minor=alt (1 <= maxAC). Added verbatim.
		p := newPlugin(2)
		v := contrastTestVariant([]string{"0/0", "0/0", "0/1", "0/0"})
		if _, err := p.Process(v); err != nil {
			t.Fatal(err)
		}
		// nals = [ctrlRef, ctrlAlt, caseRef, caseAlt] = [4,0,3,1]
		if got := p.enrichNals; got != [4]int{4, 0, 3, 1} {
			t.Fatalf("enrichNals = %v, want [4 0 3 1]", got)
		}
	})

	t.Run("ref_is_minor_columns_swapped", func(t *testing.T) {
		// Make REF the rarer allele: ctrl 1/1,1/1 -> alt=4; case 1/1,0/1 ->
		// ref=1, alt=3. ref total (0+1)=1 <= maxAC, alt total (4+3)=7. minor=ref,
		// so columns swapped: enrich[0]+=ctrlAlt, [1]+=ctrlRef, [2]+=caseAlt,
		// [3]+=caseRef.
		p := newPlugin(2)
		v := contrastTestVariant([]string{"1/1", "1/1", "1/1", "0/1"})
		if _, err := p.Process(v); err != nil {
			t.Fatal(err)
		}
		// nals = [ctrlRef=0, ctrlAlt=4, caseRef=1, caseAlt=3].
		// swapped -> [ctrlAlt, ctrlRef, caseAlt, caseRef] = [4,0,3,1].
		if got := p.enrichNals; got != [4]int{4, 0, 3, 1} {
			t.Fatalf("enrichNals = %v, want [4 0 3 1]", got)
		}
	})

	t.Run("excluded_when_minor_exceeds_threshold", func(t *testing.T) {
		// alt total = 3 > maxAC(1): the record contributes nothing.
		p := newPlugin(1)
		v := contrastTestVariant([]string{"0/1", "0/0", "0/1", "0/1"})
		if _, err := p.Process(v); err != nil {
			t.Fatal(err)
		}
		if got := p.enrichNals; got != [4]int{0, 0, 0, 0} {
			t.Fatalf("enrichNals = %v, want all-zero (excluded)", got)
		}
	})

	t.Run("disabled_when_maxAC_unset", func(t *testing.T) {
		p := newPlugin(2)
		p.maxACSet = false
		v := contrastTestVariant([]string{"0/0", "0/0", "0/1", "0/0"})
		if _, err := p.Process(v); err != nil {
			t.Fatal(err)
		}
		if got := p.enrichNals; got != [4]int{0, 0, 0, 0} {
			t.Fatalf("enrichNals = %v, want all-zero when -f unset", got)
		}
	})
}

// trioTestVariant builds a biallelic record with per-sample GTs for a single or
// multiple trios laid out as child,father,mother,child,father,mother,...
func trioTestVariant(ref, alt string, gts []string) *vcf.Variant {
	v := &vcf.Variant{Chrom: "chr1", Pos: 100, Ref: ref, Alt: []string{alt}, Format: []string{"GT"}}
	for i, gt := range gts {
		v.Samples = append(v.Samples, vcf.Sample{Name: "X" + string(rune('1'+i)), Data: map[string]string{"GT": gt}})
	}
	return v
}

func TestUnitTrioStatsAltTriosDeferred(t *testing.T) {
	// Two trios sharing the same ALT allele as an untransmitted singleton:
	//   trio A: child 0/0, father 0/1, mother 0/0 -> ALT is a parent singleton
	//   trio B: child 0/0, father 0/1, mother 0/0 -> ALT is a parent singleton
	// The ALT allele therefore appears in 2 trios at this site.
	hdr := &vcf.Header{Samples: []string{"cA", "fA", "mA", "cB", "fB", "mB"}}
	trios := []trioStatsTrio{{child: 0, father: 1, mother: 2}, {child: 3, father: 4, mother: 5}}
	v := trioTestVariant("A", "T", []string{"0/0", "0/1", "0/0", "0/0", "0/1", "0/0"})

	run := func(maxAlt int) (int, int) {
		p := &trioStatsPlugin{hdr: hdr, trios: trios, maxAltTrios: maxAlt, out: &bytes.Buffer{}}
		flt := &trioStatsFilter{label: "all", stats: make([]trioStatsCounters, len(trios))}
		p.filters = []*trioStatsFilter{flt}
		if err := p.processOne(v, flt); err != nil {
			t.Fatal(err)
		}
		var sing int
		for _, s := range flt.stats {
			sing += int(s.nsingleton)
		}
		return sing, int(flt.stats[0].nsingleton)
	}

	// -a 0 (unlimited / default): both trios count their singleton.
	if total, _ := run(0); total != 2 {
		t.Fatalf("-a 0: total singletons = %d, want 2", total)
	}
	// -a 1: the ALT allele is present in 2 trios > 1, so NO singleton is counted.
	if total, _ := run(1); total != 0 {
		t.Fatalf("-a 1: total singletons = %d, want 0 (allele in 2 trios > 1)", total)
	}
	// -a 2: 2 trios <= 2, so both singletons are counted (deferred then applied).
	if total, _ := run(2); total != 2 {
		t.Fatalf("-a 2: total singletons = %d, want 2", total)
	}
}

// TestUnitIndelStatsPEDNoADRobust is the fix-on-port reproducer for the
// upstream abort documented in docs/UPSTREAM_BUGS.md: indel-stats.c's PED mode
// calls error("Incorrect GT allele") and exits 255 on an indel-bearing trio VCF
// that lacks FORMAT/AD (nad1==0 trips the als>=nad1 guard). Our port skips the
// AD-derived contributions and still produces a clean report. This runs through
// RunPlugin with NO upstream binary.
func TestUnitIndelStatsPEDNoADRobust(t *testing.T) {
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "noad.vcf")
	pedPath := filepath.Join(dir, "trios.ped")
	// A het DNM insertion in CHILD1 (parents 0/0), GT:GQ only, no FORMAT/AD.
	vcfBody := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=100000>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHILD1\tFATHER1\tMOTHER1\n" +
		"chr1\t100\t.\tA\tAT\t50\tPASS\t.\tGT:GQ\t0/1:40\t0/0:40\t0/0:40\n"
	if err := os.WriteFile(vcfPath, []byte(vcfBody), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(pedPath, []byte("FAM1\tCHILD1\tFATHER1\tMOTHER1\t1\t0\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	var out, errBuf bytes.Buffer
	err := RunPlugin(PluginOptions{
		Name:         "indel-stats",
		Args:         []string{"-p", pedPath},
		InputFile:    vcfPath,
		OutputFormat: OutputVCF,
	}, &out, &errBuf)
	if err != nil {
		t.Fatalf("expected a clean report on an AD-less PED indel VCF, got error: %v", err)
	}
	got := out.String()
	if !bytes.Contains(out.Bytes(), []byte("SN0\t")) {
		t.Fatalf("report missing SN0 line:\n%s", got)
	}
	// The de-novo insertion is counted genotype-wise even without AD: npass_gt
	// (column 5) and one insertion (column 6) must be non-zero. SN columns:
	// SN0 <nsmpl/ntrio> <nsites> <npass> <npass_gt> <nins> <ndel> ...
	for _, line := range bytes.Split(out.Bytes(), []byte("\n")) {
		if bytes.HasPrefix(line, []byte("SN0\t")) {
			cols := bytes.Split(line, []byte("\t"))
			if len(cols) < 7 {
				t.Fatalf("SN0 line malformed: %q", line)
			}
			if string(cols[1]) != "1" { // one trio
				t.Fatalf("SN0 trio count = %q, want 1", cols[1])
			}
			if string(cols[4]) == "0" || string(cols[5]) == "0" {
				t.Fatalf("expected non-zero de-novo genotype/insertion counts, got SN0 = %q", line)
			}
		}
	}
}

func TestUnitStatsReportWriter(t *testing.T) {
	t.Run("stdout_passthrough", func(t *testing.T) {
		var buf bytes.Buffer
		w, closeFn, err := statsReportWriter("", &buf)
		if err != nil {
			t.Fatal(err)
		}
		if w != &buf {
			t.Fatal("empty -o should return the stdout writer unchanged")
		}
		w.Write([]byte("hello"))
		if err := closeFn(); err != nil {
			t.Fatal(err)
		}
		if buf.String() != "hello" {
			t.Fatalf("got %q", buf.String())
		}
	})

	t.Run("dash_passthrough", func(t *testing.T) {
		var buf bytes.Buffer
		w, _, err := statsReportWriter("-", &buf)
		if err != nil {
			t.Fatal(err)
		}
		if w != &buf {
			t.Fatal(`"-" should return the stdout writer unchanged`)
		}
	})

	t.Run("file_target", func(t *testing.T) {
		path := filepath.Join(t.TempDir(), "report.txt")
		var buf bytes.Buffer
		w, closeFn, err := statsReportWriter(path, &buf)
		if err != nil {
			t.Fatal(err)
		}
		w.Write([]byte("report-bytes"))
		if err := closeFn(); err != nil {
			t.Fatal(err)
		}
		if buf.Len() != 0 {
			t.Fatalf("stdout should be untouched, got %q", buf.String())
		}
		got, err := os.ReadFile(path)
		if err != nil {
			t.Fatal(err)
		}
		if string(got) != "report-bytes" {
			t.Fatalf("file = %q", got)
		}
	})
}
