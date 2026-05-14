package mosdepth

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
	"github.com/yassineS/bio_ai_experiment/tools/tabix/pkg/tabix"
)

// makeBAM constructs an in-memory BGZF-wrapped BAM containing the given
// header + records. Returns the raw bytes ready to be handed to Run.
func makeBAM(t *testing.T, refs []sam.Reference, recs []*sam.Record) []byte {
	t.Helper()
	hdr := &sam.Header{Refs: refs}
	for _, r := range refs {
		hdr.Lines = append(hdr.Lines, sam.HeaderLine{
			Tag: "SQ",
			Fields: []sam.HeaderField{
				{Tag: "SN", Value: r.Name},
				{Tag: "LN", Value: intToStr(int(r.Length))},
			},
		})
	}
	var buf bytes.Buffer
	bw := sam.NewBAMWriter(&buf)
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, rec := range recs {
		if err := bw.Write(rec); err != nil {
			t.Fatalf("Write rec: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes()
}

func intToStr(v int) string {
	// Avoid strconv dependency cycle in test helpers.
	if v == 0 {
		return "0"
	}
	neg := v < 0
	if neg {
		v = -v
	}
	digits := []byte{}
	for v > 0 {
		digits = append([]byte{byte('0' + v%10)}, digits...)
		v /= 10
	}
	if neg {
		return "-" + string(digits)
	}
	return string(digits)
}

func mustParseCigar(t *testing.T, s string) sam.Cigar {
	t.Helper()
	c, err := sam.ParseCigar(s)
	if err != nil {
		t.Fatalf("ParseCigar(%q): %v", s, err)
	}
	return c
}

// fixtureRefs returns a small chr1/chr2 reference set used by most tests.
func fixtureRefs() []sam.Reference {
	return []sam.Reference{
		{Name: "chr1", Length: 50},
		{Name: "chr2", Length: 20},
	}
}

// fixtureRecords mirrors depthSAM in samtools/depth_test.go but uses tiny
// references so we can hand-compute expectations.
//
// Reads (POS is 1-based; CIGAR shown):
//
//	r1: chr1:10 5M       — covers 0-based 9..13
//	r2: chr1:12 5M       — covers 0-based 11..15
//	r3: chr1:20 1S3M1S   — covers 0-based 19..21 (clips don't consume ref)
//	r4: chr1:30 2M2I2M   — covers 0-based 29..30 and 31..32
//	r5: chr2:5 3M        — covers 0-based 4..6 on chr2
//	r6: chr1:40 5M, MAPQ=0   (for MAPQ filter)
//	r7: chr1:1 5M, FLAG=0x4 (unmapped; dropped by default exclude flag)
func fixtureRecords(t *testing.T) []*sam.Record {
	mk := func(name, ref string, pos int32, cigar string, mapq uint8, flag uint16) *sam.Record {
		return &sam.Record{
			QName: name, RName: ref, Pos: pos, MapQ: mapq, Flag: flag,
			Cigar: mustParseCigar(t, cigar),
			Seq:   strings.Repeat("A", mustParseCigar(t, cigar).QueryLength()),
		}
	}
	return []*sam.Record{
		mk("r1", "chr1", 10, "5M", 60, 0),
		mk("r2", "chr1", 12, "5M", 60, 0),
		mk("r3", "chr1", 20, "1S3M1S", 60, 0),
		mk("r4", "chr1", 30, "2M2I2M", 60, 0),
		mk("r5", "chr2", 5, "3M", 60, 0),
		mk("r6", "chr1", 40, "5M", 0, 0),
		mk("r7", "chr1", 1, "5M", 60, sam.FlagUnmapped),
	}
}

// readGz decompresses a gzip-wrapped (BGZF is gzip-compatible) file into
// its lines.
func readGz(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip.NewReader: %v", err)
	}
	gr.Multistream(true)
	data, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("read gz: %v", err)
	}
	return strings.Split(strings.TrimRight(string(data), "\n"), "\n")
}

func readLines(t *testing.T, path string) []string {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return strings.Split(strings.TrimRight(string(data), "\n"), "\n")
}

// TestRunPerBase verifies the per-base output matches a hand-computed table
// for the fixture BAM.
func TestRunPerBase(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ExcludeFlag: DefaultExcludeFlag, MinMAPQ: 0}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readGz(t, prefix+".per-base.bed.gz")
	// Look for some hand-computed runs on chr1:
	//   0..9 depth 0; 9..11 depth 1; 11..14 depth 2; 14..16 depth 1; ...
	wantRuns := map[string]bool{
		"chr1\t0\t9\t0":   true,
		"chr1\t9\t11\t1":  true,
		"chr1\t11\t14\t2": true,
		"chr1\t14\t16\t1": true,
		"chr1\t19\t22\t1": true,
		"chr1\t29\t33\t1": true, // r4 2M2I2M — the two M runs are adjacent on the reference (I consumes query only), so the per-base writer fuses them into one run
		"chr1\t39\t44\t1": true, // r6 still kept (mapq 0 allowed because MinMAPQ==0)
		"chr2\t4\t7\t1":   true,
	}
	for _, ln := range lines {
		delete(wantRuns, ln)
	}
	if len(wantRuns) > 0 {
		t.Errorf("missing per-base runs: %v\nfull output:\n%s", wantRuns, strings.Join(lines, "\n"))
	}
}

// TestRunNoPerBaseSuppressesOutput verifies --no-per-base prevents the file
// from being created.
func TestRunNoPerBaseSuppressesOutput(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ExcludeFlag: DefaultExcludeFlag, NoPerBase: true}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if _, err := os.Stat(prefix + ".per-base.bed.gz"); !os.IsNotExist(err) {
		t.Errorf("per-base.bed.gz should not exist (got err=%v)", err)
	}
	// Summary still emitted.
	if _, err := os.Stat(prefix + ".mosdepth.summary.txt"); err != nil {
		t.Errorf("summary.txt should exist: %v", err)
	}
}

// TestRunRegionsBED uses a 3-column BED to summarise per-region mean depth.
func TestRunRegionsBED(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	bedPath := filepath.Join(dir, "regions.bed")
	if err := os.WriteFile(bedPath, []byte("chr1\t9\t16\tROI1\nchr1\t19\t22\tROI2\n"), 0644); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	prefix := filepath.Join(dir, "out")
	opts := Options{
		Prefix:      prefix,
		ByBED:       bedPath,
		ExcludeFlag: DefaultExcludeFlag,
	}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if _, err := os.Stat(prefix + ".per-base.bed.gz"); !os.IsNotExist(err) {
		t.Errorf("per-base.bed.gz should NOT exist when --by is set (err=%v)", err)
	}
	lines := readGz(t, prefix+".regions.bed.gz")
	// ROI1 = chr1:9..16, depths per pos: 1,1,2,2,2,1,1 → sum=10, mean=10/7≈1.43.
	// ROI2 = chr1:19..22, depths: 1,1,1 → sum=3, mean=1.00.
	wantPrefix := map[string]string{
		"chr1\t9\t16":  "ROI1\t1.43",
		"chr1\t19\t22": "ROI2\t1.00",
	}
	for _, ln := range lines {
		parts := strings.SplitN(ln, "\t", 3)
		if len(parts) < 3 {
			continue
		}
		key := parts[0] + "\t" + parts[1] + "\t" + parts[2][:strings.Index(parts[2], "\t")]
		// shape: chrom\tstart\tend ; extras are after the third tab.
		// Actually start/end are first two fields. Recompute.
		_ = key
	}
	// Simpler: just check the literal expected lines exist.
	expected := map[string]bool{
		"chr1\t9\t16\tROI1\t1.43":  true,
		"chr1\t19\t22\tROI2\t1.00": true,
	}
	_ = wantPrefix
	for _, ln := range lines {
		delete(expected, ln)
	}
	if len(expected) > 0 {
		t.Errorf("missing region rows: %v\ngot lines:\n%s", expected, strings.Join(lines, "\n"))
	}
}

// TestRunRegionsWindow tests `-b 10` window summarisation.
func TestRunRegionsWindow(t *testing.T) {
	dir := t.TempDir()
	refs := []sam.Reference{{Name: "chr1", Length: 30}}
	recs := []*sam.Record{
		{QName: "r1", RName: "chr1", Pos: 1, Cigar: mustParseCigar(t, "10M"), MapQ: 60},
		{QName: "r2", RName: "chr1", Pos: 5, Cigar: mustParseCigar(t, "10M"), MapQ: 60},
	}
	bam := makeBAM(t, refs, recs)
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ByWindow: 10, ExcludeFlag: DefaultExcludeFlag}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readGz(t, prefix+".regions.bed.gz")
	// Windows: [0,10), [10,20), [20,30)
	// Depths at positions:
	//   r1 covers 0..9 (1)
	//   r2 covers 4..13 (1)
	// So:
	//   0..3 depth 1, 4..9 depth 2, 10..13 depth 1, 14..29 depth 0.
	// Window 0..10: sum = 1*4 + 2*6 = 16; mean = 1.60.
	// Window 10..20: sum = 1*4 + 0*6 = 4; mean = 0.40.
	// Window 20..30: sum = 0; mean = 0.00.
	want := map[string]bool{
		"chr1\t0\t10\t1.60":  true,
		"chr1\t10\t20\t0.40": true,
		"chr1\t20\t30\t0.00": true,
	}
	for _, ln := range lines {
		delete(want, ln)
	}
	if len(want) > 0 {
		t.Errorf("missing window rows: %v\ngot:\n%s", want, strings.Join(lines, "\n"))
	}
}

// TestRunThresholds exercises -T 1,5 on a region.
func TestRunThresholds(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	bedPath := filepath.Join(dir, "regions.bed")
	if err := os.WriteFile(bedPath, []byte("chr1\t9\t16\tROI1\n"), 0644); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	prefix := filepath.Join(dir, "out")
	opts := Options{
		Prefix:      prefix,
		ByBED:       bedPath,
		Thresholds:  []int{1, 2, 5},
		ExcludeFlag: DefaultExcludeFlag,
	}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readGz(t, prefix+".thresholds.bed.gz")
	// Header + one row.
	if len(lines) != 2 {
		t.Fatalf("thresholds lines: got %d, want 2:\n%s", len(lines), strings.Join(lines, "\n"))
	}
	header := lines[0]
	if !strings.HasPrefix(header, "#chrom\tstart\tend\tregion\t1X\t2X\t5X") {
		t.Errorf("bad header: %q", header)
	}
	// Row: chr1\t9\t16\tROI1\t<>=1>\t<>=2>\t<>=5>
	// Depths on chr1:9..16 = 1,1,2,2,2,1,1 → ≥1=7, ≥2=3, ≥5=0.
	row := lines[1]
	if !strings.HasPrefix(row, "chr1\t9\t16\tROI1\t7\t3\t0") {
		t.Errorf("bad row: %q", row)
	}
}

// TestRunMAPQFilter drops r6 (MAPQ 0) when MinMAPQ=10.
func TestRunMAPQFilter(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ExcludeFlag: DefaultExcludeFlag, MinMAPQ: 10}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readGz(t, prefix+".per-base.bed.gz")
	// chr1:39..43 (r6 region) should now be depth 0.
	for _, ln := range lines {
		if strings.HasPrefix(ln, "chr1\t39\t44\t") && !strings.HasSuffix(ln, "\t0") {
			t.Errorf("r6 region should be depth 0 after MAPQ filter: %q", ln)
		}
	}
	// Also confirm there is no positive-depth run covering 39..44.
	for _, ln := range lines {
		if strings.HasPrefix(ln, "chr1\t39") || strings.HasPrefix(ln, "chr1\t40") {
			if !strings.HasSuffix(ln, "\t0") {
				t.Errorf("unexpected non-zero run: %q", ln)
			}
		}
	}
}

// TestRunIncludeFlagFilter: include-flag requires ALL bits — set to
// FlagReverse (0x10) and only reads with that bit survive. None of our
// fixture reads have it set, so depth everywhere should be 0.
func TestRunIncludeFlagFilter(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	prefix := filepath.Join(dir, "out")
	opts := Options{
		Prefix:      prefix,
		ExcludeFlag: DefaultExcludeFlag,
		IncludeFlag: sam.FlagReverse,
	}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readGz(t, prefix+".per-base.bed.gz")
	for _, ln := range lines {
		f := strings.Split(ln, "\t")
		if f[3] != "0" {
			t.Errorf("expected all-zero depth with FlagReverse include-only, got %q", ln)
		}
	}
}

// TestRunFastModeIgnoresCIGAR — under fast mode, 2M2I2M should fill the
// whole 4-base ref span (no gap at the I).
func TestRunFastModeIgnoresCIGAR(t *testing.T) {
	dir := t.TempDir()
	refs := []sam.Reference{{Name: "chr1", Length: 20}}
	recs := []*sam.Record{
		{QName: "r1", RName: "chr1", Pos: 5, Cigar: mustParseCigar(t, "2M2D2M"), MapQ: 60},
	}
	bam := makeBAM(t, refs, recs)
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, FastMode: true, ExcludeFlag: DefaultExcludeFlag}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readGz(t, prefix+".per-base.bed.gz")
	// Whole 6-base run covered: 4..10 at depth 1.
	want := "chr1\t4\t10\t1"
	found := false
	for _, ln := range lines {
		if ln == want {
			found = true
		}
	}
	if !found {
		t.Errorf("expected %q in per-base output, got:\n%s", want, strings.Join(lines, "\n"))
	}
}

// TestRunChromRestrict — --chrom chr2 should only emit chr2 records.
func TestRunChromRestrict(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ExcludeFlag: DefaultExcludeFlag, Chrom: "chr2"}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readGz(t, prefix+".per-base.bed.gz")
	for _, ln := range lines {
		if !strings.HasPrefix(ln, "chr2\t") {
			t.Errorf("non-chr2 line emitted: %q", ln)
		}
	}
	// Summary must only mention chr2 + total.
	sumLines := readLines(t, prefix+".mosdepth.summary.txt")
	if len(sumLines) != 3 { // header + chr2 + total
		t.Errorf("summary lines: got %d, want 3:\n%s", len(sumLines), strings.Join(sumLines, "\n"))
	}
}

// TestRunDistributionCDF inspects the global distribution file.
func TestRunDistributionCDF(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ExcludeFlag: DefaultExcludeFlag}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readLines(t, prefix+".mosdepth.global.dist.txt")
	// Must contain rows starting with "total" listing depths down to 0.
	hasTotal0 := false
	for _, ln := range lines {
		if strings.HasPrefix(ln, "total\t0\t") {
			hasTotal0 = true
			// Proportion at depth >= 0 must be 1.00 (every base counted).
			if !strings.HasSuffix(ln, "\t1.00") {
				t.Errorf("total>=0 should be 1.00, got %q", ln)
			}
		}
	}
	if !hasTotal0 {
		t.Errorf("missing 'total\\t0' row in dist file:\n%s", strings.Join(lines, "\n"))
	}
}

// TestRunSummaryPerChrom inspects per-chrom totals.
func TestRunSummaryPerChrom(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ExcludeFlag: DefaultExcludeFlag}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := readLines(t, prefix+".mosdepth.summary.txt")
	// Header + chr1 + chr2 + total.
	if len(lines) != 4 {
		t.Fatalf("summary lines: got %d, want 4:\n%s", len(lines), strings.Join(lines, "\n"))
	}
	// chr2 row: length=20, bases = 3 (one 3M r5), mean=3/20=0.15.
	wantChr2 := "chr2\t20\t3\t0.15\t0\t1"
	if lines[2] != wantChr2 {
		t.Errorf("chr2 summary row: got %q, want %q", lines[2], wantChr2)
	}
}

// TestRunD4Rejected ensures D4Output triggers a clear error.
func TestRunD4Rejected(t *testing.T) {
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	err := Run(bytes.NewReader(bam), Options{Prefix: "/tmp/m", D4Output: true})
	if err == nil || !strings.Contains(err.Error(), "D4") {
		t.Errorf("expected D4 error, got %v", err)
	}
}

// TestRunByConflict surfaces an error when both ByBED and ByWindow are set.
func TestRunByConflict(t *testing.T) {
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	err := Run(bytes.NewReader(bam), Options{Prefix: "/tmp/m", ByBED: "x.bed", ByWindow: 100})
	if err == nil || !strings.Contains(err.Error(), "cannot specify both") {
		t.Errorf("expected by-conflict error, got %v", err)
	}
}

// TestRunEmptyPrefix rejects a blank prefix.
func TestRunEmptyPrefix(t *testing.T) {
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	err := Run(bytes.NewReader(bam), Options{})
	if err == nil || !strings.Contains(err.Error(), "empty output prefix") {
		t.Errorf("expected empty-prefix error, got %v", err)
	}
}

// TestRunRegionsTbiReadable confirms the regions .tbi round-trips through
// tabix.QueryBytes on a single chromosome.
func TestRunRegionsTbiReadable(t *testing.T) {
	dir := t.TempDir()
	refs := []sam.Reference{{Name: "chr1", Length: 30}}
	recs := []*sam.Record{
		{QName: "r1", RName: "chr1", Pos: 1, Cigar: mustParseCigar(t, "10M"), MapQ: 60},
	}
	bam := makeBAM(t, refs, recs)
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ByWindow: 10, ExcludeFlag: DefaultExcludeFlag}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	dataPath := prefix + ".regions.bed.gz"
	idx, err := tabix.ReadFile(dataPath + ".tbi")
	if err != nil {
		t.Fatalf("ReadFile tbi: %v", err)
	}
	rows, err := idx.QueryBytes(dataPath, "chr1", 0, 15)
	if err != nil {
		t.Fatalf("QueryBytes: %v", err)
	}
	if len(rows) == 0 {
		t.Errorf("expected at least one row from regions tbi query")
	}
}

// TestRunTbiReadable rebuilds the TBI for the per-base file and queries
// chr1:10-15 via tabix to confirm the index round-trips.
func TestRunTbiReadable(t *testing.T) {
	dir := t.TempDir()
	bam := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	prefix := filepath.Join(dir, "out")
	opts := Options{Prefix: prefix, ExcludeFlag: DefaultExcludeFlag}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	dataPath := prefix + ".per-base.bed.gz"
	tbiPath := dataPath + ".tbi"
	if _, err := os.Stat(tbiPath); err != nil {
		t.Fatalf(".tbi missing: %v", err)
	}
	idx, err := tabix.ReadFile(tbiPath)
	if err != nil {
		t.Fatalf("tabix.ReadFile: %v", err)
	}
	rows, err := idx.QueryBytes(dataPath, "chr1", 10, 15)
	if err != nil {
		t.Fatalf("tabix QueryBytes: %v", err)
	}
	if len(rows) == 0 {
		t.Errorf("tabix query returned no rows")
	}
}

// TestRunByArg_NumberAndPath covers parseByArg's two branches.
func TestParseByArg(t *testing.T) {
	w, p, err := ParseByArg("100")
	if err != nil || w != 100 || p != "" {
		t.Errorf("ParseByArg(100): got w=%d p=%q err=%v", w, p, err)
	}
	w, p, err = ParseByArg("")
	if err != nil || w != 0 || p != "" {
		t.Errorf("ParseByArg empty: got w=%d p=%q err=%v", w, p, err)
	}
	_, _, err = ParseByArg("-5")
	if err == nil {
		t.Errorf("ParseByArg(-5): expected error")
	}
	_, _, err = ParseByArg("/nope/does/not/exist.bed")
	if err == nil {
		t.Errorf("ParseByArg with missing file: expected error")
	}
	// Real file path.
	dir := t.TempDir()
	bedPath := filepath.Join(dir, "x.bed")
	os.WriteFile(bedPath, []byte("chr1\t0\t1\n"), 0644)
	w, p, err = ParseByArg(bedPath)
	if err != nil || w != 0 || p != bedPath {
		t.Errorf("ParseByArg(bed): got w=%d p=%q err=%v", w, p, err)
	}
}

// TestParseReadGroups covers the OPS prefix mode and de-trimming.
func TestParseReadGroups(t *testing.T) {
	got := ParseReadGroups(" a, b ,c,")
	if len(got) != 3 || got[0] != "a" || got[1] != "b" || got[2] != "c" {
		t.Errorf("ParseReadGroups: got %v", got)
	}
	got = ParseReadGroups("")
	if got != nil {
		t.Errorf("expected nil for empty input, got %v", got)
	}
}

// TestKeepRecord_FragLen exercises min/max fragment length filtering.
func TestKeepRecord_FragLen(t *testing.T) {
	cases := []struct {
		name string
		tlen int32
		min  int
		max  int
		want bool
	}{
		{"in-range", 200, 100, 500, true},
		{"below-min", 50, 100, 500, false},
		{"above-max", 600, 100, 500, false},
		{"negative-abs-keeps", -200, 100, 500, true},
		{"only-min", 200, 100, 0, true},
		{"only-max-fail", 600, 0, 500, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			rec := &sam.Record{
				QName: "r", RName: "chr1", Pos: 1, MapQ: 60,
				TLen: tc.tlen,
			}
			opts := Options{MinFragLen: tc.min, MaxFragLen: tc.max}
			if got := keepRecord(rec, opts); got != tc.want {
				t.Errorf("got %v want %v (tlen=%d min=%d max=%d)", got, tc.want, tc.tlen, tc.min, tc.max)
			}
		})
	}
}

// TestKeepRecord_ReadGroupRG checks the RG aux filter.
func TestKeepRecord_ReadGroupRG(t *testing.T) {
	rec := &sam.Record{
		QName: "r", RName: "chr1", Pos: 1, MapQ: 60,
		Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg2"}},
	}
	if !keepRecord(rec, Options{ReadGroups: []string{"rg1", "rg2"}}) {
		t.Errorf("expected RG match")
	}
	if keepRecord(rec, Options{ReadGroups: []string{"rg3"}}) {
		t.Errorf("expected RG mismatch")
	}
	// OPS mode: switches to OPS aux tag.
	recOps := &sam.Record{
		QName: "r", RName: "chr1", Pos: 1, MapQ: 60,
		Aux: []sam.Aux{{Tag: "OPS", Type: 'Z', Value: "X"}},
	}
	if !keepRecord(recOps, Options{ReadGroups: []string{"OPS:X", "Y"}}) {
		t.Errorf("expected OPS match")
	}
	if keepRecord(recOps, Options{ReadGroups: []string{"OPS:Z"}}) {
		t.Errorf("expected OPS mismatch")
	}
}

// TestOpenAndRun reads a BAM from disk end-to-end.
func TestOpenAndRun(t *testing.T) {
	dir := t.TempDir()
	bamBytes := makeBAM(t, fixtureRefs(), fixtureRecords(t))
	bamPath := filepath.Join(dir, "in.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0644); err != nil {
		t.Fatalf("write bam: %v", err)
	}
	prefix := filepath.Join(dir, "out")
	if err := OpenAndRun(bamPath, Options{Prefix: prefix, ExcludeFlag: DefaultExcludeFlag}); err != nil {
		t.Fatalf("OpenAndRun: %v", err)
	}
	if _, err := os.Stat(prefix + ".per-base.bed.gz"); err != nil {
		t.Errorf("per-base missing: %v", err)
	}
}

// TestOpenAndRun_MissingFile surfaces an open error.
func TestOpenAndRun_MissingFile(t *testing.T) {
	err := OpenAndRun("/no/such/file.bam", Options{Prefix: "/tmp/x"})
	if err == nil {
		t.Errorf("expected error opening missing file")
	}
}
