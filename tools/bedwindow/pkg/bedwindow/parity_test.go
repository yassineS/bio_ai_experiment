package bedwindow

// Parity tests for `bedtools window`.
//
// The upstream bedtools test corpus ships no `window/` subdirectory, so the
// spec-driven cases below encode the documented behaviour, and a second set
// (TestParity_Window_Live*) asserts byte-for-byte equality against the LIVE
// upstream bedtools binary for the order-sensitive and column-preserving paths
// that the fixed bugs touch:
//
//   - per-A B-hit order follows upstream's UCSC bin traversal (finest level
//     first, bin number ascending, then file order) — NOT plain file order;
//   - the default window (-w) is 1000 bp;
//   - BED12 (and any wide) B records are echoed verbatim, not truncated.
//
// Window adds the slop to A (not B), so -l grows the window upstream of A and
// -r downstream; the asymmetric cases below reflect that.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

// --- spec-driven cases (binary-free) ---------------------------------------

// TestParity_NoExpansionDefault — default writer is A<TAB>B for each overlap.
func TestParity_NoExpansionDefault(t *testing.T) {
	a := strings.NewReader("chr1\t0\t10\n")
	b := strings.NewReader("chr1\t5\t15\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1", n)
	}
	if !strings.Contains(out.String(), "chr1\t0\t10\tchr1\t5\t15") {
		t.Errorf("expected A<TAB>B in output, got %q", out.String())
	}
}

// TestParity_WindowExpandsAWindow — symmetric -w 100 expansion of A's window
// pulls a non-touching B into the overlap set.
func TestParity_WindowExpandsAWindow(t *testing.T) {
	a := strings.NewReader("chr1\t300\t310\n")
	b := strings.NewReader("chr1\t200\t210\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{Left: 100, Right: 100})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1 (A window expanded to [200,410))", n)
	}
}

// TestParity_CountMode — -c emits A<TAB>count, one row per A even with no hit.
func TestParity_CountMode(t *testing.T) {
	a := strings.NewReader("chr1\t0\t10\nchr2\t0\t10\n")
	b := strings.NewReader("chr1\t5\t6\nchr1\t8\t9\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{Count: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	got := out.String()
	if !strings.Contains(got, "chr1\t0\t10\t2\n") {
		t.Errorf("expected chr1 count=2 line, got:\n%s", got)
	}
	if !strings.Contains(got, "chr2\t0\t10\t0\n") {
		t.Errorf("expected chr2 count=0 line, got:\n%s", got)
	}
}

// TestParity_InvertMode — -v emits only A records with NO B overlap.
func TestParity_InvertMode(t *testing.T) {
	a := strings.NewReader("chr1\t0\t10\nchr2\t0\t10\n")
	b := strings.NewReader("chr1\t5\t6\n")
	var out bytes.Buffer
	if _, err := Window(a, b, &out, Options{Invert: true}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	got := strings.TrimSpace(out.String())
	if !strings.Contains(got, "chr2\t0\t10") {
		t.Errorf("expected chr2 in invert output, got %q", got)
	}
	if strings.Contains(got, "chr1\t0\t10") {
		t.Errorf("chr1 has overlap; should NOT appear in invert output: %q", got)
	}
}

// TestParity_AsymmetricExpansion — -l 0 -r 100 grows only the downstream side
// of A's window, so a B upstream of A still doesn't overlap.
func TestParity_AsymmetricExpansion(t *testing.T) {
	a := strings.NewReader("chr1\t300\t310\n")
	b := strings.NewReader("chr1\t100\t150\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{Left: 0, Right: 100})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 0 {
		t.Errorf("n = %d, want 0 (downstream-only expansion shouldn't reach upstream B)", n)
	}
}

// TestParity_LeftSlopReachesUpstreamB — -l 250 -r 0 grows A's window upstream,
// pulling in a B that lies before A. The low end clips at 0.
func TestParity_LeftSlopReachesUpstreamB(t *testing.T) {
	a := strings.NewReader("chr1\t300\t310\n")
	b := strings.NewReader("chr1\t100\t150\n")
	var out bytes.Buffer
	n, err := Window(a, b, &out, Options{Left: 250, Right: 0})
	if err != nil {
		t.Fatalf("Window: %v", err)
	}
	if n != 1 {
		t.Errorf("n = %d, want 1 (A window [50,310) reaches B at [100,150))", n)
	}
}

// --- live upstream binary parity ------------------------------------------

var (
	upstreamOnce sync.Once
	upstreamPath string
	upstreamErr  error
)

func upstreamBedtools(t *testing.T) string {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamOnce.Do(func() {
		_, file, _, _ := runtime.Caller(0)
		dir := filepath.Dir(file)
		var root string
		for {
			if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
				root = dir
				break
			}
			parent := filepath.Dir(dir)
			if parent == dir {
				upstreamErr = os.ErrNotExist
				return
			}
			dir = parent
		}
		btDir := filepath.Join(root, "reference_code", "bedtools")
		bin := filepath.Join(btDir, "bin", "bedtools")
		if _, statErr := os.Stat(bin); statErr != nil {
			cmd := exec.Command("make", "-j", "4")
			cmd.Dir = btDir
			if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
				upstreamErr = &exec.ExitError{Stderr: out}
				return
			}
		}
		upstreamPath = bin
	})
	if upstreamErr != nil || upstreamPath == "" {
		t.Skipf("upstream bedtools unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools && "+
			"(cd reference_code/bedtools && make -j\"$(nproc)\")", upstreamErr)
	}
	return upstreamPath
}

func writeTemp(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

func upstreamWindow(t *testing.T, bin, aFile, bFile string, flags ...string) []byte {
	t.Helper()
	args := append([]string{"window", "-a", aFile, "-b", bFile}, flags...)
	var out, errBuf bytes.Buffer
	cmd := exec.Command(bin, args...)
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream window %v: %v\nstderr: %s", args, err, errBuf.String())
	}
	return out.Bytes()
}

// TestParity_Window_Live_HitOrderBins asserts the per-A B-hit order matches
// upstream's UCSC bin traversal, on B records that fall in different bin levels
// in a non-bin file order.
func TestParity_Window_Live_HitOrderBins(t *testing.T) {
	bin := upstreamBedtools(t)
	aContent := "chr1\t100000\t110000\ta1\t10\t+\n"
	// fine1/fine2/fineDup share the finest bin; wide1 and wider are coarser.
	// File order deliberately interleaves levels.
	bContent := "chr1\t100000\t100100\tfine1\t5\t+\n" +
		"chr1\t50000\t200000\twide1\t5\t+\n" +
		"chr1\t105000\t105100\tfine2\t5\t+\n" +
		"chr1\t100000\t100100\tfineDup\t5\t+\n" +
		"chr1\t16000\t600000\twider\t5\t+\n"
	aFile := writeTemp(t, "a.bed", aContent)
	bFile := writeTemp(t, "b.bed", bContent)

	want := upstreamWindow(t, bin, aFile, bFile, "-w", "50000")
	var got bytes.Buffer
	if _, err := Window(strings.NewReader(aContent), strings.NewReader(bContent), &got,
		Options{Left: 50000, Right: 50000}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("bin hit-order mismatch.\nwant (upstream):\n%s\ngot:\n%s", want, got.String())
	}
}

// TestParity_Window_Live_DefaultWindow asserts the default window (no -w/-l/-r)
// is 1000 bp, matching upstream.
func TestParity_Window_Live_DefaultWindow(t *testing.T) {
	bin := upstreamBedtools(t)
	// B is 800 bp upstream of A: in range for the 1000 bp default, out of range
	// for a 500 bp window.
	aContent := "chr1\t2000\t2100\ta1\t0\t+\n"
	bContent := "chr1\t1100\t1200\tb1\t0\t+\n"
	aFile := writeTemp(t, "a.bed", aContent)
	bFile := writeTemp(t, "b.bed", bContent)

	want := upstreamWindow(t, bin, aFile, bFile) // no flags -> default 1000
	var got bytes.Buffer
	// Default-window resolution lives in the CLI; the library default for the
	// "no window" call is Left=Right=1000 to match.
	if _, err := Window(strings.NewReader(aContent), strings.NewReader(bContent), &got,
		Options{Left: 1000, Right: 1000}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("default-window mismatch.\nwant (upstream):\n%s\ngot:\n%s", want, got.String())
	}
	if len(want) == 0 {
		t.Fatal("expected at least one hit at the 1000bp default window")
	}
}

// TestParity_Window_Live_BED12B asserts a BED12 B record is echoed verbatim
// (all 12 columns, including block columns), not truncated to 6 columns.
func TestParity_Window_Live_BED12B(t *testing.T) {
	bin := upstreamBedtools(t)
	aContent := "chr1\t1000\t2000\ta1\t10\t+\n"
	bContent := "chr1\t1500\t1800\tbgene\t60\t+\t1520\t1780\t255,0,0\t2\t100,50,\t0,250,\n"
	aFile := writeTemp(t, "a.bed", aContent)
	bFile := writeTemp(t, "b.bed", bContent)

	// Upstream `bedtools window` has no -wa/-wb: the default already prints the
	// full A and B records, so the verbatim BED12 B column is exercised by the
	// default mode.
	want := upstreamWindow(t, bin, aFile, bFile, "-w", "0")
	var got bytes.Buffer
	if _, err := Window(strings.NewReader(aContent), strings.NewReader(bContent), &got,
		Options{Left: 0, Right: 0}); err != nil {
		t.Fatalf("Window: %v", err)
	}
	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("BED12-B mismatch.\nwant (upstream):\n%s\ngot:\n%s", want, got.String())
	}
	if !strings.Contains(string(want), "100,50,") {
		t.Fatal("sanity: upstream output should contain BED12 block sizes")
	}
}

// TestParity_Window_Live_Strand asserts the -sm/-Sm strand filters match
// upstream across a mixed-strand B set.
func TestParity_Window_Live_Strand(t *testing.T) {
	bin := upstreamBedtools(t)
	aContent := "chr1\t1000\t2000\ta\t0\t+\n"
	bContent := "chr1\t1500\t1600\tbplus\t0\t+\n" +
		"chr1\t1700\t1800\tbminus\t0\t-\n"
	aFile := writeTemp(t, "a.bed", aContent)
	bFile := writeTemp(t, "b.bed", bContent)

	cases := []struct {
		flag string
		opts Options
	}{
		{"-sm", Options{Left: 0, Right: 0, StrandSpec: true}},
		{"-Sm", Options{Left: 0, Right: 0, InverseStrand: true}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.flag, func(t *testing.T) {
			want := upstreamWindow(t, bin, aFile, bFile, "-w", "0", tc.flag)
			var got bytes.Buffer
			if _, err := Window(strings.NewReader(aContent), strings.NewReader(bContent), &got, tc.opts); err != nil {
				t.Fatalf("Window: %v", err)
			}
			if !bytes.Equal(got.Bytes(), want) {
				t.Fatalf("%s mismatch.\nwant (upstream):\n%s\ngot:\n%s", tc.flag, want, got.String())
			}
		})
	}
}
