package prinseq

// Tests for the `--qual_noscale` graph-data knob (upstream
// prinseq-lite.pl:281-283, 989-993, 4777). The flag is a graph-data
// concern only: when set, upstream's $scale flips from 1 to 0, which
//   1. writes `"scale":0` (instead of 1) into the emitted .gd JSON, and
//   2. suppresses the relative (100-bin) `quals` table entirely while
//      leaving the absolute per-position `quala` table untouched.
//
// It has no effect on the stats/filter paths (where there is no scaling
// of quality data to a fixed bin count). These tests pin both the unit
// behaviour and, when a Perl interpreter is available, byte-for-byte
// parity against the live upstream oracle.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// TestQualNoscaleSuppressesRelativeTable verifies the in-memory
// collector contract: with --qual_noscale the relative `Quals` table is
// empty while the absolute `Quala` table is fully populated, and the
// emitted JSON carries "scale":0 with no "quals" block. The default
// (scaled) run is asserted as the contrasting control.
func TestQualNoscaleSuppressesRelativeTable(t *testing.T) {
	fastqPath := filepath.Join("..", "..", "testdata", "parity", "graphdata_example1.fastq")

	collect := func(noscale bool) *GraphData {
		f, err := os.Open(fastqPath)
		if err != nil {
			t.Fatalf("open fastq: %v", err)
		}
		defer f.Close()
		opts := DefaultGraphDataOptions()
		opts.Filename1 = "6578616d706c65312e6661737471"
		opts.QualNoscale = noscale
		g, err := CollectGraphData(f, true, opts)
		if err != nil {
			t.Fatalf("collect (noscale=%v): %v", noscale, err)
		}
		return g
	}

	// Default: relative table populated.
	def := collect(false)
	if len(def.Quals) == 0 {
		t.Fatalf("default run: relative Quals table is empty, expected it populated")
	}
	if len(def.Quala) == 0 {
		t.Fatalf("default run: absolute Quala table is empty")
	}

	// --qual_noscale: relative table suppressed, absolute kept.
	ns := collect(true)
	if len(ns.Quals) != 0 {
		t.Fatalf("qual_noscale: relative Quals table has %d entries, expected 0", len(ns.Quals))
	}
	if len(ns.Quala) == 0 {
		t.Fatalf("qual_noscale: absolute Quala table is empty, expected it populated")
	}
	// The absolute table must be identical regardless of the flag.
	if len(ns.Quala) != len(def.Quala) {
		t.Fatalf("qual_noscale changed Quala size: %d vs default %d", len(ns.Quala), len(def.Quala))
	}

	// Emitted JSON: scale field flips, quals block disappears.
	var bufDef, bufNS bytes.Buffer
	if err := def.EmitGD(&bufDef, GDHeader{}); err != nil {
		t.Fatalf("emit default: %v", err)
	}
	if err := ns.EmitGD(&bufNS, GDHeader{}); err != nil {
		t.Fatalf("emit qual_noscale: %v", err)
	}
	if !bytes.Contains(bufDef.Bytes(), []byte(`"scale":1`)) {
		t.Errorf("default emit missing \"scale\":1")
	}
	if !bytes.Contains(bufNS.Bytes(), []byte(`"scale":0`)) {
		t.Errorf("qual_noscale emit missing \"scale\":0")
	}
	if !bytes.Contains(bufDef.Bytes(), []byte(`"quals":`)) {
		t.Errorf("default emit missing \"quals\" block")
	}
	if bytes.Contains(bufNS.Bytes(), []byte(`"quals":`)) {
		t.Errorf("qual_noscale emit unexpectedly contains \"quals\" block")
	}
	// The absolute table is emitted as "qualsbin" (binned per-position)
	// plus the per-read "qualsmean" histogram; both must survive in
	// both modes since --qual_noscale only suppresses the relative
	// "quals" block.
	for _, key := range [][]byte{[]byte(`"qualsbin":`), []byte(`"qualsmean":`)} {
		if !bytes.Contains(bufNS.Bytes(), key) {
			t.Errorf("qual_noscale emit missing %s block", key)
		}
		if !bytes.Contains(bufDef.Bytes(), key) {
			t.Errorf("default emit missing %s block", key)
		}
	}
}

// upstreamQualNoscaleInputName is the fixed basename the live-oracle
// helper writes the FASTQ under. Upstream hex-encodes this basename into
// the .gd "filename1" field, so the Go emit must use the same hex
// (qualNoscaleFilename1Hex) for the parity comparison to pass.
const upstreamQualNoscaleInputName = "qns_input.fastq"

// qualNoscaleFilename1Hex is the upstream convertStringToInt encoding
// (lowercase hex of each byte) of upstreamQualNoscaleInputName.
func qualNoscaleFilename1Hex() string {
	var b strings.Builder
	for i := 0; i < len(upstreamQualNoscaleInputName); i++ {
		// %x of a byte, matching prinseq-lite.pl:4853.
		const hexdigits = "0123456789abcdef"
		c := upstreamQualNoscaleInputName[i]
		b.WriteByte(hexdigits[c>>4])
		b.WriteByte(hexdigits[c&0xf])
	}
	return b.String()
}

// runUpstreamPrinseqQualNoscale runs the vendored upstream
// prinseq-lite.pl over the given FASTQ with --graph_data, optionally
// adding --qual_noscale, and returns the JSON body of the produced .gd
// file (the two leading #-comment lines stripped). The input is copied
// to a fixed basename so the embedded "filename1" hex is predictable. It
// t.Skips only when the Perl interpreter or the submodule script is
// unavailable; any other failure is fatal.
func runUpstreamPrinseqQualNoscale(t *testing.T, fastqPath string, noscale bool) []byte {
	t.Helper()
	perl, err := exec.LookPath("perl")
	if err != nil {
		t.Skipf("perl not available: %v", err)
	}
	script := filepath.Join("..", "..", "..", "..", "reference_code", "prinseq", "prinseq-lite.pl")
	if _, err := os.Stat(script); err != nil {
		t.Skipf("upstream prinseq-lite.pl not available (submodule not initialised): %v", err)
	}

	tmp := t.TempDir()
	src, err := os.ReadFile(fastqPath)
	if err != nil {
		t.Fatalf("read fastq: %v", err)
	}
	inPath := filepath.Join(tmp, upstreamQualNoscaleInputName)
	if err := os.WriteFile(inPath, src, 0o600); err != nil {
		t.Fatalf("write fastq copy: %v", err)
	}
	gdPath := filepath.Join(tmp, "out.gd")
	args := []string{script, "-fastq", inPath, "-graph_data", gdPath}
	if noscale {
		args = append(args, "-qual_noscale")
	}
	cmd := exec.Command(perl, args...)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream prinseq-lite.pl failed: %v\nstderr:\n%s", err, stderr.String())
	}

	raw, err := os.ReadFile(gdPath)
	if err != nil {
		t.Fatalf("read upstream .gd: %v", err)
	}
	// Strip leading #-comment lines to expose the JSON body.
	idx := 0
	for idx < len(raw) && raw[idx] == '#' {
		nl := bytes.IndexByte(raw[idx:], '\n')
		if nl < 0 {
			idx = len(raw)
			break
		}
		idx += nl + 1
	}
	return raw[idx:]
}

// TestQualNoscaleLiveParity drives the live upstream oracle with and
// without --qual_noscale and asserts the Go port's emitted graph-data
// matches it structurally (up to the usual numeric tolerance) for both
// settings. This proves the flag is fully ported end-to-end.
func TestQualNoscaleLiveParity(t *testing.T) {
	fastqPath := filepath.Join("..", "..", "testdata", "parity", "graphdata_example1.fastq")

	for _, noscale := range []bool{false, true} {
		noscale := noscale
		name := "scaled"
		if noscale {
			name = "qual_noscale"
		}
		t.Run(name, func(t *testing.T) {
			upstreamBody := runUpstreamPrinseqQualNoscale(t, fastqPath, noscale)

			f, err := os.Open(fastqPath)
			if err != nil {
				t.Fatalf("open fastq: %v", err)
			}
			defer f.Close()
			opts := DefaultGraphDataOptions()
			// Upstream embeds filename1 as the hex-encoded basename of
			// the input it was given (the fixed copy name).
			opts.Filename1 = qualNoscaleFilename1Hex()
			opts.QualNoscale = noscale
			g, err := CollectGraphData(f, true, opts)
			if err != nil {
				t.Fatalf("collect: %v", err)
			}
			var buf bytes.Buffer
			if err := g.EmitGD(&buf, GDHeader{}); err != nil {
				t.Fatalf("emit: %v", err)
			}

			upstream := normaliseGD(t, upstreamBody)
			ours := normaliseGD(t, buf.Bytes())

			// Sanity: confirm the oracle actually toggled scale, so a
			// silently-no-op flag can't pass this test.
			um, ok := upstream.(map[string]any)
			if !ok {
				t.Fatalf("upstream .gd is not a JSON object")
			}
			wantScale := 1.0
			if noscale {
				wantScale = 0.0
			}
			if got, ok := um["scale"]; !ok || got != wantScale {
				t.Fatalf("upstream scale = %v (present=%v), want %v", got, ok, wantScale)
			}
			if _, hasQuals := um["quals"]; hasQuals == noscale {
				t.Fatalf("upstream quals presence=%v but qual_noscale=%v (expected absence iff noscale)", hasQuals, noscale)
			}

			diffs := gdDiff(upstream, ours, "", 1e-3)
			if len(diffs) > 0 {
				sort.Strings(diffs)
				max := len(diffs)
				if max > 60 {
					max = 60
				}
				t.Fatalf("graph-data divergence with qual_noscale=%v (%d entries):\n%s",
					noscale, len(diffs), strings.Join(diffs[:max], "\n"))
			}
		})
	}
}
