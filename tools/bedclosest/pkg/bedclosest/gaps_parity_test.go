package bedclosest

// Live-upstream parity tests for the bedclosest directional flags
// (-iu/-id/-fu/-fd). They build the real upstream `bedtools` binary from the
// vendored submodule and compare its `closest` output, byte for byte, against
// this port's Closest function. They t.Fatalf (never t.Skip) so a missing or
// unbuildable submodule is a hard failure, matching the project's parity-rig
// policy.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

var (
	upstreamBedtoolsGapsOnce sync.Once
	upstreamBedtoolsGapsPath string
	upstreamBedtoolsGapsErr  error
)

// gapsRepoRoot walks up from this test file to the module root (the dir holding
// go.mod).
func gapsRepoRoot(t *testing.T) string {
	t.Helper()
	_, file, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed")
	}
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not locate repo root above %s", file)
		}
		dir = parent
	}
}

// upstreamBedtoolsGaps builds (once) and returns the path to the upstream
// `bedtools` binary. Uniquely named so this suite's builder is independent of
// sibling packages' builders.
func upstreamBedtoolsGaps(t *testing.T) string {
	t.Helper()
	upstreamBedtoolsGapsOnce.Do(func() {
		root := gapsRepoRoot(t)
		dir := filepath.Join(root, "reference_code", "bedtools")
		if _, err := os.Stat(filepath.Join(dir, "Makefile")); err != nil {
			upstreamBedtoolsGapsErr = err
			return
		}
		bin := filepath.Join(dir, "bin", "bedtools")
		if _, err := os.Stat(bin); err != nil {
			cmd := exec.Command("make", "-j", "4")
			cmd.Dir = dir
			if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
				upstreamBedtoolsGapsErr = buildErr
				t.Logf("bedtools build output:\n%s", out)
				return
			}
		}
		upstreamBedtoolsGapsPath = bin
	})
	if upstreamBedtoolsGapsErr != nil {
		t.Fatalf("upstream bedtools unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools && "+
			"(cd reference_code/bedtools && make -j\"$(nproc)\")", upstreamBedtoolsGapsErr)
	}
	return upstreamBedtoolsGapsPath
}

// runUpstreamClosest writes A/B to temp files and captures the upstream
// `bedtools closest` stdout for the given flags.
func runUpstreamClosest(t *testing.T, bt, aData, bData string, flags ...string) []byte {
	t.Helper()
	dir := t.TempDir()
	aPath := filepath.Join(dir, "a.bed")
	bPath := filepath.Join(dir, "b.bed")
	if err := os.WriteFile(aPath, []byte(aData), 0o644); err != nil {
		t.Fatalf("write a: %v", err)
	}
	if err := os.WriteFile(bPath, []byte(bData), 0o644); err != nil {
		t.Fatalf("write b: %v", err)
	}
	args := append([]string{"closest", "-a", aPath, "-b", bPath}, flags...)
	cmd := exec.Command(bt, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream closest %v failed: %v\nstderr:\n%s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// runOurClosest drives the port's Closest with the supplied options.
func runOurClosest(t *testing.T, aData, bData string, opts Options) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := Closest(bytes.NewReader([]byte(aData)), bytes.NewReader([]byte(bData)), &out, opts); err != nil {
		t.Fatalf("Closest failed: %v", err)
	}
	return out.Bytes()
}

// TestGapsParity_ClosestDirectional asserts byte-for-byte parity for the
// -iu/-id/-fu/-fd directional flags across the three -D modes, overlap vs no
// overlap, both A strands, and tie scenarios.
func TestGapsParity_ClosestDirectional(t *testing.T) {
	bt := upstreamBedtoolsGaps(t)

	const (
		// A on the forward strand and a reverse-strand twin.
		aFwd = "chr1\t100\t200\ta1\t.\t+\n"
		aRev = "chr1\t100\t200\ta1\t.\t-\n"
		// B with an upstream, an overlapping, and a downstream feature.
		bOvl = "chr1\t50\t60\tup\t.\t+\nchr1\t150\t160\tovl\t.\t+\nchr1\t300\t400\tdown\t.\t+\n"
		// B with only upstream + downstream (no overlap).
		bNoOvl = "chr1\t50\t60\tup\t.\t+\nchr1\t300\t400\tdown\t.\t+\n"
		// Reverse-strand B twins for the -D b mode.
		bNoOvlRev = "chr1\t50\t60\tup\t.\t-\nchr1\t300\t400\tdown\t.\t-\n"
		// Equidistant upstream ties plus two downstream features.
		aTie = "chr1\t1000\t1100\tq\t.\t+\n"
		bTie = "chr1\t800\t900\tu1\t.\t+\nchr1\t800\t900\tu1b\t.\t+\n" +
			"chr1\t1200\t1300\td1\t.\t+\nchr1\t1250\t1350\td2\t.\t+\n"
	)

	modeFor := func(d string) DistanceMode {
		switch d {
		case "a":
			return DistanceA
		case "b":
			return DistanceB
		default:
			return DistanceRef
		}
	}

	type tc struct {
		name   string
		a, b   string
		dmode  string
		iu, id bool
		fu, fd bool
		tie    TieBreak
	}
	cases := []tc{
		{"ref_iu_ovl", aFwd, bOvl, "ref", true, false, false, false, TieAll},
		{"ref_id_ovl", aFwd, bOvl, "ref", false, true, false, false, TieAll},
		{"ref_fu_ovl", aFwd, bOvl, "ref", false, false, true, false, TieAll},
		{"ref_fd_ovl", aFwd, bOvl, "ref", false, false, false, true, TieAll},
		{"ref_iu_noovl", aFwd, bNoOvl, "ref", true, false, false, false, TieAll},
		{"ref_id_noovl", aFwd, bNoOvl, "ref", false, true, false, false, TieAll},
		{"ref_fu_noovl", aFwd, bNoOvl, "ref", false, false, true, false, TieAll},
		{"ref_fd_noovl", aFwd, bNoOvl, "ref", false, false, false, true, TieAll},
		{"a_iu_rev", aRev, bNoOvl, "a", true, false, false, false, TieAll},
		{"a_id_rev", aRev, bNoOvl, "a", false, true, false, false, TieAll},
		{"a_fu_rev", aRev, bNoOvl, "a", false, false, true, false, TieAll},
		{"a_fd_rev", aRev, bNoOvl, "a", false, false, false, true, TieAll},
		{"b_iu_fwd", aFwd, bNoOvl, "b", true, false, false, false, TieAll},
		{"b_fd_fwd", aFwd, bNoOvl, "b", false, false, false, true, TieAll},
		{"b_iu_rev", aFwd, bNoOvlRev, "b", true, false, false, false, TieAll},
		{"b_fu_rev", aFwd, bNoOvlRev, "b", false, false, true, false, TieAll},
		{"tie_fu_all", aTie, bTie, "ref", false, false, true, false, TieAll},
		{"tie_fu_first", aTie, bTie, "ref", false, false, true, false, TieFirst},
		{"tie_fu_last", aTie, bTie, "ref", false, false, true, false, TieLast},
		{"tie_id_all", aTie, bTie, "ref", false, true, false, false, TieAll},
		{"tie_fd_all", aTie, bTie, "ref", false, false, false, true, TieAll},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			flags := []string{"-D", c.dmode}
			if c.iu {
				flags = append(flags, "-iu")
			}
			if c.id {
				flags = append(flags, "-id")
			}
			if c.fu {
				flags = append(flags, "-fu")
			}
			if c.fd {
				flags = append(flags, "-fd")
			}
			switch c.tie {
			case TieFirst:
				flags = append(flags, "-t", "first")
			case TieLast:
				flags = append(flags, "-t", "last")
			}
			want := runUpstreamClosest(t, bt, c.a, c.b, flags...)

			opts := Options{
				PrintDistance:    true,
				DistanceMode:     modeFor(c.dmode),
				TieBreak:         c.tie,
				IgnoreUpstream:   c.iu,
				IgnoreDownstream: c.id,
				ForceUpstream:    c.fu,
				ForceDownstream:  c.fd,
			}
			got := runOurClosest(t, c.a, c.b, opts)
			if !bytes.Equal(got, want) {
				t.Fatalf("mismatch %s (%v)\nupstream:\n%s\nours:\n%s", c.name, flags, want, got)
			}
		})
	}
}
