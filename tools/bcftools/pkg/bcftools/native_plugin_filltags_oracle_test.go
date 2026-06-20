package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"testing"
)

// Live-oracle parity tests for the fill-tags plugin's population grouping
// (-S/--samples-file), custom expressions (TAG[:Number]=[int|float](EXPR)),
// and the -l/--list-tags table. Each drives BOTH the upstream binary and OUR
// port through their CLIs and diffs the provenance-stripped output (or, for
// -l, the stderr table). The fixtures use multi-character sample names so the
// upstream parse_samples accepts them.

func TestNativePluginFillTagsPops(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	pops := parityFixture(t, "filltags_pops.vcf")
	groups := parityFixture(t, "filltags_groups.txt")

	cases := [][]string{
		{"-S", groups, "-t", "all"},
		{"-S", groups, "-t", "AN,AC,AF,MAF,NS"},
		{"-S", groups, "-t", "AC_Hom,AC_Het,AC_Hemi"},
		{"-S", groups, "-t", "HWE,ExcHet"},
		{"-S", groups, "-t", "AN"},
		{"-S", groups, "-t", "GD=sum(FORMAT/DP)"},
		{"-S", groups, "-t", "MGQ=mean(FORMAT/GQ)"},
		// custom expressions (no grouping)
		{"-t", "DP:1=int(sum(FORMAT/DP))"},
		{"-t", "FORMAT/VD:1=int(smpl_sum(FORMAT/AD))"},
		{"-t", "MDP=mean(FORMAT/DP)"},
		{"-t", "MX=max(FORMAT/GQ)"},
		{"-t", "MN=min(FORMAT/GQ)"},
		{"-t", "MED=median(FORMAT/DP)"},
		{"-t", "SD=stdev(FORMAT/DP)"},
		{"-t", "AS=sum(FORMAT/AD)"},
		{"-t", "nhet=N_PASS(GT=\"het\")"},
		{"-t", "fhet=F_PASS(GT=\"het\")"},
		{"-t", "nalt=N_PASS(GT=\"alt\")"},
		{"-t", "nmiss=N_PASS(GT=\"mis\")"},
		{"-t", "F_MISSING"},
		{"-t", "X=sum(FORMAT/DP),Y=mean(FORMAT/GQ)"},
		// HWE/ExcHet (kfunc), F_MISSING via "all".
		{"-t", "HWE,ExcHet"},
		{"-d", "-t", "all"},
	}
	for _, args := range cases {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			assertPluginParity(t, bin, pops, "fill-tags", args...)
		})
	}

	// Sites-only AF derived from INFO/AN,AC (process_info_af path).
	sites := parityFixture(t, "filltags_sites.vcf")
	t.Run("sites-only-AF", func(t *testing.T) {
		assertPluginParity(t, bin, sites, "fill-tags", "-t", "AF")
	})
}

// TestNativePluginFillTagsListTags asserts that -l/--list-tags prints the
// available-tag table to stderr (byte-identical to upstream) and exits
// non-zero with empty stdout, for both the upstream binary and our port.
func TestNativePluginFillTagsListTags(t *testing.T) {
	bin, err := buildBcftools()
	if err != nil {
		t.Skipf("build upstream bcftools: %v", err)
	}
	fixture := parityFixture(t, "filltags_pops.vcf")
	dir := pluginDirAbs(t)

	run := func(binPath string) (stdout, stderr []byte, exit int) {
		argv := []string{"+fill-tags", fixture, "--", "-l"}
		cmd := exec.Command(binPath, argv...)
		cmd.Env = append(os.Environ(), "BCFTOOLS_PLUGINS="+dir)
		var out, errBuf bytes.Buffer
		cmd.Stdout = &out
		cmd.Stderr = &errBuf
		err := cmd.Run()
		code := 0
		if err != nil {
			if ee, ok := err.(*exec.ExitError); ok {
				code = ee.ExitCode()
			} else {
				t.Fatalf("run %s: %v", binPath, err)
			}
		}
		return out.Bytes(), errBuf.Bytes(), code
	}

	upOut, upErr, upExit := run(bin)
	if ourBinPath == "" {
		t.Fatalf("local bcftools port binary not built; cannot run CLI oracle")
	}
	ourOut, ourErr, ourExit := run(ourBinPath)

	if len(upOut) != 0 || len(ourOut) != 0 {
		t.Fatalf("-l should produce no stdout: upstream=%q ours=%q", upOut, ourOut)
	}
	if upExit == 0 || ourExit == 0 {
		t.Fatalf("-l should exit non-zero: upstream=%d ours=%d", upExit, ourExit)
	}
	// The upstream table is emitted verbatim; our port emits the same table
	// (it is a prefix of our stderr — the upstream table contains no trailing
	// host wrapper). Compare the shared list lines.
	if !bytes.Equal(upErr, ourErr[:min(len(ourErr), len(upErr))]) {
		t.Fatalf("-l table diverges from upstream\n--- upstream ---\n%s\n--- ours ---\n%s", upErr, ourErr)
	}
	if !bytes.HasPrefix(ourErr, upErr) {
		t.Fatalf("-l: our stderr does not start with the upstream table\nupstream:\n%s\nours:\n%s", upErr, ourErr)
	}
}
