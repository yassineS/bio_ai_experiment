package bedclosest

// Live-upstream parity tests mirroring the upstream
// reference_code/bedtools/test/closest/test-closest.sh suite. They build the
// real upstream `bedtools` binary once (via the shared sync.Once builder in
// gaps_parity_test.go) and compare its `closest` output byte-for-byte against
// this port. They t.Fatalf (never t.Skip) so a missing/unbuildable submodule
// is a hard failure, matching the project's parity-rig policy.
//
// Inputs are the vendored upstream fixtures under
// reference_code/bedtools/test/closest/. Each case names the upstream test it
// reproduces and the exact flags it exercises.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// closestFixture reads an upstream fixture file from the vendored test dir.
func closestFixture(t *testing.T, name string) []byte {
	t.Helper()
	root := gapsRepoRoot(t)
	path := filepath.Join(root, "reference_code", "bedtools", "test", "closest", name)
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read fixture %s: %v", name, err)
	}
	return data
}

// runUpstreamClosestFiles runs upstream `bedtools closest` with the given -a
// file, one or more -b files, and trailing flags, returning its stdout.
func runUpstreamClosestFiles(t *testing.T, bt string, aFile string, bFiles []string, flags ...string) []byte {
	t.Helper()
	root := gapsRepoRoot(t)
	base := filepath.Join(root, "reference_code", "bedtools", "test", "closest")
	args := []string{"closest", "-a", filepath.Join(base, aFile), "-b"}
	for _, b := range bFiles {
		args = append(args, filepath.Join(base, b))
	}
	args = append(args, flags...)
	return runUpstreamArgs(t, bt, args)
}

// runOurClosestFiles drives this port's ClosestMulti with the named fixtures.
func runOurClosestFiles(t *testing.T, aFile string, bFiles []string, opts Options) []byte {
	t.Helper()
	a := closestFixture(t, aFile)
	readers := make([]io.Reader, len(bFiles))
	for i, b := range bFiles {
		readers[i] = bytes.NewReader(closestFixture(t, b))
	}
	var out bytes.Buffer
	if _, err := ClosestMulti(bytes.NewReader(a), readers, &out, opts); err != nil {
		t.Fatalf("ClosestMulti failed: %v", err)
	}
	return out.Bytes()
}

// closestParityCase couples an upstream command (a/b fixtures + flags) with the
// equivalent Options for this port.
type closestParityCase struct {
	name  string
	a     string
	b     []string
	flags []string
	opts  Options
}

// TestParity_Closest_Suite reproduces the upstream test-closest.sh cases that
// exercise the core selection, distance, tie, multi-database, strand, and
// directional behaviours, asserting byte-for-byte parity with the live binary.
func TestParity_Closest_Suite(t *testing.T) {
	bt := upstreamBedtoolsGaps(t)

	cases := []closestParityCase{
		// t1..t4: 1bp/0bp off-by-one with -d (unsigned).
		{"t1", "a.bed", []string{"b.bed"}, []string{"-d"}, Options{ReportDistance: true}},
		{"t2", "b.bed", []string{"a.bed"}, []string{"-d"}, Options{ReportDistance: true}},
		{"t3", "a.bed", []string{"b-one-bp-closer.bed"}, []string{"-d"}, Options{ReportDistance: true}},
		{"t4", "b-one-bp-closer.bed", []string{"a.bed"}, []string{"-d"}, Options{ReportDistance: true}},
		// t5/t6: names with -d, then -N.
		{"t5", "a.names.bed", []string{"b.names.bed"}, []string{"-d"}, Options{ReportDistance: true}},
		{"t6", "a.names.bed", []string{"b.names.bed"}, []string{"-d", "-N"}, Options{ReportDistance: true, DifferentNames: true}},
		// t7/t8: -s with no same-strand B (null BED6 shape), -S overlapping hit.
		{"t7", "strand-test-a.bed", []string{"strand-test-b.bed"}, []string{"-s"}, Options{SameStrand: true}},
		{"t8", "strand-test-a.bed", []string{"strand-test-b.bed"}, []string{"-S"}, Options{OppositeStrand: true}},
		// t9..t11: tie modes, with BED3 null on chr2.
		{"t9", "close-a.bed", []string{"close-b.bed"}, nil, Options{TieBreak: TieAll}},
		{"t10", "close-a.bed", []string{"close-b.bed"}, []string{"-t", "first"}, Options{TieBreak: TieFirst}},
		{"t11", "close-a.bed", []string{"close-b.bed"}, []string{"-t", "last"}, Options{TieBreak: TieLast}},
		// t13..t15: multiple databases, -names, -filenames.
		{"t13", "mq1.bed", []string{"mdb1.bed", "mdb2.bed", "mdb3.bed"}, nil, Options{}},
		{"t14", "mq1.bed", []string{"mdb1.bed", "mdb2.bed", "mdb3.bed"}, []string{"-names", "a", "b", "c"},
			Options{DBLabels: []string{"a", "b", "c"}}},
		{"t15", "mq1.bed", []string{"mdb1.bed", "mdb2.bed", "mdb3.bed"}, []string{"-filenames"},
			Options{DBLabels: []string{
				absFixture("mdb1.bed"), absFixture("mdb2.bed"), absFixture("mdb3.bed")}}},
		// t16: -mdb all.
		{"t16", "mq1.bed", []string{"mdb1.bed", "mdb2.bed", "mdb3.bed"}, []string{"-mdb", "all"},
			Options{MultiDBMode: MultiDBAll}},
		// t17..t19: 2 dbs, tie modes with -mdb all.
		{"t17", "mq1.bed", []string{"mdb1.bed", "mdb2.bed"}, []string{"-t", "all"}, Options{TieBreak: TieAll}},
		{"t18", "mq1.bed", []string{"mdb1.bed", "mdb2.bed"}, []string{"-mdb", "all", "-t", "first"},
			Options{MultiDBMode: MultiDBAll, TieBreak: TieFirst}},
		{"t19", "mq1.bed", []string{"mdb1.bed", "mdb2.bed"}, []string{"-mdb", "all", "-t", "last"},
			Options{MultiDBMode: MultiDBAll, TieBreak: TieLast}},
		// t20/t21: strand on a single db.
		{"t20", "mq1.bed", []string{"mdb1.bed"}, []string{"-s"}, Options{SameStrand: true}},
		{"t21", "mq1.bed", []string{"mdb1.bed"}, []string{"-S"}, Options{OppositeStrand: true}},
		// t22/t23: 2 dbs, tie all, -mdb all, strand.
		{"t22", "mq1.bed", []string{"mdb1.bed", "mdb2.bed"}, []string{"-t", "all", "-mdb", "all", "-s"},
			Options{TieBreak: TieAll, MultiDBMode: MultiDBAll, SameStrand: true}},
		{"t23", "mq1.bed", []string{"mdb1.bed", "mdb2.bed"}, []string{"-t", "all", "-mdb", "all", "-S"},
			Options{TieBreak: TieAll, MultiDBMode: MultiDBAll, OppositeStrand: true}},
		// t48..t50: non-overlapping ties, with -s/-S.
		{"t48", "a2.bed", []string{"b2.bed"}, nil, Options{}},
		{"t49", "a2.bed", []string{"b2.bed"}, []string{"-s"}, Options{SameStrand: true}},
		{"t50", "a2.bed", []string{"b2.bed"}, []string{"-S"}, Options{OppositeStrand: true}},
		// t52/t53: stranded sweep cache correctness.
		{"t52", "strand-test-c.bed", []string{"strand-test-d.bed"}, []string{"-s"}, Options{SameStrand: true}},
		{"t53", "strand-test-c.bed", []string{"strand-test-d.bed"}, []string{"-S"}, Options{OppositeStrand: true}},
		// t54/t55: -iu/-id with -D ref.
		{"t54", "d.bed", []string{"d_iu.bed"}, []string{"-D", "ref", "-iu"},
			Options{ReportDistance: true, DistanceMode: DistanceSignedRef, IgnoreUpstream: true}},
		{"t55", "d.bed", []string{"d_id.bed"}, []string{"-D", "ref", "-id"},
			Options{ReportDistance: true, DistanceMode: DistanceSignedRef, IgnoreDownstream: true}},
		// t56..t58: single-db ties and -iu/-id.
		{"t56", "bug157_a.bed", []string{"bug157_b.bed"}, nil, Options{}},
		{"t57", "bug157_a.bed", []string{"bug157_b.bed"}, []string{"-D", "ref", "-iu"},
			Options{ReportDistance: true, DistanceMode: DistanceSignedRef, IgnoreUpstream: true}},
		{"t58", "bug157_a.bed", []string{"bug157_b.bed"}, []string{"-D", "ref", "-id"},
			Options{ReportDistance: true, DistanceMode: DistanceSignedRef, IgnoreDownstream: true}},
		// t60..t67: bug281 cache-purge for -s/-S.
		{"t60", "bug281_a.medium.bed", []string{"bug281_b.bed"}, []string{"-s"}, Options{SameStrand: true}},
		{"t61", "bug281_a.medium.bed", []string{"bug281_b.bed"}, []string{"-S"}, Options{OppositeStrand: true}},
		{"t62", "bug281_a.flip.medium.bed", []string{"bug281_b.bed"}, []string{"-s"}, Options{SameStrand: true}},
		{"t63", "bug281_a.flip.medium.bed", []string{"bug281_b.bed"}, []string{"-S"}, Options{OppositeStrand: true}},
		{"t64", "bug281_a.medium.bed", []string{"bug281_b.flip.bed"}, []string{"-s"}, Options{SameStrand: true}},
		{"t65", "bug281_a.medium.bed", []string{"bug281_b.flip.bed"}, []string{"-S"}, Options{OppositeStrand: true}},
		{"t66", "bug281_a.flip.medium.bed", []string{"bug281_b.flip.bed"}, []string{"-s"}, Options{SameStrand: true}},
		{"t67", "bug281_a.flip.medium.bed", []string{"bug281_b.flip.bed"}, []string{"-S"}, Options{OppositeStrand: true}},
		// t68: multidb null shape (BED4+).
		{"t68", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-names", "b", "c"},
			Options{DBLabels: []string{"b", "c"}}},
		// t69..t71: -mdb all/each with -filenames -d on mixed overlap.
		{"t69", "dmr.bed", []string{"islands.bed", "tfbs.bed", "shores.bed"},
			[]string{"-filenames", "-d", "-mdb", "all"},
			Options{ReportDistance: true, MultiDBMode: MultiDBAll, DBLabels: []string{
				absFixture("islands.bed"), absFixture("tfbs.bed"), absFixture("shores.bed")}}},
		{"t70", "dmr.bed", []string{"islands.2.bed", "tfbs.bed", "shores.bed"},
			[]string{"-filenames", "-d", "-mdb", "all"},
			Options{ReportDistance: true, MultiDBMode: MultiDBAll, DBLabels: []string{
				absFixture("islands.2.bed"), absFixture("tfbs.bed"), absFixture("shores.bed")}}},
		{"t71", "dmr.bed", []string{"islands.2.bed", "tfbs.bed", "shores.bed"},
			[]string{"-filenames", "-d", "-mdb", "each"},
			Options{ReportDistance: true, MultiDBMode: MultiDBEach, DBLabels: []string{
				absFixture("islands.2.bed"), absFixture("tfbs.bed"), absFixture("shores.bed")}}},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			want := runUpstreamClosestFiles(t, bt, c.a, c.b, c.flags...)
			got := runOurClosestFiles(t, c.a, c.b, c.opts)
			if !bytes.Equal(got, want) {
				t.Fatalf("mismatch %s (flags %v)\nupstream:\n%s\nours:\n%s", c.name, c.flags, want, got)
			}
		})
	}
}

// runUpstreamKClosest / runOurKClosest run the kclosest sub-suite fixtures
// (which live in the kclosest/ subdirectory).
func runUpstreamKClosestFiles(t *testing.T, bt, aFile string, bFiles []string, flags ...string) []byte {
	t.Helper()
	root := gapsRepoRoot(t)
	base := filepath.Join(root, "reference_code", "bedtools", "test", "closest", "kclosest")
	args := []string{"closest", "-a", filepath.Join(base, aFile), "-b"}
	for _, b := range bFiles {
		args = append(args, filepath.Join(base, b))
	}
	args = append(args, flags...)
	return runUpstreamArgs(t, bt, args)
}

func kFixture(t *testing.T, name string) []byte {
	t.Helper()
	root := gapsRepoRoot(t)
	data, err := os.ReadFile(filepath.Join(root, "reference_code", "bedtools", "test", "closest", "kclosest", name))
	if err != nil {
		t.Fatalf("read kfixture %s: %v", name, err)
	}
	return data
}

func runOurKClosestFiles(t *testing.T, aFile string, bFiles []string, opts Options) []byte {
	t.Helper()
	a := kFixture(t, aFile)
	readers := make([]io.Reader, len(bFiles))
	for i, b := range bFiles {
		readers[i] = bytes.NewReader(kFixture(t, b))
	}
	var out bytes.Buffer
	if _, err := ClosestMulti(bytes.NewReader(a), readers, &out, opts); err != nil {
		t.Fatalf("ClosestMulti failed: %v", err)
	}
	return out.Bytes()
}

// TestParity_Closest_KClosest reproduces the upstream kclosest sub-suite,
// exercising -k with the all/first/last tie modes, -io, and -iu/-id under the
// -D ref/a/b modes (the same matrix as kclosest/test-kclosest.sh).
func TestParity_Closest_KClosest(t *testing.T) {
	bt := upstreamBedtoolsGaps(t)

	sref := DistanceMode(DistanceSignedRef)
	type kc struct {
		name  string
		a     string
		b     []string
		flags []string
		opts  Options
	}
	cases := []kc{
		{"k3", "q1.bed", []string{"d1.bed"}, []string{"-k", "3"}, Options{KClosest: 3}},
		{"k5", "q1.bed", []string{"d1.bed"}, []string{"-k", "5"}, Options{KClosest: 5}},
		{"k6", "q1.bed", []string{"d1.bed"}, []string{"-k", "6"}, Options{KClosest: 6}},
		{"k7", "q1.bed", []string{"d1.bed"}, []string{"-k", "7"}, Options{KClosest: 7}},
		{"k4_first", "q1.bed", []string{"d1.bed"}, []string{"-k", "4", "-t", "first"}, Options{KClosest: 4, TieBreak: TieFirst}},
		{"k4_last", "q1.bed", []string{"d1.bed"}, []string{"-k", "4", "-t", "last"}, Options{KClosest: 4, TieBreak: TieLast}},
		{"d2_k4", "q1.bed", []string{"d2.bed"}, []string{"-k", "4"}, Options{KClosest: 4}},
		{"d2_k4_first", "q1.bed", []string{"d2.bed"}, []string{"-k", "4", "-t", "first"}, Options{KClosest: 4, TieBreak: TieFirst}},
		{"d2_k4_first_io", "q1.bed", []string{"d2.bed"}, []string{"-k", "4", "-t", "first", "-io"}, Options{KClosest: 4, TieBreak: TieFirst, IgnoreOverlaps: true}},
		{"d2_k7", "q1.bed", []string{"d2.bed"}, []string{"-k", "7"}, Options{KClosest: 7}},
		{"d2_k7_first", "q1.bed", []string{"d2.bed"}, []string{"-k", "7", "-t", "first"}, Options{KClosest: 7, TieBreak: TieFirst}},
		{"d2_k7_last", "q1.bed", []string{"d2.bed"}, []string{"-k", "7", "-t", "last"}, Options{KClosest: 7, TieBreak: TieLast}},
		{"d2_k7_id", "q1.bed", []string{"d2.bed"}, []string{"-k", "7", "-id", "-D", "ref"}, Options{KClosest: 7, IgnoreDownstream: true, ReportDistance: true, DistanceMode: sref}},
		{"d2_k7_id_first", "q1.bed", []string{"d2.bed"}, []string{"-k", "7", "-id", "-D", "ref", "-t", "first"}, Options{KClosest: 7, IgnoreDownstream: true, ReportDistance: true, DistanceMode: sref, TieBreak: TieFirst}},
		{"d3_k15_iu", "q1.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "ref", "-iu"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: sref, IgnoreUpstream: true}},
		{"d3_k15_io", "q1.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "ref", "-io"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: sref, IgnoreOverlaps: true}},
		{"d3_k15_id", "q1.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "ref", "-id"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: sref, IgnoreDownstream: true}},
		{"d3_k15_Da_iu", "q2.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "a", "-iu"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: DistanceA, IgnoreUpstream: true}},
		{"d3_k15_Da_id", "q2.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "a", "-id"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: DistanceA, IgnoreDownstream: true}},
		{"d3_k15_Db_iu", "q2.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "b", "-iu"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: DistanceB, IgnoreUpstream: true}},
		{"d3_k15_Db_id", "q2.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "b", "-id"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: DistanceB, IgnoreDownstream: true}},
		{"d3_k15_Da_s_iu", "q2.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "a", "-s", "-iu"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: DistanceA, SameStrand: true, IgnoreUpstream: true}},
		{"d3_k15_Da_s_id", "q2.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "a", "-s", "-id"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: DistanceA, SameStrand: true, IgnoreDownstream: true}},
		{"d3_k15_Db_s_iu", "q2.bed", []string{"d3.bed"}, []string{"-k", "15", "-D", "b", "-s", "-iu"}, Options{KClosest: 15, ReportDistance: true, DistanceMode: DistanceB, SameStrand: true, IgnoreUpstream: true}},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			want := runUpstreamKClosestFiles(t, bt, c.a, c.b, c.flags...)
			got := runOurKClosestFiles(t, c.a, c.b, c.opts)
			if !bytes.Equal(got, want) {
				t.Fatalf("mismatch %s (flags %v)\nupstream:\n%s\nours:\n%s", c.name, c.flags, want, got)
			}
		})
	}
}

// absFixture returns the absolute path upstream prints in its -filenames column
// (it echoes the path exactly as given on the command line). The test passes
// absolute paths to -b, so the label is the absolute path too.
func absFixture(name string) string {
	root, _ := os.Getwd()
	// gapsRepoRoot is the reliable anchor; getwd points at the package dir.
	// Walk up to the module root.
	dir := root
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			break
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			break
		}
		dir = parent
	}
	return filepath.Join(dir, "reference_code", "bedtools", "test", "closest", name)
}
