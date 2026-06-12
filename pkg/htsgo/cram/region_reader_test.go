package cram

import (
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// multiContainerCRAM encodes the 569-read, coordinate-sorted mpileup.1.sam
// fixture (single reference "17", with a matching reference FASTA) into a
// CRAM v3 file whose slices are capped at seqsPerSlice records, forcing
// several containers. It returns the CRAM path and the reference FASTA path.
// It uses the vendored upstream samtools binary as the encoder so the
// RegionReader is exercised against a real samtools-written CRAM; a build
// failure is a hard error per the project's testing rules.
func multiContainerCRAM(t *testing.T, seqsPerSlice int) (cramPath, refPath string) {
	t.Helper()
	samtools := upstreamSamtoolsCram(t)
	srcSAM := filepath.Join(samtoolsTestDir, "dat/mpileup.1.sam")
	refPath = filepath.Join(samtoolsTestDir, "dat/mpileup.ref.fa")
	for _, p := range []string{srcSAM, refPath} {
		if _, err := os.Stat(p); err != nil {
			t.Fatalf("fixture %s missing: %v", p, err)
		}
	}
	cramPath = filepath.Join(t.TempDir(), "multi.cram")
	runUpstream(t, samtools, "view", "-C", "-T", refPath,
		"--output-fmt-option", "seqs_per_slice="+strconv.Itoa(seqsPerSlice),
		"-o", cramPath, srcSAM)
	return cramPath, refPath
}

// oracleRecords decodes the whole CRAM with the sequential RecordReader,
// attaching the reference FASTA, and returns every record. It is the
// self-consistent oracle the RegionReader is checked against: filtering this
// full set in Go must reproduce exactly what a seek-based Query returns.
func oracleRecords(t *testing.T, cramPath, refPath string) (*sam.Header, []*sam.Record) {
	t.Helper()
	rr, err := OpenRecords(cramPath)
	if err != nil {
		t.Fatalf("OpenRecords: %v", err)
	}
	defer rr.Close()
	if err := rr.SetReferenceFASTA(refPath); err != nil {
		t.Fatalf("SetReferenceFASTA: %v", err)
	}
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	return rr.Header(), recs
}

// filterOracle returns the oracle records overlapping reg, applying the same
// coordinate-region rule the RegionReader implements.
func filterOracle(recs []*sam.Record, reg region.ResolvedRegion, regionName string) []*sam.Record {
	var out []*sam.Record
	for _, rec := range recs {
		if regionOverlap(rec, reg, regionName) {
			out = append(out, rec)
		}
	}
	return out
}

// recordKey is a stable identity for a reconstructed record used to compare
// the seek-based and oracle result sets order-independently.
func recordKey(r *sam.Record) string {
	return r.QName + "\x00" + strconv.Itoa(int(r.Flag)) + "\x00" + r.RName + "\x00" +
		strconv.Itoa(int(r.Pos)) + "\x00" + r.Cigar.String() + "\x00" + r.Seq
}

func keysOf(recs []*sam.Record) []string {
	keys := make([]string, len(recs))
	for i, r := range recs {
		keys[i] = recordKey(r)
	}
	sort.Strings(keys)
	return keys
}

// TestRegionReaderMatchesOracle proves the seek-based RegionReader returns
// exactly the records overlapping a region, cross-checked against a
// whole-file RecordReader decode filtered in Go — a self-consistent oracle
// that needs no network. It sweeps several regions, including an empty
// region and regions that span multiple containers, plus an open-ended
// region covering the whole reference.
func TestRegionReaderMatchesOracle(t *testing.T) {
	cramPath, refPath := multiContainerCRAM(t, 50)

	craiPath := filepath.Join(filepath.Dir(cramPath), "multi.cram.crai")
	if err := CreateCRAI(cramPath, craiPath); err != nil {
		t.Fatalf("CreateCRAI: %v", err)
	}
	idx, err := OpenCRAI(craiPath)
	if err != nil {
		t.Fatalf("OpenCRAI: %v", err)
	}
	if len(idx.Entries) < 2 {
		t.Fatalf("expected a multi-container fixture, got %d crai entries", len(idx.Entries))
	}
	distinct := map[int64]struct{}{}
	for _, e := range idx.Entries {
		distinct[e.ContainerOffset] = struct{}{}
	}
	if len(distinct) < 2 {
		t.Fatalf("expected multiple containers, got %d distinct offsets", len(distinct))
	}
	t.Logf("fixture has %d crai entries across %d containers", len(idx.Entries), len(distinct))

	oracleHdr, oracleRecs := oracleRecords(t, cramPath, refPath)

	regions := []string{
		"17:1-200",       // first container only.
		"17:1000-2000",   // spans several middle containers.
		"17:1",           // single base near the start.
		"17",             // whole reference (open-ended).
		"17:4000-4200",   // tail of the reference.
		"17:99000-99999", // empty: beyond every read.
	}

	for _, spec := range regions {
		t.Run(spec, func(t *testing.T) {
			f, err := os.Open(cramPath)
			if err != nil {
				t.Fatalf("open: %v", err)
			}
			defer f.Close()
			rg, err := NewRegionReader(f, idx)
			if err != nil {
				t.Fatalf("NewRegionReader: %v", err)
			}
			if err := rg.SetReferenceFASTA(refPath); err != nil {
				t.Fatalf("SetReferenceFASTA: %v", err)
			}
			defer rg.Close()

			resolved, _, err := region.ResolveRegions([]string{spec},
				func(n string) int { return rg.Header().RefIndex(n) })
			if err != nil {
				t.Fatalf("ResolveRegions: %v", err)
			}
			if len(resolved) != 1 {
				t.Fatalf("region %q resolved to %d regions, want 1", spec, len(resolved))
			}
			got, err := rg.Query(resolved[0])
			if err != nil {
				t.Fatalf("Query: %v", err)
			}

			regionName := oracleHdr.Refs[resolved[0].RefID].Name
			want := filterOracle(oracleRecs, resolved[0], regionName)

			gk, wk := keysOf(got), keysOf(want)
			if len(gk) != len(wk) {
				t.Fatalf("region %q: Query returned %d records, oracle %d", spec, len(gk), len(wk))
			}
			for i := range gk {
				if gk[i] != wk[i] {
					t.Fatalf("region %q record %d mismatch:\n got=%q\nwant=%q", spec, i, gk[i], wk[i])
				}
			}
		})
	}
}

// runUpstream runs the upstream samtools binary with args, failing the test
// (never skipping) on error.
func runUpstream(t *testing.T, bin string, args ...string) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream samtools %v: %v\n%s", args, err, out)
	}
}
