package samtools

import (
	"bytes"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestViewCRAMToBAMRawAuxByteIdentical drives the actual samtools-view CRAM→BAM
// path and asserts its BAM output is byte-identical to an oracle that decodes
// the same CRAM eagerly (rec.Aux) and re-encodes through the BAM writer. View
// engages the memory-lean raw-aux passthrough (rec.RawAux) for this
// configuration — BAM output, no aux-touching filter — so equality proves the
// CLI wiring is byte-exact, not just the cram package in isolation.
func TestViewCRAMToBAMRawAuxByteIdentical(t *testing.T) {
	data := buildMultiContainerCRAMBytes(t, 25000)

	// Oracle: eager decode → BAM, with no raw-aux passthrough.
	oracle := func() []byte {
		rr, err := cram.NewRecordReader(bytes.NewReader(data))
		if err != nil {
			t.Fatalf("NewRecordReader: %v", err)
		}
		var out bytes.Buffer
		bw := sam.NewBAMWriter(&out)
		if err := bw.WriteHeader(rr.Header()); err != nil {
			t.Fatalf("WriteHeader: %v", err)
		}
		for {
			rec, rerr := rr.Read()
			if rerr != nil {
				break
			}
			if err := bw.Write(rec); err != nil {
				t.Fatalf("BAM Write: %v", err)
			}
		}
		if err := bw.Close(); err != nil {
			t.Fatalf("BAM Close: %v", err)
		}
		return out.Bytes()
	}()

	for _, threads := range []int{0, 4} {
		var out bytes.Buffer
		opts := ViewOptions{OutputBAM: true, WithHeader: true, Threads: threads}
		if _, err := View(bytes.NewReader(data), &out, opts); err != nil {
			t.Fatalf("View(threads=%d): %v", threads, err)
		}
		if !bytes.Equal(out.Bytes(), oracle) {
			t.Fatalf("view -b CRAM→BAM (threads=%d, %d B) differs from eager oracle (%d B)",
				threads, out.Len(), len(oracle))
		}
	}
}

// TestRawAuxBAMSinkEligiblePredicate pins the view opt-in predicate: the
// passthrough engages only for a BAM sink with no aux-touching filter, and stays
// OFF for SAM/CRAM output, count mode, and RG / -d / -D filters.
func TestRawAuxBAMSinkEligiblePredicate(t *testing.T) {
	cases := []struct {
		name string
		opts ViewOptions
		want bool
	}{
		{"bam plain", ViewOptions{OutputBAM: true}, true},
		{"bam uncompressed", ViewOptions{Uncompressed: true}, true},
		{"sam text", ViewOptions{}, false},
		{"cram out", ViewOptions{OutputCRAM: true}, false},
		{"bam but cram wins", ViewOptions{OutputBAM: true, OutputCRAM: true}, false},
		{"count", ViewOptions{OutputBAM: true, Count: true}, false},
		{"rg filter", ViewOptions{OutputBAM: true, ReadGroup: "rg1"}, false},
		{"rg set filter", ViewOptions{OutputBAM: true, ReadGroupSet: map[string]struct{}{"rg1": {}}}, false},
		{"tag filter", ViewOptions{OutputBAM: true, TagFilters: []TagFilter{{Tag: "NH", ExistsOnly: true}}}, false},
		// Non-aux filters do not block the passthrough.
		{"flag filter ok", ViewOptions{OutputBAM: true, ExcludeFlags: 0x4}, true},
		{"mapq filter ok", ViewOptions{OutputBAM: true, MinMAPQ: 20}, true},
		{"subsample ok", ViewOptions{OutputBAM: true, Subsample: 0.5}, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if got := rawAuxBAMSinkEligible(&c.opts); got != c.want {
				t.Errorf("rawAuxBAMSinkEligible(%s) = %v, want %v", c.name, got, c.want)
			}
		})
	}
}
