package mosdepth

import (
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

func TestParseQuantize(t *testing.T) {
	cases := []struct {
		spec    string
		want    []int
		wantErr bool
	}{
		{"", nil, false},
		{"0:1:4", []int{0, 1, 4}, false},
		{":1:4", []int{0, 1, 4}, false},
		{"1:4:", []int{1, 4, maxQuantizeBound}, false},
		{"4:1:0", []int{0, 1, 4}, false}, // sorted ascending
		{"5", nil, true},                 // needs >= 2 boundaries
		{"a:b", nil, true},
		{"1::4", nil, true}, // empty middle segment
	}
	for _, tc := range cases {
		got, err := ParseQuantize(tc.spec)
		if tc.wantErr {
			if err == nil {
				t.Errorf("ParseQuantize(%q): want error, got %v", tc.spec, got)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseQuantize(%q): unexpected error %v", tc.spec, err)
			continue
		}
		if !reflect.DeepEqual(got, tc.want) {
			t.Errorf("ParseQuantize(%q) = %v, want %v", tc.spec, got, tc.want)
		}
	}
}

func TestQuantizeLabels(t *testing.T) {
	// Default labels.
	got := quantizeLabels([]int{0, 1, 4})
	want := []string{"0:1", "1:4"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("labels = %v, want %v", got, want)
	}
	// Open-ended top bin renders as "N:inf".
	got = quantizeLabels([]int{1, 4, maxQuantizeBound})
	want = []string{"1:4", "4:inf"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("labels = %v, want %v", got, want)
	}
	// Environment override.
	t.Setenv("MOSDEPTH_Q0", "NONE")
	t.Setenv("MOSDEPTH_Q1", "LOW")
	got = quantizeLabels([]int{0, 1, 4})
	want = []string{"NONE", "LOW"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("labels with env = %v, want %v", got, want)
	}
}

func TestQuantizeBin(t *testing.T) {
	quants := []int{1, 4, 10}
	cases := []struct {
		d    int
		want int
	}{
		{0, -1},  // below quants[0]
		{1, 0},   // exactly the low boundary -> bin 0
		{3, 0},   // within [1,4)
		{4, 1},   // boundary -> bin 1
		{9, 1},   // within [4,10)
		{10, -1}, // == finite top boundary: no half-open bin covers it
		{11, -1}, // above quants[high]
	}
	for _, tc := range cases {
		if got := quantizeBin(tc.d, quants); got != tc.want {
			t.Errorf("quantizeBin(%d, %v) = %d, want %d", tc.d, quants, got, tc.want)
		}
	}

	// Open-ended top bin (trailing ':' appends the sentinel): any realistic
	// depth >= the last finite boundary lands in the open bin, not -1.
	open := []int{1, 4, maxQuantizeBound}
	if got := quantizeBin(4, open); got != 1 {
		t.Errorf("quantizeBin(4, open) = %d, want 1 (open top bin)", got)
	}
	if got := quantizeBin(1000000, open); got != 1 {
		t.Errorf("quantizeBin(1000000, open) = %d, want 1 (open top bin)", got)
	}
}

// TestAddFragment covers the full-fragment span arithmetic.
func TestAddFragment(t *testing.T) {
	a := newCovAccum(1000)
	// Read at POS=101 (1-based) -> 0-based start 100; mate at PNext=201 ->
	// 0-based 200; |TLEN|=150 -> fragment covers [100, 250).
	a.addFragment(&sam.Record{Pos: 101, PNext: 201, TLen: 150})
	var maxEnd int
	a.emitRuns(func(start, end int, depth int32) {
		if depth == 1 {
			if start != 100 {
				t.Errorf("fragment start = %d, want 100", start)
			}
			if end != 250 {
				t.Errorf("fragment end = %d, want 250", end)
			}
		}
		if end > maxEnd {
			maxEnd = end
		}
	})
	if maxEnd != 1000 {
		t.Errorf("trailing run end = %d, want refLen 1000", maxEnd)
	}
}
