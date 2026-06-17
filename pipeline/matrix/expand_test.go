package matrix

import "testing"

// TestExpandSingleFlags verifies the expander emits a baseline plus one entry
// per single flag value, and that boolean and value-bearing flags are handled.
func TestExpandSingleFlags(t *testing.T) {
	spec := ExpandSpec{
		Tool:     "samtools",
		Input:    InputBAM,
		BaseArgs: []string{"{bam}"},
		Flags: []Flag{
			{Name: "-c", Bool: true},
			{Name: "-q", Values: []string{"10", "30"}},
		},
	}
	got := spec.Expand()
	// baseline + (-c) + (-q 10) + (-q 30) = 4 entries.
	if len(got) != 4 {
		t.Fatalf("got %d entries, want 4: %+v", len(got), got)
	}
	if got[0].Name != "base" {
		t.Errorf("first entry name = %q, want base", got[0].Name)
	}
	// Every entry must carry the base args.
	for _, e := range got {
		if len(e.Args) == 0 || e.Args[0] != "{bam}" {
			t.Errorf("entry %q missing base arg: %v", e.Name, e.Args)
		}
	}
}

// TestExpandCombos verifies curated combos are appended verbatim and that the
// expander does NOT produce the power set.
func TestExpandCombos(t *testing.T) {
	spec := ExpandSpec{
		Tool: "x",
		Flags: []Flag{
			{Name: "-a", Bool: true},
			{Name: "-b", Bool: true},
			{Name: "-c", Bool: true},
		},
		Combos: []Combo{
			{Name: "ab", Flags: []string{"-a", "-b"}},
		},
	}
	got := spec.Expand()
	// baseline + 3 single flags + 1 combo = 5; power set would be 2^3 = 8.
	if len(got) != 5 {
		t.Fatalf("got %d entries, want 5 (no power set): %+v", len(got), names(got))
	}
	var sawCombo bool
	for _, e := range got {
		if e.Name == "combo_ab" {
			sawCombo = true
			if len(e.Args) != 2 || e.Args[0] != "-a" || e.Args[1] != "-b" {
				t.Errorf("combo args = %v, want [-a -b]", e.Args)
			}
		}
	}
	if !sawCombo {
		t.Errorf("combo_ab not found in %v", names(got))
	}
}

// TestCrossProduct verifies the bounded cross-product helper.
func TestCrossProduct(t *testing.T) {
	combos := CrossProduct(
		Flag{Name: "-x", Values: []string{"1", "2"}},
		Flag{Name: "-y", Values: []string{"a"}},
	)
	if len(combos) != 2 { // 2*1
		t.Fatalf("got %d combos, want 2: %+v", len(combos), combos)
	}
}

// TestRegistryFilter verifies tool filtering.
func TestRegistryFilter(t *testing.T) {
	r := &Registry{}
	r.Add(
		Entry{Tool: "a", Name: "n1"},
		Entry{Tool: "b", Name: "n2"},
	)
	if got := r.FilterTools(map[string]bool{"a": true}); len(got) != 1 || got[0].Tool != "a" {
		t.Fatalf("filter returned %+v", got)
	}
	if got := r.FilterTools(nil); len(got) != 2 {
		t.Fatalf("nil filter returned %d, want 2", len(got))
	}
}

func names(es []Entry) []string {
	out := make([]string, len(es))
	for i, e := range es {
		out[i] = e.Name
	}
	return out
}
