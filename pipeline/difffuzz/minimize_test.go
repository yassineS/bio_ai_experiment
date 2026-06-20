package difffuzz

import (
	"bytes"
	"testing"
)

// TestMinimizeShrinksToTrigger feeds a large input where the divergence is
// triggered by the presence of a single marker line. The minimizer should
// shrink it down to (essentially) just that line.
func TestMinimizeShrinksToTrigger(t *testing.T) {
	var big []byte
	for i := 0; i < 200; i++ {
		big = append(big, []byte("filler line that does not matter\n")...)
	}
	// Inject the trigger somewhere in the middle.
	trigger := []byte("TRIGGER\n")
	idx := len(big) / 2
	withTrigger := append(append(append([]byte{}, big[:idx]...), trigger...), big[idx:]...)

	pred := func(b []byte) bool { return bytes.Contains(b, []byte("TRIGGER")) }

	got := Minimize(withTrigger, pred, 5000)
	if !pred(got) {
		t.Fatalf("minimized input no longer triggers: %q", got)
	}
	if len(got) >= len(withTrigger)/2 {
		t.Fatalf("minimization barely shrank: %d -> %d", len(withTrigger), len(got))
	}
	// It should be close to just the trigger token.
	if len(got) > len("TRIGGER")+8 {
		t.Logf("minimized to %d bytes: %q (trigger is %d)", len(got), got, len("TRIGGER"))
	}
}

// TestMinimizeByteLevel exercises the byte-granularity shrink on a single-line
// input where only a substring matters.
func TestMinimizeByteLevel(t *testing.T) {
	in := []byte("aaaaaaaaaaNEEDLEbbbbbbbbbb")
	pred := func(b []byte) bool { return bytes.Contains(b, []byte("NEEDLE")) }
	got := Minimize(in, pred, 5000)
	if !pred(got) {
		t.Fatalf("lost trigger: %q", got)
	}
	if !bytes.Equal(got, []byte("NEEDLE")) {
		t.Logf("minimized to %q (ideal NEEDLE)", got)
	}
	if len(got) > len("NEEDLE")+4 {
		t.Fatalf("byte minimization too weak: %d bytes %q", len(got), got)
	}
}

// TestMinimizeNoTriggerReturnsOriginal ensures Minimize returns the original
// input untouched when the predicate never fires on it.
func TestMinimizeNoTriggerReturnsOriginal(t *testing.T) {
	in := []byte("hello world\n")
	pred := func(b []byte) bool { return false }
	got := Minimize(in, pred, 100)
	if !bytes.Equal(got, in) {
		t.Fatalf("expected original returned, got %q", got)
	}
}

// TestMinimizeRespectsStepBudget ensures a tiny step budget terminates and
// still returns a predicate-satisfying result.
func TestMinimizeRespectsStepBudget(t *testing.T) {
	in := bytes.Repeat([]byte("x\n"), 500)
	in = append(in, []byte("KEY\n")...)
	pred := func(b []byte) bool { return bytes.Contains(b, []byte("KEY")) }
	got := Minimize(in, pred, 3) // absurdly small budget
	if !pred(got) {
		t.Fatalf("budget-limited minimize broke the trigger: %q", got)
	}
}
