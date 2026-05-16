package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

// TestRename_NoPrefixFasta verifies the default (no-prefix) behaviour:
// records are renumbered with bare integers as names.
func TestRename_NoPrefixFasta(t *testing.T) {
	in := ">seqA\nACGT\n>seqB\nGGG\n>seqC\nTT\n"
	want := ">1\nACGT\n>2\nGGG\n>3\nTT\n"
	var out bytes.Buffer
	if err := Rename(strings.NewReader(in), &out, ""); err != nil {
		t.Fatalf("Rename: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("Rename no-prefix:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestRename_PrefixFasta verifies the prefix variant.
func TestRename_PrefixFasta(t *testing.T) {
	in := ">a\nACGT\n>b\nGG\n"
	want := ">PX1\nACGT\n>PX2\nGG\n"
	var out bytes.Buffer
	if err := Rename(strings.NewReader(in), &out, "PX"); err != nil {
		t.Fatalf("Rename: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("Rename prefix:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestRename_PairsShareIndex covers the "pair detection" path: two
// adjacent records whose names compare equal modulo a "/<digit>"
// suffix share the same numeric index.
func TestRename_PairsShareIndex(t *testing.T) {
	in := strings.Join([]string{
		"@r/1", "ACGT", "+", "IIII",
		"@r/2", "TGCA", "+", "IIII",
		"@s/1", "AAAA", "+", "####",
		"@s/2", "TTTT", "+", "$$$$",
		"",
	}, "\n")
	want := strings.Join([]string{
		"@S_1", "ACGT", "+", "IIII",
		"@S_1", "TGCA", "+", "IIII",
		"@S_2", "AAAA", "+", "####",
		"@S_2", "TTTT", "+", "$$$$",
		"",
	}, "\n")
	var out bytes.Buffer
	if err := Rename(strings.NewReader(in), &out, "S_"); err != nil {
		t.Fatalf("Rename: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("Rename pairs:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestRename_SingletonsAdvance covers the path where adjacent records
// are NOT a pair: each gets its own index.
func TestRename_SingletonsAdvance(t *testing.T) {
	in := strings.Join([]string{
		">a", "AAA",
		">b", "BBB",
		">c", "CCC",
		"",
	}, "\n")
	want := strings.Join([]string{
		">X1", "AAA",
		">X2", "BBB",
		">X3", "CCC",
		"",
	}, "\n")
	var out bytes.Buffer
	if err := Rename(strings.NewReader(in), &out, "X"); err != nil {
		t.Fatalf("Rename: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("Rename singletons:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestRename_StickyComment_Quirk pins the upstream cpy_kstr leak
// (seqtk.c:1210) that we reproduce on purpose: a record without a
// comment that follows one with a comment inherits the previous
// comment text. See rename.go's header comment for the full rationale.
func TestRename_StickyComment_Quirk(t *testing.T) {
	in := strings.Join([]string{
		">a", "AAAA",
		">b with-comment", "GGGG",
		">c", "TTTT",
		"",
	}, "\n")
	// "a" is emitted with empty comment (sticky="" so no leak).
	// "b" is emitted with its own "with-comment" and updates sticky.
	// "c" is then emitted carrying the leaked "with-comment" comment.
	want := strings.Join([]string{
		">1", "AAAA",
		">2 with-comment", "GGGG",
		">3 with-comment", "TTTT",
		"",
	}, "\n")
	var out bytes.Buffer
	if err := Rename(strings.NewReader(in), &out, ""); err != nil {
		t.Fatalf("Rename: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("Rename sticky-comment quirk:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestRename_EmptyInput verifies an empty stream produces no output.
func TestRename_EmptyInput(t *testing.T) {
	var out bytes.Buffer
	if err := Rename(strings.NewReader(""), &out, "X"); err != nil {
		t.Fatalf("Rename empty: %v", err)
	}
	if out.Len() != 0 {
		t.Errorf("Rename empty produced %q", out.String())
	}
}

// TestCommentOf is a unit test for the description-splitting helper.
func TestCommentOf(t *testing.T) {
	cases := []struct {
		desc, id, want string
	}{
		{"seqA", "seqA", ""},
		{"seqA description", "seqA", "description"},
		{"seqA  two-space then text", "seqA", "two-space then text"},
		{"seqA\tdesc", "seqA", "desc"},
		{"seqA  ", "seqA", ""},
		{"", "", ""},
	}
	for _, tc := range cases {
		if got := commentOf(tc.desc, tc.id); got != tc.want {
			t.Errorf("commentOf(%q, %q) = %q, want %q", tc.desc, tc.id, got, tc.want)
		}
	}
}
