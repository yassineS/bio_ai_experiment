package mosdepth

import "sort"

// stdSortEvents sorts a covEvent slice by pos ascending using the standard
// library. Hoisted into its own file so callers needing a different sort
// implementation can swap it without touching coverage.go.
func stdSortEvents(s []covEvent) {
	sort.Slice(s, func(i, j int) bool { return s[i].pos < s[j].pos })
}
