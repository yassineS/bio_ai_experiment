//go:build linux

package bench

// maxRSSToKiB normalises rusage.Maxrss to KiB. On Linux ru_maxrss is already
// reported in KiB (kibibytes), so the value passes through unchanged — keeping
// the historical Linux behaviour byte-identical.
func maxRSSToKiB(maxrss int64) int64 { return int64(maxrss) }
