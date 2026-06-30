//go:build !linux

package realbench

// maxRSSToKiB normalises rusage.Maxrss to KiB. On darwin (and the BSDs)
// ru_maxrss is reported in BYTES, not KiB as on Linux, so dividing by 1024
// converts to KiB and the rss_kb field carries the same unit on every host
// (matching pipeline/bench's portable normalisation).
func maxRSSToKiB(maxrss int64) int64 { return maxrss / 1024 }
