//go:build linux

package realbench

// maxRSSToKiB normalises rusage.Maxrss to KiB. On Linux ru_maxrss is already
// reported in KiB (kibibytes), so the value passes through unchanged — keeping
// the Linux measurement host's numbers identical to realparity's.
func maxRSSToKiB(maxrss int64) int64 { return maxrss }
