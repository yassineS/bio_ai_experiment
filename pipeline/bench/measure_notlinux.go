//go:build !linux

package bench

// maxRSSToKiB normalises rusage.Maxrss to KiB. On darwin (and the BSDs)
// ru_maxrss is reported in BYTES, not KiB as on Linux — e.g. a 100 MB
// allocation reports ~119,881,728. Dividing by 1024 converts to KiB so the
// shared mb() (KiB/1024) and the manuscript figures see the same unit as on a
// Linux measurement host. This non-Linux file covers darwin and the BSDs, all
// of which report ru_maxrss in bytes; Linux (the KiB pass-through) lives in
// measure_linux.go.
func maxRSSToKiB(maxrss int64) int64 { return int64(maxrss) / 1024 }
