//go:build !linux

package samtools

// tview_tty_other.go is the non-Linux stub for the interactive `-d C` viewer's
// terminal layer. The raw-mode termios / TIOCGWINSZ ioctl code in
// tview_tty_linux.go is Linux-specific (it drives the kernel's TCGETS/TCSETS
// directly via syscall, with no third-party dependency), so on every other GOOS
// the interactive viewer is unavailable and openInteractiveTerminal reports a
// clear, build-green error. The pure-Go control loop, key map, and frame
// renderer (tview_interactive.go) still compile and are still unit-tested on
// all platforms.

import (
	"errors"
	"os"
)

// errInteractiveUnsupported is returned on non-Linux platforms where the
// raw-mode terminal layer is not implemented.
var errInteractiveUnsupported = errors.New("samtools tview: interactive mode (-d C) requires Linux; use -d T or -d H")

// openInteractiveTerminal is the non-Linux stub: it always reports that the
// interactive viewer is unsupported on this platform.
func openInteractiveTerminal(_, _ *os.File) (tvTerminal, error) {
	return nil, errInteractiveUnsupported
}
