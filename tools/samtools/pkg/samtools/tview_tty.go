package samtools

// tview_tty.go holds the OS-independent surface of the interactive `-d C`
// viewer's terminal layer: the shared error returned when stdin/stdout are not
// TTYs, and the exported entry point the CLI calls. The actual raw-mode termios
// / ioctl code is platform-specific (tview_tty_linux.go for Linux; a stub in
// tview_tty_other.go for every other GOOS), wired in through
// openInteractiveTerminal.

import (
	"errors"
	"os"
)

// errTTYRequired is returned by openInteractiveTerminal when stdin or stdout is
// not a terminal. The interactive viewer cannot run over a pipe, and the
// text/HTML backends (`-d T` / `-d H`) are the pipeline-friendly alternatives.
var errTTYRequired = errors.New("samtools tview: interactive mode (-d C) requires a terminal; pipe to -d T or -d H instead")

// RunTviewInteractiveStdio runs the interactive `-d C` viewer against the real
// stdin/stdout. It opens a raw-mode terminal (failing with a clear message if
// either stream is not a TTY), guarantees the terminal is restored on return —
// including on panic — and drives the control loop. It is the single entry the
// CLI calls for the curses-equivalent mode.
func RunTviewInteractiveStdio(opts TviewOptions) error {
	term, err := openInteractiveTerminal(os.Stdin, os.Stdout)
	if err != nil {
		return err
	}
	// Restore the terminal no matter how the loop exits (normal return, error,
	// or panic): Close resets the saved termios.
	defer term.Close()
	return RunTviewInteractive(term, opts)
}
