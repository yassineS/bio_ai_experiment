//go:build linux

package samtools

// tview_tty_linux.go is the small, OS-specific terminal layer for the
// interactive `-d C` viewer on Linux. It is deliberately the ONLY file that
// touches the terminal: it puts stdin into raw (cbreak/no-echo) mode via the
// TCGETS/TCSETS termios ioctls, queries the window size via TIOCGWINSZ, reads
// key events (decoding ESC sequences for arrows/Home), and restores the
// original termios on Close. Everything above it (the key map, the state
// machine, the frame renderer) is pure Go in tview_interactive.go and is
// unit-tested without a TTY; this file is the tiny untestable core.
//
// No third-party deps: the termios struct and ioctl numbers are declared
// directly and driven through syscall.Syscall(SYS_IOCTL, ...), the same
// approach golang.org/x/term uses internally — but kept in-tree per this
// project's stdlib-only rule.

import (
	"io"
	"os"
	"syscall"
	"unsafe"
)

// Linux/amd64 ioctl request numbers for the termios and window-size calls.
// These are the generic asm-generic values used by all mainstream Linux
// architectures this project targets.
const (
	ioctlTCGETS     = 0x5401
	ioctlTCSETS     = 0x5402
	ioctlTIOCGWINSZ = 0x5413
)

// termios mirrors struct termios from <termios.h> (Linux). Only the fields the
// raw-mode transform touches are named precisely; the layout matches the kernel
// ABI so TCGETS/TCSETS round-trip correctly.
type termios struct {
	Iflag  uint32
	Oflag  uint32
	Cflag  uint32
	Lflag  uint32
	Line   uint8
	Cc     [19]uint8
	Ispeed uint32
	Ospeed uint32
}

// Local-flag (Lflag) and input-flag (Iflag) bits used to enter raw mode. Values
// are the Linux constants; we keep them local to avoid depending on the
// platform-tagged names in the syscall package.
const (
	tcECHO   = 0x00000008
	tcICANON = 0x00000002
	tcISIG   = 0x00000001
	tcIEXTEN = 0x00008000
	tcIXON   = 0x00000400
	tcICRNL  = 0x00000100
	tcVMIN   = 6
	tcVTIME  = 5
)

// winsize mirrors struct winsize for the TIOCGWINSZ ioctl.
type winsize struct {
	Row    uint16
	Col    uint16
	Xpixel uint16
	Ypixel uint16
}

// linuxTerminal is the raw-mode terminal used by the interactive viewer. It
// holds the input/output files and the saved original termios for restore.
type linuxTerminal struct {
	in     *os.File
	out    *os.File
	saved  termios
	rawSet bool
}

// tcGet reads the termios for fd via the TCGETS ioctl.
func tcGet(fd uintptr) (termios, error) {
	var t termios
	_, _, errno := syscall.Syscall(syscall.SYS_IOCTL, fd, ioctlTCGETS, uintptr(unsafe.Pointer(&t)))
	if errno != 0 {
		return t, errno
	}
	return t, nil
}

// tcSet writes t as the termios for fd via the TCSETS ioctl.
func tcSet(fd uintptr, t *termios) error {
	_, _, errno := syscall.Syscall(syscall.SYS_IOCTL, fd, ioctlTCSETS, uintptr(unsafe.Pointer(t)))
	if errno != 0 {
		return errno
	}
	return nil
}

// getWinsize queries the terminal window size for fd via TIOCGWINSZ. It returns
// (cols, rows) and ok=false when the ioctl fails (e.g. fd is not a TTY).
func getWinsize(fd uintptr) (cols, rows int, ok bool) {
	var ws winsize
	_, _, errno := syscall.Syscall(syscall.SYS_IOCTL, fd, ioctlTIOCGWINSZ, uintptr(unsafe.Pointer(&ws)))
	if errno != 0 || ws.Col == 0 || ws.Row == 0 {
		return 0, 0, false
	}
	return int(ws.Col), int(ws.Row), true
}

// isTTY reports whether f is a terminal (its TCGETS ioctl succeeds).
func isTTY(f *os.File) bool {
	if f == nil {
		return false
	}
	_, err := tcGet(f.Fd())
	return err == nil
}

// newLinuxTerminal puts in into raw (cbreak, no-echo) mode and returns a
// terminal whose Close restores the original settings. Both in and out must be
// TTYs; the caller checks that first.
func newLinuxTerminal(in, out *os.File) (*linuxTerminal, error) {
	orig, err := tcGet(in.Fd())
	if err != nil {
		return nil, err
	}
	raw := orig
	// cbreak/no-echo: disable canonical mode, echo, signals, extended input,
	// CR->NL translation, and software flow control; deliver each byte as it
	// arrives (VMIN=1, VTIME=0).
	raw.Lflag &^= tcECHO | tcICANON | tcISIG | tcIEXTEN
	raw.Iflag &^= tcIXON | tcICRNL
	raw.Cc[tcVMIN] = 1
	raw.Cc[tcVTIME] = 0
	if err := tcSet(in.Fd(), &raw); err != nil {
		return nil, err
	}
	return &linuxTerminal{in: in, out: out, saved: orig, rawSet: true}, nil
}

// Size returns the terminal (cols, rows), falling back to 80x40 when the ioctl
// fails (matching the non-interactive backend's 80-column default).
func (t *linuxTerminal) Size() (int, int) {
	if cols, rows, ok := getWinsize(t.out.Fd()); ok {
		return cols, rows
	}
	return DefaultTviewWidth, 40
}

// Write writes p to the terminal output.
func (t *linuxTerminal) Write(p []byte) error {
	_, err := t.out.Write(p)
	return err
}

// ReadKey reads one key event: a single byte for an ordinary key, or a full
// ANSI escape sequence (ESC '[' final / ESC 'O' final, or ESC [ 1 ~) for
// arrow/Home keys. A lone ESC (no following bytes) is returned as a one-byte
// event. It returns io.EOF at end of input.
func (t *linuxTerminal) ReadKey() ([]byte, error) {
	b, err := t.readByte()
	if err != nil {
		return nil, err
	}
	if b != keyEsc {
		return []byte{b}, nil
	}
	// Escape sequence: try to read '[' or 'O' then the final byte(s). Reads use
	// the same blocking fd; a bare ESC (user pressed Escape) is unusual to
	// distinguish without a timeout, so we attempt one more byte and treat a
	// read error as a lone ESC.
	b2, err := t.readByte()
	if err != nil {
		return []byte{keyEsc}, nil
	}
	if b2 != '[' && b2 != 'O' {
		// Not a recognised CSI/SS3 introducer; surface ESC alone.
		return []byte{keyEsc}, nil
	}
	b3, err := t.readByte()
	if err != nil {
		return []byte{keyEsc, b2}, nil
	}
	// ESC [ 1 ~ (Home on some terminals): read the trailing '~'.
	if b2 == '[' && b3 == '1' {
		b4, err := t.readByte()
		if err == nil && b4 == '~' {
			return []byte{keyEsc, b2, b3, b4}, nil
		}
		return []byte{keyEsc, b2, b3}, nil
	}
	return []byte{keyEsc, b2, b3}, nil
}

// readByte reads a single byte from the input, mapping a zero-length read to
// io.EOF.
func (t *linuxTerminal) readByte() (byte, error) {
	var buf [1]byte
	n, err := t.in.Read(buf[:])
	if err != nil {
		return 0, err
	}
	if n == 0 {
		return 0, io.EOF
	}
	return buf[0], nil
}

// ReadLine temporarily returns the terminal to a line-buffered echo mode, reads
// a line for the goto prompt, then restores raw mode. It echoes typed bytes and
// handles backspace so the prompt is usable.
func (t *linuxTerminal) ReadLine(prompt string) (string, error) {
	_ = t.Write([]byte("\r\n" + prompt))
	var line []byte
	for {
		b, err := t.readByte()
		if err == io.EOF {
			break
		}
		if err != nil {
			return string(line), err
		}
		switch b {
		case '\r', '\n':
			return string(line), nil
		case keyBackspace, keyCtrlH:
			if len(line) > 0 {
				line = line[:len(line)-1]
				_ = t.Write([]byte("\b \b"))
			}
		case keyCtrlC, keyEsc:
			return "", nil
		default:
			line = append(line, b)
			_ = t.Write([]byte{b}) // echo
		}
	}
	return string(line), nil
}

// Close restores the original termios settings.
func (t *linuxTerminal) Close() error {
	if !t.rawSet {
		return nil
	}
	t.rawSet = false
	return tcSet(t.in.Fd(), &t.saved)
}

// openInteractiveTerminal verifies in and out are TTYs and returns a raw-mode
// terminal. It returns a non-nil error (ttyRequiredErr) when either side is not
// a TTY, so the CLI can print a clear message and exit without garbling a pipe.
func openInteractiveTerminal(in, out *os.File) (tvTerminal, error) {
	if !isTTY(in) || !isTTY(out) {
		return nil, errTTYRequired
	}
	return newLinuxTerminal(in, out)
}
