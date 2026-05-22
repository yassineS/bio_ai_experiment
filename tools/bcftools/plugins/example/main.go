// Command example is a minimal reference bcftools plugin. It demonstrates the
// subprocess plugin protocol documented in docs/PLUGIN_PROTOCOL.md: it reads
// an uncompressed VCF from stdin, writes an uncompressed VCF to stdout, and
// uses its argv for plugin-specific options.
//
// The transform is deliberately trivial: it resets the FILTER column of every
// data record to "PASS" (i.e. it "clears" soft filters) and reports the number
// of records seen on stderr. With `--about` it prints a one-line description
// and exits, which lets `bcftools plugin -lv` probe it.
//
// It is shipped only as a test fixture and usage example; bcftools does not
// bundle a library of plugins.
package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

const about = "example: reset FILTER to PASS on every record (reference plugin)"

func main() {
	for _, a := range os.Args[1:] {
		if a == "--about" {
			fmt.Println(about)
			return
		}
	}

	in := bufio.NewReader(os.Stdin)
	out := bufio.NewWriter(os.Stdout)
	defer out.Flush()

	records := 0
	for {
		line, err := in.ReadString('\n')
		if len(line) > 0 {
			trimmed := strings.TrimRight(line, "\n")
			if strings.HasPrefix(trimmed, "#") || trimmed == "" {
				// Header or blank line: pass through unchanged.
				out.WriteString(line)
			} else {
				records++
				fields := strings.Split(trimmed, "\t")
				if len(fields) > 6 {
					fields[6] = "PASS" // FILTER is column 7 (0-based index 6).
				}
				out.WriteString(strings.Join(fields, "\t"))
				out.WriteByte('\n')
			}
		}
		if err != nil {
			break
		}
	}
	fmt.Fprintf(os.Stderr, "example: processed %d record(s)\n", records)
}
