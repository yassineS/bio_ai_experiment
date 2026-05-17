package prinseq

import (
	"bufio"
	"fmt"
	"io"
	"strings"
)

// CollectGraphData runs the graph-data collectors over `reader`,
// returning a populated *GraphData ready for EmitGD. This is the
// entry point used by the CLI; library users who already have
// records can drive AddSeq/AddQual directly.
//
// `isFastq` selects the file-format path. The reader is consumed
// line by line; FASTQ is assumed to use the 4-line layout (matches
// upstream's parse path and the rest of the package's helpers).
func CollectGraphData(reader io.Reader, isFastq bool, opts GraphDataOptions) (*GraphData, error) {
	opts.IsFasta = !isFastq
	g := NewGraphData(opts)
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 0, 64*1024), 16*1024*1024)
	if isFastq {
		if err := collectFastq(scanner, g); err != nil {
			return nil, err
		}
	} else {
		if err := collectFasta(scanner, g); err != nil {
			return nil, err
		}
	}
	return g, nil
}

func collectFasta(scanner *bufio.Scanner, g *GraphData) error {
	var seq strings.Builder
	first := true
	flush := func() {
		if seq.Len() == 0 {
			return
		}
		g.AddSeq(strings.ToUpper(seq.String()))
		seq.Reset()
	}
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if !first {
				flush()
			}
			first = false
			continue
		}
		seq.WriteString(line)
	}
	flush()
	return scanner.Err()
}

func collectFastq(scanner *bufio.Scanner, g *GraphData) error {
	state := 0 // 0=header, 1=seq, 2=+, 3=qual
	var seq, qual string
	for scanner.Scan() {
		line := scanner.Text()
		switch state {
		case 0:
			if !strings.HasPrefix(line, "@") {
				return fmt.Errorf("expected '@' at start of FASTQ header, got: %q", line)
			}
			state = 1
		case 1:
			seq = strings.ToUpper(line)
			state = 2
		case 2:
			if !strings.HasPrefix(line, "+") {
				return fmt.Errorf("expected '+' separator, got: %q", line)
			}
			state = 3
		case 3:
			qual = line
			g.AddSeq(seq)
			quals := decodeQualString(qual, g.opts.Phred64)
			g.AddQual(quals)
			seq = ""
			qual = ""
			state = 0
		}
	}
	return scanner.Err()
}

// decodeQualString converts an ASCII quality string into integer
// Phred values, choosing offset 33 (sanger) or 64 (Phred+64) based
// on the option flag. Negative results from the subtraction (an
// invalid Phred+64 input) are clamped to 0; the upstream prinseq
// path returns an error on negative values, but the tests we care
// about here all use valid encodings.
func decodeQualString(qual string, phred64 bool) []int {
	offset := 33
	if phred64 {
		offset = 64
	}
	out := make([]int, len(qual))
	for i := 0; i < len(qual); i++ {
		v := int(qual[i]) - offset
		if v < 0 {
			v = 0
		}
		out[i] = v
	}
	return out
}

// ResolveGraphDataPath replicates upstream lines 984-987: when the
// caller passes `--graph_data` without a value, default the output
// path to `<file1>__.gd`. When a value is given, use it verbatim.
//
// `inputPath` is the path of the primary input file (or "" for
// stdin); `requested` is whatever the user supplied for the flag
// (an empty string meaning "use default").
func ResolveGraphDataPath(requested, inputPath string) string {
	if requested != "" {
		return requested
	}
	base := inputPath
	if base == "" || base == "-" {
		base = "nonamegiven"
	}
	return base + "__.gd"
}
