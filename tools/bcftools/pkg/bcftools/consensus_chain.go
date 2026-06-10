package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
)

// openChainFile creates the chain output file. Upstream writes the chain as
// plain (uncompressed) text via fopen, so this does the same.
func openChainFile(path string) (io.WriteCloser, error) {
	return os.Create(path)
}

// Chain accumulates the liftover blocks that map reference coordinates to
// consensus coordinates for a single chromosome/sequence, then renders them
// in UCSC chain format (https://genome.ucsc.edu/goldenPath/help/chain.html).
//
// It mirrors the chain_t structure and the init_chain / push_chain_gap /
// print_chain routines in upstream bcftools consensus.c. A new Chain is
// created per output sequence; chainID is a process-wide auto-incrementing
// identifier shared across sequences.
type Chain struct {
	// oriPos is the 0-based reference start position of the sequence
	// region (0 for whole-contig output). It is the chain block origin.
	oriPos int

	// blockLengths holds the length of each ungapped alignment block.
	// refGaps/altGaps hold the gap on the reference / alternate sequence
	// between this block and the next (slices are parallel).
	blockLengths []int
	refGaps      []int
	altGaps      []int

	// refLastBlockOri / altLastBlockOri track the reference / alternate
	// origin of the block currently being extended.
	refLastBlockOri int
	altLastBlockOri int
}

// NewChain allocates a Chain whose first block begins at the 0-based
// reference origin refOriPos. It mirrors upstream's init_chain.
func NewChain(refOriPos int) *Chain {
	return &Chain{
		oriPos:          refOriPos,
		refLastBlockOri: refOriPos,
		altLastBlockOri: refOriPos,
	}
}

// pushGap records one indel as a gap between two ungapped blocks. The
// arguments are 0-based: refStart/refLen describe the gap on the reference
// sequence and altStart/altLen the gap on the consensus (alternate)
// sequence. It mirrors upstream's push_chain_gap, including the back-to-back
// merge of variants that abut the previous block.
func (c *Chain) pushGap(refStart, refLen, altStart, altLen int) {
	num := len(c.blockLengths)
	if num > 0 && refStart <= c.refLastBlockOri {
		// This variant is back-to-back with the previous one: extend the
		// previous gap rather than opening a new (zero-length) block.
		c.refLastBlockOri = refStart + refLen
		c.altLastBlockOri = altStart + altLen
		c.refGaps[num-1] += refLen
		c.altGaps[num-1] += altLen
		return
	}
	// Close the current ungapped block and open a new gap.
	c.blockLengths = append(c.blockLengths, refStart-c.refLastBlockOri)
	c.refGaps = append(c.refGaps, refLen)
	c.altGaps = append(c.altGaps, altLen)
	c.refLastBlockOri = refStart + refLen
	c.altLastBlockOri = altStart + altLen
}

// print writes the chain in UCSC chain format to w for a sequence named
// chrom whose original reference length is faLength. It mirrors upstream's
// print_chain, including the trailing blank line. chainID is the running
// identifier that the caller pre-increments per sequence.
func (c *Chain) print(w io.Writer, chrom string, faLength, chainID int) error {
	refEndPos := faLength + c.oriPos
	lastBlockSize := refEndPos - c.refLastBlockOri
	altEndPos := c.altLastBlockOri + lastBlockSize

	score := 0
	for _, bl := range c.blockLengths {
		score += bl
	}
	score += lastBlockSize

	if _, err := fmt.Fprintf(w, "chain %d %s %d + %d %d %s %d + %d %d %d\n",
		score, chrom, refEndPos, c.oriPos, refEndPos,
		chrom, altEndPos, c.oriPos, altEndPos, chainID); err != nil {
		return err
	}
	for i := range c.blockLengths {
		if _, err := fmt.Fprintf(w, "%d %d %d\n", c.blockLengths[i], c.refGaps[i], c.altGaps[i]); err != nil {
			return err
		}
	}
	_, err := fmt.Fprintf(w, "%d\n\n", lastBlockSize)
	return err
}

// chainWriter wraps the chain output file with a buffered writer and a
// running chain identifier. The zero value is not usable; use newChainWriter.
type chainWriter struct {
	bw      *bufio.Writer
	closer  io.Closer
	chainID int
}

// writeChain renders one sequence's chain, pre-incrementing the running
// identifier exactly as upstream does (++args->chain_id).
func (cw *chainWriter) writeChain(c *Chain, chrom string, faLength int) error {
	cw.chainID++
	return c.print(cw.bw, chrom, faLength, cw.chainID)
}

// flush flushes the buffered writer.
func (cw *chainWriter) flush() error { return cw.bw.Flush() }
