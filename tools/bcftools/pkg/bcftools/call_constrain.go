// Package bcftools — call_constrain.go.
//
// Foundation for `bcftools call -C alleles -T sites.tsv`: the
// constrained-alleles path that restricts the multiallelic caller's EM
// search to a user-supplied allele set per record.
//
// Faithful port of mcall_constrain_alleles in
// reference_code/bcftools/mcall.c (~lines 1274-1413) plus the
// targets-file parser that vcfcall.c builds via regidx (tgt_parse,
// ~lines 415-459).
//
// The format of the sites file is one record per line:
//
//	CHROM\tPOS\tREF,ALT[,ALT...]
//
// Lines beginning with '#' are skipped as comments. Both REF and at
// least one ALT must be present; upstream's tgt_parse rejects shorter
// rows with the message wired in errConstrainShortRow.
//
// The runtime contract:
//
//  1. Only records whose (CHROM, POS) match a sites entry are passed
//     through the caller. Records that don't match are dropped — the
//     "missed line" insertion behaviour gated by upstream's `-i` flag
//     is intentionally not implemented (it is itself gated by `-C
//     alleles` and is an optional add-on; the residual tracker in
//     docs/PARITY_ROADMAP.md mentions it explicitly).
//
//  2. For matching records the REF/ALT set is rewritten to the
//     sites-file values. New alleles that don't appear in the
//     mpileup record borrow their PLs from the unseen-allele (<*>)
//     column; alleles that do appear keep their existing PLs through
//     the pl_map.
//
//  3. INFO/QS is reprojected onto the new allele indexes (with zeros
//     filling alleles that weren't present), and any FORMAT/AD field
//     is re-indexed identically. The trimmed record is then handed to
//     mcallSite() for normal EM / GT calling.
//
// Companion CLI wiring (call.go CallOptions.Constrain / ConstrainSites
// and main.go runCall -C flag) is the follow-up slice — this file
// supplies the parser + projection helpers in standalone form.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// errConstrainShortRow mirrors upstream's exact message text so error
// channels stay parity-comparable when callers grep for it.
const errConstrainShortRow = "Unable to parse the -T file; expected CHROM\\tPOS\\tREF,ALT with -C alleles"

// CallConstrain selects the user-supplied-alleles constraint for the
// multiallelic caller. CallConstrainNone is the no-constraint default.
type CallConstrain int

const (
	CallConstrainNone CallConstrain = iota
	CallConstrainAlleles
	CallConstrainTrio
)

type constrainSite struct {
	chrom   string
	pos     int
	alleles []string
	used    bool // for `-i/--insert-missed`: flips true once mpileup
	// has surfaced this site (so the flush-missed
	// pass at end-of-stream doesn't re-emit it).
}

type ConstrainAlleles struct {
	byKey map[string]*constrainSite
	// order tracks insertion order so the -i flush can emit missed
	// sites in the same order they appeared in the input file.
	order []*constrainSite
}

func LoadConstrainAlleles(path string) (*ConstrainAlleles, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("bcftools call -C alleles: %w", err)
	}
	defer in.Close()
	return parseConstrainAlleles(in)
}

func parseConstrainAlleles(r io.Reader) (*ConstrainAlleles, error) {
	out := &ConstrainAlleles{byKey: map[string]*constrainSite{}}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 4*1024*1024)
	for sc.Scan() {
		line := sc.Text()
		if line == "" || line[0] == '#' {
			continue
		}
		fields := splitWhitespace3(line)
		if len(fields) < 3 {
			return nil, fmt.Errorf("bcftools call -C alleles: %s but found instead:\n\t%s", errConstrainShortRow, line)
		}
		pos, err := strconv.Atoi(fields[1])
		if err != nil || pos == 0 {
			return nil, fmt.Errorf("bcftools call -C alleles: could not parse tab line, expected 1-based coordinate: %s", line)
		}
		alleleField := strings.SplitN(fields[2], "\t", 2)[0]
		alleleField = strings.SplitN(alleleField, " ", 2)[0]
		alleles := strings.Split(alleleField, ",")
		if len(alleles) < 2 {
			return nil, fmt.Errorf("bcftools call -C alleles: %s but found instead:\n\t%s", errConstrainShortRow, line)
		}
		key := fields[0] + "\x00" + fields[1]
		if _, exists := out.byKey[key]; exists {
			continue
		}
		site := &constrainSite{chrom: fields[0], pos: pos, alleles: alleles}
		out.byKey[key] = site
		out.order = append(out.order, site)
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("bcftools call -C alleles: %w", err)
	}
	return out, nil
}

func splitWhitespace3(s string) []string {
	out := make([]string, 0, 3)
	i := 0
	for len(out) < 2 {
		for i < len(s) && (s[i] == ' ' || s[i] == '\t') {
			i++
		}
		if i >= len(s) {
			break
		}
		start := i
		for i < len(s) && s[i] != ' ' && s[i] != '\t' {
			i++
		}
		out = append(out, s[start:i])
	}
	for i < len(s) && (s[i] == ' ' || s[i] == '\t') {
		i++
	}
	if i < len(s) {
		out = append(out, s[i:])
	}
	return out
}

func (c *ConstrainAlleles) Lookup(chrom string, pos int) *constrainSite {
	if c == nil {
		return nil
	}
	return c.byKey[chrom+"\x00"+strconv.Itoa(pos)]
}

func applyConstrainAlleles(v *vcf.Variant, site *constrainSite, stderrW io.Writer) (ok bool, fatal bool) {
	if site == nil {
		return false, false
	}
	origAlts := v.Alt
	if len(origAlts) == 1 && origAlts[0] == "." {
		origAlts = nil
	}
	oriAlleles := make([]string, 0, 1+len(origAlts))
	oriAlleles = append(oriAlleles, v.Ref)
	oriAlleles = append(oriAlleles, origAlts...)
	nOri := len(oriAlleles)

	oriUnseen := -1
	for i := 1; i < nOri; i++ {
		if oriAlleles[i] == "<*>" || oriAlleles[i] == "<X>" {
			oriUnseen = i
			break
		}
	}

	if site.alleles[0] != v.Ref {
		if stderrW != nil {
			fmt.Fprintf(stderrW, "The reference alleles are not compatible at %s:%d .. %s vs %s\n",
				v.Chrom, v.Pos, site.alleles[0], v.Ref)
		}
		return false, true
	}

	nNew := len(site.alleles)
	hasNew := false
	alsMap := make([]int, 0, nNew+1)
	newAlleles := make([]string, 0, nNew+1)

	alsMap = append(alsMap, 0)
	newAlleles = append(newAlleles, site.alleles[0])

	for i := 1; i < nNew; i++ {
		newAlleles = append(newAlleles, site.alleles[i])
		found := -1
		for j := 1; j < nOri; j++ {
			if oriAlleles[j] == site.alleles[i] {
				found = j
				break
			}
		}
		if found >= 0 {
			alsMap = append(alsMap, found)
		} else {
			if oriUnseen >= 0 {
				alsMap = append(alsMap, oriUnseen)
			} else {
				alsMap = append(alsMap, nOri-1)
			}
			hasNew = true
		}
	}
	if oriUnseen >= 0 {
		alsMap = append(alsMap, oriUnseen)
		newAlleles = append(newAlleles, oriAlleles[oriUnseen])
	}

	nalsNew := len(newAlleles)

	if !hasNew && nalsNew == nOri {
		identity := true
		for i := 0; i < nOri; i++ {
			if alsMap[i] != i || oriAlleles[i] != newAlleles[i] {
				identity = false
				break
			}
		}
		if identity {
			return true, false
		}
	}

	v.Ref = newAlleles[0]
	if nalsNew > 1 {
		v.Alt = append([]string(nil), newAlleles[1:]...)
	} else {
		v.Alt = []string{"."}
	}

	nplsNew := nalsNew * (nalsNew + 1) / 2
	plMap := make([]int, nplsNew)
	k := 0
	for i := 0; i < nalsNew; i++ {
		for j := 0; j <= i; j++ {
			a, b := alsMap[i], alsMap[j]
			if a > b {
				plMap[k] = a*(a+1)/2 + b
			} else {
				plMap[k] = b*(b+1)/2 + a
			}
			k++
		}
	}

	nplsOri := nOri * (nOri + 1) / 2
	for i := range v.Samples {
		s := &v.Samples[i]
		if s.Data == nil {
			continue
		}
		raw, hasPL := s.Data["PL"]
		if !hasPL {
			continue
		}
		oriPL, ok := decodePLInts(raw, nplsOri)
		if !ok {
			continue
		}
		for len(oriPL) < nplsOri {
			oriPL = append(oriPL, plMissing)
		}
		newPL := make([]string, nplsNew)
		for kk := 0; kk < nplsNew; kk++ {
			idx := plMap[kk]
			val := plMissing
			if idx < len(oriPL) {
				val = oriPL[idx]
			}
			if val == plMissing && oriUnseen >= 0 {
				newIa, newIb := pairFromPLIndex(kk)
				oa := alsMap[newIa]
				ob := alsMap[newIb]
				kOri := bcfAlleles2Gt(oa, oriUnseen)
				if kOri >= len(oriPL) || oriPL[kOri] == plMissing {
					kOri = bcfAlleles2Gt(ob, oriUnseen)
				}
				if kOri >= len(oriPL) || oriPL[kOri] == plMissing {
					kOri = bcfAlleles2Gt(oriUnseen, oriUnseen)
				}
				if kOri < len(oriPL) && oriPL[kOri] != plMissing {
					val = oriPL[kOri]
				}
			}
			if val == plMissing {
				newPL[kk] = "."
			} else {
				newPL[kk] = strconv.Itoa(val)
			}
		}
		s.Data["PL"] = strings.Join(newPL, ",")
	}

	if qsRaw, ok := v.Info["QS"]; ok {
		qsVals := parseFloatList(qsRaw)
		newQS := make([]string, nalsNew)
		for i := 0; i < nalsNew; i++ {
			oi := alsMap[i]
			if oi < len(qsVals) {
				newQS[i] = formatFloat32G(qsVals[oi])
			} else {
				newQS[i] = "0"
			}
		}
		v.Info["QS"] = strings.Join(newQS, ",")
	}

	for i := range v.Samples {
		s := &v.Samples[i]
		if s.Data == nil {
			continue
		}
		adRaw, ok := s.Data["AD"]
		if !ok {
			continue
		}
		adVals := parseIntList(adRaw)
		newAD := make([]string, nalsNew)
		for k := 0; k < nalsNew; k++ {
			oi := alsMap[k]
			if oi < len(adVals) {
				newAD[k] = strconv.Itoa(adVals[oi])
			} else {
				newAD[k] = "0"
			}
		}
		s.Data["AD"] = strings.Join(newAD, ",")
	}

	return true, false
}

func pairFromPLIndex(k int) (a, b int) {
	a = 0
	for (a+1)*(a+2)/2 <= k {
		a++
	}
	b = k - a*(a+1)/2
	return a, b
}

func bcfAlleles2Gt(a, b int) int {
	if a < b {
		a, b = b, a
	}
	return a*(a+1)/2 + b
}
