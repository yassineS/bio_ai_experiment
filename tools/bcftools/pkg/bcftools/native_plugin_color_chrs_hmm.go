// HMM machinery and emission/segmentation logic for the color-chrs plugin,
// porting the relevant parts of htslib/HMM.c (hmm_init, hmm_set_tprob,
// _set_tprob, hmm_run_viterbi) and color-chrs.c (init_hmm_*, set_observed_prob_*,
// flush_viterbi). Only the Viterbi path is needed; fwd/bwd and Baum-Welch are
// not used by this plugin.
//
// Matrices are stored row-major as flat []float64 with MAT(m,n,i,j)=m[n*i+j],
// representing the transition j->i, exactly as the C MAT macro and the C code's
// loop order, so the floating-point accumulation is bit-identical.
package bcftools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// ccNTProb is the number of precomputed transition matrices, matching the
// hmm_init(..., 10000) call in color-chrs.c.
const ccNTProb = 10000

// colorChrsState carries the running HMM state and accumulators, mirroring
// args_t in color-chrs.c plus the embedded hmm_t fields color-chrs uses.
type colorChrsState struct {
	hdr    *vcf.Header
	mode   int
	pij    float64
	pgtErr float64

	imother, ifather, ichild int
	isample, jsample         int

	nstates   int
	tprobArr  []float64 // ccNTProb precomputed nstates*nstates matrices, flattened
	currTprob []float64 // scratch for the active per-step matrix
	tmp       []float64 // scratch for matrix multiply
	hapSwitch [8][8]int // trio switch classification (mother/father)

	sites []uint32  // per-site positions for the current chromosome
	eprob []float64 // per-site emission probs, nstates per site

	nsites      int
	nhetMother  int
	nhetFather  int
	prevRID     string
	wroteHeader bool

	fp io.Writer
}

// mat returns MAT(m,n,i,j) = m[n*i+j].
func ccMat(m []float64, n, i, j int) float64 { return m[n*i+j] }

// setMat sets MAT(m,n,i,j).
func ccSetMat(m []float64, n, i, j int, v float64) { m[n*i+j] = v }

// multiplyMatrix computes dst = a*b for n x n matrices, mirroring
// multiply_matrix() in HMM.c (including the a==dst || b==dst aliasing guard via
// the scratch buffer).
func ccMultiplyMatrix(n int, a, b, dst, tmp []float64) {
	out := dst
	aliased := sameSlice(a, dst) || sameSlice(b, dst)
	if aliased {
		out = tmp
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			val := 0.0
			for k := 0; k < n; k++ {
				val += ccMat(a, n, i, k) * ccMat(b, n, k, j)
			}
			ccSetMat(out, n, i, j, val)
		}
	}
	if aliased {
		copy(dst, out[:n*n])
	}
}

// sameSlice reports whether two slices share the same backing array start,
// emulating the C pointer-equality aliasing check in multiply_matrix.
func sameSlice(a, b []float64) bool {
	if len(a) == 0 || len(b) == 0 {
		return false
	}
	return &a[0] == &b[0]
}

// setTprobArr fills tprobArr with the initial matrix followed by its successive
// matrix powers, mirroring hmm_set_tprob() with ntprob=ccNTProb.
func (st *colorChrsState) setTprobArr(tprob []float64) {
	n := st.nstates
	nmat := n * n
	st.tprobArr = make([]float64, nmat*ccNTProb)
	st.currTprob = make([]float64, nmat)
	st.tmp = make([]float64, nmat)
	copy(st.tprobArr[:nmat], tprob)
	for i := 1; i < ccNTProb; i++ {
		ccMultiplyMatrix(n, st.tprobArr[:nmat], st.tprobArr[(i-1)*nmat:i*nmat], st.tprobArr[i*nmat:(i+1)*nmat], st.tmp)
	}
}

// setTprob loads currTprob for a given position gap, mirroring _set_tprob():
// pick the (posDiff % ntprob)-th precomputed matrix, then matrix-power-jump by
// the last precomputed matrix for each full ntprob block.
func (st *colorChrsState) setTprob(posDiff int) {
	n := st.nstates
	nmat := n * n
	idx := posDiff % ccNTProb
	copy(st.currTprob, st.tprobArr[idx*nmat:(idx+1)*nmat])
	jumps := posDiff / ccNTProb
	last := st.tprobArr[(ccNTProb-1)*nmat : ccNTProb*nmat]
	for i := 0; i < jumps; i++ {
		ccMultiplyMatrix(n, last, st.currTprob, st.currTprob, st.tmp)
	}
}

// initHMMTrio builds the 8-state trio transition matrix, mirroring
// init_hmm_trio() in color-chrs.c.
func (st *colorChrsState) initHMMTrio() {
	st.nstates = 8
	n := st.nstates
	tprob := make([]float64, n*n)

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			st.hapSwitch[i][j] = 0
		}
	}
	st.hapSwitch[trioAD][trioAC] = swFather
	st.hapSwitch[trioBC][trioAC] = swMother
	st.hapSwitch[trioBD][trioAC] = swMother | swFather
	st.hapSwitch[trioAC][trioAD] = swFather
	st.hapSwitch[trioBC][trioAD] = swMother | swFather
	st.hapSwitch[trioBD][trioAD] = swMother
	st.hapSwitch[trioAC][trioBC] = swMother
	st.hapSwitch[trioAD][trioBC] = swMother | swFather
	st.hapSwitch[trioBD][trioBC] = swFather
	st.hapSwitch[trioAC][trioBD] = swMother | swFather
	st.hapSwitch[trioAD][trioBD] = swMother
	st.hapSwitch[trioBC][trioBD] = swFather

	st.hapSwitch[trioDA][trioCA] = swFather
	st.hapSwitch[trioCB][trioCA] = swMother
	st.hapSwitch[trioDB][trioCA] = swMother | swFather
	st.hapSwitch[trioCA][trioDA] = swFather
	st.hapSwitch[trioCB][trioDA] = swMother | swFather
	st.hapSwitch[trioDB][trioDA] = swMother
	st.hapSwitch[trioCA][trioCB] = swMother
	st.hapSwitch[trioDA][trioCB] = swMother | swFather
	st.hapSwitch[trioDB][trioCB] = swFather
	st.hapSwitch[trioCA][trioDB] = swMother | swFather
	st.hapSwitch[trioDA][trioDB] = swMother
	st.hapSwitch[trioCB][trioDB] = swFather

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if st.hapSwitch[i][j] == 0 {
				ccSetMat(tprob, n, i, j, 0)
			} else {
				val := 1.0
				if st.hapSwitch[i][j]&swMother != 0 {
					val *= st.pij
				}
				if st.hapSwitch[i][j]&swFather != 0 {
					val *= st.pij
				}
				ccSetMat(tprob, n, i, j, val)
			}
		}
	}
	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			if i != j {
				sum += ccMat(tprob, n, i, j)
			}
		}
		ccSetMat(tprob, n, i, i, 1-sum)
	}
	st.setTprobArr(tprob)
}

// initHMMUnrelated builds the 7-state unrelated transition matrix, mirroring
// init_hmm_unrelated() in color-chrs.c.
func (st *colorChrsState) initHMMUnrelated() {
	st.nstates = 7
	n := st.nstates
	tprob := make([]float64, n*n)

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			ccSetMat(tprob, n, i, j, st.pij)
		}
	}
	ccSetMat(tprob, n, unrl0101, unrlXXXX, st.pij*st.pij)
	ccSetMat(tprob, n, unrl0110, unrlXXXX, st.pij*st.pij)
	ccSetMat(tprob, n, unrlX0x0, unrl0x0x, st.pij*st.pij)
	ccSetMat(tprob, n, unrl0110, unrl0x0x, st.pij*st.pij)
	ccSetMat(tprob, n, unrlX00x, unrl0xx0, st.pij*st.pij)
	ccSetMat(tprob, n, unrl0101, unrl0xx0, st.pij*st.pij)
	ccSetMat(tprob, n, unrl0101, unrlX00x, st.pij*st.pij)
	ccSetMat(tprob, n, unrl0110, unrlX0x0, st.pij*st.pij)
	ccSetMat(tprob, n, unrl0110, unrl0101, st.pij*st.pij)

	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			ccSetMat(tprob, n, i, j, ccMat(tprob, n, j, i))
		}
	}
	for i := 0; i < n; i++ {
		sum := 0.0
		for j := 0; j < n; j++ {
			if i != j {
				sum += ccMat(tprob, n, i, j)
			}
		}
		ccSetMat(tprob, n, i, i, 1-sum)
	}
	st.setTprobArr(tprob)
}

// probShared mirrors prob_shared() in color-chrs.c.
func (st *colorChrsState) probShared(a, b int) float64 {
	if a == b {
		return 1 - st.pgtErr
	}
	return st.pgtErr
}

// probNotShared mirrors prob_not_shared() in color-chrs.c (af is fixed at 0.5).
func (st *colorChrsState) probNotShared(af float64, a, b int) float64 {
	if a != b {
		return af * (1 - af)
	}
	if a == 0 {
		return (1 - af) * (1 - af)
	}
	return af * af
}

// ccPhasedAlleles returns the two biallelic alleles of sample i and whether the
// genotype is usable: present (no missing) and phased (the second allele carries
// the phase bit, matching bcf_gt_is_phased). It returns ok=false otherwise.
func ccPhasedAlleles(v *vcf.Variant, i int) (a, b int, ok bool) {
	gt, parsed := sampleGT(v, i)
	if !parsed || len(gt.alleles) != 2 {
		return 0, 0, false
	}
	if gt.alleles[0] == missingAllele || gt.alleles[1] == missingAllele {
		return 0, 0, false
	}
	// Upstream: return if NEITHER allele is phased. In our model phasing is
	// encoded on the separator preceding the second allele (gt.phased[1]).
	if !gt.phased[1] {
		return 0, 0, false
	}
	return gt.alleles[0], gt.alleles[1], true
}

// appendSite grows the per-chromosome site/eprob arrays by one and returns the
// emission-probability slice for the new site.
func (st *colorChrsState) appendSite(pos uint32) []float64 {
	st.nsites++
	st.sites = append(st.sites, pos)
	base := (st.nsites - 1) * st.nstates
	for len(st.eprob) < base+st.nstates {
		st.eprob = append(st.eprob, 0)
	}
	return st.eprob[base : base+st.nstates]
}

// setObservedProbTrio mirrors set_observed_prob_trio() in color-chrs.c.
func (st *colorChrsState) setObservedProbTrio(v *vcf.Variant) {
	a, b, ok := ccPhasedAlleles(v, st.imother)
	if !ok {
		return
	}
	c, d, ok := ccPhasedAlleles(v, st.ifather)
	if !ok {
		return
	}
	e, f, ok := ccPhasedAlleles(v, st.ichild)
	if !ok {
		return
	}

	mother := (1 << a) | (1 << b)
	father := (1 << c) | (1 << d)
	child := (1 << e) | (1 << f)
	if mother&child == 0 || father&child == 0 {
		return // Mendelian-inconsistent site, skip
	}
	if a != b {
		st.nhetMother++
	}
	if c != d {
		st.nhetFather++
	}

	prob := st.appendSite(uint32(v.Pos - 1))
	prob[trioAC] = st.probShared(e, a) * st.probShared(f, c)
	prob[trioAD] = st.probShared(e, a) * st.probShared(f, d)
	prob[trioBC] = st.probShared(e, b) * st.probShared(f, c)
	prob[trioBD] = st.probShared(e, b) * st.probShared(f, d)
	prob[trioCA] = st.probShared(e, c) * st.probShared(f, a)
	prob[trioDA] = st.probShared(e, d) * st.probShared(f, a)
	prob[trioCB] = st.probShared(e, c) * st.probShared(f, b)
	prob[trioDB] = st.probShared(e, d) * st.probShared(f, b)
}

// setObservedProbUnrelated mirrors set_observed_prob_unrelated() in
// color-chrs.c.
func (st *colorChrsState) setObservedProbUnrelated(v *vcf.Variant) {
	const af = 0.5
	a, b, ok := ccPhasedAlleles(v, st.isample)
	if !ok {
		return
	}
	c, d, ok := ccPhasedAlleles(v, st.jsample)
	if !ok {
		return
	}

	prob := st.appendSite(uint32(v.Pos - 1))
	prob[unrlXXXX] = st.probNotShared(af, a, c) * st.probNotShared(af, a, d) * st.probNotShared(af, b, c) * st.probNotShared(af, b, d)
	prob[unrl0x0x] = st.probShared(a, c) * st.probNotShared(af, b, d)
	prob[unrl0xx0] = st.probShared(a, d) * st.probNotShared(af, b, c)
	prob[unrlX00x] = st.probShared(b, c) * st.probNotShared(af, a, d)
	prob[unrlX0x0] = st.probShared(b, d) * st.probNotShared(af, a, c)
	prob[unrl0101] = st.probShared(a, c) * st.probShared(b, d)
	prob[unrl0110] = st.probShared(a, d) * st.probShared(b, c)
}

// runViterbi runs the Viterbi decode over the current chromosome's sites,
// returning the per-site decoded state path, mirroring hmm_run_viterbi() in
// HMM.c (uniform initial state probs, no snapshots, no set_tprob hook).
func (st *colorChrsState) runViterbi() []int {
	n := st.nsites
	nstates := st.nstates
	if n == 0 {
		return nil
	}
	// vpath[i*nstates + j] holds the best predecessor of state j at step i; after
	// trace-back vpath[i*nstates] holds the decoded state at step i.
	vpath := make([]int, n*nstates)
	vprob := make([]float64, nstates)
	vprobTmp := make([]float64, nstates)
	for i := 0; i < nstates; i++ {
		vprob[i] = 1.0 / float64(nstates) // hmm_init_states with probs==NULL
	}
	prevPos := st.sites[0]

	for i := 0; i < n; i++ {
		eprob := st.eprob[i*nstates : (i+1)*nstates]
		var posDiff int
		if st.sites[i] != prevPos {
			posDiff = int(st.sites[i] - prevPos - 1)
		}
		st.setTprob(posDiff)
		prevPos = st.sites[i]

		vnorm := 0.0
		for j := 0; j < nstates; j++ {
			vmax := 0.0
			kVmax := 0
			for k := 0; k < nstates; k++ {
				pval := vprob[k] * ccMat(st.currTprob, nstates, j, k)
				if vmax < pval {
					vmax = pval
					kVmax = k
				}
			}
			vpath[i*nstates+j] = kVmax
			vprobTmp[j] = vmax * eprob[j]
			vnorm += vprobTmp[j]
		}
		for j := 0; j < nstates; j++ {
			vprobTmp[j] /= vnorm
		}
		vprob, vprobTmp = vprobTmp, vprob
	}

	// Most likely final state.
	iptr := 0
	for i := 1; i < nstates; i++ {
		if vprob[iptr] < vprob[i] {
			iptr = i
		}
	}
	// Trace back, reusing vpath[i*nstates] for the decoded state.
	for i := n - 1; i >= 0; i-- {
		iptrPrev := vpath[i*nstates+iptr]
		vpath[i*nstates] = iptr
		iptr = iptrPrev
	}
	path := make([]int, n)
	for i := 0; i < n; i++ {
		path[i] = vpath[i*nstates]
	}
	return path
}

// ccHapSwitchAt returns hap_switch[i][j] using C's flat 8x8 row-major indexing,
// so that j==-1 reads hap_switch[i-1][7] exactly as the C plugin does on the
// first segment (prev_state==-1). The slot before the array (i==0, j==-1) is 0,
// matching the static BSS layout. Out-of-range flat indices yield 0.
func ccHapSwitchAt(hs *[8][8]int, i, j int) int {
	flat := i*8 + j
	if flat < 0 || flat >= 64 {
		return 0
	}
	return hs[flat/8][flat%8]
}

// flushViterbi decodes the current chromosome, writes its SG segment lines and
// the two SW switch-count lines to the .dat file, then resets the per-chromosome
// accumulators, mirroring flush_viterbi() in color-chrs.c.
func (st *colorChrsState) flushViterbi() {
	var s1, s2, s3 string
	if st.mode == ccUnrl {
		s1 = st.hdr.Samples[st.isample]
		s2 = st.hdr.Samples[st.jsample]
	} else {
		s1 = st.hdr.Samples[st.imother]
		s3 = st.hdr.Samples[st.ifather]
		s2 = st.hdr.Samples[st.ichild]
	}

	if !st.wroteHeader {
		fmt.Fprintf(st.fp, "# SG, shared segment\t[2]Chromosome\t[3]Start\t[4]End\t[5]%s:1\t[6]%s:2\n", s2, s2)
		fmt.Fprintf(st.fp, "# SW, number of switches\t[3]Sample\t[4]Chromosome\t[5]nHets\t[5]nSwitches\t[6]switch rate\n")
		st.wroteHeader = true
	}

	path := st.runViterbi()
	iprev := -1
	prevState := -1
	nswitchMother := 0
	nswitchFather := 0
	chr := st.prevRID
	for i := 0; i < st.nsites; i++ {
		state := path[i]
		if state != prevState || i+1 == st.nsites {
			var start, end uint32 = 1, 1
			if iprev >= 0 {
				start = st.sites[iprev] + 1
			}
			if i > 0 {
				end = st.sites[i-1]
			}
			if st.mode == ccUnrl {
				switch prevState {
				case unrl0x0x:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:1\t-\n", chr, start, end, s1)
				case unrl0xx0:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t-\t%s:1\n", chr, start, end, s1)
				case unrlX00x:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:2\t-\n", chr, start, end, s1)
				case unrlX0x0:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t-\t%s:2\n", chr, start, end, s1)
				case unrl0101:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr, start, end, s1, s1)
				case unrl0110:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr, start, end, s1, s1)
				}
			} else {
				switch prevState {
				case trioAC:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:1\t%s:1\n", chr, start, end, s1, s3)
				case trioAD:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr, start, end, s1, s3)
				case trioBC:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr, start, end, s1, s3)
				case trioBD:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:2\t%s:2\n", chr, start, end, s1, s3)
				case trioCA:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:1\t%s:1\n", chr, start, end, s3, s1)
				case trioDA:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr, start, end, s3, s1)
				case trioCB:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr, start, end, s3, s1)
				case trioDB:
					fmt.Fprintf(st.fp, "SG\t%s\t%d\t%d\t%s:2\t%s:2\n", chr, start, end, s3, s1)
				}
				// Upstream indexes hap_switch[state][prev_state] even on the very
				// first segment, when prev_state==-1. In C that reads the int
				// immediately before hap_switch[state][0], i.e. hap_switch[state-1][7]
				// (and hap_switch[-1][7]==0 for state==0, the BSS slot before the
				// array). We reproduce that exact aliasing so the switch counts match
				// byte-for-byte; in practice column 7 (TRIO_DB) of every row is 0, so
				// this contributes no spurious switch.
				sw := ccHapSwitchAt(&st.hapSwitch, state, prevState)
				if sw&swMother != 0 {
					nswitchMother++
				}
				if sw&swFather != 0 {
					nswitchFather++
				}
			}
			iprev = i - 1
		}
		prevState = state
	}

	var mrate, frate float64
	if st.nhetMother > 1 {
		mrate = float64(float32(nswitchMother) / float32(st.nhetMother-1))
	}
	if st.nhetFather > 1 {
		frate = float64(float32(nswitchFather) / float32(st.nhetFather-1))
	}
	// Upstream writes both SW lines unconditionally. In unrelated mode s3 is the
	// NULL father pointer; glibc's printf("%s", NULL) renders it as "(null)", so
	// we reproduce that token for byte parity.
	s3out := s3
	if st.mode == ccUnrl {
		s3out = "(null)"
	}
	fmt.Fprintf(st.fp, "SW\t%s\t%s\t%d\t%d\t%f\n", s1, chr, st.nhetMother, nswitchMother, mrate)
	fmt.Fprintf(st.fp, "SW\t%s\t%s\t%d\t%d\t%f\n", s3out, chr, st.nhetFather, nswitchFather, frate)

	st.nsites = 0
	st.sites = st.sites[:0]
	st.eprob = st.eprob[:0]
	st.nhetFather = 0
	st.nhetMother = 0
}
