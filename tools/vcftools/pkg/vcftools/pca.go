// Package vcftools PCA implementation.
//
// Mirrors upstream's output_PCA / output_PCA_SNP_loadings in
// reference_code/vcftools/src/cpp/variant_file_output.cpp:4871-5246,
// following the Patterson, Price & Reich (2006) recipe:
//
//  1. Build a per-individual matrix M of centred (and optionally
//     variance-normalised) ALT-allele counts at every biallelic,
//     non-monomorphic, diploid kept site.
//  2. Compute the Genomic Relatedness Matrix X = (1/n) * M * Mᵀ.
//  3. Eigendecompose X.
//  4. Sort eigenpairs by |eigenvalue| descending.
//  5. Write `<prefix>.pca` with the per-sample loadings on each
//     principal component.
//  6. Optionally, project raw genotypes back onto the first K
//     eigenvectors and write per-site SNP loadings to
//     `<prefix>.pca.loadings`.
//
// Wave-19 fix-on-port: upstream's `output_PCA` only appends to M[i]
// when the i-th kept individual has a non-missing genotype, but
// advances the per-site index regardless. With any missing data the
// per-individual M vectors become jagged and the inner triple loop
// at variant_file_output.cpp:4988-4991 reads off the end of the
// shortest vector. This port mirrors upstream's "no imputation"
// stance but enforces it correctly by SKIPPING any site that has a
// missing GT in any kept individual — see `docs/UPSTREAM_BUGS.md`
// (fix-on-port section).
//
// Numerical strategy: X is real symmetric so we eigendecompose via
// `gonum.org/v1/gonum/mat` (SymEigen) rather than the general-purpose
// LAPACK dgeev upstream calls. SymEigen guarantees real eigenvalues
// and orthonormal eigenvectors. Eigenvector signs are arbitrary
// (both LAPACK and gonum); parity tests compare with sign-tolerance.
package vcftools

import (
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
	"gonum.org/v1/gonum/mat"
)

// pcaRunner accumulates centred / optionally normalised genotype values
// for every kept-and-non-missing site, then computes the eigendecomposition
// of the resulting GRM and writes `<prefix>.pca` (and optionally
// `<prefix>.pca.loadings`).
type pcaRunner struct {
	indvNames    []string
	useNorm      bool
	writePCA     bool // true if --pca or --pca-no-norm was set
	snpLoadingsK int  // 0 = no .pca.loadings output
	streamOut    bool // if true, .pca goes to stdout (mirroring --stdout)

	// M[i] is the centred (and optionally normalised) genotype value for
	// kept individual i at every accepted site. Always rectangular
	// (N_indv × N_sites_accepted) because we skip sites with any
	// missing kept-indv genotype — see the fix-on-port note above.
	M [][]float64

	// rawGenos[i][s] is the raw ALT-allele count (0/1/2) for indv i at
	// accepted site s. Only populated when snpLoadingsK > 0.
	rawGenos [][]float64

	// chrom / pos identify each accepted site, in input order. Only
	// populated when snpLoadingsK > 0 (the .pca output itself doesn't
	// need per-site coords).
	chrom []string
	pos   []int

	// skippedMonomorphic / skippedMissing / skippedNonDiploid count
	// per-reason drops. Reported via stderr at writeOutput time as a
	// debug aid (upstream prints no such summary).
	skippedMonomorphic int
	skippedMissing     int
}

// newPCARunner returns a runner sized for `indvNames` kept samples.
// `useNorm` toggles per-SNP variance normalisation (default true,
// disabled via --pca-no-norm). `snpLoadingsK` is the number of leading
// principal components to project genotypes onto for the
// `<prefix>.pca.loadings` output; 0 disables that output entirely
// (matching upstream where `--pca-snp-loadings` is its own opt-in flag).
func newPCARunner(indvNames []string, writePCA, useNorm bool, snpLoadingsK int, streamOut bool) *pcaRunner {
	n := len(indvNames)
	r := &pcaRunner{
		indvNames:    append([]string(nil), indvNames...),
		useNorm:      useNorm,
		writePCA:     writePCA,
		snpLoadingsK: snpLoadingsK,
		streamOut:    streamOut,
		M:            make([][]float64, n),
	}
	if snpLoadingsK > 0 {
		r.rawGenos = make([][]float64, n)
	}
	return r
}

// pcaEpsilon is the threshold below which an allele frequency is treated
// as monomorphic (matches upstream's `numeric_limits<double>::epsilon()`
// check at variant_file_output.cpp:4948).
const pcaEpsilon = 2.2204460492503131e-16

// addVariant ingests one filtered variant. Mirrors upstream's per-site
// branch in `output_PCA` (variant_file_output.cpp:4928-4972). Returns an
// error for inputs upstream itself errors on (non-biallelic, non-diploid);
// silently skips monomorphic sites (upstream `continue` at line 4949) and
// sites with any missing kept-indv genotype (fix-on-port — see file
// header).
func (p *pcaRunner) addVariant(v *vcf.Variant) error {
	if p == nil {
		return nil
	}
	// Biallelic check — upstream LOG.error at variant_file_output.cpp:4939.
	if len(v.Alt) != 1 {
		return fmt.Errorf("PCA only works for biallelic sites (site %s:%d has %d ALT alleles)", v.Chrom, v.Pos, len(v.Alt))
	}
	// Collect per-individual allele counts. Upstream's
	// `get_indv_GENOTYPE_ids` returns (-1,-1) for missing; the sum-of-two-codes
	// representation `x = geno_id.first + geno_id.second` is -1 if missing.
	// We replicate the value (`raw`) and gate the site on missingness
	// across ALL kept individuals.
	nIndv := len(p.indvNames)
	if nIndv == 0 {
		return fmt.Errorf("PCA requires at least one kept individual")
	}
	raw := make([]float64, nIndv)
	var altCount, totalCount float64
	for i := 0; i < nIndv; i++ {
		gt := getGT(v, i)
		a, b, miss := parseDiploidGT(gt)
		if miss {
			// Wave-19 fix-on-port: drop the entire site rather than
			// hit upstream's jagged-M[i] bug. See file header.
			p.skippedMissing++
			return nil
		}
		// Non-diploid would be detected here too — parseDiploidGT
		// returns missing=true for a haploid call, so the missing
		// branch above also catches that case. Upstream
		// LOG.errors on non-diploid (line 4943); the missing
		// branch's silent skip is the fix-on-port equivalent.
		if a < 0 || b < 0 {
			p.skippedMissing++
			return nil
		}
		// Biallelic site means allele indices ∈ {0,1}; anything
		// higher indicates a header/data mismatch the prior
		// biallelic check should have caught — surface it.
		if a > 1 || b > 1 {
			return fmt.Errorf("PCA: non-biallelic allele index at %s:%d sample %s (alleles %d/%d)", v.Chrom, v.Pos, p.indvNames[i], a, b)
		}
		raw[i] = float64(a + b)
		altCount += raw[i]
		totalCount += 2
	}
	if totalCount == 0 {
		// Should be unreachable — we'd have returned in the missing
		// branch above — but guard against a divide-by-zero.
		return nil
	}
	freq := altCount / totalCount
	if freq <= pcaEpsilon || freq >= 1.0-pcaEpsilon {
		p.skippedMonomorphic++
		return nil
	}
	mu := 2.0 * freq
	div := 1.0 / math.Sqrt(freq*(1.0-freq))
	for i := 0; i < nIndv; i++ {
		val := raw[i] - mu
		if p.useNorm {
			val *= div
		}
		p.M[i] = append(p.M[i], val)
	}
	if p.snpLoadingsK > 0 {
		for i := 0; i < nIndv; i++ {
			p.rawGenos[i] = append(p.rawGenos[i], raw[i])
		}
		p.chrom = append(p.chrom, v.Chrom)
		p.pos = append(p.pos, v.Pos)
	}
	return nil
}

// computeAndWrite finalises the PCA computation and writes outputs.
// Errors out if N_sites <= N_indv (upstream constraint at
// variant_file_output.cpp:4975-4976).
func (p *pcaRunner) computeAndWrite(prefix string) error {
	if p == nil {
		return nil
	}
	n := len(p.indvNames)
	if n == 0 {
		return fmt.Errorf("PCA: no kept individuals")
	}
	if len(p.M[0]) <= n {
		return fmt.Errorf("PCA computation requires that there are more sites than individuals (got %d sites, %d individuals)", len(p.M[0]), n)
	}
	nSites := len(p.M[0])

	// Build the symmetric N×N GRM. X[i][j] = (1/n_sites) * Σ_s M[i][s]*M[j][s].
	// Compute the upper triangle, then symmetrise — same shape as upstream
	// variant_file_output.cpp:4988-5001.
	xData := make([]float64, n*n)
	for i := 0; i < n; i++ {
		for j := i; j < n; j++ {
			var sum float64
			for s := 0; s < nSites; s++ {
				sum += p.M[i][s] * p.M[j][s]
			}
			sum /= float64(nSites)
			xData[i*n+j] = sum
			xData[j*n+i] = sum
		}
	}

	// Symmetric eigendecomposition via gonum.
	sym := mat.NewSymDense(n, xData)
	var es mat.EigenSym
	if ok := es.Factorize(sym, true); !ok {
		return fmt.Errorf("PCA: eigendecomposition did not converge")
	}
	// gonum returns ascending order; we need descending by magnitude
	// (upstream's dgeev_sort at dgeev.cpp:17-65 sorts by |Er|² + |Ei|²
	// which for our symmetric input reduces to |Er|).
	rawVals := es.Values(nil)
	var vecs mat.Dense
	es.VectorsTo(&vecs)
	eigVals, eigVecs := sortEigenByMagnitudeDesc(rawVals, &vecs)

	// Apply a canonical sign convention so output is deterministic across
	// LAPACK implementations: the first |component| > tolerance in each
	// eigenvector must be positive. Upstream does NOT do this (signs depend
	// on the installed LAPACK), so parity tests compare with sign tolerance.
	canonicaliseEigenvectorSigns(eigVecs)

	if p.writePCA {
		if err := p.writePCAFile(prefix, eigVals, eigVecs); err != nil {
			return err
		}
	}
	if p.snpLoadingsK > 0 {
		if err := p.writeSNPLoadings(prefix, eigVecs); err != nil {
			return err
		}
	}
	return nil
}

// writePCAFile emits `<prefix>.pca` (or stdout if streamOut). Format
// mirrors upstream variant_file_output.cpp:5017-5034 exactly:
//
//	INDV\tEIG_0\tEIG_1\t...\tEIG_{n-1}
//	EIGENVALUE\t<v0>\t<v1>\t...
//	<indv_0>\t<vec[0,0]>\t<vec[0,1]>\t...
//	...
func (p *pcaRunner) writePCAFile(prefix string, vals []float64, vecs *mat.Dense) error {
	var b strings.Builder
	n := len(p.indvNames)
	b.WriteString("INDV")
	for j := 0; j < n; j++ {
		b.WriteString("\tEIG_")
		b.WriteString(strconv.Itoa(j))
	}
	b.WriteByte('\n')
	b.WriteString("EIGENVALUE")
	for j := 0; j < n; j++ {
		b.WriteByte('\t')
		b.WriteString(formatPCAFloat(vals[j]))
	}
	b.WriteByte('\n')
	for i := 0; i < n; i++ {
		b.WriteString(p.indvNames[i])
		for j := 0; j < n; j++ {
			b.WriteByte('\t')
			b.WriteString(formatPCAFloat(vecs.At(i, j)))
		}
		b.WriteByte('\n')
	}

	var w io.Writer
	if p.streamOut {
		w = os.Stdout
	} else {
		f, err := iohelper.OpenWriter(prefix + ".pca")
		if err != nil {
			return fmt.Errorf("opening PCA output: %w", err)
		}
		defer f.Close()
		w = f
	}
	if _, err := w.Write([]byte(b.String())); err != nil {
		return fmt.Errorf("writing PCA output: %w", err)
	}
	return nil
}

// writeSNPLoadings emits `<prefix>.pca.loadings`. Mirrors upstream
// variant_file_output.cpp:5217-5241. Per accepted site:
//
//	gamma_k = Σ_i (g_i * Evecs[i][k]) / Σ_i Evecs[i][k]²
//
// where g_i is the raw (uncentred, unnormalised) ALT-allele count for
// individual i. The denominator is constant across sites for each k
// (because we drop sites with any missing GT, every site sums over the
// full eigenvector), but we compute it per-site to match upstream's loop
// body verbatim.
func (p *pcaRunner) writeSNPLoadings(prefix string, vecs *mat.Dense) error {
	var b strings.Builder
	k := p.snpLoadingsK
	b.WriteString("CHROM\tPOS")
	for j := 0; j < k; j++ {
		b.WriteString("\tGAMMA_")
		b.WriteString(strconv.Itoa(j))
	}
	b.WriteByte('\n')

	n := len(p.indvNames)
	nSites := len(p.chrom)
	for s := 0; s < nSites; s++ {
		gamma := make([]float64, k)
		aSum := make([]float64, k)
		for i := 0; i < n; i++ {
			x := p.rawGenos[i][s]
			for kk := 0; kk < k; kk++ {
				v := vecs.At(i, kk)
				gamma[kk] += x * v
				aSum[kk] += v * v
			}
		}
		b.WriteString(p.chrom[s])
		b.WriteByte('\t')
		b.WriteString(strconv.Itoa(p.pos[s]))
		for kk := 0; kk < k; kk++ {
			b.WriteByte('\t')
			if aSum[kk] == 0 {
				b.WriteString("0")
			} else {
				b.WriteString(formatPCAFloat(gamma[kk] / aSum[kk]))
			}
		}
		b.WriteByte('\n')
	}

	f, err := iohelper.OpenWriter(prefix + ".pca.loadings")
	if err != nil {
		return fmt.Errorf("opening PCA loadings output: %w", err)
	}
	defer f.Close()
	if _, err := f.Write([]byte(b.String())); err != nil {
		return fmt.Errorf("writing PCA loadings: %w", err)
	}
	return nil
}

// sortEigenByMagnitudeDesc returns eigenvalues sorted descending by |v|
// alongside the corresponding columns of `vecs` reordered to match.
// gonum's SymEigen yields ascending real eigenvalues; this is the
// upstream-compatible reorder. Returns a fresh slice for vals and a
// fresh Dense for vecs to avoid aliasing surprises.
func sortEigenByMagnitudeDesc(vals []float64, vecs *mat.Dense) ([]float64, *mat.Dense) {
	n := len(vals)
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(a, b int) bool {
		return math.Abs(vals[idx[a]]) > math.Abs(vals[idx[b]])
	})
	outVals := make([]float64, n)
	outVecs := mat.NewDense(n, n, nil)
	for k, oldCol := range idx {
		outVals[k] = vals[oldCol]
		for r := 0; r < n; r++ {
			outVecs.Set(r, k, vecs.At(r, oldCol))
		}
	}
	return outVals, outVecs
}

// canonicaliseEigenvectorSigns flips each eigenvector (column of vecs)
// so its first |component| > 1e-12 is positive. Removes the LAPACK
// implementation-dependent sign ambiguity for deterministic output.
func canonicaliseEigenvectorSigns(vecs *mat.Dense) {
	r, c := vecs.Dims()
	const tol = 1e-12
	for j := 0; j < c; j++ {
		flip := false
		for i := 0; i < r; i++ {
			v := vecs.At(i, j)
			if math.Abs(v) > tol {
				flip = v < 0
				break
			}
		}
		if !flip {
			continue
		}
		for i := 0; i < r; i++ {
			vecs.Set(i, j, -vecs.At(i, j))
		}
	}
}

// formatPCAFloat renders one PCA matrix entry with the same precision
// upstream emits via C++ `ostream <<` default — 6 significant digits.
// Uses Go's %.6g to keep the small-magnitude-exponent form
// (e.g. "-3.753e-16" for ~zero eigenvalues of the centred GRM).
func formatPCAFloat(x float64) string {
	if x == 0 {
		return "0"
	}
	return strconv.FormatFloat(x, 'g', 6, 64)
}
