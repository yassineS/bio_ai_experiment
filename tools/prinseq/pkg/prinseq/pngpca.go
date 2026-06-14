package prinseq

// pngpca.go ports the dinucleotide-odds PCA scatter plots
// (_pm.png microbial, _pv.png viral) from prinseq-graphs.pl
// (createPCAPlot / convertToPCAValues, lines 1213-1425). Upstream
// runs the user's 10-value odds-ratio row against a fixed reference
// matrix of microbial/viral metagenomes through Statistics::PCA, then
// scatters the transformed first two principal components.
//
// PCA is genuinely linear-algebra-heavy (it is an eigendecomposition
// of the covariance matrix), which is exactly the sanctioned use of
// `gonum.org/v1/gonum` per CLAUDE.md ("eigendecomp ... reuse is
// allowed"). We use gonum's stat.PrincipalComponents. The exact
// numeric scores cannot match Statistics::PCA byte-for-byte
// (different centring/scaling and sign conventions), and PNG
// byte-identity is out of scope regardless; what matters is that the
// plot renders the full reference set plus the user point with valid,
// decodable output.

import (
	"image/color"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
)

// dinucMIC is the reference microbial dinucleotide odds-ratio matrix
// (prinseq-graphs.pl $DINUCODDS_MIC, 20 rows × 10 cols).
var dinucMIC = [][]float64{
	{1.13127323, 0.853587195, 0.911041047, 1.104520778, 1.065586428, 1.021434164, 0.999734139, 1.063684014, 1.078035184, 0.733596552},
	{1.173267344, 0.840539337, 0.919534602, 1.068050141, 1.062394214, 1.051999071, 0.96770576, 1.035511729, 1.095600433, 0.72328141},
	{1.172939786, 0.84567902, 0.911836259, 1.106288994, 1.05351787, 1.026143368, 1.002308358, 1.066319771, 1.094918797, 0.710733535},
	{1.073527689, 0.850290918, 0.978455025, 1.080882178, 1.111174765, 1.010754115, 0.895668707, 1.072980666, 1.079304608, 0.754057386},
	{1.08807747, 0.837444678, 0.95824965, 1.097310298, 1.118897971, 1.030863881, 0.886827263, 1.072349394, 1.07406322, 0.733440096},
	{1.071685485, 0.861055813, 0.966566865, 1.090268118, 1.112945761, 1.012538936, 0.909535491, 1.063745603, 1.071156598, 0.755770377},
	{1.142698587, 0.867936867, 1.000612099, 0.977934257, 1.111801746, 1.018318601, 0.788556794, 0.987763594, 1.184649653, 0.784776176},
	{1.134560074, 0.876651844, 0.998190253, 0.995723123, 1.128448077, 1.014172324, 0.781776188, 0.971020602, 1.182411449, 0.786449476},
	{1.180029632, 0.787899325, 1.01316945, 0.932268406, 1.077837263, 1.211699678, 0.612128817, 1.033036699, 1.157314398, 0.74940288},
	{1.160925546, 0.788308899, 1.003702496, 0.965371236, 1.076051693, 1.188304271, 0.641536444, 1.070331188, 1.124067192, 0.740126813},
	{1.173873006, 0.790118011, 1.014718833, 0.937979878, 1.07453725, 1.207167373, 0.622279064, 1.046150047, 1.145627707, 0.742212886},
	{1.128383111, 0.870541389, 0.987269741, 0.98353238, 1.115643879, 1.040107028, 0.774505865, 1.010896432, 1.164757274, 0.775254395},
	{1.15297511, 0.853883985, 0.956393231, 1.000027661, 1.139915472, 1.01355294, 0.838843622, 1.015553125, 1.216219741, 0.70447264},
	{1.148264236, 0.852123859, 0.974568293, 0.985455546, 1.13192373, 1.015879393, 0.828987111, 1.016820786, 1.216647853, 0.71634006},
	{1.12933788, 0.831777975, 1.005434367, 0.991081409, 1.126146895, 1.07421504, 0.69343913, 1.054032466, 1.14809591, 0.728541157},
	{1.124157235, 0.828112691, 1.022348424, 0.983822386, 1.143028487, 1.081830005, 0.672594435, 1.05685982, 1.149537403, 0.684432106},
	{1.128029586, 0.841853305, 1.00983936, 0.967179139, 1.122524003, 1.094555807, 0.659238308, 1.061578854, 1.1243601, 0.740148171},
	{1.093521636, 0.855071052, 0.929160818, 1.203773691, 1.178257185, 0.881341255, 1.078305505, 1.051988532, 1.169143967, 0.555057308},
	{1.073737278, 0.877396537, 0.968017446, 1.124155374, 1.166244435, 0.909044208, 0.999147578, 1.071098934, 1.120156138, 0.607444953},
	{1.092150184, 0.863407008, 0.927040387, 1.185387013, 1.171670826, 0.882276859, 1.083058605, 1.048379554, 1.168635365, 0.580337997},
}

// dinucVIR is the reference viral dinucleotide odds-ratio matrix
// (prinseq-graphs.pl $DINUCODDS_VIR, 22 rows × 10 cols).
var dinucVIR = [][]float64{
	{1.086940308, 0.98976932, 1.034167044, 0.880024041, 1.070421277, 0.990687084, 0.890945575, 1.069957074, 0.92465631, 0.803973303},
	{1.101064857, 0.986812783, 1.038299155, 0.896162618, 1.081652847, 0.976365237, 0.867445186, 1.06727283, 0.94688543, 0.768007295},
	{1.071548411, 0.912204166, 1.196914981, 0.80628184, 1.294201511, 1.148517794, 0.269295791, 1.033948026, 0.895951033, 0.623192149},
	{1.090253719, 0.907428629, 1.203991784, 0.786359294, 1.281499107, 1.145421568, 0.235974709, 1.033437274, 0.899580091, 0.631699771},
	{1.075864745, 1.003413074, 1.01872902, 0.897841689, 0.980373171, 1.05854979, 0.934262259, 1.052477953, 0.88145851, 0.889239724},
	{1.101890467, 1.030028291, 1.019912674, 0.84191395, 1.0015174, 1.069546264, 0.900151602, 0.996269395, 0.889195343, 0.904039022},
	{1.152417359, 0.855028574, 0.91164793, 1.017415486, 1.114163672, 1.128353311, 0.846355573, 0.916745489, 1.206820475, 0.811014651},
	{1.142454218, 0.8635465, 0.923406967, 1.026242747, 1.134445058, 1.131747833, 0.79793368, 0.920767641, 1.179468556, 0.799770057},
	{1.124462747, 0.873556143, 0.945627041, 1.013755408, 1.159866153, 1.096259526, 0.757315047, 0.972924919, 1.105562567, 0.772731886},
	{1.143826972, 0.866968779, 0.995740249, 0.945859278, 1.109590621, 1.089305083, 0.76048874, 0.971561388, 1.157101408, 0.792923027},
	{1.131900141, 0.82776996, 0.996204924, 0.999433455, 1.024692372, 1.071176333, 0.921026216, 1.088936699, 1.054010776, 0.773498892},
	{1.042180476, 0.930180412, 1.019242897, 0.98909997, 1.006666828, 1.046708539, 0.959492164, 1.011183418, 1.055168776, 0.937433818},
	{1.086515695, 0.985345815, 0.930914307, 0.969581792, 1.043010232, 1.087463712, 0.939482285, 0.990551965, 0.954752469, 0.893972874},
	{1.096657826, 0.950117614, 0.936195529, 0.965619788, 1.114975275, 1.077011195, 0.843153131, 0.989128406, 1.043790912, 0.840634731},
	{1.158030995, 0.935307365, 0.874812261, 1.056236525, 1.117171274, 0.937484692, 1.057442372, 0.970079538, 1.174848738, 0.725071711},
	{1.15591506, 0.93000227, 0.883538923, 1.0567652, 1.095730954, 0.944489906, 1.074229471, 0.983993745, 1.156051409, 0.726688465},
	{1.205726473, 0.924439339, 1.049457756, 0.805718412, 0.975472778, 1.07581991, 0.726992211, 1.075025787, 0.8704929, 0.726672843},
	{1.188544681, 0.95239611, 1.049066985, 0.790031334, 1.038632598, 1.056749787, 0.665197397, 1.057566244, 0.862429061, 0.708982398},
	{1.063631482, 0.925593715, 1.014869316, 0.944904401, 1.119690731, 1.325971834, 0.273781451, 0.943347677, 1.06438014, 0.920825904},
	{1.077560287, 0.911888545, 1.044147857, 0.927758054, 1.058535939, 1.296838544, 0.421514996, 0.945722451, 1.128317986, 0.926419928},
	{1.163753415, 0.989905668, 0.893599328, 0.955641844, 1.176047687, 0.941559156, 0.950641089, 0.959741692, 1.100815282, 0.72491925},
	{1.139253929, 0.946297517, 0.922096125, 1.024801537, 1.205206793, 0.968818717, 0.915801342, 0.971626058, 1.107569276, 0.627623404},
}

// pcaPointLabels mirrors the upstream $DATA_VIR / $DATA_MIC group
// numbers (column index 1) used as the scatter labels. Index 0 is the
// user input. Colours: user is red, references blue.
func pcaLabels(typ string) []string {
	if typ == "v" {
		return []string{
			"1", "1", "2", "2", "1", "1", "3", "3", "4", "4", "5", "5",
			"6", "6", "7", "7", "8", "8", "9", "9", "10", "10", "U",
		}
	}
	return []string{
		"1", "1", "1", "2", "2", "2", "3", "3", "4", "4", "4", "3",
		"5", "5", "3", "3", "3", "6", "7", "6", "U",
	}
}

// drawPCAPlot computes the two leading principal components of the
// reference matrix (microbial or viral) with the user's odds-ratio
// row appended, then scatters them. Returns the rendered canvas.
func drawPCAPlot(userRow []float64, typ string) *canvas {
	var ref [][]float64
	if typ == "v" {
		ref = dinucVIR
	} else {
		ref = dinucMIC
	}

	// Build the data matrix: reference rows + user row (last).
	rows := len(ref) + 1
	cols := 10
	data := make([]float64, 0, rows*cols)
	for _, r := range ref {
		data = append(data, r[:cols]...)
	}
	// Pad/truncate the user row to 10 columns.
	urow := make([]float64, cols)
	copy(urow, userRow)
	data = append(data, urow...)

	m := mat.NewDense(rows, cols, data)

	var pc stat.PC
	ok := pc.PrincipalComponents(m, nil)
	xs := make([]float64, rows)
	ys := make([]float64, rows)
	if ok {
		var proj mat.Dense
		var vecs mat.Dense
		pc.VectorsTo(&vecs)
		proj.Mul(m, vecs.Slice(0, cols, 0, 2))
		for i := 0; i < rows; i++ {
			xs[i] = proj.At(i, 0)
			ys[i] = proj.At(i, 1)
		}
	}

	return scatterPCA(xs, ys, pcaLabels(typ), typ)
}

// scatterPCA renders the PCA scatter (geometry follows
// createPCAPlot). Reference points are blue, the user point red.
func scatterPCA(xs, ys []float64, labels []string, typ string) *canvas {
	const (
		size   = 5
		offset = 20
		left   = 25
		bottom = 15
		height = 500
		space  = 10
	)
	top := 20
	if typ == "v" {
		top = 35
	}
	width := left + offset*2 + height + 2*space
	totalH := top + bottom + offset*2 + height + 2*space
	c := newCanvas(width, totalH)

	plotX := left + offset
	plotY := top + offset
	c.fillRect(plotX, plotY, height+2*space, height+2*space, colBackplot)

	// data ranges
	xmin, xmax := minMax(xs)
	ymin, ymax := minMax(ys)
	xrange := xmax - xmin
	yrange := ymax - ymin
	if xrange == 0 {
		xrange = 1
	}
	if yrange == 0 {
		yrange = 1
	}

	// zero crosshairs
	zx := plotX + space + int(absF(xmin)/xrange*float64(height))
	zy := plotY + space + int(absF(ymax)/yrange*float64(height))
	c.vLine(zx, plotY, plotY+height+2*space, colHelpline)
	c.hLine(plotX, plotX+2*space+height, zy, colHelpline)

	// axis ticks + labels
	c.vLine(plotX+space, plotY+height+2*space, plotY+height+2*space+3, colTick)
	c.vLine(plotX+space+height, plotY+height+2*space, plotY+height+2*space+3, colTick)
	c.textCentered(plotX+space, plotY+height+2*space+4, format2(xmin), colTick)
	c.textCentered(plotX+space+height, plotY+height+2*space+4, format2(xmax), colTick)
	c.textRight(plotX-5, plotY+space-3, format2(ymax), colTick)
	c.textRight(plotX-5, plotY+height+space-3, format2(ymin), colTick)

	// axis labels
	c.textCentered(plotX+height/2+space, plotY+height+14+2*space, "1st Principal Component Score", colLabel)
	c.textVertical(offset-4, plotY+height/2+textWidth("2nd Principal Component Score")/2, "2nd Principal Component Score", colLabel)

	// type badge
	badge := "M"
	if typ == "v" {
		badge = "V"
	}
	c.fillCircle(offset/2+textWidth(badge)/2, offset-5, 10, color.RGBA{0, 0, 0, 128})
	c.text(offset/2, offset-3, badge, color.RGBA{255, 255, 255, 255})

	// dots
	refCol := color.RGBA{127, 127, 255, 255}
	userCol := color.RGBA{255, 127, 127, 255}
	n := len(xs)
	for i := 0; i < n; i++ {
		col := refCol
		if i == n-1 {
			col = userCol
		}
		px := plotX + space + int((xs[i]+absF(xmin))/xrange*float64(height))
		py := plotY + space + int((ys[i]+absF(ymin))/yrange*float64(height))
		c.fillCircle(px, py, size, col)
		if i < len(labels) {
			c.text(px+size+1, py+size*2-3, labels[i], colLabel)
		}
	}
	return c
}

func minMax(xs []float64) (mn, mx float64) {
	if len(xs) == 0 {
		return 0, 1
	}
	mn, mx = xs[0], xs[0]
	for _, v := range xs {
		if v < mn {
			mn = v
		}
		if v > mx {
			mx = v
		}
	}
	return mn, mx
}

func absF(f float64) float64 {
	if f < 0 {
		return -f
	}
	return f
}
