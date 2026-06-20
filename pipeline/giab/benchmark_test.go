package giab

import (
	"strings"
	"testing"
)

// happySample is a captured hap.py <prefix>.summary.csv. It has the columns the
// parser keys on and both ALL and PASS rows for SNP and INDEL.
const happySample = `Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score
INDEL,ALL,525469,521228,4241,1041083,2932,515436,0.991929,0.994424,0.495096,0.993175
INDEL,PASS,525469,505901,19568,950136,1142,442507,0.962760,0.997751,0.465734,0.979943
SNP,ALL,3365063,3354208,10855,4119016,4533,750329,0.996774,0.998652,0.182163,0.997712
SNP,PASS,3365063,3350001,15062,3835093,1675,480902,0.995524,0.999500,0.125393,0.997508
`

func TestParseHappySummary_PrefersPASS(t *testing.T) {
	ms, err := ParseHappySummary(strings.NewReader(happySample))
	if err != nil {
		t.Fatalf("ParseHappySummary: %v", err)
	}
	byType := map[string]BenchMetrics{}
	for _, m := range ms {
		byType[m.VarType] = m
	}
	snp, ok := byType["SNP"]
	if !ok {
		t.Fatal("no SNP metrics")
	}
	// Must be the PASS row (recall 0.995524), not the ALL row (0.996774).
	if !approx(snp.Recall, 0.995524) {
		t.Fatalf("SNP recall: want PASS row 0.995524, got %v", snp.Recall)
	}
	if !approx(snp.Precision, 0.999500) {
		t.Fatalf("SNP precision: %v", snp.Precision)
	}
	if !approx(snp.F1, 0.997508) {
		t.Fatalf("SNP F1: %v", snp.F1)
	}
	if snp.TruthTP != 3350001 || snp.QueryFP != 1675 {
		t.Fatalf("SNP counts: TP=%d FP=%d", snp.TruthTP, snp.QueryFP)
	}
	indel := byType["INDEL"]
	if !approx(indel.Recall, 0.962760) {
		t.Fatalf("INDEL recall (PASS): %v", indel.Recall)
	}
}

func TestParseHappySummary_MissingColumn(t *testing.T) {
	bad := "Type,Filter\nSNP,PASS\n"
	if _, err := ParseHappySummary(strings.NewReader(bad)); err == nil {
		t.Fatal("expected error for missing metric columns")
	}
}

func TestParseHappySummary_FallbackToALL(t *testing.T) {
	// Only ALL rows present -> parser should still return them.
	onlyAll := `Type,Filter,TRUTH.TP,TRUTH.FN,QUERY.FP,METRIC.Recall,METRIC.Precision,METRIC.F1_Score
SNP,ALL,100,5,2,0.95,0.98,0.965
`
	ms, err := ParseHappySummary(strings.NewReader(onlyAll))
	if err != nil {
		t.Fatalf("ParseHappySummary: %v", err)
	}
	if len(ms) != 1 || !approx(ms[0].Recall, 0.95) {
		t.Fatalf("fallback to ALL failed: %+v", ms)
	}
}

// vcfevalSample is a captured RTG vcfeval summary.txt (fixed-width, dashes
// separator, threshold rows, final aggregate "None" row).
const vcfevalSample = `Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
----------------------------------------------------------------------------------------------------------
   17.000            3349000        3349500       1500      16000     0.9996       0.9952     0.9974
    0.000            3354208        3355000       4533      10855     0.9987       0.9968     0.9977
     None            3354208        3355000       4533      10855     0.9987       0.9968     0.9977
`

func TestParseVcfevalSummary(t *testing.T) {
	ms, err := ParseVcfevalSummary(strings.NewReader(vcfevalSample))
	if err != nil {
		t.Fatalf("ParseVcfevalSummary: %v", err)
	}
	if len(ms) != 1 {
		t.Fatalf("expected 1 aggregate metric, got %d", len(ms))
	}
	m := ms[0]
	if !approx(m.Precision, 0.9987) {
		t.Fatalf("precision: %v", m.Precision)
	}
	if !approx(m.Recall, 0.9968) {
		t.Fatalf("recall (Sensitivity): %v", m.Recall)
	}
	if !approx(m.F1, 0.9977) {
		t.Fatalf("F1 (F-measure): %v", m.F1)
	}
	if m.QueryFP != 4533 || m.TruthFN != 10855 {
		t.Fatalf("counts: FP=%d FN=%d", m.QueryFP, m.TruthFN)
	}
}

func TestParseVcfevalSummary_NoData(t *testing.T) {
	if _, err := ParseVcfevalSummary(strings.NewReader("")); err == nil {
		t.Fatal("expected error on empty summary")
	}
}

func approx(a, b float64) bool {
	d := a - b
	if d < 0 {
		d = -d
	}
	return d < 1e-6
}
