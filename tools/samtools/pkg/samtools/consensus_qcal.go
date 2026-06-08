// Package samtools — samtools consensus quality-calibration tables.
//
// Upstream bam_consensus.c bundles six static_qcal[] presets (FLAT,
// HiFi, HiSeq, ONT R10.4 super, ONT R10.4 duplex, Ultima) each with
// three 101-element arrays — smap (substitution), umap (undercall),
// omap (overcall) — that remap the input Phred quality score before
// the bayesian posterior calculation. Selection happens via
// `-t/--qual-calibration` (file or `:preset` shorthand) or
// `-X/--config` (preset that also overrides several bayesian knobs).
//
// This file ports the six presets verbatim and provides
// loadQcalFile for the file shape (`QUAL <q> <smap> <umap> <omap>`
// per line, comments with `#`, ascending q). selectQcalPreset
// maps an option string to a preset table.
package samtools

import (
	"bufio"
	"fmt"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// qcalTables holds the three per-quality remap arrays.
type qcalTables struct {
	S [101]int // substitution
	U [101]int // undercall (DEL)
	O [101]int // overcall (INS)
}

// QcalTables is the exported alias used by consensus CLI runners
// (the unexported type's fields would otherwise be inaccessible
// from the cmd package).
type QcalTables = qcalTables

// ConfigPreset is the exported alias for the bayesian bundle that
// `-X NAME` selects.
type ConfigPreset = configPreset

// LoadQcalFile is the public wrapper around the unexported
// loadQcalFile.
func LoadQcalFile(path string) (*QcalTables, error) { return loadQcalFile(path) }

// qcalFlat is the identity FLAT preset (every quality maps to
// itself, capped at 99). Matches upstream static_qcal[QCAL_FLAT].
var qcalFlat = func() qcalTables {
	var t qcalTables
	for i := 0; i < 101; i++ {
		v := i
		if v > 99 {
			v = 99
		}
		t.S[i] = v
		t.U[i] = v
		t.O[i] = v
	}
	return t
}()

// qcalHiFi is upstream static_qcal[QCAL_HIFI] (bam_consensus.c:486-520).
var qcalHiFi = qcalTables{
	S: [101]int{
		10, 11, 11, 12, 13, 14, 15, 16, 18, 19,
		20, 21, 22, 23, 24, 25, 27, 28, 29, 30,
		31, 32, 33, 33, 34, 35, 36, 36, 37, 38,
		38, 39, 39, 40, 40, 41, 41, 41, 41, 42,
		42, 42, 42, 43, 43, 43, 43, 43, 43, 43,
		44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
		44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
		44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
		44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
		44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
	},
	U: [101]int{
		4, 4, 4, 4, 5, 6, 6, 7, 8, 9,
		10, 11, 11, 12, 13, 14, 15, 15, 16, 17,
		18, 19, 19, 20, 20, 21, 22, 23, 23, 24,
		25, 25, 25, 26, 26, 26, 27, 27, 28, 28,
		28, 28, 27, 27, 27, 28, 28, 28, 28, 27,
		27, 27, 27, 27, 27, 27, 27, 27, 27, 27,
		27, 27, 26, 26, 25, 26, 26, 27, 27, 27,
		26, 26, 26, 26, 26, 26, 26, 26, 27, 27,
		28, 29, 28, 28, 28, 27, 27, 27, 27, 27,
		27, 28, 28, 30, 30, 30, 30, 30, 30, 30,
	},
	O: [101]int{
		8, 8, 8, 8, 9, 10, 11, 12, 13, 14,
		15, 15, 16, 17, 18, 19, 19, 20, 20, 21,
		21, 22, 22, 23, 23, 23, 24, 24, 24, 25,
		25, 25, 25, 25, 25, 26, 26, 26, 26, 27,
		27, 27, 27, 27, 27, 28, 28, 28, 28, 28,
		29, 29, 29, 29, 29, 29, 30, 30, 30, 30,
		30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
		30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
		30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
		30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
	},
}

// qcalHiSeq is upstream static_qcal[QCAL_HISEQ] (bam_consensus.c:522-556).
var qcalHiSeq = qcalTables{
	S: [101]int{
		2, 2, 2, 3, 3, 4, 5, 5, 6, 7,
		8, 9, 10, 11, 11, 12, 13, 14, 15, 16,
		17, 17, 18, 19, 20, 21, 22, 22, 23, 24,
		25, 26, 27, 28, 28, 29, 30, 31, 32, 33,
		34, 34, 35, 36, 37, 38, 39, 39, 40, 41,
		42, 43, 44, 45, 45, 46, 47, 48, 49, 50,
		51, 51, 52, 53, 54, 55, 56, 56, 57, 58,
		59, 60, 61, 62, 62, 63, 64, 65, 66, 67,
		68, 68, 69, 70, 71, 72, 73, 73, 74, 75,
		76, 77, 78, 79, 79, 80, 81, 82, 83, 84,
	},
	U: [101]int{
		1, 2, 3, 4, 5, 7, 8, 9, 10, 11,
		13, 14, 15, 16, 17, 19, 20, 21, 22, 23,
		25, 26, 27, 28, 29, 31, 32, 33, 34, 35,
		37, 38, 39, 40, 41, 43, 44, 45, 46, 47,
		49, 50, 51, 52, 53, 55, 56, 57, 58, 59,
		61, 62, 63, 64, 65, 67, 68, 69, 70, 71,
		73, 74, 75, 76, 77, 79, 80, 81, 82, 83,
		85, 86, 87, 88, 89, 91, 92, 93, 94, 95,
		97, 98, 99, 100, 101, 103, 104, 105, 106, 107,
		109, 110, 111, 112, 113, 115, 116, 117, 118, 119,
	},
	O: [101]int{
		1, 2, 3, 4, 5, 7, 8, 9, 10, 11,
		13, 14, 15, 16, 17, 19, 20, 21, 22, 23,
		25, 26, 27, 28, 29, 31, 32, 33, 34, 35,
		37, 38, 39, 40, 41, 43, 44, 45, 46, 47,
		49, 50, 51, 52, 53, 55, 56, 57, 58, 59,
		61, 62, 63, 64, 65, 67, 68, 69, 70, 71,
		73, 74, 75, 76, 77, 79, 80, 81, 82, 83,
		85, 86, 87, 88, 89, 91, 92, 93, 94, 95,
		97, 98, 99, 100, 101, 103, 104, 105, 106, 107,
		109, 110, 111, 112, 113, 115, 116, 117, 118, 119,
	},
}

// qcalONTSup is upstream static_qcal[QCAL_ONT_R10_4_SUP] (bam_consensus.c:557-591).
var qcalONTSup = qcalTables{
	S: [101]int{
		0, 2, 2, 2, 3, 4, 4, 5, 6, 7,
		7, 8, 9, 12, 13, 14, 15, 15, 16, 17,
		18, 19, 20, 22, 24, 25, 26, 27, 28, 29,
		30, 31, 33, 34, 36, 37, 38, 38, 39, 39,
		40, 40, 40, 40, 40, 40, 40, 41, 40, 40,
		41, 41, 40, 40, 40, 40, 41, 40, 40, 40,
		40, 41, 41, 40, 40, 41, 40, 40, 39, 41,
		40, 41, 40, 40, 41, 41, 41, 40, 40, 40,
		40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
		40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
	},
	U: [101]int{
		0, 2, 2, 2, 3, 4, 5, 6, 7, 8,
		8, 9, 9, 10, 10, 10, 11, 12, 12, 13,
		13, 13, 14, 14, 15, 16, 16, 17, 18, 18,
		19, 19, 20, 21, 22, 23, 24, 25, 25, 25,
		25, 25, 25, 25, 25, 25, 26, 26, 26, 26,
		26, 26, 26, 26, 27, 27, 27, 27, 27, 27,
		27, 27, 27, 27, 27, 27, 27, 28, 28, 28,
		28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
		28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
		28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
	},
	O: [101]int{
		0, 4, 6, 6, 6, 7, 7, 8, 9, 9,
		9, 10, 10, 11, 11, 12, 12, 13, 13, 14,
		15, 15, 15, 16, 16, 17, 17, 18, 18, 19,
		19, 20, 20, 21, 22, 22, 23, 23, 24, 24,
		24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
		24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
		24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
		24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
		24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
		24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
	},
}

// qcalONTDup is upstream static_qcal[QCAL_ONT_R10_4_DUP] — the
// upstream comment notes it's currently a copy of HiFi
// (bam_consensus.c:592-626).
var qcalONTDup = qcalHiFi

// qcalUltima is upstream static_qcal[QCAL_ULTIMA] (bam_consensus.c:627-661).
var qcalUltima = qcalTables{
	S: [101]int{
		2, 2, 3, 4, 5, 6, 6, 7, 8, 9,
		10, 10, 11, 12, 13, 14, 14, 15, 16, 17,
		18, 18, 19, 21, 22, 23, 23, 24, 25, 26,
		27, 27, 28, 29, 30, 31, 31, 32, 33, 34,
		35, 35, 36, 37, 38, 39, 39, 40, 42, 43,
		44, 44, 45, 46, 47, 48, 48, 49, 50, 51,
		52, 52, 53, 54, 55, 56, 56, 57, 58, 59,
		60, 60, 61, 63, 64, 65, 65, 66, 67, 68,
		69, 69, 70, 71, 72, 73, 73, 74, 75, 76,
		77, 77, 78, 79, 80, 81, 81, 82, 84, 85,
	},
	U: [101]int{
		1, 1, 2, 2, 3, 3, 4, 4, 4, 4,
		5, 5, 6, 6, 7, 7, 8, 8, 9, 10,
		10, 10, 11, 12, 13, 13, 13, 14, 15, 16,
		16, 16, 17, 18, 18, 19, 19, 20, 20, 21,
		21, 22, 22, 22, 22, 23, 23, 24, 24, 25,
		25, 25, 25, 25, 25, 25, 26, 26, 26, 26,
		26, 26, 27, 27, 27, 27, 27, 27, 27, 27,
		27, 28, 28, 28, 28, 28, 28, 28, 28, 28,
		28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
		28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
	},
	O: [101]int{
		1, 1, 2, 2, 3, 3, 4, 4, 4, 4,
		5, 5, 6, 6, 7, 7, 8, 8, 9, 10,
		10, 10, 11, 12, 13, 13, 13, 14, 15, 16,
		16, 16, 17, 18, 18, 19, 19, 20, 20, 21,
		21, 22, 22, 22, 22, 23, 23, 24, 24, 25,
		25, 25, 25, 25, 25, 25, 26, 26, 26, 26,
		26, 26, 27, 27, 27, 27, 27, 27, 27, 27,
		27, 28, 28, 28, 28, 28, 28, 28, 28, 28,
		28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
		28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
	},
}

// selectQcalPreset returns the named preset table, or nil when the
// name is unknown. Matches the `:NAME` shorthand upstream load_qcal
// accepts on the `-t` option.
func selectQcalPreset(name string) *qcalTables {
	switch strings.ToLower(name) {
	case ":flat", "flat":
		t := qcalFlat
		return &t
	case ":hifi", "hifi":
		t := qcalHiFi
		return &t
	case ":hiseq", "hiseq":
		t := qcalHiSeq
		return &t
	case ":r10.4_sup", "r10.4_sup":
		t := qcalONTSup
		return &t
	case ":r10.4_dup", "r10.4_dup":
		t := qcalONTDup
		return &t
	case ":ultima", "ultima":
		t := qcalUltima
		return &t
	}
	return nil
}

// loadQcalFile ports upstream load_qcal (bam_consensus.c:672) for
// non-preset paths. The file is `QUAL <q> <smap> <umap> <omap>` per
// line; `#` lines are comments; quality values must ascend; gaps
// between defined quality points are filled by carrying the
// previous mapping forward (linear hold, not interpolation —
// matching upstream's `while (v > last_qual) { ...smap[last+1] =
// smap[last]; ... }` loop).
func loadQcalFile(path string) (*qcalTables, error) {
	t := qcalFlat
	if path == "" || strings.EqualFold(path, ":flat") {
		return &t, nil
	}
	if pre := selectQcalPreset(path); pre != nil {
		return pre, nil
	}
	r, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, fmt.Errorf("samtools consensus: open --qual-calibration %s: %w", path, err)
	}
	defer r.Close()
	sc := bufio.NewScanner(r)
	last := 0
	maxQ := 0
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || line[0] == '#' {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) != 5 || fields[0] != "QUAL" {
			return nil, fmt.Errorf("samtools consensus: bad qcal line %q (want `QUAL <q> <s> <u> <o>`)", line)
		}
		q, e1 := strconv.Atoi(fields[1])
		sv, e2 := strconv.Atoi(fields[2])
		uv, e3 := strconv.Atoi(fields[3])
		ov, e4 := strconv.Atoi(fields[4])
		if e1 != nil || e2 != nil || e3 != nil || e4 != nil {
			return nil, fmt.Errorf("samtools consensus: bad qcal line %q", line)
		}
		for q > last && last+1 < 101 {
			t.S[last+1] = t.S[last]
			t.U[last+1] = t.U[last]
			t.O[last+1] = t.O[last]
			last++
		}
		if q >= 0 && q < 100 {
			t.S[q] = sv
			t.U[q] = uv
			t.O[q] = ov
		}
		if q < maxQ {
			return nil, fmt.Errorf("samtools consensus: qcal file is not in ascending quality order")
		}
		maxQ = q
		last = q
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	for i := maxQ + 1; i < 101; i++ {
		t.S[i] = t.S[maxQ]
		t.U[i] = t.U[maxQ]
		t.O[i] = t.O[maxQ]
	}
	return &t, nil
}

// configPreset captures the bayesian-knob bundle that upstream's
// `-X NAME` applies (bam_consensus.c:3171-3229). Returns nil for
// unrecognised names. The caller MUST also load the matching qcal
// preset; the field naming mirrors ConsensusOptions exactly.
type configPreset struct {
	QCal          *qcalTables
	Mode          ConsensusMode
	HomopolyFix   float64
	HomopolyRedux float64
	LowMQual      int
	ScaleMQ       float64
	HetScale      float64
	HasHomoFix    bool
	HasHomoRedux  bool
	HasLowMQ      bool
	HasScaleMQ    bool
	HasHetScale   bool
}

// SelectConsensusConfig returns the upstream `-X NAME` parameter
// bundle, or nil when name is empty / unrecognised.
func SelectConsensusConfig(name string) *configPreset {
	switch strings.ToLower(name) {
	case "hifi":
		t := qcalHiFi
		return &configPreset{
			QCal: &t, Mode: ConsensusModeBayesian,
			HomopolyFix: 0.3, HasHomoFix: true,
			HomopolyRedux: 0.01, HasHomoRedux: true,
			LowMQual: 5, HasLowMQ: true,
			ScaleMQ: 1.5, HasScaleMQ: true,
			HetScale: 0.37, HasHetScale: true,
		}
	case "hiseq":
		t := qcalHiSeq
		return &configPreset{
			QCal: &t, Mode: ConsensusModeBayesian,
			HomopolyRedux: 0.01, HasHomoRedux: true,
		}
	case "r10.4_sup":
		t := qcalONTSup
		return &configPreset{
			QCal: &t, Mode: ConsensusModeBayesian,
			HomopolyFix: 0.3, HasHomoFix: true,
			HomopolyRedux: 0.01, HasHomoRedux: true,
			LowMQual: 5, HasLowMQ: true,
			ScaleMQ: 1.5, HasScaleMQ: true,
			HetScale: 0.37, HasHetScale: true,
		}
	case "r10.4_dup":
		t := qcalONTDup
		return &configPreset{
			QCal: &t, Mode: ConsensusModeBayesian,
			HomopolyFix: 0.3, HasHomoFix: true,
			HomopolyRedux: 0.01, HasHomoRedux: true,
			LowMQual: 5, HasLowMQ: true,
			ScaleMQ: 1.5, HasScaleMQ: true,
			HetScale: 0.37, HasHetScale: true,
		}
	case "ultima":
		t := qcalUltima
		return &configPreset{
			QCal: &t, Mode: ConsensusModeBayesian,
			HomopolyFix: 0.3, HasHomoFix: true,
			HomopolyRedux: 0.01, HasHomoRedux: true,
			HetScale: 0.37, HasHetScale: true,
			ScaleMQ: 2, HasScaleMQ: true,
			LowMQual: 10, HasLowMQ: true,
		}
	}
	return nil
}
