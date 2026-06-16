// Native registration of the upstream `+mendelian2` plugin
// (plugins/mendelian2.c). The plugin form shares its entire engine with the
// existing `bcftools mendelian2` subcommand: rather than re-deriving any
// Mendelian logic, this file parses the plugin's own argv and delegates to the
// shared Mendelian2 / Mendelian2File engine, so `bcftools +mendelian2 ...`
// produces exactly what `bcftools mendelian2 ...` (and therefore upstream's
// `+mendelian2`) produces.
//
// mendelian2 is a run()-style plugin (its options precede the input file with
// no `--` separator). It is registered as a fullPlugin so runNativePlugin hands
// it the whole invocation. The supported options mirror the subcommand:
//
//	-p/--pfm [1X:|2X:]P,F,M    single-trio shortcut
//	-P/--ped FILE             PED file
//	-m/--mode c|[adeEgmMS]     output mode (default c)
//	-i/--include / -e/--exclude EXPR  record-level filter (no per-sample yet)
//	-o/--output FILE          output file (default stdout)
//	-O/--output-type u|b|v|z   output container
//	--rules ASSEMBLY / --rules-file FILE   inheritance rules
//	-v/--verbosity / --no-version          accepted, no effect on parity
//
// Region/target selection (-r/-R/-t/-T) is supported: the shared
// regionTargetFilter is parsed by the framework and applied to the input
// records before the Mendelian accounting. -o writes the result to a file, and
// -W/--write-index indexes it (a CSI by default, a TBI for -W=tbi on VCF.gz)
// when the chosen mode emits VCF/BCF output; the text-counts mode (-m c) emits
// no VCF and is left unindexed, matching upstream. The
// --regions-overlap/--targets-overlap tuning knobs are still reported as a clean
// unsupported Init error.
package bcftools

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("mendelian2", func() NativePlugin { return &mendelian2Plugin{} }) }

// mendelian2Plugin adapts the `+mendelian2` plugin form onto the shared
// Mendelian2 engine. It carries the shared region/target selection; RunFull does
// all the work.
type mendelian2Plugin struct {
	rt regionTargetFilter
}

// SetRegionTarget records the shared -r/-R/-t/-T selection the framework parsed
// out of mendelian2's argv; it is forwarded to the Mendelian2 engine.
func (p *mendelian2Plugin) SetRegionTarget(f regionTargetFilter) { p.rt = f }

// Name returns the plugin name.
func (p *mendelian2Plugin) Name() string { return "mendelian2" }

// RegionTargetCaps opts mendelian2 into the shared -r/-R/-t/-T region/target
// filter, forwarded to the Mendelian2 engine via SetRegionTarget.
func (p *mendelian2Plugin) RegionTargetCaps() regionTargetCaps { return allRegionTargetCaps }

// About returns the one-line description, matching mendelian2.c about().
func (p *mendelian2Plugin) About() string {
	return "Count Mendelian consistent / inconsistent genotypes.\n"
}

// RunStyle reports that mendelian2 is a run()-style plugin (options precede the
// input file, no `--` separator).
func (p *mendelian2Plugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of mendelian2's value-taking flags consumes
// the following CLI token, used by the host to split the input-file positional.
func (p *mendelian2Plugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-p", "--pfm", "-P", "--ped", "-m", "--mode",
		"-i", "--include", "-e", "--exclude", "-o", "--output",
		"-O", "--output-type", "--rules", "--rules-file",
		"-r", "--regions", "-R", "--regions-file",
		"-t", "--targets", "-T", "--targets-file",
		"--regions-overlap", "--targets-overlap", "-v", "--verbosity":
		return true
	}
	return false
}

// The NativePlugin lifecycle methods are unused (RunFull owns the invocation),
// but are required to satisfy the interface for registration.

// Init is unused: RunFull bypasses the per-record pipeline.
func (p *mendelian2Plugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return hdr, nil
}

// Process is unused.
func (p *mendelian2Plugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy is unused.
func (p *mendelian2Plugin) Destroy() error { return nil }

// RunFull parses the plugin argv and delegates to the shared Mendelian2File
// engine, which owns input reading and output writing (text count summary
// and/or filtered VCF/BCF). This is the faithful reuse path: no Mendelian logic
// is duplicated here.
func (p *mendelian2Plugin) RunFull(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	var (
		pfm, ped    string
		modeStr     = "c"
		includeExpr string
		excludeExpr string
		outputType  = "v"
		outputFile  string
		rules       string
		rulesFile   string
		writeIndex  = writeIndexOff
	)
	args := opts.Args
	for i := 0; i < len(args); i++ {
		a := args[i]
		val := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("mendelian2: %s requires a value", a)
			}
			i++
			return args[i], nil
		}
		var err error
		switch a {
		case "-p", "--pfm":
			pfm, err = val()
		case "-P", "--ped":
			ped, err = val()
		case "-m", "--mode":
			modeStr, err = val()
		case "-i", "--include":
			includeExpr, err = val()
		case "-e", "--exclude":
			excludeExpr, err = val()
		case "-o", "--output":
			outputFile, err = val()
		case "-O", "--output-type":
			outputType, err = val()
		case "--rules":
			rules, err = val()
		case "--rules-file":
			rulesFile, err = val()
		case "-v", "--verbosity":
			_, err = val() // accepted, no effect on parity
		case "--no-version":
			// accepted, no effect
		case "--regions-overlap", "--targets-overlap":
			return fmt.Errorf("mendelian2: %s is not supported by the native plugin", a)
		default:
			if sel, handled, werr := parseWriteIndexArg(a); handled {
				if werr != nil {
					return fmt.Errorf("mendelian2: %w", werr)
				}
				writeIndex = sel
				continue
			}
			// Attached -O<x> / -o<path> getopt forms (e.g. `-Oz`, `-oout.vcf`).
			if strings.HasPrefix(a, "-O") && len(a) > 2 {
				outputType = a[2:]
				continue
			}
			if strings.HasPrefix(a, "-o") && len(a) > 2 {
				outputFile = a[2:]
				continue
			}
			return fmt.Errorf("mendelian2: unsupported option %q", a)
		}
		if err != nil {
			return err
		}
	}

	if pfm == "" && ped == "" {
		return fmt.Errorf("mendelian2: missing the -p or -P option")
	}
	if pfm != "" && ped != "" {
		return fmt.Errorf("mendelian2: -p/--pfm and -P/--ped are mutually exclusive")
	}
	if rules != "" && rulesFile != "" {
		return fmt.Errorf("mendelian2: --rules and --rules-file are mutually exclusive")
	}

	modeBits, err := ParseMendelian2Mode(modeStr)
	if err != nil {
		return fmt.Errorf("mendelian2: %w", err)
	}
	format, err := ParseOutputFormat(outputType)
	if err != nil {
		return fmt.Errorf("mendelian2: %w", err)
	}

	var ruleSet *MendelianRules
	switch {
	case rulesFile != "":
		rs, err := LoadMendelianRulesFile(rulesFile)
		if err != nil {
			return fmt.Errorf("mendelian2: %w", err)
		}
		ruleSet = rs
	case rules != "":
		rs, err := LoadMendelianRulesByName(rules)
		if err != nil {
			return fmt.Errorf("mendelian2: %w", err)
		}
		ruleSet = rs
	}

	m2 := Mendelian2Options{
		PEDFile:       ped,
		Mode:          modeBits,
		IncludeExpr:   includeExpr,
		ExcludeExpr:   excludeExpr,
		OutputFormat:  format,
		CompressLevel: opts.CompressLevel,
		Threads:       opts.Threads,
		Rules:         ruleSet,
		RegionTarget:  p.rt,
	}
	if pfm != "" {
		parsed, perr := ParseMendelian2PFM(pfm)
		if perr != nil {
			return fmt.Errorf("mendelian2: %w", perr)
		}
		m2.PFM = &parsed
	}

	// When -o names a file, write there directly; otherwise stream to the host
	// stdout. -W indexes the file, but only when the chosen mode actually
	// produces VCF/BCF output (any non-count mode); the text-counts mode
	// (-m c) writes no VCF and is left unindexed, matching upstream.
	dst := out
	var outFile *os.File
	toFile := outputFile != "" && outputFile != "-"
	producesVCF := modeBits&^Mendelian2Count != 0
	if toFile {
		f, ferr := os.Create(outputFile)
		if ferr != nil {
			return fmt.Errorf("mendelian2: cannot write to %q: %w", outputFile, ferr)
		}
		outFile = f
		dst = f
	} else if writeIndex != writeIndexOff && producesVCF {
		// Indexing stdout is impossible — reproduce upstream's error.
		return fmt.Errorf("mendelian2: failed to initialise index for -")
	}

	input := opts.InputFile
	if input == "" {
		input = "-"
	}
	if _, err = Mendelian2File(input, dst, m2); err != nil {
		if outFile != nil {
			_ = outFile.Close()
		}
		return err
	}
	if outFile != nil {
		if err := outFile.Close(); err != nil {
			return err
		}
		if producesVCF {
			if err := writeIndexFor(outputFile, format, writeIndex); err != nil {
				return fmt.Errorf("mendelian2: %w", err)
			}
		}
	}
	return nil
}
