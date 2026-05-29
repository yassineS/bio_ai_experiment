// CLI runner for the `bcftools plugin` subcommand and the `bcftools +<name>`
// shorthand. Unlike upstream bcftools, which loads plugins as native shared
// objects via dlopen, this port runs a plugin as an ordinary child process
// and streams VCF text over its standard streams. See docs/PLUGIN_PROTOCOL.md.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const pluginUsage = `bcftools plugin - run a user plugin over a VCF/BCF.

Usage:
  bcftools plugin <name>  [plugin-opts] [-- <input>]
  bcftools +<name>        [plugin-opts] [-- <input>]
  bcftools plugin -l

A plugin is an executable found in one of the colon-separated directories
named by the BCFTOOLS_PLUGINS environment variable. The host streams the
input VCF as uncompressed text to the plugin's stdin and reads the plugin's
stdout back as VCF; the plugin is therefore "a filter from VCF on stdin to
VCF on stdout". Plugin options after the name are passed verbatim as the
plugin's argv. Everything after a literal -- is the host input file (and
optional regions). See docs/PLUGIN_PROTOCOL.md for the contract.

Options (host side; must precede the plugin name or --):
  -l, --list-plugins         List discoverable plugins and exit.
  -lv                        List plugins verbosely (path + --about line).
  -o, --output PATH          Output file (default stdout).
  -O, --output-type {v|z}    v=VCF (default), z=VCF.gz. (b/u need a BCF writer.)
  -r, --regions chr[:b-e]    Restrict input records fed to the plugin.
  -R, --regions-file PATH    BED-like regions file.
      --compression-level N  gzip level for -O z output.
  -?, --help                 Show this help.
      --version              Show version.
`

// runPlugin implements `bcftools plugin <name> ...`. The plugin name and the
// plugin's own arguments are NOT parsed as host flags: host flags are only
// recognised before the plugin name (or before a `--`). pluginName, when
// non-empty, is the name supplied via the `+<name>` shorthand, in which case
// args holds everything after `+<name>`.
func runPlugin(args []string, pluginName string) int {
	// Split host options from the plugin name + plugin arguments. The first
	// non-flag token is the plugin name (unless supplied via `+name`); a
	// literal `--` terminates plugin args and introduces the host input.
	var hostArgs []string
	var pluginArgs []string
	var inputArgs []string
	{
		rest := args
		if pluginName == "" {
			// Find the plugin name: the first token that is not a host flag.
			i := 0
			for i < len(rest) {
				tok := rest[i]
				if tok == "--" {
					break
				}
				if len(tok) > 0 && tok[0] == '-' {
					hostArgs = append(hostArgs, tok)
					// -o/-O/-r/-R/-l take a value as the next token.
					if needsValue(tok) && i+1 < len(rest) {
						i++
						hostArgs = append(hostArgs, rest[i])
					}
					i++
					continue
				}
				// First bare word is the plugin name.
				pluginName = tok
				rest = rest[i+1:]
				goto haveName
			}
			// No plugin name found (only host flags / a `--`): keep rest as-is.
			rest = rest[len(hostArgs):]
		}
	haveName:
		// Everything up to a `--` is plugin argv; after `--` is host input.
		for i := 0; i < len(rest); i++ {
			if rest[i] == "--" {
				inputArgs = rest[i+1:]
				break
			}
			pluginArgs = append(pluginArgs, rest[i])
		}
	}

	fs := flag.NewFlagSet("bcftools plugin", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		listPlugins   bool
		listVerbose   bool
		outputPath    string
		outputType    string
		regions       string
		regionsFile   string
		compressLevel int
		showHelp      bool
		showVer       bool
	)
	cliflag.BoolVar(fs, &listPlugins, "l", "list-plugins", false, "List plugins")
	fs.BoolVar(&listVerbose, "lv", false, "List plugins verbosely")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&compressLevel, "compression-level", -1, "gzip level")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, hostArgs); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, pluginUsage)
		return 2
	}
	if showHelp {
		fmt.Print(pluginUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	if listPlugins || listVerbose {
		verbose := listVerbose
		plugins := bcftools.ListPlugins(verbose)
		if len(plugins) == 0 {
			fmt.Fprintf(os.Stderr, "bcftools plugin: no plugins found in %s\n", "BCFTOOLS_PLUGINS")
			return 0
		}
		fmt.Print(bcftools.FormatPluginList(plugins, verbose))
		return 0
	}

	if pluginName == "" {
		fmt.Fprintln(os.Stderr, "bcftools plugin: missing plugin name")
		fmt.Fprint(os.Stderr, pluginUsage)
		return 2
	}

	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	opts := bcftools.PluginOptions{
		Name:          pluginName,
		Args:          pluginArgs,
		OutputFormat:  format,
		CompressLevel: compressLevel,
		RegionsFile:   regionsFile,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if len(inputArgs) > 0 {
		opts.InputFile = inputArgs[0]
		opts.Regions = append(opts.Regions, inputArgs[1:]...)
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools plugin: %v\n", err)
		return 1
	}
	defer out.Close()

	if err := bcftools.RunPlugin(opts, out, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools plugin: %v\n", err)
		return 1
	}
	return 0
}

// needsValue reports whether a host plugin flag consumes the following token
// as its value (used while splitting host flags from plugin arguments).
func needsValue(flagTok string) bool {
	switch flagTok {
	case "-o", "--output", "-O", "--output-type",
		"-r", "--regions", "-R", "--regions-file",
		"--compression-level":
		return true
	}
	// `-l` is the boolean --list-plugins here, so it takes no value.
	return false
}
