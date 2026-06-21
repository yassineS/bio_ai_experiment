// Command flagcompat is described in main.go. This file holds the
// source-scanning extractor that derives the set of CLI flag names a ported
// tool registers, directly from its Go source via the go/ast parser.
package main

import (
	"go/ast"
	"go/parser"
	"go/token"
	"io/fs"
	"path/filepath"
	"strings"
)

// flagRegistrars maps a flag-registration helper name to the zero-based
// positions of its short and long name string-literal arguments. A position of
// -1 means "that form is not present in this signature".
//
// Two registration families are recognised:
//
//   - The project's pkg/cliflag helpers, e.g.
//     cliflag.StringVar(fs, &dest, "s", "long", def, usage) — short at arg 2,
//     long at arg 3. cliflag.Var(fs, &val, "s", "long", usage) has the same
//     short/long positions.
//   - The standard library's flag package, used directly on a *flag.FlagSet
//     (fs.StringVar(&dest, "name", def, usage)) or the default CommandLine
//     (flag.Bool("name", def, usage), flag.Func("name", usage, fn)). These
//     register a single name, treated as short when one character long and as a
//     long name otherwise.
type flagRegistrar struct {
	shortArg int // index of the short-name literal, or -1
	longArg  int // index of the long-name literal, or -1
	oneArg   int // index of a single-name literal (stdlib flag), or -1
}

// cliflagFuncs are the pkg/cliflag helpers that take (fs, ptr, short, long, …).
var cliflagFuncs = map[string]flagRegistrar{
	"StringVar":   {shortArg: 2, longArg: 3, oneArg: -1},
	"IntVar":      {shortArg: 2, longArg: 3, oneArg: -1},
	"Int64Var":    {shortArg: 2, longArg: 3, oneArg: -1},
	"Uint64Var":   {shortArg: 2, longArg: 3, oneArg: -1},
	"Uint64":      {shortArg: 2, longArg: 3, oneArg: -1},
	"Float64Var":  {shortArg: 2, longArg: 3, oneArg: -1},
	"BoolVar":     {shortArg: 2, longArg: 3, oneArg: -1},
	"DurationVar": {shortArg: 2, longArg: 3, oneArg: -1},
	"Var":         {shortArg: 2, longArg: 3, oneArg: -1},
}

// stdlibVarFuncs are the stdlib XxxVar forms: f.StringVar(&dest, "name", …) —
// the name literal is at arg 1.
var stdlibVarFuncs = map[string]flagRegistrar{
	"StringVar":   {shortArg: -1, longArg: -1, oneArg: 1},
	"IntVar":      {shortArg: -1, longArg: -1, oneArg: 1},
	"Int64Var":    {shortArg: -1, longArg: -1, oneArg: 1},
	"UintVar":     {shortArg: -1, longArg: -1, oneArg: 1},
	"Uint64Var":   {shortArg: -1, longArg: -1, oneArg: 1},
	"Float64Var":  {shortArg: -1, longArg: -1, oneArg: 1},
	"BoolVar":     {shortArg: -1, longArg: -1, oneArg: 1},
	"DurationVar": {shortArg: -1, longArg: -1, oneArg: 1},
	"Var":         {shortArg: -1, longArg: -1, oneArg: 1},
	"TextVar":     {shortArg: -1, longArg: -1, oneArg: 1},
}

// stdlibPtrFuncs are the stdlib pointer-returning forms: flag.Bool("name", …),
// flag.Func("name", usage, fn) — the name literal is at arg 0.
var stdlibPtrFuncs = map[string]flagRegistrar{
	"String":   {shortArg: -1, longArg: -1, oneArg: 0},
	"Int":      {shortArg: -1, longArg: -1, oneArg: 0},
	"Int64":    {shortArg: -1, longArg: -1, oneArg: 0},
	"Uint":     {shortArg: -1, longArg: -1, oneArg: 0},
	"Uint64":   {shortArg: -1, longArg: -1, oneArg: 0},
	"Float64":  {shortArg: -1, longArg: -1, oneArg: 0},
	"Bool":     {shortArg: -1, longArg: -1, oneArg: 0},
	"Duration": {shortArg: -1, longArg: -1, oneArg: 0},
	"Func":     {shortArg: -1, longArg: -1, oneArg: 0},
	"BoolFunc": {shortArg: -1, longArg: -1, oneArg: 0},
}

// extractOurFlags walks every non-test .go file under dir and returns the set
// of distinct flag names the port registers (both short and long forms),
// normalised to their bare token (no leading dashes). It recognises the
// pkg/cliflag helpers and the stdlib flag package registration calls.
func extractOurFlags(dir string) (map[string]struct{}, error) {
	found := map[string]struct{}{}
	walkErr := filepath.WalkDir(dir, func(path string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() {
			return nil
		}
		if !strings.HasSuffix(path, ".go") || strings.HasSuffix(path, "_test.go") {
			return nil
		}
		fset := token.NewFileSet()
		file, perr := parser.ParseFile(fset, path, nil, 0)
		if perr != nil {
			return perr
		}
		ast.Inspect(file, func(n ast.Node) bool {
			switch node := n.(type) {
			case *ast.CallExpr:
				harvestRegistration(node, found)
			case *ast.SwitchStmt:
				harvestSwitchCases(node, found)
			case *ast.CompositeLit:
				harvestFlagMapKeys(node, found)
			}
			return true
		})
		return nil
	})
	if walkErr != nil {
		return nil, walkErr
	}
	return found, nil
}

// harvestRegistration records the flag names from a single pkg/cliflag or
// stdlib flag registration call, if the call is one.
func harvestRegistration(call *ast.CallExpr, found map[string]struct{}) {
	sel, ok := call.Fun.(*ast.SelectorExpr)
	if !ok {
		return
	}
	pkgIdent, ok := sel.X.(*ast.Ident)
	if !ok {
		return
	}
	fn := sel.Sel.Name
	var reg flagRegistrar
	var matched bool
	switch {
	case pkgIdent.Name == "cliflag":
		reg, matched = cliflagFuncs[fn]
	default:
		// Any other receiver (flag, fs, flags, f, …): treat as the stdlib flag
		// package. Disambiguate XxxVar (name at arg 1) from the
		// pointer-returning Xxx (name at arg 0).
		if reg, matched = stdlibVarFuncs[fn]; !matched {
			reg, matched = stdlibPtrFuncs[fn]
		}
	}
	if !matched {
		return
	}
	for _, idx := range []int{reg.shortArg, reg.longArg, reg.oneArg} {
		if idx < 0 || idx >= len(call.Args) {
			continue
		}
		if name, ok := stringLit(call.Args[idx]); ok && name != "" {
			found[normalizeFlag(name)] = struct{}{}
		}
	}
}

// argSwitchTags are the conventional identifier names a hand-rolled argument
// loop switches over (e.g. `switch arg {` / `switch name {`). A few bedtools
// ports parse flags this way instead of registering them with the flag
// package; their case clauses are the authoritative accepted-flag list.
var argSwitchTags = map[string]bool{
	"arg": true, "args": true, "a": true, "name": true,
	"flag": true, "opt": true, "optname": true, "o": true, "f": true,
}

// harvestSwitchCases records flag names from a `switch <argvar> { case "-x": … }`
// dispatch over command-line arguments. It only considers switches whose tag is
// a conventional argument identifier (argSwitchTags) AND that contain at least
// one dash-prefixed case literal, so ordinary switches over unrelated strings
// are ignored. Case literals are normalised by stripping leading dashes.
func harvestSwitchCases(sw *ast.SwitchStmt, found map[string]struct{}) {
	ident, ok := sw.Tag.(*ast.Ident)
	if !ok || !argSwitchTags[ident.Name] {
		return
	}
	var literals []string
	hasDashed := false
	for _, stmt := range sw.Body.List {
		cc, ok := stmt.(*ast.CaseClause)
		if !ok {
			continue
		}
		for _, expr := range cc.List {
			if s, ok := stringLit(expr); ok && s != "" {
				literals = append(literals, s)
				if len(s) > 0 && s[0] == '-' {
					hasDashed = true
				}
			}
		}
	}
	if !hasDashed {
		return
	}
	for _, s := range literals {
		if n := normalizeFlag(s); n != "" {
			found[n] = struct{}{}
		}
	}
}

// harvestFlagMapKeys records flag names from a flag-binding map literal such as
// `map[string]*bool{"wa": &o.writeA, …}`. Some ports keep their boolean and
// value flags in a `map[string]*T` keyed by the bare (dash-stripped) flag name.
// Only string-keyed maps whose value type is a pointer are considered, which is
// the flag-binding idiom and avoids harvesting unrelated lookup tables.
func harvestFlagMapKeys(cl *ast.CompositeLit, found map[string]struct{}) {
	mt, ok := cl.Type.(*ast.MapType)
	if !ok {
		return
	}
	if id, ok := mt.Key.(*ast.Ident); !ok || id.Name != "string" {
		return
	}
	if _, ok := mt.Value.(*ast.StarExpr); !ok {
		return
	}
	for _, elt := range cl.Elts {
		kv, ok := elt.(*ast.KeyValueExpr)
		if !ok {
			continue
		}
		if s, ok := stringLit(kv.Key); ok && s != "" {
			if n := normalizeFlag(s); n != "" {
				found[n] = struct{}{}
			}
		}
	}
}

// normalizeFlag strips leading dashes from a flag token so that "-i", "--input"
// and "input" all compare equal.
func normalizeFlag(s string) string {
	return strings.TrimLeft(s, "-")
}

// stringLit returns the unquoted value of a basic string-literal expression.
func stringLit(e ast.Expr) (string, bool) {
	lit, ok := e.(*ast.BasicLit)
	if !ok || lit.Kind != token.STRING {
		return "", false
	}
	// lit.Value includes the surrounding quotes; strip them. Use a tolerant
	// unquote that handles the common double-quoted case.
	v := lit.Value
	if len(v) >= 2 && (v[0] == '"' || v[0] == '`') {
		return v[1 : len(v)-1], true
	}
	return "", false
}
