# Go Implementation Guide

This document explains the approach, structure, and best practices for reimplementing bioinformatics tools in Go.

## Project Structure

```
bio_ai_experiment/
├── pkg/
│   └── bioformats/          # Shared format libraries
│       ├── fasta/           # FASTA format
│       ├── fastq/           # FASTQ format
│       ├── vcf/             # VCF format
│       ├── bed/             # BED format
│       └── README.md        # Format library docs
├── tools/
│   ├── seqtk/              # Example tool implementation
│   │   ├── cmd/seqtk/      # CLI entry point
│   │   ├── pkg/seqtk/      # Core library
│   │   ├── tests/          # Integration tests
│   │   ├── testdata/       # Test data
│   │   ├── docs/           # Tool-specific docs
│   │   └── README.md       # Tool documentation
│   └── [other tools]/
├── docs/
│   ├── tools/              # Tool analyses
│   └── GOLANG_GUIDE.md     # This file
└── go.mod                  # Go module definition
```

## Design Principles

### 1. Shared Libraries First

Create reusable format parsers in `pkg/htsgo/`:

- Reduces code duplication
- Ensures consistency
- Makes testing easier
- Simplifies tool implementation

Example:

```go
// Don't reimplement FASTA parsing in each tool
// Use the shared library:
import "github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"

reader := fasta.NewReader(file)
```

### 2. CLI Best Practices

Follow standard CLI conventions:

```bash
tool <command> [options] <input>
```

Features:

- Subcommands for different operations
- `-h` flag for help
- `-o` flag for output file (default: stdout)
- Clear error messages
- Exit codes (0 = success, 1 = error)

Example:

```go
fs := flag.NewFlagSet("command", flag.ExitOnError)
output := fs.String("o", "", "output file (default: stdout)")
quality := fs.Int("q", 20, "quality threshold")

fs.Usage = func() {
    fmt.Fprintf(os.Stderr, "Usage: tool command [options]\n")
}

fs.Parse(os.Args[2:])
```

### 3. Standard Format Support

Always support standard formats:

- FASTA/FASTQ for sequences
- VCF for variants
- BAM/SAM for alignments
- BED for regions
- GFF for annotations

Use available Go libraries or the shared bioformats package.

### 4. Error Handling

Provide clear, actionable error messages:

```go
if err != nil {
    fmt.Fprintf(os.Stderr, "Error reading file %s: %v\n", filename, err)
    os.Exit(1)
}
```

### 5. Testing

Aim for >80% code coverage:

```go
func TestFeature(t *testing.T) {
    // Arrange
    input := "test data"
    
    // Act
    result, err := ProcessData(input)
    
    // Assert
    if err != nil {
        t.Fatalf("Unexpected error: %v", err)
    }
    if result != expected {
        t.Errorf("Expected %v, got %v", expected, result)
    }
}
```

### 6. Documentation

Three levels of documentation:

1. **godoc Comments** - For API documentation

   ```go
   // CalculateStats computes sequence statistics.
   // Returns error if the input is invalid.
   func CalculateStats(data []byte) (*Stats, error)
   ```

2. **README.md** - User-facing documentation
   - Installation instructions
   - Usage examples
   - Command reference
   - Performance notes

3. **Parity notes** - the per-tool `README.md` plus the authoritative gap
   list in [`PARITY_ROADMAP.md`](PARITY_ROADMAP.md)
   - Original tool comparison
   - Implementation notes
   - Remaining feature gaps

## Implementation Workflow

### 1. Analyze Original Tool

Capture the analysis in the tool's `README.md` and the parity gap list:

- What does the tool do?
- What formats does it use?
- What are its strengths/weaknesses?
- What features are most important?

### 2. Design Go Implementation

Plan the structure:

```
tools/[tool]/
├── cmd/[tool]/main.go      # CLI
├── pkg/[tool]/[tool].go    # Library
├── pkg/[tool]/[tool]_test.go
└── README.md
```

### 3. Implement Shared Libraries

Before implementing the tool:

- Check if required format parsers exist
- Create/update shared libraries as needed
- Test thoroughly

### 4. Implement Core Functionality

Start with library (`pkg/[tool]/`):

```go
package toolname

// Core types
type Record struct {
    // ...
}

// Core functions
func Process(input io.Reader, output io.Writer) error {
    // ...
}
```

### 5. Add CLI Interface

Implement in `cmd/[tool]/main.go`:

```go
package main

import "github.com/yassineS/bio_ai_experiment/tools/[tool]/pkg/[tool]"

func main() {
    // Parse flags
    // Call library functions
    // Handle errors
}
```

### 6. Write Tests

Test coverage checklist:

- [ ] Happy path tests
- [ ] Error cases
- [ ] Edge cases (empty input, large files)
- [ ] Format validation
- [ ] Performance tests (optional)

### 7. Document

Create comprehensive documentation:

- [ ] godoc comments on all exported items
- [ ] README.md with examples and parity notes
- [ ] Update the parity gap list (`PARITY_ROADMAP.md`)
- [ ] Update main README if needed

## Code Style

### Go Idioms

1. **Error Handling**

   ```go
   // Good
   if err != nil {
       return err
   }
   
   // Avoid
   if err == nil {
       // ... lots of code
   }
   ```

2. **Naming**

   ```go
   // Good
   type FastqReader struct { }
   func NewReader() *FastqReader { }
   
   // Avoid
   type fastq_reader struct { }
   func new_reader() *fastq_reader { }
   ```

3. **Interface Design**

   ```go
   // Accept interfaces
   func Process(r io.Reader) error { }
   
   // Return concrete types
   func NewReader(r io.Reader) *Reader { }
   ```

### Project Conventions

1. **Package Names**: Lowercase, no underscores
   - `fasta` not `FASTA` or `fasta_parser`

2. **File Names**: Lowercase, underscores for multi-word
   - `fasta.go`, `fasta_test.go`

3. **Import Grouping**:

   ```go
   import (
       // Standard library
       "fmt"
       "io"
       
       // External
       "github.com/pkg/errors"
       
       // Internal
       "github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
   )
   ```

## Testing Strategy

### Unit Tests

Test individual functions:

```go
func TestParseRecord(t *testing.T) {
    input := ">seq1\nACGT\n"
    reader := strings.NewReader(input)
    
    record, err := ParseRecord(reader)
    if err != nil {
        t.Fatalf("Unexpected error: %v", err)
    }
    
    if record.ID != "seq1" {
        t.Errorf("Expected ID 'seq1', got '%s'", record.ID)
    }
}
```

### Integration Tests

Test full workflows:

```go
func TestEndToEnd(t *testing.T) {
    // Create temp input file
    input := createTestFile(t)
    defer os.Remove(input)
    
    // Run command
    output := runCommand(t, "comp", input)
    
    // Verify output
    if !strings.Contains(output, "3 sequences") {
        t.Errorf("Unexpected output: %s", output)
    }
}
```

### Table-Driven Tests

For multiple test cases:

```go
func TestValidation(t *testing.T) {
    tests := []struct {
        name    string
        input   string
        wantErr bool
    }{
        {"valid", "ACGT", false},
        {"invalid", "XYZ", true},
        {"empty", "", true},
    }
    
    for _, tt := range tests {
        t.Run(tt.name, func(t *testing.T) {
            err := Validate(tt.input)
            if (err != nil) != tt.wantErr {
                t.Errorf("Validate() error = %v, wantErr %v", err, tt.wantErr)
            }
        })
    }
}
```

## Performance Considerations

### 1. Streaming I/O

Always stream for large files:

```go
// Good - streaming
reader := bufio.NewReader(file)
for {
    line, err := reader.ReadString('\n')
    if err == io.EOF {
        break
    }
    // Process line
}

// Avoid - loading entire file
data, _ := ioutil.ReadFile(filename)
```

### 2. Buffer Sizes

Use appropriate buffer sizes:

```go
scanner := bufio.NewScanner(file)
buf := make([]byte, 0, 64*1024)          // 64KB initial
scanner.Buffer(buf, 10*1024*1024)        // 10MB max
```

### 3. Avoid Allocations

Reuse buffers when possible:

```go
// Good - reuse buffer
var buf bytes.Buffer
for _, record := range records {
    buf.Reset()
    buf.WriteString(record.Sequence)
    // Process buf
}

// Avoid - allocate each time
for _, record := range records {
    s := string(record.Sequence)  // Allocation!
    // Process s
}
```

### 4. Profiling

Use Go's built-in profiling:

```bash
go test -cpuprofile=cpu.prof -bench=.
go tool pprof cpu.prof
```

## Common Patterns

### Reading a File

```go
func ProcessFile(filename string) error {
    file, err := os.Open(filename)
    if err != nil {
        return fmt.Errorf("open %s: %w", filename, err)
    }
    defer file.Close()
    
    reader := fasta.NewReader(file)
    for {
        record, err := reader.Read()
        if err == io.EOF {
            break
        }
        if err != nil {
            return fmt.Errorf("read record: %w", err)
        }
        // Process record
    }
    return nil
}
```

### Writing a File

```go
func WriteRecords(filename string, records []*Record) error {
    file, err := os.Create(filename)
    if err != nil {
        return fmt.Errorf("create %s: %w", filename, err)
    }
    defer file.Close()
    
    writer := fasta.NewWriter(file, 80)
    defer writer.Flush()
    
    return writer.WriteAll(records)
}
```

### Command Pattern

```go
func mainCommand() {
    if len(os.Args) < 2 {
        printUsage()
        os.Exit(1)
    }
    
    command := os.Args[1]
    switch command {
    case "comp":
        compCommand()
    case "convert":
        convertCommand()
    default:
        fmt.Fprintf(os.Stderr, "Unknown command: %s\n", command)
        os.Exit(1)
    }
}
```

## Checklist for New Tools

Before considering a tool complete:

### Code

- [ ] Core functionality implemented
- [ ] CLI interface implemented
- [ ] Error handling comprehensive
- [ ] Uses shared bioformats libraries
- [ ] Follows Go best practices
- [ ] No external dependencies (if possible)

### Testing

- [ ] Unit tests written
- [ ] Test coverage >80%
- [ ] Edge cases tested
- [ ] Integration tests (if applicable)
- [ ] All tests passing

### Documentation

- [ ] godoc comments on exports
- [ ] README.md with examples
- [ ] Tool analysis document
- [ ] Usage examples
- [ ] Migration guide (if replacing existing tool)

### Quality

- [ ] `go vet` passes
- [ ] `gofmt` applied
- [ ] No compiler warnings
- [ ] Performance acceptable
- [ ] Memory usage reasonable

## Examples

See `tools/seqtk/` for a complete reference implementation demonstrating:

- Shared library usage
- CLI design
- Testing strategy
- Documentation
- Best practices

## Resources

- [Effective Go](https://golang.org/doc/effective_go)
- [Go Code Review Comments](https://github.com/golang/go/wiki/CodeReviewComments)
- [Standard Go Project Layout](https://github.com/golang-standards/project-layout)

## Getting Help

- Check existing implementations in `tools/`
- Review `pkg/htsgo/` for format handling
- Read each tool's `README.md` and the gap list in `docs/PARITY_ROADMAP.md`
- Open an issue for questions
