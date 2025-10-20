# Bioinformatics Tools - Go Implementations

This directory contains Go reimplementations of popular bioinformatics tools. Each tool is designed to be performant, well-documented, and thoroughly tested.

## Available Tools

Currently, this directory is empty. Tools will be added as they are recoded from their original implementations.

### Planned Tools

The following tools are being considered for implementation (subject to change based on analysis):

1. Sequence alignment tools
2. Quality control tools
3. Format converters
4. Statistical analysis tools
5. Visualization tools
6. Assembly tools
7. Variant calling tools
8. Annotation tools
9. Read mappers
10. De novo assemblers

See [analysis/](../analysis/) directory for detailed assessments of each tool.

## Tool Structure

Each tool follows this standard structure:

```
[tool-name]/
├── cmd/
│   └── [tool-name]/
│       └── main.go           # CLI entry point
├── pkg/
│   └── [tool-name]/
│       ├── core.go           # Core functionality
│       ├── io.go             # Input/output handling
│       ├── process.go        # Processing logic
│       └── utils.go          # Utility functions
├── tests/
│   ├── unit/                 # Unit tests
│   ├── integration/          # Integration tests
│   ├── edge_cases/           # Edge case tests
│   └── benchmarks/           # Performance benchmarks
├── testdata/
│   ├── input/                # Test input files
│   └── expected/             # Expected output files
├── docs/
│   ├── API.md                # API documentation
│   ├── USAGE.md              # Usage guide
│   ├── EXAMPLES.md           # Code examples
│   └── MIGRATION.md          # Migration from original
├── go.mod
├── go.sum
└── README.md
```

## Building Tools

### Prerequisites

- Go 1.21 or later
- Git

### Building a Specific Tool

```bash
cd [tool-name]
go build ./cmd/[tool-name]
```

### Building All Tools

```bash
# From the tools directory
for dir in */; do
    (cd "$dir" && go build ./...)
done
```

### Installing a Tool

```bash
cd [tool-name]
go install ./cmd/[tool-name]
```

## Testing Tools

### Running Tests

```bash
cd [tool-name]

# Run all tests
go test ./...

# Run with coverage
go test -cover ./...

# Run specific test suite
go test ./tests/unit/...
go test ./tests/integration/...
go test ./tests/edge_cases/...

# Run benchmarks
go test -bench=. ./tests/benchmarks/...
```

### Test Coverage

All tools aim for >80% test coverage. Check coverage with:

```bash
go test -coverprofile=coverage.out ./...
go tool cover -html=coverage.out
```

## Documentation

Each tool includes comprehensive documentation:

- **README.md**: Tool overview and quick start
- **docs/API.md**: Complete API reference
- **docs/USAGE.md**: Detailed usage guide with examples
- **docs/EXAMPLES.md**: Code examples and common patterns
- **docs/MIGRATION.md**: Guide for users of the original tool

## Performance

All tools are benchmarked against their original implementations. Performance results are documented in each tool's README.md.

### Running Benchmarks

```bash
cd [tool-name]
go test -bench=. -benchmem ./tests/benchmarks/
```

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines on contributing to tool development.

### Adding a New Tool

1. Create tool structure using the template above
2. Implement core functionality
3. Add comprehensive tests
4. Document thoroughly
5. Submit pull request

## Quality Standards

All tools must meet these standards:

- ✓ Functional parity with original tool
- ✓ >80% test coverage
- ✓ Complete documentation
- ✓ Pass all linting checks (`go vet`, `gofmt`)
- ✓ Performance benchmarks completed
- ✓ Edge cases handled and tested

## License

All tools in this directory are licensed under the Apache License 2.0, the same as the parent project. See [../LICENSE](../LICENSE) for details.

## Support

For questions or issues with specific tools, please open an issue on GitHub with the tool name in the title.
