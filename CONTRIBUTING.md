# Contributing to Bio AI Experiment

Thank you for your interest in contributing to this project! This guide will help you get started.

## Code of Conduct

This project adheres to a code of conduct that promotes a welcoming and inclusive environment. By participating, you are expected to uphold this code.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates. When creating a bug report, include:

- **Clear title and description**
- **Steps to reproduce** the issue
- **Expected behavior** vs. actual behavior
- **Environment details** (OS, Go version, etc.)
- **Code samples** or test cases if applicable

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. When creating an enhancement suggestion, include:

- **Clear title and description**
- **Rationale** for the enhancement
- **Possible implementation** approach
- **Examples** of how it would be used

### Pull Requests

1. **Fork the repository** and create your branch from `main`
2. **Make your changes** following our coding standards
3. **Add tests** for new functionality
4. **Update documentation** as needed
5. **Ensure tests pass** (`go test ./...`)
6. **Run linters** (`go vet`, `gofmt`)
7. **Submit the pull request**

## Development Setup

### Prerequisites

- Go 1.21 or later
- Git
- golangci-lint (recommended)

### Setting Up Your Development Environment

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/bio_ai_experiment.git
cd bio_ai_experiment

# Add upstream remote
git remote add upstream https://github.com/yassineS/bio_ai_experiment.git

# Create a branch for your work
git checkout -b feature/your-feature-name
```

### Building

```bash
# Build all tools
cd tools/[tool-name]
go build ./...

# Build specific package
go build ./pkg/...
```

### Testing

```bash
# Run all tests
go test ./...

# Run tests with coverage
go test -cover ./...

# Run tests verbosely
go test -v ./...

# Run benchmarks
go test -bench=. ./...
```

### Code Quality

```bash
# Format code
gofmt -w .

# Run go vet
go vet ./...

# Run staticcheck (if installed)
staticcheck ./...

# Run golangci-lint (if installed)
golangci-lint run
```

## Coding Standards

### Go Code Style

- Follow [Effective Go](https://golang.org/doc/effective_go.html)
- Use `gofmt` for formatting
- Write idiomatic Go code
- Keep functions small and focused
- Use meaningful variable names

### Documentation

- Document all exported functions, types, and packages
- Use complete sentences in comments
- Provide examples for complex functionality
- Keep documentation up-to-date with code changes

Example:
```go
// ProcessSequence processes a DNA sequence and returns its reverse complement.
//
// The function validates the input sequence and returns an error if it contains
// invalid characters. Valid characters are A, T, C, G (case-insensitive).
//
// Example:
//   result, err := ProcessSequence("ATCG")
//   if err != nil {
//       log.Fatal(err)
//   }
//   fmt.Println(result) // Output: CGAT
func ProcessSequence(sequence string) (string, error) {
    // Implementation
}
```

### Testing

- Write tests for all new functionality
- Include edge cases in tests
- Use table-driven tests for multiple scenarios
- Aim for >80% code coverage
- Write benchmarks for performance-critical code

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb in the imperative mood
- Keep the first line under 50 characters
- Add detailed description if needed

Good examples:
```
Add reverse complement function for DNA sequences
Fix edge case in FASTA parser for empty files
Update documentation for ProcessSequence function
```

## Project Structure

Understanding the project structure will help you navigate the codebase:

```
bio_ai_experiment/
├── .github/
│   └── agents/          # Agent configuration files
├── tools/               # Recoded bioinformatics tools
│   └── [tool-name]/
│       ├── cmd/         # Command-line interface
│       ├── pkg/         # Library code
│       ├── tests/       # Tests
│       └── docs/        # Tool-specific docs
├── analysis/            # Tool analysis reports
├── mcp-servers/         # MCP server implementations
└── docs/                # Project documentation
```

## Working with Agents

This project uses AI agents to assist with various tasks. When contributing:

1. **Respect agent guidelines** in `.github/agents/`
2. **Follow established patterns** from agent work
3. **Document agent interactions** when relevant
4. **Provide feedback** on agent effectiveness

## Review Process

1. **Automated checks** run on all PRs (tests, linting)
2. **Code review** by maintainers
3. **Address feedback** and update PR
4. **Merge** once approved and checks pass

## Getting Help

- **Issues**: Ask questions by opening an issue
- **Discussions**: Use GitHub Discussions for broader topics
- **Documentation**: Check existing docs first

## Recognition

Contributors are recognized in:
- Git commit history
- Release notes (for significant contributions)
- Project documentation

## License

By contributing, you agree that your contributions will be licensed under the Apache License 2.0.
