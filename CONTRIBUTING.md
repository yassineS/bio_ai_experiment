# Contributing to Bio AI Experiment

Thank you for your interest in contributing to this project! This guide will help you get started.

## Code of Conduct

This project adheres to a code of conduct that promotes a welcoming and inclusive environment. By participating, you are expected to uphold this code.

## How Can I Contribute?

### Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates. We have a structured bug report template that will guide you through providing all necessary information.

**To report a bug:**

1. Go to [Issues](https://github.com/yassineS/bio_ai_experiment/issues/new/choose)
2. Select "Bug Report" template
3. Fill in all required fields
4. Submit the issue

The bug report template will ask for:

- Component affected
- Clear description of the bug
- Expected vs. actual behavior
- Steps to reproduce
- Environment details
- Sample data (if applicable)

### Suggesting Enhancements

Enhancement suggestions are tracked as GitHub issues. We have a feature request template to help structure your suggestion.

**To suggest a feature:**

1. Go to [Issues](https://github.com/yassineS/bio_ai_experiment/issues/new/choose)
2. Select "Feature Request" template
3. Fill in all required fields
4. Submit the issue

The feature request template will ask for:

- Problem statement
- Proposed solution
- Use cases
- Priority and impact
- Implementation ideas (optional)

### Requesting Tool Analysis

If you'd like to suggest a bioinformatics tool for analysis and potential recoding:

**To request tool analysis:**

1. Go to [Issues](https://github.com/yassineS/bio_ai_experiment/issues/new/choose)
2. Select "Tool Analysis Request" template
3. Fill in tool details and justification
4. Submit the issue

### Starting Discussions

For questions, ideas, or general discussions that don't require tracking as issues:

**Use GitHub Discussions:**

- 💬 [Q&A](https://github.com/yassineS/bio_ai_experiment/discussions/categories/q-a) - Ask questions
- 💡 [Ideas](https://github.com/yassineS/bio_ai_experiment/discussions/categories/ideas) - Brainstorm features
- 🚀 [Show and Tell](https://github.com/yassineS/bio_ai_experiment/discussions/categories/show-and-tell) - Share your work
- 📣 [Announcements](https://github.com/yassineS/bio_ai_experiment/discussions/categories/announcements) - Project updates

See [DISCUSSIONS_SETUP.md](.github/DISCUSSIONS_SETUP.md) for more details.

### Pull Requests

We have a comprehensive pull request template to ensure all necessary information is provided.

**To submit a pull request:**

1. **Fork the repository** and create your branch from `main`
2. **Make your changes** following our coding standards
3. **Add tests** for new functionality
4. **Update documentation** as needed
5. **Ensure tests pass** (`go test ./...`)
6. **Run linters** (`go vet`, `gofmt`)
7. **Create pull request** using our template
8. **Fill in all sections** of the PR template

**The PR template includes:**

- Description of changes
- Type of change (bug fix, feature, etc.)
- Testing performed
- Code quality checklist
- Documentation updates
- Security considerations

See our [Pull Request Template](.github/PULL_REQUEST_TEMPLATE.md) for the complete checklist.

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

This repository is a **single Go module** (`github.com/yassineS/bio_ai_experiment`);
there is no per-tool `go.mod`. Run all Go commands from the repo root.

```bash
# Build everything
go build ./...

# Build a single tool's binary
go build ./tools/seqtk/cmd/seqtk

# Build the shared packages
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
├── go.mod               # single Go module; no per-tool go.mod
├── .github/
│   └── agents/          # Agent role descriptions
├── pkg/                 # shared libraries (htsgo, cliflag, ...)
├── tools/               # Recoded bioinformatics tools
│   └── <tool>/
│       ├── cmd/<tool>/main.go   # CLI entry point
│       ├── pkg/<tool>/          # tool logic + inline *_test.go
│       └── README.md            # per-tool usage + parity notes
├── test/                # validation + paper: nextflow, manuscript, figures, scripts
├── mcp-servers/         # MCP server implementations (descoped)
└── docs/                # Project documentation
```

Tests live **inline** next to the code as `*_test.go` under
`tools/<tool>/pkg/<tool>/` — there is no separate `tests/` or per-tool
`docs/` subtree (older docs that describe one are stale). For the current
tool-by-tool status see [`PROJECT_STATUS.md`](PROJECT_STATUS.md).

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

- **Questions?** Ask in [Discussions Q&A](https://github.com/yassineS/bio_ai_experiment/discussions/categories/q-a)
- **Ideas?** Share in [Discussions Ideas](https://github.com/yassineS/bio_ai_experiment/discussions/categories/ideas)
- **Found a bug?** [Open an issue](https://github.com/yassineS/bio_ai_experiment/issues/new/choose) using the Bug Report template
- **Documentation?** Check existing docs first in [docs/](docs/) and [test/](test/)

## Project Organization

### Issue Templates

We provide structured issue templates to help you provide all necessary information:

- **Bug Report** - For reporting bugs and unexpected behavior
- **Feature Request** - For suggesting new features or enhancements
- **Tool Analysis Request** - For requesting analysis of bioinformatics tools

See [.github/ISSUE_TEMPLATE/](.github/ISSUE_TEMPLATE/) for all templates.

### Project Boards

We use GitHub Projects to organize and track work. See [PROJECT_BOARDS.md](.github/PROJECT_BOARDS.md) for details on:

- Tool Development Board
- Feature/Enhancement Board
- Bug Tracking Board
- Research & Analysis Board

### GitHub Discussions

For community conversations, Q&A, and brainstorming, use GitHub Discussions. See [DISCUSSIONS_SETUP.md](.github/DISCUSSIONS_SETUP.md) for:

- Discussion categories
- How to participate
- Community guidelines

## Recognition

Contributors are recognized in:

- Git commit history
- Release notes (for significant contributions)
- Project documentation

## License

By contributing, you agree that your contributions will be licensed under the Apache License 2.0.
