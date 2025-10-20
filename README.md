# Bio AI Experiment

Deploying agents to recode and modernise bioinformatics tools.

## Overview

This repository is an experimental project that leverages AI agents and automated tools to improve the bioinformatics software ecosystem. The primary focus is to identify, analyze, recode, test, and document popular bioinformatics tools, making them more robust, performant, and easier to use.

## Project Goals

The main objectives of this repository are to use AI agents to:

1. **Build a Comprehensive Tool Database**
   - Compile a list of the top 100 most used bioinformatics tools
   - Document each tool's name, link, citation, and download statistics
   - Create a structured database for easy reference

2. **Identify Areas for Improvement**
   - Analyze tools for poor code quality, performance issues, and documentation gaps
   - Identify edge cases and use cases that are not well-handled
   - Document findings systematically

3. **Recode Tools in GoLang**
   - Reimplement selected tools using Go for improved performance and maintainability
   - Maintain compatibility with original tool functionality
   - Leverage Go's concurrency and performance features

4. **Comprehensive Testing**
   - Provide test data for each tool
   - Write unit tests for main functionalities
   - Write unit tests for edge cases
   - Ensure robust test coverage

5. **Documentation**
   - Document all code comprehensively
   - Create user guides and API documentation
   - Document known issues and workarounds
   - Maintain up-to-date documentation throughout iterations

6. **Iterative Improvement**
   - Repeat the analysis-recode-test-document cycle until tools are robust
   - Continuously refine based on findings and feedback

7. **MCP Integration**
   - Write Model Context Protocol (MCP) servers for each tool
   - Simplify integration of bioinformatics tools with Large Language Models
   - Enable easier access to bioinformatics capabilities through AI interfaces

## Ultimate Goals

- **Improve Usability**: Make bioinformatics tools more accessible and easier to use
- **Enhance Documentation**: Provide clear, comprehensive documentation for all tools
- **Boost Performance**: Leverage modern programming practices and languages for better performance
- **Document AI Agent Utility**: Track and document the effectiveness (or lack thereof) of coding agents in this process

## Repository Structure

```
bio_ai_experiment/
├── .github/
│   └── agents/          # Agent configuration files
├── tools/               # Directory for recoded tools
│   └── [tool-name]/
│       ├── src/         # Go source code
│       ├── tests/       # Unit tests
│       ├── testdata/    # Test data
│       └── docs/        # Tool-specific documentation
├── analysis/            # Tool analysis and findings
├── mcp-servers/         # MCP server implementations
└── docs/                # General documentation

```

## Getting Started

### Prerequisites

- Go 1.21 or later
- Git
- Basic understanding of bioinformatics tools

### Installation

```bash
git clone https://github.com/yassineS/bio_ai_experiment.git
cd bio_ai_experiment
```

### Building Tools

Individual tools can be built from their respective directories:

```bash
cd tools/[tool-name]
go build ./...
```

### Running Tests

```bash
cd tools/[tool-name]
go test ./...
```

## Contributing

Contributions are welcome! This project is designed to be collaborative and benefits from diverse perspectives.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Contribution Guidelines

- Follow Go best practices and conventions
- Include tests for new functionality
- Update documentation as needed
- Ensure all tests pass before submitting PR

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- The bioinformatics community for creating the original tools
- AI/LLM technologies that make this experiment possible
- All contributors to this project

## Project Status

This is an active experimental project. Progress and findings will be documented regularly.

## Contact

For questions, suggestions, or collaboration opportunities, please open an issue in this repository.
