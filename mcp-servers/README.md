# MCP Servers for Bioinformatics Tools

This directory contains Model Context Protocol (MCP) server implementations for bioinformatics tools. These servers enable seamless integration of bioinformatics capabilities with Large Language Models.

## What is MCP?

The Model Context Protocol (MCP) is a standard protocol for connecting AI assistants to external tools and data sources. MCP servers expose tool functionality in a way that LLMs can easily discover and use.

## Available MCP Servers

Currently, this directory is empty. MCP servers will be added as tools are implemented and tested.

### Planned Servers

MCP servers will be created for each recoded bioinformatics tool, enabling LLM-powered interfaces to:

1. Sequence analysis
2. Quality control
3. Format conversion
4. Statistical analysis
5. Data visualization
6. Genome assembly
7. Variant detection
8. Annotation
9. Read mapping
10. And more...

## MCP Server Structure

Each MCP server follows this standard structure:

```
[tool-name]-mcp/
├── src/
│   ├── index.ts           # Main MCP server entry
│   ├── tools.ts           # Tool implementations
│   ├── resources.ts       # Resource handlers
│   └── prompts.ts         # Prompt templates (optional)
├── tests/
│   ├── tools.test.ts      # Tool tests
│   └── integration.test.ts # Integration tests
├── docs/
│   ├── README.md          # Server documentation
│   └── INTEGRATION.md     # Integration guide
├── package.json
├── tsconfig.json
└── .env.example           # Example configuration
```

## Using MCP Servers

### Prerequisites

- Node.js 18+ and npm
- Go binaries for the tools (built from ../tools/)

### Setting Up a Server

```bash
cd [tool-name]-mcp

# Install dependencies
npm install

# Build the server
npm run build

# Run tests
npm test
```

### Running a Server

```bash
# Start the server
npm start

# Or run in development mode
npm run dev
```

### Configuration

Each server can be configured via environment variables. Copy `.env.example` to `.env` and adjust settings:

```bash
cp .env.example .env
# Edit .env with your settings
```

## Integrating with LLMs

### Claude Desktop Integration

Add the server to your Claude Desktop configuration (`~/Library/Application Support/Claude/claude_desktop_config.json` on macOS):

```json
{
  "mcpServers": {
    "[tool-name]": {
      "command": "node",
      "args": [
        "/path/to/bio_ai_experiment/mcp-servers/[tool-name]-mcp/dist/index.js"
      ],
      "env": {
        "TOOL_PATH": "/path/to/bio_ai_experiment/tools/[tool-name]/[tool-name]"
      }
    }
  }
}
```

### Other MCP Clients

MCP servers can be used with any MCP-compatible client. See each server's `docs/INTEGRATION.md` for specific integration guides.

## Available Tools

Each MCP server exposes tools that can be called by LLMs. Tools are described with schemas that specify:

- Tool name and description
- Input parameters and types
- Output format
- Examples

### Example Tool

```typescript
{
  name: "analyze_sequence",
  description: "Analyze a DNA sequence for quality metrics",
  inputSchema: {
    type: "object",
    properties: {
      sequence: {
        type: "string",
        description: "DNA sequence (ATCG)",
      },
      options: {
        type: "object",
        properties: {
          quality_threshold: {
            type: "number",
            description: "Minimum quality score (0-100)",
          },
        },
      },
    },
    required: ["sequence"],
  },
}
```

## Development

### Creating a New MCP Server

1. Use the template structure above
2. Implement tool handlers in `src/tools.ts`
3. Add tests in `tests/`
4. Document in `docs/README.md`
5. Test with an MCP client

### Testing

```bash
# Run unit tests
npm test

# Run tests in watch mode
npm run test:watch

# Run integration tests
npm run test:integration
```

### Debugging

Enable debug logging:

```bash
DEBUG=mcp:* npm start
```

## Best Practices

### Error Handling

- Return clear error messages
- Include suggestions for fixing errors
- Don't expose sensitive information
- Log errors appropriately

### Performance

- Stream large outputs when possible
- Implement timeouts for long operations
- Cache results where appropriate
- Monitor resource usage

### Security

- Validate all inputs
- Sanitize file paths
- Implement rate limiting
- Use secure authentication
- Audit access logs

## Documentation

Each MCP server includes:

1. **README.md**: Server overview, installation, configuration
2. **INTEGRATION.md**: Integration guides for different clients
3. **API documentation**: Complete tool and resource reference

## Testing with MCP Inspector

Use the MCP Inspector for interactive testing:

```bash
npx @modelcontextprotocol/inspector node dist/index.js
```

This opens a web interface for testing tools and inspecting server capabilities.

## Common Issues

### Server Won't Start

- Check that Node.js version is 18+
- Ensure all dependencies are installed (`npm install`)
- Verify the tool binary path is correct
- Check environment variables

### Tools Not Appearing

- Verify the server started successfully
- Check the MCP client configuration
- Review server logs for errors
- Ensure tools are properly exported

### Performance Issues

- Check tool binary performance
- Monitor resource usage
- Review timeout settings
- Consider implementing caching

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) for guidelines on contributing MCP server implementations.

### Adding a New Server

1. Create server structure from template
2. Implement tool handlers
3. Add comprehensive tests
4. Document thoroughly
5. Test with multiple MCP clients
6. Submit pull request

## Quality Standards

All MCP servers must meet:

- ✓ Complete tool documentation
- ✓ Comprehensive tests (unit + integration)
- ✓ Error handling for all edge cases
- ✓ Performance optimization
- ✓ Security best practices
- ✓ Integration guide

## Resources

- [MCP Documentation](https://modelcontextprotocol.io/)
- [MCP TypeScript SDK](https://github.com/modelcontextprotocol/typescript-sdk)
- [MCP Specification](https://spec.modelcontextprotocol.io/)

## License

All MCP servers in this directory are licensed under the Apache License 2.0. See [../LICENSE](../LICENSE) for details.

## Support

For questions or issues with MCP servers, please open an issue on GitHub with the server name in the title.
