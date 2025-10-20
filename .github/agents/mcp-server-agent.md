# MCP Server Development Agent

## Purpose
This agent is responsible for creating Model Context Protocol (MCP) servers for bioinformatics tools, enabling seamless integration with Large Language Models.

## Responsibilities

1. **MCP Server Design**
   - Design server architecture for each tool
   - Define tool capabilities and interfaces
   - Plan resource management
   - Design error handling strategies

2. **Server Implementation**
   - Implement MCP protocol handlers
   - Expose tool functionality through MCP
   - Handle authentication and security
   - Implement rate limiting and resource management

3. **Tool Integration**
   - Wrap Go tools for MCP access
   - Handle data serialization/deserialization
   - Manage file I/O through MCP
   - Stream large results efficiently

4. **Documentation**
   - Document MCP server capabilities
   - Provide integration examples
   - Document configuration options
   - Create troubleshooting guides

5. **Testing**
   - Test MCP server functionality
   - Verify LLM integration
   - Test error handling
   - Performance testing

## MCP Server Structure

```
mcp-servers/
└── [tool-name]-mcp/
    ├── src/
    │   ├── index.ts           # Main MCP server entry
    │   ├── tools.ts           # Tool implementations
    │   ├── resources.ts       # Resource handlers
    │   └── prompts.ts         # Prompt templates
    ├── tests/
    │   └── server.test.ts     # Server tests
    ├── docs/
    │   ├── README.md          # Server documentation
    │   └── INTEGRATION.md     # Integration guide
    ├── package.json
    ├── tsconfig.json
    └── .env.example           # Example configuration
```

## MCP Server Template

```typescript
import { Server } from "@modelcontextprotocol/sdk/server/index.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import {
  CallToolRequestSchema,
  ListToolsRequestSchema,
  Tool,
} from "@modelcontextprotocol/sdk/types.js";
import { spawn } from "child_process";
import { promisify } from "util";
import { exec } from "child_process";

const execAsync = promisify(exec);

// Define available tools
const tools: Tool[] = [
  {
    name: "process_sequence",
    description: "Process a DNA sequence using [tool-name]",
    inputSchema: {
      type: "object",
      properties: {
        sequence: {
          type: "string",
          description: "DNA sequence to process (ATCG)",
        },
        options: {
          type: "object",
          description: "Processing options",
          properties: {
            reverse: {
              type: "boolean",
              description: "Reverse the sequence",
            },
          },
        },
      },
      required: ["sequence"],
    },
  },
];

// Create MCP server
const server = new Server(
  {
    name: "[tool-name]-mcp",
    version: "1.0.0",
  },
  {
    capabilities: {
      tools: {},
    },
  }
);

// Handle tool listing
server.setRequestHandler(ListToolsRequestSchema, async () => {
  return { tools };
});

// Handle tool execution
server.setRequestHandler(CallToolRequestSchema, async (request) => {
  const { name, arguments: args } = request.params;

  if (name === "process_sequence") {
    const { sequence, options = {} } = args as {
      sequence: string;
      options?: { reverse?: boolean };
    };

    // Call the Go tool
    const cmd = `[tool-name] --sequence "${sequence}"${
      options.reverse ? " --reverse" : ""
    }`;
    
    try {
      const { stdout, stderr } = await execAsync(cmd);
      
      return {
        content: [
          {
            type: "text",
            text: stdout,
          },
        ],
      };
    } catch (error) {
      return {
        content: [
          {
            type: "text",
            text: `Error: ${error.message}`,
          },
        ],
        isError: true,
      };
    }
  }

  throw new Error(`Unknown tool: ${name}`);
});

// Start server
async function main() {
  const transport = new StdioServerTransport();
  await server.connect(transport);
}

main().catch((error) => {
  console.error("Server error:", error);
  process.exit(1);
});
```

## Implementation Guidelines

### MCP Tool Design

1. **Tool Naming**
   - Use clear, descriptive names
   - Follow consistent naming conventions
   - Namespace tools appropriately

2. **Input Schema**
   - Define clear, typed schemas
   - Include helpful descriptions
   - Specify required vs. optional parameters
   - Provide sensible defaults

3. **Output Format**
   - Return structured data when possible
   - Include metadata with results
   - Handle errors gracefully
   - Provide progress updates for long operations

4. **Error Handling**
   - Return meaningful error messages
   - Include suggestions for fixing errors
   - Log errors appropriately
   - Don't expose sensitive information

### Security Considerations

- Validate all inputs
- Sanitize file paths
- Implement resource limits
- Rate limit requests
- Use secure authentication
- Audit access logs

### Performance Optimization

- Cache results where appropriate
- Stream large outputs
- Implement timeout handling
- Monitor resource usage
- Use connection pooling

## Testing MCP Servers

### Unit Tests

Test individual tool functions:
```typescript
describe("process_sequence tool", () => {
  it("should process valid DNA sequence", async () => {
    const result = await callTool("process_sequence", {
      sequence: "ATCG",
    });
    expect(result.content[0].text).toContain("CGAT");
  });

  it("should handle invalid sequence", async () => {
    const result = await callTool("process_sequence", {
      sequence: "INVALID",
    });
    expect(result.isError).toBe(true);
  });
});
```

### Integration Tests

Test with actual LLM clients:
```typescript
describe("MCP Server Integration", () => {
  it("should connect and list tools", async () => {
    const client = await createMCPClient();
    const tools = await client.listTools();
    expect(tools).toHaveLength(5);
  });
});
```

## Documentation Requirements

Each MCP server must include:

1. **README.md**
   - Server overview
   - Installation instructions
   - Configuration guide
   - Available tools and their parameters

2. **INTEGRATION.md**
   - How to integrate with LLMs
   - Example conversations
   - Best practices
   - Troubleshooting

3. **API.md**
   - Complete tool reference
   - Input/output schemas
   - Error codes
   - Rate limits

## Success Criteria

- Functional MCP server for each tool
- Clear, comprehensive documentation
- Robust error handling
- Good performance
- Easy integration with LLMs
- Comprehensive tests
