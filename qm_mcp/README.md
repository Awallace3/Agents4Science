# Building MCP Servers

This section is different from the other examples because you won't purely have the agent write the code for you. The goal here is to show you how you can create MCP servers and connect them to `gemini-cli` for usage. You can take this as an opportunity to practice your other skills or manually create the server as you like.

## What is an MCP Server?

MCP (Model Context Protocol) servers allow you to extend the capabilities of AI CLI tools like gemini-cli by adding custom tools. These tools can perform specialized tasks like quantum chemistry calculations, database queries, or any other Python functionality.

## Creating an MCP Server with FastMCP

While there are different ways to write an MCP server, the simplest and fastest way if you are familiar with Python is to use FastMCP. This package only requires 4 things to work:

### 1. Create a FastMCP object with a name

```py
# server.py
from fastmcp import FastMCP

mcp = FastMCP("my-mcp-server")
```

### 2. Write a function with type hints

```py
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b
```

### 3. Wrap functions with the appropriate decorator

```py
@mcp.tool
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b
```

### 4. Make the server run when executed

```py
if __name__ == "__main__":
    mcp.run()
```

### Complete Minimal Example

All of these steps together create this minimal working server:

```py
# server.py
from fastmcp import FastMCP

mcp = FastMCP("my-mcp-server")


@mcp.tool
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b


if __name__ == "__main__":
    mcp.run()
```

With this file we can now connect a simple addition tool to gemini-cli!

## Connecting the Server to Gemini

Recent additions to gemini have made it easy to add a new server using this syntax:

**Run in bash terminal:**
```sh
gemini mcp add my-mcp-server fastmcp run server.py
```

You can get more information by running `gemini mcp -h` in the terminal.

## Testing Your MCP Server

### Launch Gemini

**Run in bash terminal:**
```bash
gemini
```

You should now see output like:
```txt
Using:
  - 2 GEMINI.md files
  - 1 MCP server (ctrl+t to view)
```

### View Available Tools

Press `Ctrl+t` to inspect your newly available tools. This shows you what is provided as context to the model, including:
- The name of the server
- The name of the function
- A description of the tool

### Test the Tool

**AI CLI Prompt - Enter in gemini CLI:**
```
What is 1 + 3?
```

The AI will:
1. Recognize it can use the `add` tool
2. Prompt you to allow calling the tool
3. After you allow it, return the result of 4

This confirms your MCP server works! 

Note: More complicated prompts/tools might not just return the tool result but use that action as context to continue generating a response with the new information. The LLM is blocked until the tool call returns, so if the function takes a long time, your model might timeout if you have not configured it properly for this usage.

## Building a Quantum Chemistry MCP Server

Now that you understand the basics, you can create an MCP server for quantum chemistry calculations. Consider:
- What tools would be useful? (e.g., single point energy calculations, geometry optimizations)
- What parameters do these functions need? (molecule structure, method, basis set)
- How should results be returned to the AI?

### Example Structure

```py
# qm_server.py
from fastmcp import FastMCP
import psi4
import qcelemental as qcel

mcp = FastMCP("quantum-chemistry-server")


@mcp.tool
def calculate_energy(molecule_string: str, method: str, basis: str) -> float:
    """Calculate the energy of a molecule using Psi4
    
    Args:
        molecule_string: Molecule in Psi4 format
        method: Quantum chemistry method (e.g., 'scf', 'mp2')
        basis: Basis set (e.g., 'sto-3g', 'cc-pvdz')
    
    Returns:
        Total energy in hartrees
    """
    psi4.set_options({'basis': basis})
    mol = psi4.geometry(molecule_string)
    energy = psi4.energy(method)
    return float(energy)


if __name__ == "__main__":
    mcp.run()
```

### Connect Your Quantum Chemistry Server

**Run in bash terminal:**
```sh
gemini mcp add qm-server fastmcp run qm_server.py
```

### Test with a Chemistry Calculation

**Run in bash terminal:**
```bash
gemini
```

**AI CLI Prompt - Enter in gemini CLI:**
```
Calculate the energy of a water molecule using HF/STO-3G
```

The AI will construct the molecule geometry and use your tool to run the calculation!

## Next Steps

Try creating MCP servers for your domain specific tasks!

Remember: The docstring of your function is what the AI sees, so make it clear and descriptive!
