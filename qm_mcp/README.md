# Building MCP Servers
This section is a little different than the other examples because you won't
purely have the agent write the code for you. The goal here is to show you how
you can create MCP servers and connect them to `gemini-cli` for usage. You
could take this as an opportunity to practice your other skills or just manually
create the server as you like.

## FastMCP
While there are different ways to write an MCP server, the simplest and 
fastest way if you are familiar with Python is to use FastMCP. This
package only requires 4 things to work:
1. Create a FastMCP object with a name:
```py
# server.py
from fastmcp import FastMCP

mcp = FastMCP("my-mcp-server")
```
2. Write a function with type hints:
```py
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b
```
3. Wrap any functions you want to be a tool/resource with the appropriate decorator:
```py
@mcp.tool
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b
```
4. Have the server run when executed
```py
if __name__ == "__main__":
    mcp.run()
```
All of these steps together should give you the following minimal example:
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

# gemini MCP setup
Recent additions to gemini have made it even easier to add a new server by
using the following syntax `gemini mcp add <name> <commandOrUrl>`. You can get
more information if you run `geminii mcp -h` in the terminal, but for our
specific example you can just execute the following:
```sh
gemini mcp add my-mcp-server fastmcp run server.py
```
Re-launching `gemini` in the terminal, you should now see lines like
```txt
Using:
  - 2 GEMINI.md files
  - 1 MCP server (ctrl+t to view)
```
and pressing ctrl+t will allow us to inspect our newly available tool for
adding. Importantly, this will tell us what is provided as context to the
model that includes the name of the server, the name of the function, and a
description on the tool. To test that it works successfully, just ask gemini
to add two numbers like 1+3. It will prompt you to allow calling our tool,
`add`, and after allowing, provide back a result of 4, meaning our server
works! More complicated prompts/tools might not just return the tool result
but actually use that action as more context to continue reasoning/building a
response with the new information. Note the LLM is blocked until the tool call
returns, so if the function is very long, your model might timeout if you have
not configured it properly.


# Quantum Chemistry MCP Example
TODO
