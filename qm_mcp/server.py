# server.py
from fastmcp import FastMCP

mcp = FastMCP("my-mcp-server")


@mcp.tool
def add(a: int, b: int) -> int:
    """Add two numbers"""
    return a + b


if __name__ == "__main__":
    mcp.run()
