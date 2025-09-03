# AI4Science Workshop Gemini-CLI Agents

This project contains three Gemini-CLI agents designed for an AI4Science workshop. Each agent is tailored for a specific scientific task.

## 1. Scientific Plotting Agent

This agent quickly plots scientific data from CSV or PKL files.

**Functionality:**
- Reads data from a CSV or PKL file.
- Analyzes the data to understand the columns.
- Generates a plot using Matplotlib.
- Saves the plot as an image file.

**Example Prompts:**
- "Plot the data in `results.csv`. The first column is the x-axis, and the second column is the y-axis."
- "I have a pandas DataFrame in `data.pkl`. Create a scatter plot of the 'temperature' vs. 'pressure' columns."

## 2. Equation Implementation Agent

This agent translates equations from scientific papers (in PDF format) into Python code.

**Functionality:**
- Reads a PDF of a scientific paper.
- Identifies equations and the surrounding text.
- Translates the equations into a Python function.
- Writes the Python code to a file.

**Example Prompts:**
- "Read the paper `paper.pdf` and implement the equation for the frobulator in Python."
- "I need a Python implementation of the main equation in `theory.pdf`."

## 3. Quantum Chemistry Agent

This agent connects to a server to run quantum chemistry calculations.

**Functionality:**
- Connects to a Molecular Computation Platform (MCP) server.
- Estimates errors and compute times for quantum chemistry calculations.
- Runs the calculations on the server.
- Returns the results to the user.

**Example Prompts:**
- "Connect to the MCP server and calculate the ground state energy of a water molecule."
- "Estimate the compute time for a DFT calculation on a benzene molecule."
