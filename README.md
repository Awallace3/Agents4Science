# AI4Science Workshop Gemini-CLI Agents

This project provides three example agents for a workshop on using AI in science. Each agent is designed to be used with the `gemini-cli` tool.

## How to Use

To use an agent, you will invoke the `gemini-cli` with the `-C` flag to specify the agent's directory, followed by your prompt.

```bash
gemini
```

Before running the plotting agent, make sure to install the required dependencies:

```bash
pip install -r requirements.txt
```

## Agents

### 1. Scientific Plotting Agent

This agent creates plots from data in CSV or PKL files.

**Command:**
```bash
gemini 
"scientific_plotting_agent: Plot the data in 'data.csv' with the first column as the x-axis and the second as the y-axis."
```

### 2. Equation Implementation Agent

This agent reads a scientific paper in PDF format and translates the equations into Python code.

**Command:**
```bash
gemini
```

### 3. Quantum Chemistry Agent

This agent simulates connecting to a server to run quantum chemistry calculations.

**Command:**
```bash
gemini
```
