# AI4Science Workshop: Gemini-CLI Agents

This project provides three example agents for a workshop on using AI in science. Each agent is designed to be used with the `gemini-cli` tool.


To use an agent, you will need to install npm, which we encourage you to install using Node Version Manager [(nvm)](https://github.com/nvm-sh/nvm?tab=readme-ov-file#installing-and-updating). The following can be used to install nvm and npm:
```bash
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
. ~/.bashrc
nvm install --lts
nvm use --lts
```

Then we can install gemini-cli through

```bash
npm install -g @google/gemini-cli
```

Alternatively on MAC, you can install gemini-cli through brew
```bash
brew install gemini-cli
```

More installation details can be found at the official gemini-cli page [here](https://github.com/google-gemini/gemini-cli)

## How to Use

Launch gemini and start prompting
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
