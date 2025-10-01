# AI4Science Workshop: Opencode/Gemini-CLI

This project provides three example agents for a workshop on using AI in science. Each agent is designed to be used with the `gemini-cli` or `opencode` tool.


# Installation
To use either [gemini-cli](https://github.com/google-gemini/gemini-cli) or [opencode](https://opencode.ai/), you will need to install npm, which we encourage you to install using Node Version Manager [(nvm)](https://github.com/nvm-sh/nvm?tab=readme-ov-file#installing-and-updating). The following can be used to install nvm and then npm:
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

Or opencode via
```bash
npm install -g opencode-ai
```

Other installation details can be found at the official gemini-cli page [here](https://github.com/google-gemini/gemini-cli) or opencode page [here](https://opencode.ai/).

# Usage
## 1. Scientific Plotting Example

Plotting figures is an essential task for scientists, so this example shows how
you can use agentic command line tools to interpret data and plot data to
reproduce literature results. The example specifically highlights reproducing
violin plots from a recent
[work](https://chemrxiv.org/engage/chemrxiv/article-details/67fe885f6e70d6fb2e033804)
comparing levels of symmetry adapted perturbation theory methods in quantum
chemistry. You don't have to understand the chemistry or physics to engage in
this example for it is supposed to be general enough to really stand in for any
data you might encounter.

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
