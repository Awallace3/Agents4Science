# AI4Science Workshop: Opencode/Gemini-CLI

![header](extras/Header.png)

This repository provides information for both tutorials apart of the IDEaS
One-Day Tutorial on AI4Science. Some general links are provided directly below
for the google colab notebooks and the details for Tutorial 1 with three
separate applications are included below. Note, Tutorial 1 is designed to work
with `gemini-cli` or `opencode` tool (see installation details below).

## General Links

- Workshop [website](https://sites.gatech.edu/ai4science-tutorial/)
- Tutorial 1 Slides (after workshop)
- Tutorial 2 Slides (after workshop): 
    - notebook links (day of workshop)

# Tutorial 1
## Installation

Create a conda environment for python packages
```py
conda env create -f environment.yml
conda activate p4_qcml
```

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

### Optional
We recommend installing and using VSCode for these demos to facilitate the
integration between a terminal (can open with ctrl+\` or cmd+\`). When
launching gemini or opencode the first time, allow vscode to install the cli
extension for integrated support to see file changes easier from the models.

## Usage
## 1. Scientific Plotting and Data Extension Example

Plotting figures is an essential task for scientists, so this example shows how
you can use agentic command line tools to interpret data and plot data to
reproduce literature results. The example specifically highlights reproducing
violin plots from a recent
[work](https://chemrxiv.org/engage/chemrxiv/article-details/67fe885f6e70d6fb2e033804)
comparing levels of symmetry adapted perturbation theory methods in quantum
chemistry. You don't have to understand the chemistry or physics to engage in
this example for it is supposed to be general enough to really stand in for any
data you might encounter.

This agent interprets and plots from data in CSV or PKL files. 

### 2. Equation Implementation Example

This agent reads a scientific paper in PDF format and translates the equations
into Python code using an agentic test-driven development approach

### 3. Creating an MCP Server

Goes through creating and adding an MCP server to gemini-cli

