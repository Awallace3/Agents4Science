# Prompts
Leveraging our `./GEMINI.md` file, the agent already has context for our tasks;
hence, the prompt can be quite minimal. Think of AGENTS.md, CLAUDE.md, and
GEMINI.md files as re-usable context for different tasks. Certain CLI tools,
like opencode and claude-code, support creating agents in dot files that can
save this type of and expose only certain tools to the agent. gemini-cli does
not support this directly yet, so we are demonstrating it explicitly with a
file here. When using agents, the AGENTS.md file might be used to provide
general context of the project as determined by a `/init` call.

## Levels of SAPT II Violin Plot
```txt
Compare SAPT0 versus SAPT2+3(CCD)
```

```txt
The units of SAPT methods are in hartree and need to be converted to kcal/mol
```

```txt
Make the violins of the error
```

## Want to extend the dataset?
What if we wanted to take on of the dissociation curves and add on more
datapoints? Let's select out a specific curve, determine the displacement 
vector to add data to this curve, and extend the results for SAPT0 adz.
```txt

```
