# Intro

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

# Workflow

First, read through `./GEMINI.md` to get a general understanding of the context
we are providing to the model. Key aspects are defining the role, workflow,
output format, and software stack (Dependencies). Then launch `gemini` in a terminal
in this current directory with this `README.md` open as well. Copy and paste
the prompts in the txt blocks below to work through this example the first time
and then try to modify prompts and/or use different models with opencode.


## Levels of SAPT II Violin Plot

Leveraging our `./GEMINI.md` file, the agent already has context for our tasks;
hence, the prompt can be quite minimal. Think of AGENTS.md, CLAUDE.md, and
GEMINI.md files as re-usable context for different tasks. Certain CLI tools,
like opencode and claude-code, support creating agents in dot files that can
save this type of and expose only certain tools to the agent. gemini-cli does
not support this directly yet, so we are demonstrating it explicitly with a
file here. When using agents, the AGENTS.md file might be used to provide
general context of the project as determined by a slash command, `/init`, in
either gemini or opencode. Because this is an empty project, we will skip the
`/init` call. Go ahead and start with the following prompt:
```txt
Compare SAPT0 versus SAPT2+3(CCD)
```

With this minimal prompt and context file, you should have created a violin
plot; however, you might notice the numbers do not agree with the paper...

```txt
The units of SAPT methods are in hartree and need to be converted to kcal/mol
```

Now the violin plots should be comparable to Figure 1 of the Levels of SAPT II
paper. Sometimes you might need to explicitly tell gemini to create violin plots
and can be done with the prompt below (skip if it worked already):
```txt
Make the violins of the error
```

## Want to extend the dataset?

What if we wanted to take on of the dissociation curves and add on more
datapoints? Let's select out a specific curve, determine the displacement
vector to add data to this curve, and extend the results for SAPT0 adz. We will
provide a code snippet to select a specific curve; however, you could iterate
through processing the data without this call. We would still recommend
converting data to a qcelemental.models.Molecule object to ensure a rudimentary
validity check on the system before the quantum calculations.

```txt
# Task
I want to select out a specific dissociation curve from the dataset through the python code snippet below in extended_pes.py. Using the df_subset, create a displacement_geometry function that 1) identifies which fragment is being adjusted between iloc[0] and iloc[1]. 2) computes a displacement vector for the shifted atoms. 3) use the displacement vector to create an arbitrary dimer potential energy surface with the minimum contact distance of starting_distance and increments of increment up to the minimum contact distance of ending_distance. Compute the minimum contact distance for each geometry and keep these values for plotting. The distance for geometry will be in bohr, so convert angstrom distances to bohr before adjusting the geometries. Extend the df_subset dataframe with these systems, setting all other columns to have nan values. Then use psi4 to compute SAPT0/aug-cc-pVDZ energies for each system using qcel_molecule.to_string('psi4') over these new rows with nans in the SAPT0 adz columns. Save the updated df_subset to 'extended_water_dimer.pkl'. 

# Final Outputs
Plot the SAPT0 PES over these points with total (black), elst (red), exch (green), ind (blue) and disp (orange). Have the original datapoints plotted with a X and the new datapoints with a O with all energies converted to kcalmol. Save the plot to water_pes.png.

Create a gif of the water PES to visualize the waters being pulled apart. Display the minimum contact distance while pulling these apart.

# Code Snippet
'''
import pandas as pd
import qcelemental as qcel
from qm_tools_aw import tools


def isolate_curve(df: pd.DataFrame, system_id_contains: str = "01_Water-Water"):
    """
    Isolate a subset of the DataFrame for a specific system_id.
    """
    df_subset = df[df["system_id"].str.contains(system_id_contains)].copy()
    print(df_subset)

    def qcel_mols(row) -> qcel.models.Molecule:
        """
        Convert the row to a qcel molecule
        """
        atomic_numbers = [
            row["atomic_numbers"][row["monAs"]],
            row["atomic_numbers"][row["monBs"]],
        ]
        coords = [row["coordinates"][row["monAs"]], row["coordinates"][row["monBs"]]]
        cm = [
            [row["monA_charge"], row["monA_multiplicity"]],
            [row["monB_charge"], row["monB_multiplicity"]],
        ]
        return tools.convert_pos_carts_to_mol(atomic_numbers, coords, cm)

    df_subset["qcel_molecule"] = df_subset.apply(qcel_mols, axis=1)
    return df_subset
'''
```

This should take a water-water dimer curve from the dataset, compute SAPT0
interaction energies energies, and plot the results. To restore state to the
beginning, run `git restore .` to clear all outputs to begin again with different
prompting and/or models! 

# Opencode
Opencode works very similar to gemini; however, you will want to authenticate with
GitHub Copilot before using it by running
```sh
opencode auth login
```
Then select GitHub Copilot
```txt
┌  Add credential
│
◆  Select provider

│  Search:
│  ○ opencode zen
│  ○ Anthropic
│  ● GitHub Copilot
│  ○ OpenAI
│  ○ Google
│  ○ OpenRouter
│  ○ Vercel AI Gateway
│  ...
│  ↑/↓ to select • Enter: confirm • Type: to search
└
```

Next launch `opencode` that should open a prompt very similar to gemini-cli.
Select a model with a slash command (`/models`) and pressing tab or space to
open a selection for models. I recommend selecting the `Claude Sonnet 4.5
(Preview) GitHub Copilot` model to use one of the best models as of 2025-10-03.
If you get an API error, you likely need to enable the model you selected from
copilot under your GitHub Copilot Settings
[here](https://github.com/settings/copilot/features).
