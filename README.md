# easyMD
A simple tool for creating and running MD simulations, using Jupyter and OpenMM.
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Description

This Jupyter notebook will walk you through setting up an MD simulation using the PDB. It can be run as-is, or configured easily using the OpenMM python API.

This notebook is not meant to be a black box: where useful, I package up some extra code under the `utils/` directory, but I try to leave as much of the code as possible visible to the user for easy editing.

The simulations can be run locally, or queued to a cluster from within the notebook.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Installation

The notebooks and relevant code are packaged into this repository. To install, first clone this github repo using the following command:

`git clone https://github.com/degrado-lab/easyMD.git`

Then, enter the directory.

A conda environment is provided, in the form of a .yml file. To create it, use the following command:

`conda env create -f easyMD.yml`

Once created, enter the environment using:

`conda activate easyMD`

Make a copy of the `prep_and_run_template.ipynb` notebook for your analysis, and get started! For now, all notebooks must stay in the original directory, to reference local code.

### Tips for Wynton:
On Wynton, you must first load the conda and conda-stage modules:

`module load CBI miniconda3`

`module load CBI conda-stage`

Create the environment using `conda env create`, as above.

This can be painfully slow on Wynton, but in my experience, shouldn't take more than 20 minutes for this environment.

Then, set up local staging using the following command (this makes things run faster on the cluster!):

`conda activate easyMD`

`conda-stage --auto-stage=enable`

Finally, deactivate and activate once more before running the notebook, to make sure it's staged correctly:

`conda deactivate`

`conda activate easyMD`

## Usage
This code requires a GPU to run efficiently. You should install and run the notebook on a machine with a GPU, like a personal computer or the gpu development nodes on Wynton.

### Setting up simulations:
This notebook guides you through prepping and running a simulation - that means user input will be required throughout, to specify the input files and the specifics of preparing and running your system.

Not all cells will need to be run. Under each section, where a cell is preceded by "OPTION", please only run the cell that matches your needs (e.g. Are you loading a local PDB file? Or from the database?)

### Running simulations:
If on Wynton, please don't run your long production simulations on the development nodes, since those are shared. Instead, use the appropriate cell to queue your simulation, then check back later and do analysis. (This feature is WIP).

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

- [Nicholas Freitas](https://github.com/njf042)
- [Project Repository](https://github.com/degrado-lab/easyMD)

## Work in Progress:
- SUBMISSION ON WYNTON: Currently, the job submission script is for SLURM, and does not yet work on Wynton's SGE queue.
- SMALL MOLECULES:      Will eventually use OpenFF to allow for parameterization for small-molecules. This is an important and missing feature.

## References
- nglview help: https://projects.volkamerlab.org/teachopencadd/talktorials/T017_advanced_nglview_usage.html
- md traj analysis examples: https://github.com/mdtraj/mdtraj/tree/master/examples
- OpenMM cookbook: https://openmm.github.io/openmm-cookbook/latest/cookbook
- PDBFixer Guide: https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html
- Example of eventual OpenFF implementation: https://github.com/openforcefield/openff-toolkit/blob/stable/examples/toolkit_showcase/toolkit_showcase.ipynb

