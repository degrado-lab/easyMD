# easyMD
A simple tool for creating and running MD simulations, using Jupyter and OpenMM.
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Description

EasyMD is a software package that allows for preparing and running MD simulations. It is packaged as an Apptainer for easy installation.

A provided Jupyter notebook will walk you through setting up an MD simulation from a PDB file. It can be run as-is, or configured easily using the OpenMM 
python API.

The simulations can be run locally, or queued to a cluster from within the notebook.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Installation

First, download and extract the latest [release file](https://github.com/degrado-lab/easyMD/releases/latest), easyMD_vX.X.X.tar.gz, onto your computer or compute cluster.

Extract the directory with `tar xzf easyMD_vX.X.X.tar.gz` 

Once extracted, you will see three files:
- easyMD_setup.sh: this script configures your Wynton profile to run properly, and dowloads the virtual machine containing all the software.
- prep_and_run_template.ipynb: A template Jupyter notebook to prepare and begin an MD simulation.

On Wynton, run the setup script with the following commands:

`chmod u+x easyMD_setup.sh`

`./easyMD_setup.sh`

EasyMD runs as a virtual machine, with all software pre-installed. This is managed by the software Apptainer, which Wynton already has installed.
If you wish to run this on your personal computer, install Apptainer as described here for [Linux](https://apptainer.org/docs/admin/main/installation.html#install-ubuntu-packages) or [Mac](https://apptainer.org/docs/admin/main/installation.html#mac).

## Usage
This code requires a GPU to run efficiently. You should install and run the notebook on a machine with a GPU, like a personal computer or the gpu development nodes on Wynton.

Launching this virtual machine will load up a Jupyter server that you can use to prepare a simulation with the provided notebook. When you're finished, the notebook will let you launch queue the simulations on Wynton.

### Setting up simulations with the notebook:
To prepare a simulation using the provided Jupyter notebook, first make a copy of it with a memorable name.

Next, we will start the Jupyter server with the following command:
`./easymd_vX.X.X.sif notebook`

This will launch the server. In the output you will see the address of the server, which looks like `http://localhost:54321/tree?token=...`. 
You may connect to this server using your web browser or in VS Code.

This notebook guides you through prepping and running a simulation - that means user input will be required throughout, to specify the input files and the specifics of preparing and running your system.

This notebook is not meant to be a black box: where useful, I package up some extra code under the `utils/` directory, but I try to leave as much of the code as possible visible to the user for easy editing.

Not all cells will need to be run. Under each section, where a cell is preceded by "OPTION", please only run the cell that matches your needs (e.g. Are you loading a local PDB file? Or from the database?)

Note: If you have an error about the software not recognizing `nvidia-smi`, this means the container is not seeing your devices GPU. Please restart the container, manually telling it to enable the GPU using the command: `apptainer run --nv easymd_vX.X.X.sif notebook`

### Preparing Ligands:
You can simulate ligands that are in your PDB by specifying the residue name (for example, 'BNZ') when preparing the structure. Currently there is a bug that only allows for one ligand type.

To use a ligand, your PDB file must have the correct hydrogens specified in the ligand structure (this is because, the PDB file does not store bond orders, so we cannot infer the correct number of hydrogens from only the PDB). I recommend doing this in PyMOL.

Additionally, to make sure the bond orders are correct a template SDF file for your residue must be used. By default, the notebook will try to download a matching template from the PDB. If this does not work, you should manually load a matching SDF.

### Running simulations:
At the end of the notebook is the option to queue your jobs to Wynton. Please specify a time duration that will be sufficient to finish your simulation. If your job ends early the data should be fine, just truncated. You can queue another job using the same cell to continue it, and it will pick up from the last checkpoint and add to the existing DCD file.

If on Wynton, please **don't** run your long production simulations on the development nodes with run_sim_local(), since those are shared. Instead, use the next cell to queue your simulation, then check back later and do analysis.

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

- [Nicholas Freitas](https://github.com/njf042)
- [Project Repository](https://github.com/degrado-lab/easyMD)

## Work in Progress:
- FIX LIGAND PARAMETERIZATION: This is working on-and-off for now.
- REPEAT FRAMES DURING SERIAL SIMULATIONS: When queueing multiple wynton jobs serially for the same simulation, if one is shut down before writing a checkpoint, the next will re-run those frames and append to the DCD. This may cause duplicate frames.

## References
- nglview help: https://projects.volkamerlab.org/teachopencadd/talktorials/T017_advanced_nglview_usage.html
- md traj analysis examples: https://github.com/mdtraj/mdtraj/tree/master/examples
- OpenMM cookbook: https://openmm.github.io/openmm-cookbook/latest/cookbook
- PDBFixer Guide: https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html
- Example of eventual OpenFF implementation: https://github.com/openforcefield/openff-toolkit/blob/stable/examples/toolkit_showcase/toolkit_showcase.ipynb

