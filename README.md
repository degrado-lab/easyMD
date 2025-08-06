# easyMD
A CLI tool for preparing and running molecular dynamics simulations, using OpenMM.
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Description

Never run an MD simulation before? With EasyMD, it's as simple as:
```easymd run my_protein.pdb ```

> [!WARNING]  
> EasyMD is a Work In Progress, and is currently in pre-release form. Feel free to use and offer suggestions!

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Installation
EasyMD requires packages installed via the Conda package manager. 

First, clone this repo:
```bash
git clone https://github.com/degrado-lab/easyMD.git
```

Inside the repo, create a new conda environment using the provided file:
```bash
conda env create -f environment.yml -y
```

Then, enter the environment:
```bash
conda activate easymd
```

Install EasyMD into the conda environment:
```
pip install .
```

> [!NOTE]  
> Soon, this installation will be simplified with an official PyPi release of EasyMD.

## Usage
Like all MD software, EasyMD requires a GPU to run efficiently. You should run EasyMD on a machine with a GPU, like a personal computer or the gpu development nodes on Wynton.

### Running a basic protein simulation
For an unprocessed input structure (PDB or CIF), you can prepare and run a simulation using the `easymd run` command:

```easymd run 1ohp.pdb --duration 100    # run for 100 ns```

What's going on when we run this command?

First, EasyMD prepares the protein structure for simulation. This includes handling non-standard residues or molecules, adding hydrogens, and solvating the protein. Next, it runs an energy mimimization, followed by an initial equilibration & relaxation simulation (by default, 0.5ns of NVT, then 0.5ns of NPT). These initial simulations allow the system to settle into a realistic pose and energy level for the subsequent production simulation.

Finally, EasyMD will run the production simulation. The processed output structure will be saved as a PDB file, and the position data for the simulation is saved as a separate DCD file. These can be viewed in software like PyMol by loading first the PDB, then the DCD file.

To change the output file, relaxation simulation duration, and system preparation variables, view the full list of parameters using:

```easymd run --help```

### Changing forcefields
OpenMM (the "engine" inside EasyMD) comes with many pre-packaged MD forcefields. To view them, run the command:

```easymd get-forcefields```

To run a simulation with a different forcefield, simply pass in the names of the forcefield XMLs to the ```--forcefield``` argument. You can also use forcefields downloaded to your computer, as long as they're in the OpenMM-compatible XML format.

```bash
easymd run 1ohp.cif --forcefield amber96.xml --forcefield tip3p.xml    # Use Amber96, with tip3p water
```

### Running a simulation with a ligand
Most protein forcefields are designed to recognize the typical protein amino acids, and some solvents. For any other molecule in your simulation (which I'll refer to as ligands), EasyMD will attempt to recognize it and generate forcefield parameters using the GAFF forcefield.

To do this, EasyMD needs a copy of each ligand as a separate file, with the correct bond orders and hydrogens. You can provide an SDF file of your ligand using the `-l [ligand.SDF]` flag (use an additional flag for each new ligand). The ligand in the SDF _does not_ have to be in the same pose as in your PDB - this is used as a template for bond orders and hydrogen positions.

```bash
easymd run 1ohp.cif -l 6NT.sdf    # manually adding ligand here
```

Alternatively, if your ligand is in the Protein DataBank's ligand repository, you can exclude the input SDF - EasyMD will use the ligand residue name to download an idealized SDF at runtime. If you'd like to utilize this feature, simply make sure your ligands have a residue name which matches their entry in the PDB (e.g. (6NT)[https://www.rcsb.org/ligand/6NT]).

```bash
easymd run 1ohp.cif    # ligand will be downloaded from PDB
```

> [!NOTE]  
> If you've parameterized your ligand yourself, you may provide that as a forcefield file using the ```--forcefield``` flag.
> Additionally, more general ligand forcefield generators (SMIRNOFF and ESPALOMA) will be added soon.

### Output Processing
After your simulation is run, the protein will likely be tumbling, with any ligands crossing over periodic boundaries frequently. This makes analysis difficult - so, EasyMD will automatically attempt to center and align your protein at the end of the simulation.
The centered trajectory will be saved with the suffix `_aligned.dcd`.

### Adding Hydrogens
EasyMD has a simple utility for adding hydrogens to a structure. Using custom forcefield files and provided ligand SDFs, it will use OpenMM to add hydrogens and energy-minimize their pose (without affecting heavy atoms). The command is:

```easymd reduce 1ohp.cif```

And it takes the same forcefield and ligand parameters as `easymd run`.

## Troubleshooting

Facing a CUDA error? 

```OpenMMException: There is no registered Platform called "CUDA"```

This means OpenMM was built for a different cuda version than your system. Run the following command to re-install openmm:

```conda install -c conda-forge openmm cuda-version=12```

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

- [Nicholas Freitas](https://github.com/njf042)
- [Project Repository](https://github.com/degrado-lab/easyMD)

## References
- md traj analysis examples: https://github.com/mdtraj/mdtraj/tree/master/examples
- OpenMM cookbook: https://openmm.github.io/openmm-cookbook/latest/cookbook
- PDBFixer Guide: https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html
- OpenFF implementation: https://github.com/openforcefield/openff-toolkit/blob/stable/examples/toolkit_showcase/toolkit_showcase.ipynb

