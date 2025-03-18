import logging
#from openmm import app, unit, Platform, openmm
import openmm
from pathlib import Path
from tempfile import NamedTemporaryFile

from easyMD.prep import structure as structure_prep
from easyMD.prep import ligand as ligand_prep

from easyMD.prep.data import common_residues

logger = logging.getLogger(__name__)
# This will suppress INFO warnings from openmm forcefields. This may be bad!
logging.getLogger("openmmforcefields.generators.template_generators").setLevel(logging.WARNING)

def create_modeller(system_pdb_path):
    """
    Builds an OpenMM `Modeller` from a PDB topology/positions file.
    """
    # Load the PDB file
    pdb = openmm.app.PDBFile(str(system_pdb_path))
    pdb_topology = pdb.topology
    pdb_positions = pdb.positions

    logger.info("Creating OpenMM Modeller from PDB file %s.", system_pdb_path)
    modeller = openmm.app.Modeller(pdb_topology, pdb_positions)

    return modeller


def create_simulation(system, 
                      modeller,
                      integrator,
                      integrator_args,
                      platform_name="CUDA"):
    """
    Creates a Simulation object with a given integrator and settings.
    Args:
        system: OpenMM System object
        modeller: OpenMM Modeller object
        integrator: str, name of the integrator to use
            Choose from: 'LangevinIntegrator', 'VerletIntegrator', 'NoseHooverIntegrator'
        integrator_args: dict, arguments to pass to the integrator
            These follow OpenMM's integrator constructors.
            E.g. For LangevinIntegrator: {'temperature': 300*unit.kelvin, 'friction': 1/unit.picoseconds, 'timestep': 2*unit.femtoseconds}
        platform_name: str, name of the OpenMM platform to use (default: "CUDA")
    Returns:
        simulation: OpenMM Simulation object
    """
    integrators_dict = {
        'LangevinIntegrator': openmm.LangevinIntegrator,
        'VerletIntegrator': openmm.VerletIntegrator,
        'NoseHooverIntegrator': openmm.NoseHooverIntegrator,
    }

    # The integrators use positional arguments instead of named ones; FML.
    # the arguments that each integrator takes, in order:
    integrator_args_order = {
        'LangevinIntegrator': ['temperature', 'friction', 'timestep'],
        'VerletIntegrator': ['timestep'],
        'NoseHooverIntegrator': ['temperature', 'collision_rate', 'timestep']
    }
    # We create a list of arguments to pass to the integrator constructor:
    ordered_args = [integrator_args[arg] for arg in integrator_args_order[integrator]]

    topology, positions = modeller.topology, modeller.positions
    integrator = integrators_dict[integrator](*ordered_args)
    platform = openmm.Platform.getPlatformByName(platform_name)
    simulation = openmm.app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    return simulation
