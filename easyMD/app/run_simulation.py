import openmm.app as app
# import System and Modellers:
from openmm.app import Modeller
from openmm import unit
from easyMD.prep import system as prep_system
import logging

logger = logging.getLogger(__name__)

def run_simulation(system, modeller: app.Modeller, output_file: str, duration: unit.Quantity, output_frequency: unit.Quantity = 1000*unit.femtoseconds, energy_minimize: bool = True, relax: bool = True, relax_duration: unit.Quantity = 1*unit.nanoseconds):
    """
    Run a simulation for a given duration, and output to a DCD file.
    By default, the system is energy minimized until convergence, and equilibrated for 1 ns.
    Inputs:
        simulation: OpenMM Simulation object
        modeller: OpenMM Modeller object
        output_file: str, path to output file
        duration: How long to run the simulation. A duration in OpenMM time units (e.g. 100*unit.nanoseconds)
        output_frequency: frequency of output in OpenMM time units (default: 0.1*unit.nanoseconds)
        energy_minimize: bool, whether to energy minimize the system (default: True)
        equilibrate: bool, whether to equilibrate the system (default: True)
        equilibration_duration: duration of equilibration in OpenMM time units (default: 1*unit.nanoseconds)
    """
    from sys import stdout
    from openmm import MonteCarloBarostat

    # Step 1: Create the simulation:
    simulation = prep_system.create_simulation(
        system=system,
        modeller=modeller,
        integrator='LangevinIntegrator',
        integrator_args={
            'temperature': 300*unit.kelvin,
            'friction': 1/unit.picoseconds,
            'timestep': 2*unit.femtoseconds
        }
    )

    # What's the time step?
    time_step = simulation.integrator.getStepSize()

    # How many steps to run?
    relax_steps = int(relax_duration / time_step)
    prod_steps = int(duration / time_step)
    output_frequency_steps = int(output_frequency / time_step)

    logger.info(f'Duration: {duration} \n Output Frequency: {output_frequency} ')

    if energy_minimize:
        logger.info("Minimizing energy...")
        simulation.minimizeEnergy(
            #tolerance=10 * unit.kilojoules_per_mole, #I'm not sure what units this should actually be...
            maxIterations= 0 # Until convergence
        )
    
    # Write the output file:
    simulation.reporters.append(app.StateDataReporter(stdout, output_frequency_steps, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
    
    # Equilibration is NVT, then NPT
    if relax:
        # NVT
        logger.info("Equilibrating...")
        simulation.step(int( relax_steps / 2))

        # NPT
        simulation.system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin))
        simulation.context.reinitialize(preserveState=True)
        simulation.step(int( relax_steps / 2))

    # Add the DCD reporter and run the production simulation:
    logger.info("Running production simulation...")
    simulation.reporters.append(app.DCDReporter(output_file, output_frequency_steps))
    simulation.step(prod_steps)

'''
Eventually, it would be nice for the sim to be defined by a comprehensive dictionary, like the following.
For now, we'll be a bit more rigid but practical.
{
    'EM:': {
        'type': 'minimize',
        'tolerance': 10 * unit.kilojoules_per_mole,
        'maxIterations': 0
    },
    'relax_NVT': {
        'duration': 1*unit.nanoseconds,
        'time_step': 2*unit.femtoseconds,
        'integrator': {
            'type': 'LangevinIntegrator',
            'temperature': 300*unit.kelvin,
            'friction': 1/unit.picoseconds,
        },
        'Reporter': {
            'type': 'StateDataReporter',
            'stdout': True,
            'output_frequency': 1000*unit.femtoseconds,
            'step': True,
            'potentialEnergy': True,
            'temperature': True,
            'volume': True,
            'density': True
            },
        'output_frequency': None,
    },
    'relax_NPT': {
        'duration': 1*unit.nanoseconds,
        'time_step': 2*unit.femtoseconds,
        'integrator': {
            'type': 'LangevinIntegrator',
            'temperature': 300*unit.kelvin,
            'friction': 1/unit.picoseconds,
        },
        'barostat': {
            'type': 'MonteCarloBarostat',
            'pressure': 1*unit.atmosphere,
        },
        'output_frequency': None,
        'Reporter': {
            'type': 'StateDataReporter',
            'output_frequency': 1000*unit.femtoseconds,
            },
        'Reporter': {
            'type': 'DCDReporter',
            'output_file': output_file,
            'output_frequency': 1000*unit.femtoseconds,
    },


        }
}
'''

# def run_simulation(system: app.System, output_file: str, duration: unit.Quantity, output_frequency: unit.Quantity = 1000*unit.femtoseconds, energy_minimize: bool = True, relax: bool = True, relax_duration: unit.Quantity = 1*unit.nanoseconds):
#     """
#     Run a simulation for a given duration, and output to a DCD file.
#     By default, the system is energy minimized until convergence, and equilibrated for 1 ns.
#     Inputs:
#         simulation: OpenMM Simulation object
#         output_file: str, path to output file
#         duration: How long to run the simulation. A duration in OpenMM time units (e.g. 100*unit.nanoseconds)
#         output_frequency: frequency of output in OpenMM time units (default: 0.1*unit.nanoseconds)
#         energy_minimize: bool, whether to energy minimize the system (default: True)
#         equilibrate: bool, whether to equilibrate the system (default: True)
#         equilibration_duration: duration of equilibration in OpenMM time units (default: 1*unit.nanoseconds)
#     """
#     from sys import stdout
#     from openmm import MonteCarloBarostat

#     # What's the time step?
#     time_step = simulation.integrator.getStepSize()
#     # How many steps to run?
#     relax_steps = int(relax_duration / time_step)
#     prod_steps = int(duration / time_step)
#     output_frequency_steps = int(output_frequency / time_step)

#     logger.info(f'Duration: {duration} \n Output Frequency: {output_frequency} ')

#     if energy_minimize:
#         logger.info("Minimizing energy...")
#         simulation.minimizeEnergy()
    
#     # Write the output file:
#     simulation.reporters.append(app.StateDataReporter(stdout, output_frequency_steps, step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
    
#     # Equilibration is NVT, then NPT
#     if relax:
#         # NVT
#         logger.info("Equilibrating...")
#         simulation.step(int( relax_steps / 2))

#         # NPT
#         simulation.system.addForce(MonteCarloBarostat(1*unit.atmosphere, 300*unit.kelvin))
#         simulation.context.reinitialize(preserveState=True)
#         simulation.step(int( relax_steps / 2))

#     # Add the DCD reporter and run the production simulation:
#     logger.info("Running production simulation...")
#     simulation.reporters.append(app.DCDReporter(output_file, output_frequency_steps))
#     simulation.step(prod_steps)
    