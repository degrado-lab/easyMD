def write_pdb_from_simulation(simulation, filename):
    """
    Write a PDB file from an OpenMM Simulation object.
    """
    import openmm.app as app
    # output initial topology file:
    with open(filename, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)