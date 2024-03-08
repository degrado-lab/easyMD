import nglview as nv
import mdtraj as md

def show_pdb(pdb_path, water=False, size=[1600, 800]):
    view = nv.show_file(str(pdb_path))
    view.clear_representations()
    view.add_representation('cartoon', selection='all', color='residueindex')
    #show heteroatoms larger:
    view.add_representation('ball+stick', selection='hetero and not water', aspectRatio=4.0)
    if water:
        view.add_representation('ball+stick', selection='water', opacity=0.1, aspectRatio=4.0)
    #get size as string
    size_str = [str(s) + 'px' for s in size]
    view._set_size(*size_str)
    return view


def show_traj(pdb_path, dcd_path, size=[1600, 800]):
    #show the trajectory from PDB with nglview

    traj = md.load(
                str(dcd_path),
        top=    str(pdb_path)
    )

    #center the trajectory around the protein (to stop periodic boundary artifacts)
    protein = traj.topology.select('protein')
    protein_anchor_molecules = [(atom for i, atom in enumerate(traj.topology.atoms) if i in protein)] #this is a strange way to create a list of a set of atom objects, which is used by image_molecules() to align the trajectory.
    centered_traj = traj.image_molecules(anchor_molecules=protein_anchor_molecules)

    #first let's make copy of trajectory without water:
    traj_no_water = centered_traj.atom_slice(centered_traj.topology.select('not water'))

    view = nv.show_mdtraj(traj_no_water)
    view.clear_representations()
    view.add_representation('cartoon', selection='all', color='residueindex')
    view.add_representation('ball+stick', selection='hetero and not water', aspectRatio=4.0)
    
    #get size as string
    size_str = [str(s) + 'px' for s in size]
    view._set_size(*size_str)
    return view