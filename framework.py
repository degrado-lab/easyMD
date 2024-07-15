### What do I want this process to look like? 

# command line:
# easymd --input_structure complex.pdb --output output_dir/

# Calls CLI. CLI handles input/output, and calls the main function.

# Python interface:

universe = easymd.init('complex.pdb') #what exactly does this do? Split up the complex into its components, and save as cif?
#returns an "easyMD universe", which wraps up:
# openMM system
# list of components (ligands, proteins, ions, water)
# other?

prep_params = {'pH': 7.4, 'forcefield': 'amber14', 'water_model': 'tip3p'} #default parameters...
easymd.prep(universe, prep_params) #Do prep including adding Hs, solvating, etc
#optionally, provide your own prep function. Prep_function = default...
#this function will take an input and prep params. will look like:

def prep(universe, prep_params):
    # Load file input stuff
    # Load what?

    # Run user-specified prep function
    prepped_universe = prep_function(universe, prep_params)

    # Save output
    # Save prepped system to input/prepped dir

    return universe

easymd.pre_production(universe, pre_prod_params) #Do pre-production minimization and equilibration

easymd.production(universe, prod_params) #Do production run

easymd.post_production(universe, post_prod_params) #Do post-production analysis




### Universe object?
# Universe object will contain:
class easyMD_universe:
    def __init__(self, directory):
        self.directory = directory
        self.prep = None
        self.pre_production = None
        self.production = None
        self.post_production = None

class easyMD_segment:
    def __init__(self, name, files):
        # Contains all the info for each segment of an easyMD run...
        
        self.name = name
        self.files 