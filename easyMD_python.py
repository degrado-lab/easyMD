from easyMD.initialize import init
from easyMD.prepare import prep
from easyMD.run import pre_production, production

init('1AK4', 'output/1AK4', chains=['A', 'B'], ligands=['HEM'], model_termini=True)

prep_params = {'pH': 7.0, 
               'ion_concentration_MOLAR': 0.15, 
               'padding_ANSTROM': 10}

prep('output/1AK4', {'pH': 7.0, 'ion_concentration_MOLAR': 0.15, 'padding_ANSTROM': 10})

sim_params = {"minimization_length_NS":  0,
            "equilibration_length_NS": 1,
            "simulation_length_NS":    1,
            "step_size_FS":            2,

            "temperature_K":          300, 
            "pressure_ATM":           1.0,
            "forcefield":           ["amber14-all.xml", "amber14/tip3p.xml"],
            
            "save_frame_interval_NS":  0.01, 
            "checkpoint_interval_NS":  0.01,}

pre_production('output/1AK4', sim_params)

production('output/1AK4', sim_params)

