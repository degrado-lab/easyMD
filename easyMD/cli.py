import argparse
from easyMD.initialize import init
from easyMD.prepare import prep_wrapper
from easyMD.run import pre_production_wrapper, production_wrapper

def init_cli():
    parser = argparse.ArgumentParser(description='easyMD: A set of utilities and notebooks for running OpenMM simulations.')
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')

    # Initialize
    init_parser = subparsers.add_parser('init', help='Initialize an easyMD run.')
    init_parser.add_argument('input_file', help='Input structure file.')
    init_parser.add_argument('experiment_dir', help='Directory to save the output files.')
    init_parser.add_argument('--chains', nargs='+', help='Chains to keep.', default=[])
    init_parser.add_argument('--ligands', nargs='+', help='Ligands to split out.', default=[])
    init_parser.add_argument('--ligand_files', nargs='+', help='Paths to ligand SDF files to use instead of splitting out ligands.', default=[])
    init_parser.add_argument('--model_termini', action='store_true', help='Model the termini.', default=False)
    init_parser.set_defaults(func=init)

    # Prepare
    prep_parser = subparsers.add_parser('prep', help='Prepare the input structure for simulation.')
    prep_parser.add_argument('experiment_dir', help='Directory to save the output files.')
    prep_parser.add_argument('parameter_file', help='JSON file with parameters for the preparation.')
    prep_parser.add_argument('--custom_func', type=str, help='Path to custom Python file containing the prep() function.', default=None)
    prep_parser.set_defaults(func=prep_wrapper)

    # Pre-Production (Energy Minimization and Equilibration)
    pre_prod_parser = subparsers.add_parser('pre_prod', help='Run energy minimization and equilibration.')
    pre_prod_parser.add_argument('experiment_dir', help='Directory to save the output files.')
    pre_prod_parser.add_argument('parameter_file', help='JSON file with parameters for the preparation.')
    pre_prod_parser.add_argument('--custom_func', type=str, help='Path to custom Python file containing the pre_prod() function.', default=None)
    pre_prod_parser.set_defaults(func=pre_production_wrapper)

    # Production
    prod_parser = subparsers.add_parser('prod', help='Run production simulation.')
    prod_parser.add_argument('experiment_dir', help='Directory to save the output files.')
    prod_parser.add_argument('parameter_file', help='JSON file with parameters for the preparation.')
    prod_parser.add_argument('--custom_func', type=str, help='Path to custom Python file containing the prod() function.', default=None)
    prod_parser.set_defaults(func=production_wrapper)

    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
    elif args.command == 'init':
        args.func(args.input_file, args.experiment_dir, args.chains, args.ligands, args.ligand_files, args.model_termini)
    elif args.command == 'prep':
        args.func(args.experiment_dir, args.parameter_file, args.custom_func)
    elif args.command == 'pre_prod':
        args.func(args.experiment_dir, args.parameter_file, args.custom_func)
    elif args.command == 'prod':
        args.func(args.experiment_dir, args.parameter_file, args.custom_func)

# def init_cli():
#     parser = argparse.ArgumentParser(description='Initialize an easyMD run.')
#     parser.add_argument('input_file', help='Input structure file.')
#     parser.add_argument('experiment_dir', help='Directory to save the output files.')
#     parser.add_argument('--chains', nargs='+', help='Chains to keep.', default=[])
#     parser.add_argument('--ligands', nargs='+', help='Ligands to split out.', default=[])
#     parser.add_argument('--model_termini', action='store_true', help='Model the termini.', default=False)
#     args = parser.parse_args()

#     init(args.input_file, args.experiment_dir, args.chains, args.ligands, args.model_termini)