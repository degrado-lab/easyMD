import argparse
import os

def main():
    '''This script is the entrypoint for our EasyMD apptainer.
    It will take two subcommands, 'sim' or 'notebook', and run the appropriate
    script.
    sim: takes as argument a parameter file and runs the simulation using easyMD.utils.run_sim_local()
    notebook: starts a jupyter notebook server
    '''
    parser = argparse.ArgumentParser(description='EasyMD apptainer')
    subparsers = parser.add_subparsers(required=True, dest='subcommand')
    parser_sim = subparsers.add_parser('sim', help='Run a simulation')
    parser_sim.add_argument('param_file', type=str, help='Parameter file for the simulation')
    parser_sim.set_defaults(subcommand='sim')

    parser_notebook = subparsers.add_parser('notebook', help='Start a jupyter notebook server')
    parser_notebook.set_defaults(subcommand='notebook')

    args = parser.parse_args()

    if args.subcommand == 'sim':
        from easyMD.utils import run_sim_local
        run_sim_local(args.param_file)
    elif args.subcommand == 'notebook':
        os.system('jupyter notebook --port=$JUPYTER_PORT')


if __name__ == '__main__':
    main()