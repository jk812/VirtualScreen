#!/usr/bin/env python3

import argparse

from workflow_functions import run_workflow
from workflow_io import generate_default_workflow_config

parser = argparse.ArgumentParser(
                    prog='ANNalog with docking',
                    description='This program combines the ANNalog Seq2Seq model with molecular docking to bias generation towards a target protein.',
                    )

parser.add_argument('-c', '--config', default='workflow_config.ini', help='Path to the configuration file. Default: workflow_config.ini')
parser.add_argument('--regen_default_config', action='store_true', help='Generate a new default config file (workflow_config.ini) and overwrite any existing one')

args = parser.parse_args()

if args.regen_default_config:
    generate_default_workflow_config()
else:
    print(f"Using configuration file: {args.config}")
    run_workflow(config_file_path=args.config)
