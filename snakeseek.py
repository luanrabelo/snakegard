#!/usr/bin/env python3

import os
import subprocess
import argparse
import logging
import time
import shutil
import sys
import yaml

import pandas as pd

__authors__     = ""
__license__     = ""
__version__     = "0.0.0"
__maintainers__ = ""
__email__       = "luan.rabelo@pq.itv.org"
__date__        = "2025/10/15"
__github__      = "luanrabelo/snakeseek"
__status__      = "Development"
__tool__        = "SnakeSeek"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y/%m/%d - %H:%M:%S'
)

parser = argparse.ArgumentParser(
    prog=__tool__,
    usage='%(prog)s [options]',
    #description=pipenote_ascii_art,
    epilog=f"In Development...",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=True,
    )

excel = parser.add_argument_group(
    title="Excel with RNA-Seq data information. See documentation for details.",
    description="Specify the absolute path to the Excel file. File must be in .xlsx or .xls format."
    )

excel.add_argument(
    "-ex", "--excel",
    default="example/snakeseek.xlsx",
    required=True,
    type=str,
    help="Path to the Excel file."
    )

threads = parser.add_argument_group(
    title="Threads",
    description="Number of threads to use."
    )
threads.add_argument(
    "-t", "--threads",
    default=8,
    type=int,
    help="Number of threads to use. Default is 8."
    )

test = parser.add_argument_group(
    title="Test workflow with --dry-run option",
    description="Run in test mode with example data."
    )
test.add_argument(
    "--dry-run",
    action="store_true",
    help="Run the workflow in dry-run mode without executing any commands."
    )

verbose = parser.add_argument_group(
    title="Verbose",
    description="Enable verbose output."
    )
verbose.add_argument(
    "-v", "--verbose",
    action="store_true",
    help="Enable verbose output."
    )

args = parser.parse_args()

verbose = True if args.verbose else False

#print(pipenote_ascii_art)
print(f"\n{'#'*70}\n")
print(f"Version:     {__version__}")
print(f"Status:      {__status__}")
print(f"Maintainers: {__maintainers__}")
print(f"License:     {__license__}")
print(f"GitHub:      {__github__}")
print(f"\n{'#'*70}\n")

if __name__ == "__main__":
    tool_env = "snakegard"
    
    if not shutil.which("conda"):
        logging.error("Conda is not installed or not found in PATH. Please install Conda and try again.")
        sys.exit(1)

    conda_envs = subprocess.run(
        ["conda", "env", "list"],
        capture_output=True,
        text=True
        ).stdout.splitlines()
    
    if not any(tool_env in env for env in conda_envs):
        logging.error(f"The '{tool_env}' conda environment does not exist. Please create it using the provided environment.yml file.")
        sys.exit(1)

    conda_env = os.getenv('CONDA_DEFAULT_ENV')
    if conda_env != tool_env:
        logging.warning(f"Activate the '{tool_env}' conda environment before running the {__tool__}. Current environment: {conda_env}. Run 'conda activate {tool_env}' and try again.")
        sys.exit(1)

    if not os.path.exists(args.excel) or not args.excel.endswith(('.xlsx', '.xls')) or os.path.getsize(args.excel) == 0:
        logging.error(f"Excel file {args.excel} is missing, not a valid Excel file, or empty. Please provide a valid Excel file.")
        sys.exit(1)
    
    df = pd.read_excel(args.excel, sheet_name=0, engine='openpyxl')
    if df.empty:
        logging.error(f"Excel file {args.excel} is empty. Please provide a valid Excel file with data.")
        sys.exit(1)

    with open('config.yaml', 'w') as config_file:
        config_file.write("samples:\n")
        for index, row in df.iterrows():
            config_file.write(f"  - species: {row['Species']}\n")
            config_file.write(f"    sra:\n")
            config_file.write(f"      - id: {row['SRA Run ID']}\n")
            config_file.write(f"        reference: {row['Reference ID']}\n")
            config_file.write(f"        reference_type: {row['Reference Type']}\n")
            config_file.write(f"        mapper: {row['Mapper']}\n")
    logging.info(f"Configuration file 'config.yml' generated successfully from {args.excel}.")

    # Check tmp directory
    tmp_dir = os.path.join(os.getcwd(), "tmp")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir, exist_ok=True, mode=0o755)
        logging.info(f"Temporary directory '{tmp_dir}' created.")

    os.environ['TMPDIR'] = tmp_dir
    os.environ['TEMP'] = tmp_dir
    os.environ['TMP'] = tmp_dir

    # Check and unlock if necessary
    lock_dir = os.path.join(os.getcwd(), ".snakemake", "locks")
    if os.path.exists(lock_dir):
        logging.warning("Snakemake lock detected. Attempting to unlock...")
        unlock_cmd = [
            "snakemake",
            "-d", str(os.getcwd()),
            "-s", "workflow/Snakefile",
            "--unlock"
        ]
        subprocess.run(unlock_cmd, check=False, capture_output=True)
        logging.info("Lock removed successfully, continuing with workflow execution.")

        
    cmd_args = [
        "snakemake",
        "-d", str(os.getcwd()),
        "-s", "workflow/Snakefile",
        "--cores", str(args.threads),
    ]

    if args.dry_run:
        cmd_args.append("--dry-run")
        logging.info("Running Snakemake in dry-run mode...")

    logging.info(f"Running Snakemake workflow with {args.threads} threads...")
    status = subprocess.run(cmd_args, check=True, text=True, env=os.environ.copy())
    if status.returncode != 0:
        logging.error("Failed to run Snakemake workflow. Please check the logs above for details.")
        sys.exit(1)
    else:
        shutil.rmtree(tmp_dir)
        logging.info(f"Snakemake workflow completed successfully. See 'results' directory for output files.")
