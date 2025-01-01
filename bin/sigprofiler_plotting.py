#!/usr/bin/env python3

from sigProfilerPlotting import plotSBS, plotDBS, plotID
from pathlib import Path
from glob import glob
import logging 

import argparse, sys

logging.basicConfig(filename='sigprofiler_plotting.log', filemode='w', level = logging.INFO)

def plot_spectra(matrix_dir, project, output_dir):
    matrix_dir = Path(matrix_dir)
    # SBS 
    matrix_paths = glob(f'{matrix_dir}/**/*SBS*.all', recursive=True)
    output_path = f'{output_dir}/SBS/'
    for matrix_path in matrix_paths:
        plot_type=Path(matrix_path).name.split('.')[1][3:]
        logging.info(f'plotting SBS {plot_type}')
        plotSBS(str(matrix_path), output_path, project, plot_type, percentage=False)
        
    # DBS
    matrix_paths = glob(f'{matrix_dir}/**/*DBS*.all', recursive=True)
    output_path = f'{output_dir}/DBS/'
    for matrix_path in matrix_paths:
        plot_type=Path(matrix_path).name.split('.')[1][3:]
        logging.info(f'plotting DBS {plot_type}')
        plotDBS(str(matrix_path), output_path, project, plot_type, percentage=False)

    # ID
    matrix_paths = glob(f'{matrix_dir}/**/*ID*.all', recursive=True)
    output_path = f'{output_dir}/ID/'
    for matrix_path in matrix_paths:
        plot_type=Path(matrix_path).name.split('.')[1][2:]
        logging.info(f'plotting ID {plot_type}')
        if plot_type == '96': # there's an issue with format checking for ID96
            continue
        plotID(str(matrix_path), output_path, project, plot_type, percentage=False)
   
    
def get_arguments():
    parser = argparse.ArgumentParser(description='plot SBS and DBS spectra for SNPs, or ID for Indels')
    parser.add_argument('--matrix_dir', required=True,
                        help='path to the directory that contains the matrices, ie output from SigProfilerMatrixGenerator', type = str)
    parser.add_argument('--project', required=True,
                        help='path to the concordance file', type = str)
    parser.add_argument('--output_path', '-o', required=True,
                        help='output directory, has to have a slash for sigprofiler plotting to work', type = str)
    return parser

def main(args):
    plot_spectra(args.matrix_dir, args.project, args.output_path)

if __name__ == "__main__":
    args = get_arguments().parse_args()
    sys.exit(main(args))

            
