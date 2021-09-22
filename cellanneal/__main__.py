"""This script defines the command line entry point for cell anneal and
implements a pipeline including identification of highly variable genes,
deconvolution of each bulk sample, ouput of deconvolution results and
production of a set of plots."""

import argparse
from pathlib import Path
import pandas as pd

from .pipelines import cellanneal_pipe, repeatanneal_pipe


def init_parser(parser):
    """Initialize parser arguments."""
    parser.add_argument(
        'bulk_data_path', type=str,
        help=("""Path to bulk data file, csv format with sample names as
            columns and genes as rows.""")
    )

    parser.add_argument(
        'celltype_data_path', type=str,
        help=("""Path to celltype reference data file, csv format with sample
            names as columns and genes as rows. """)
    )

    parser.add_argument(
        'output_path', type=str,
        help=("""Path to folder in which to store the results. """)
    )

    parser.add_argument(
        '--bulk_min', type=float, default=1e-5,
        help=("""Minimum relative expression in the bulk sample for a gene
            to be considered.""")
    )

    parser.add_argument(
        '--bulk_max', type=float, default=0.01,
        help=("""Maximum relative expression in the bulk sample for a gene
            to be considered.""")
    )

    parser.add_argument(
        '--disp_min', type=float, default=0.5,
        help=("""Minimum scaled dispersion (var/mean) relative to other
            genes of similar expression for a gene to be considered.""")
    )

    parser.add_argument(
        '--maxiter', type=int, default=1000,
        help=("""Number of initial values from which to run simulated
            annealing.""")
    )

    return parser


def main():
    """ cellanneal. User-friendly deconvolution of bulk RNA-Seq data.

    Deconvovles bulk RNA-Seq data into single-cell populations based
    on Spearman's correlation between bulk sample gene expression and
    single-cell based mixture gene expression. Simulated annealing is
    used to find the optimal mixture.

    Input:
            bulk_data_path
            celltype_data_path
            output_path
            bulk_min
            bulk_max
            disp_min
            maxiter

    Output:

    Example:


    """
    # get a parser object, initialise the inputs and read them into args
    my_parser = argparse.ArgumentParser(
        description=('cellanneal deconvolves bulk RNA-Seq data.'))
    args = init_parser(my_parser).parse_args()

    # grab the individual inputs for further use
    bulk_data_path = Path(args.bulk_data_path)
    celltype_data_path = Path(args.celltype_data_path)
    output_path = Path(args.output_path)
    bulk_min = args.bulk_min
    bulk_max = args.bulk_max
    disp_min = args.disp_min
    maxiter = args.maxiter

    print("\nWelcome to cellanneal!\n")

    """ 1) Import bulk and cell type data """
    print('1A. Importing bulk data ...')
    try:
        bulk_df = pd.read_csv(bulk_data_path, index_col=0, header=0)
        bulk_names = bulk_df.columns.tolist()
        print('{} bulk samples identified: {}\n'.format(len(bulk_names),
                                                        bulk_names))
    except:
        print("""Your bulk data file could not be imported.
        Please check the documentation for format requirements
        and look at the example bulk data file.""")

    print('1B. Importing celltype reference data ...')
    # import single cell based reference
    try:
        celltype_df = pd.read_csv(celltype_data_path, index_col=0, header=0)
        celltypes = celltype_df.columns.tolist()
        print('{} cell types identified: {}\n'.format(len(celltypes), celltypes))
    except:
        print("""Your celltype data file could not be imported.
        Please check the documentation for format requirements
        and look at the example celltype data file.""")

    # start pipeline
    cellanneal_pipe(
        celltype_data_path,
        celltype_df,
        bulk_data_path,
        bulk_df,
        disp_min,
        bulk_min,
        bulk_max,
        maxiter,
        output_path)


def repeat():
    """ repeatanneal. Runs cellanneal a given number of times and averages
    the results (takes the mean of cell fractions produced in each of the
    runs). Plots and gene expression predictions as known from cellanneal
    are produced for the mean result. Additional plots showing the spread
    of the repeats are also produced.

    Input:
            bulk_data_path
            celltype_data_path
            output_path
            bulk_min
            bulk_max
            disp_min
            maxiter
            N_repeat

    Output:

    Example:


    """
    # get a parser object, initialise the inputs and read them into args
    # in a first step, the same arguments as in the main function are
    # required/allowed
    my_parser = argparse.ArgumentParser(
        description=('cellanneal deconvolves bulk RNA-Seq data.'))
    my_parser = init_parser(my_parser)

    # however, we also need additional arguments for the subsampling which
    # are unique unique to the repeat function and will be added below
    my_parser.add_argument(
        '--N_repeat', type=int, default=10,
        help=("""Number of times a subset of genes is drawn and deconvolution
            is performed on that set of genes.""")
    )

    args = my_parser.parse_args()

    # grab the individual inputs for further use
    bulk_data_path = Path(args.bulk_data_path)
    celltype_data_path = Path(args.celltype_data_path)
    output_path = Path(args.output_path)
    bulk_min = args.bulk_min
    bulk_max = args.bulk_max
    disp_min = args.disp_min
    maxiter = args.maxiter
    N_repeat = args.N_repeat

    print("""\nWelcome to cellanneal!\n\n
    You have chosen to repeatedly anneal on subsets of the full gene list.""")

    """ 1) Import bulk and cell type data """
    print('1A. Importing bulk data ...')
    try:
        bulk_df = pd.read_csv(bulk_data_path, index_col=0, header=0)
        bulk_names = bulk_df.columns.tolist()
        print('{} bulk samples identified: {}\n'.format(len(bulk_names),
                                                        bulk_names))
    except:
        print("""Your bulk data file could not be imported.
        Please check the documentation for format requirements
        and look at the example bulk data file.""")

    print('1B. Importing celltype reference data ...')
    # import single cell based reference
    try:
        celltype_df = pd.read_csv(celltype_data_path, index_col=0, header=0)
        celltypes = celltype_df.columns.tolist()
        print('{} cell types identified: {}\n'.format(len(celltypes), celltypes))
    except:
        print("""Your celltype data file could not be imported.
        Please check the documentation for format requirements
        and look at the example celltype data file.""")

    # start pipeline
    repeatanneal_pipe(
        celltype_data_path,
        celltype_df,
        bulk_data_path,
        bulk_df,
        disp_min,
        bulk_min,
        bulk_max,
        maxiter,
        N_repeat,
        output_path)
