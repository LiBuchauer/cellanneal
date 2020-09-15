"""This script defines the command line entry point for cell anneal and implements
a pipeline including identification of highly variable genes, deconvolution of each
bulk sample, ouput of deconvolution results and production of a set of plots."""

import argparse
from pathlib import Path

import pandas as pd

from .general import make_gene_dictionary, deconvolve
from .plots import (plot_pies_from_df, plot_mix_heatmap,
                    plot_1D_lines, plot_repeats, plot_scatter)
from .stability import repeat_annealing
import sys


def init_parser(parser):
    """Initialize parser arguments."""
    parser.add_argument(
        'bulk_data', type=str,
        help=("Path to bulk data file, csv format with sample names as columns and genes as rows.")
    )

    parser.add_argument(
        'celltype_data', type=str,
        help=("Path to celltype reference data file, csv format with sample names as columns and genes as rows. ")
    )

    parser.add_argument(
        '--bulk_min', type=float, default=1e-5,
        help=("Minimum relative expression in the bulk sample for a gene to be considered.")
    )

    parser.add_argument(
        '--bulk_max', type=float, default=0.01,
        help=("Maximum relative expression in the bulk sample for a gene to be considered.")
    )

    parser.add_argument(
        '--disp_min', type=float, default=0.5,
        help=("Minimum scaled dispersion (var/mean) relative to other genes of similar expression for a gene to be considered.")
    )

    parser.add_argument(
        '--maxiter', type=int, default=1000,
        help=("Number of initial values from which to run simulated annealing.")
    )

    return parser



def main():
    """ cellanneal. User-friendly deconvolution of bulk RNA-Seq data.

    Deconvovles bulk RNA-Seq data into single-cell populations based
    on Spearman's correlation between bulk sample gene expression and
    single-cell based mixture gene expression. Simulated annealing is
    used to find the optimal mixture.

    Input:
            bulk_data
            celltype_data
            bulk_min
            bulk_max
            disp_min
            maxiter

    Output:

    Example:


    """
    # get a parser object, initialise the inputs and read them into args
    my_parser = argparse.ArgumentParser(description=('cellanneal deconvolves bulk RNA-Seq data.'))
    args = init_parser(my_parser).parse_args()

    # grab the individual inputs for further use
    bulk_import_path = args.bulk_data
    sc_ref_import_path = args.celltype_data
    bulk_min = args.bulk_min
    bulk_max = args.bulk_max
    disp_min = args.disp_min
    maxiter = args.maxiter

    print("\nWelcome to cellanneal Version 0.1.0!\n")

    """ 1) Import bulk and cell type data """
    print('1A. Importing bulk data ...')
    bulk_df = pd.read_csv(bulk_import_path, index_col=0, header=0)
    bulk_names = bulk_df.columns.tolist()
    print('{} bulk samples identified: {}\n'.format(len(bulk_names), bulk_names))

    print('1B. Importing celltype reference data ...')
    # import single cell based reference
    sc_ref_df = pd.read_csv(sc_ref_import_path, index_col=0, header=0)
    celltypes = sc_ref_df.columns.tolist()
    print('{} cell types identified: {}\n'.format(len(celltypes), celltypes))

    """ 2) Identify highly variable genes and genes that pass the thresholds for each bulk. """
    # produce lists of genes on which to base deconvolution
    print('2. Constructing gene sets ...')
    gene_dict = make_gene_dictionary(
                    sc_ref_df,
                    bulk_df,
                    disp_min=disp_min,
                    bulk_min=bulk_min,
                    bulk_max=bulk_max,
                    remove_mito=True)

    """ 3) Run cellanneal. """
    print('\n3. Running cellanneal ...')
    all_mix_df = deconvolve(sc_ref_df=sc_ref_df,
                     bulk_df=bulk_df,
                     maxiter=maxiter,
                     gene_dict=gene_dict,
                     no_local_search=False)

    """ 4) Write results to file."""
    print('\n4. Writing results to file in folder "results" ...')
    output_name = 'results_' + bulk_import_path
    result_path = Path('results/') / output_name
    all_mix_df.to_csv(result_path, header=True, index=True, sep=',')

    """ 5) Produce plots and save to folder"""
    print('\n5. Storing figures in folder "figures" ...')
    # plot results
    figure_path = Path('figures/')
    pie_path = figure_path / 'pies_{}.pdf'.format(bulk_import_path.split('.')[-2])
    plot_pies_from_df(all_mix_df, save_path=pie_path)

    heat_path = figure_path / 'heat_{}.pdf'.format(bulk_import_path.split('.')[-2])
    plot_mix_heatmap(all_mix_df, rownorm=False, save_path=heat_path)

    line_path = figure_path / 'lines_{}.pdf'.format(bulk_import_path.split('.')[-2])
    plot_1D_lines(all_mix_df, line_path)

    scatter_path = figure_path / 'scatter_{}.pdf'.format(bulk_import_path.split('.')[-2])
    plot_scatter(all_mix_df, bulk_df, sc_ref_df, gene_dict, save_path=scatter_path)

    print('\nAll done! :-)\n')



def repeat():
    """ Repeat anneal.


    """
    # get a parser object, initialise the inputs and read them into args
    # in a first step, the same arguments as in the main function are required/allowed
    my_parser = argparse.ArgumentParser(description=('cellanneal deconvolves bulk RNA-Seq data.'))
    my_parser = init_parser(my_parser)

    # however, we also need additional arguments for the subsampling which are unique
    # unique to the repeat function and will be added below
    my_parser.add_argument(
        '--N_repeat', type=int, default=10,
        help=("Number of times a subset of genes is drawn and deconvolution is performed on that set of genes.")
    )

    my_parser.add_argument(
        '--sample_fraction', type=float, default=0.9,
        help=("Fraction of the original gene list to include in each random sample, (0, 1].")
    )

    args = my_parser.parse_args()

    # grab the individual inputs for further use
    bulk_import_path = args.bulk_data
    sc_ref_import_path = args.celltype_data
    bulk_min = args.bulk_min
    bulk_max = args.bulk_max
    disp_min = args.disp_min
    maxiter = args.maxiter
    N_repeat = args.N_repeat
    sample_fraction = args.sample_fraction

    print("\nWelcome to cellanneal Version 0.1.0!\n\nYou have chosen to repeatedly anneal on subsets of the full gene list.")

    """ 0. Check the inputs """
    if sample_fraction < 0 or sample_fraction > 1:
        sys.exit('Error: Please enter a sample_fraction between 0 and 1.')

    """ 1) Import bulk and cell type data """
    print('1A. Importing bulk data ...')
    bulk_df = pd.read_csv(bulk_import_path, index_col=0, header=0)
    bulk_names = bulk_df.columns.tolist()
    print('{} bulk samples identified: {}\n'.format(len(bulk_names), bulk_names))

    print('1B. Importing celltype reference data ...')
    # import single cell based reference
    sc_ref_df = pd.read_csv(sc_ref_import_path, index_col=0, header=0)
    celltypes = sc_ref_df.columns.tolist()
    print('{} cell types identified: {}\n'.format(len(celltypes), celltypes))

    """ 2) Identify highly variable genes and genes that pass the thresholds for each bulk. """
    # produce lists of genes on which to base deconvolution
    print('2. Constructing base gene sets ...')
    gene_dict = make_gene_dictionary(
                    sc_ref_df,
                    bulk_df,
                    disp_min=disp_min,
                    bulk_min=bulk_min,
                    bulk_max=bulk_max,
                    remove_mito=True)

    """ 3) Run cellanneal. """
    print('\n3. Running cellanneal {} times...'.format(N_repeat))
    master_df = repeat_annealing(
                    sc_ref_df,
                    bulk_df,
                    gene_dict,
                    no_local_search=False,
                    N_repeat=N_repeat,
                    sample_fraction=sample_fraction,
                    maxiter=maxiter)

    """ 4) Write results to file."""
    print('\n4. Writing results to file in folder "results" ...')
    # all results
    output_name = 'repeat_N={}_fraction={}_'.format(N_repeat, sample_fraction) + bulk_import_path
    result_path = Path('results/') / output_name
    master_df.to_csv(result_path, header=True, index=True, sep=',')

    # mean result
    mean_name = 'mean_repeat_N={}_fraction={}_'.format(N_repeat, sample_fraction) + bulk_import_path
    result_path = Path('results/') / mean_name
    mean_df = master_df.groupby(['bulk', 'celltype']).mean()
    mean_df = mean_df.drop(['run'], axis=1)
    mean_df.to_csv(result_path, header=True, index=True, sep=',')

    """ 5) Produce plot and save to folder"""
    print('\n5. Storing figure in folder "figures" ...')
    # plot results
    figure_path = Path('figures/')
    repeat_path = figure_path / 'boxes_repeat_N={}_fraction={}_{}.pdf'.format(N_repeat, sample_fraction, bulk_import_path.split('.')[-2])
    plot_repeats(master_df, save_path=repeat_path)

    print('\nAll done! :-)\n')
