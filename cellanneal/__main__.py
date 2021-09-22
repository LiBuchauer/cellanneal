"""This script defines the command line entry point for cell anneal and
implements a pipeline including identification of highly variable genes,
deconvolution of each bulk sample, ouput of deconvolution results and
production of a set of plots."""

import argparse
from pathlib import Path
import time

import pandas as pd

from .general import make_gene_dictionary, deconvolve, calc_gene_expression
from .plots import (plot_pies_from_df, plot_mix_heatmap,
                    plot_1D_lines, plot_repeats, plot_scatter)
from .stability import repeat_annealing
import sys


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
        sc_ref_df = pd.read_csv(celltype_data_path, index_col=0, header=0)
        celltypes = sc_ref_df.columns.tolist()
        print('{} cell types identified: {}\n'.format(len(celltypes), celltypes))
    except:
        print("""Your celltype data file could not be imported.
        Please check the documentation for format requirements
        and look at the example celltype data file.""")

    # start pipeline
    cellanneal_pipe(
        celltype_data_path,
        sc_ref_df,
        bulk_data_path,
        bulk_df,
        disp_min,
        bulk_min,
        bulk_max,
        maxiter,
        output_path)


def cellanneal_pipe(
        celltype_data_path,  # path object!
        sc_ref_df,
        bulk_data_path,  # path object!
        bulk_df,
        disp_min,
        bulk_min,
        bulk_max,
        maxiter,
        output_path):  # path object!
    """ Serves as entrypoint into cellanneal pipeline for both gui and cli
    once all data and parameters have been collected. """

    """ 2) Identify highly variable genes and genes that pass the thresholds
    for each bulk. """
    # extract names
    bulk_names = bulk_df.columns.tolist()
    celltypes = sc_ref_df.columns.tolist()
    # produce lists of genes on which to base deconvolution
    print('2. Constructing gene sets ...')
    gene_dict = make_gene_dictionary(
                    sc_ref_df,
                    bulk_df,
                    disp_min=disp_min,
                    bulk_min=bulk_min,
                    bulk_max=bulk_max)

    """ 3) Run cellanneal. """
    print('\n3. Running cellanneal ...')
    all_mix_df = deconvolve(
                    sc_ref_df=sc_ref_df,
                    bulk_df=bulk_df,
                    maxiter=maxiter,
                    gene_dict=gene_dict,
                    no_local_search=False)

    """ 4) Write results to file."""
    print('\n4. Writing results to file ...')

    # make top level folder for all results from this run
    # get timestamp for labelling
    timestamp = time.asctime().replace(' ', '_').replace(':', '-')
    bulk_file_name = bulk_data_path.name
    bulk_file_ID = bulk_file_name.split(".")[0]
    top_folder_name = 'cellanneal_' + bulk_file_ID + '_' + timestamp
    top_folder_path = output_path / top_folder_name
    top_folder_path.mkdir(parents=True, exist_ok=True)

    # make subfolders for deconv, gen expr and figures
    deconv_folder_path = top_folder_path / 'deconvolution_results'
    deconv_folder_path.mkdir(parents=True, exist_ok=True)
    figure_folder_path = top_folder_path / 'figures'
    figure_folder_path.mkdir(parents=True, exist_ok=True)
    genexpr_folder_path = top_folder_path / 'genewise_comparison'
    genexpr_folder_path.mkdir(parents=True, exist_ok=True)

    # first, write the mix matrix to csv
    deconv_name = 'deconvolution_' + bulk_file_ID + '.csv'
    result_path = deconv_folder_path / deconv_name
    all_mix_df.sort_index(axis=0, inplace=True)
    all_mix_df.to_csv(result_path, header=True, index=True, sep=',')

    # next, write the actual and estimated gene expression to file,
    # this has to be done per sample as the genes are sample specific
    # a mix_df version without correlation entries is needed
    all_mix_df_no_corr = all_mix_df[celltypes]
    for sample_name in bulk_names:
        gene_comp_df = calc_gene_expression(
                            mix_vec=all_mix_df_no_corr.loc[sample_name],
                            bulk_vec=bulk_df[sample_name],
                            sc_ref_df=sc_ref_df,
                            gene_list=gene_dict[sample_name])
        # construct export path for this sample
        sample_gene_name = 'expression_' + bulk_file_ID + '_' + \
            sample_name + '.csv'
        sample_gene_path = genexpr_folder_path / sample_gene_name
        gene_comp_df.sort_index(axis=0, inplace=True)
        gene_comp_df.to_csv(sample_gene_path, header=True, index=True, sep=',')

    """ 5) Produce plots and save to folder"""
    print('\n5. Storing figures in folder "figures" ...')
    # plot results

    pie_path = figure_folder_path / 'pies_{}.pdf'.format(bulk_file_ID)
    plot_pies_from_df(all_mix_df, save_path=pie_path)

    heat_path = figure_folder_path / 'heat_{}.pdf'.format(bulk_file_ID)
    plot_mix_heatmap(all_mix_df, rownorm=False, save_path=heat_path)

    scatter_path = figure_folder_path / 'scatter_{}.pdf'.format(bulk_file_ID)
    plot_scatter(all_mix_df, bulk_df, sc_ref_df, gene_dict,
                 save_path=scatter_path)

    print('\nAll done! :-)\n')


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
        sc_ref_df = pd.read_csv(celltype_data_path, index_col=0, header=0)
        celltypes = sc_ref_df.columns.tolist()
        print('{} cell types identified: {}\n'.format(len(celltypes), celltypes))
    except:
        print("""Your celltype data file could not be imported.
        Please check the documentation for format requirements
        and look at the example celltype data file.""")

    # start pipeline
    repeatanneal_pipe(
        celltype_data_path,
        sc_ref_df,
        bulk_data_path,
        bulk_df,
        disp_min,
        bulk_min,
        bulk_max,
        maxiter,
        N_repeat,
        output_path)


def repeatanneal_pipe(
        celltype_data_path,
        sc_ref_df,
        bulk_data_path,
        bulk_df,
        disp_min,
        bulk_min,
        bulk_max,
        maxiter,
        N_repeat,
        output_path):

    """ 2) Identify highly variable genes and genes that pass the thresholds
    for each bulk. """
    # extract names
    bulk_names = bulk_df.columns.tolist()
    celltypes = sc_ref_df.columns.tolist()
    # produce lists of genes on which to base deconvolution
    print('2. Constructing base gene sets ...')
    gene_dict = make_gene_dictionary(
                    sc_ref_df,
                    bulk_df,
                    disp_min=disp_min,
                    bulk_min=bulk_min,
                    bulk_max=bulk_max)

    """ 3) Run cellanneal. """
    print('\n3. Running cellanneal {} times...'.format(N_repeat))
    parent_df = repeat_annealing(
                    sc_ref_df,
                    bulk_df,
                    gene_dict,
                    no_local_search=False,
                    N_repeat=N_repeat,
                    maxiter=maxiter)

    """ 4) Write results to file."""
    print('\n4. Writing results to file ...')

    # make top level folder for all results from this run
    # get timestamp for labelling
    timestamp = time.asctime().replace(' ', '_').replace(':', '-')
    bulk_file_name = bulk_data_path.name
    bulk_file_ID = bulk_file_name.split(".")[0]
    top_folder_name = 'repeatanneal_N={}_'.format(N_repeat) + bulk_file_ID + '_' + timestamp
    top_folder_path = output_path / top_folder_name
    top_folder_path.mkdir(parents=True, exist_ok=True)

    # make subfolders for deconv, gen expr and figures
    deconv_folder_path = top_folder_path / 'deconvolution_results'
    deconv_folder_path.mkdir(parents=True, exist_ok=True)
    figure_folder_path = top_folder_path / 'figures'
    figure_folder_path.mkdir(parents=True, exist_ok=True)
    genexpr_folder_path = top_folder_path / 'genewise_comparison'
    genexpr_folder_path.mkdir(parents=True, exist_ok=True)

    # first, write individual deconvolution results to file
    deconv_name = 'deconvolution_individual_' + bulk_file_ID + '.csv'
    result_path = deconv_folder_path / deconv_name
    parent_df.sort_index(axis=0, inplace=True)
    parent_df.to_csv(result_path, header=True, index=True, sep=',')

    # next, calculate a mean deconvolution result and store
    mean_name = 'deconvolution_mean_' + bulk_file_ID + '.csv'
    mean_path = deconv_folder_path / mean_name
    # make mean df
    mean_df_long = parent_df.groupby(['bulk', 'celltype']).mean()
    mean_df_long = mean_df_long.drop(['run'], axis=1)
    # the mean df is in long form currently, we want wide form
    mean_df = pd.pivot_table(mean_df_long, index='bulk', columns='celltype')
    mean_df.columns = mean_df.columns.droplevel(0)
    mean_df.index.name = None
    mean_df.sort_index(axis=0, inplace=True)
    mean_df.to_csv(mean_path, header=True, index=True, sep=',')

    # next, write the actual and estimated gene expression to file for the
    # mean result, this has to be done per sample as the genes are sample-
    # specific
    # a mix_df version without correlation entries is needed
    mean_df_no_corr = mean_df[celltypes]
    for sample_name in bulk_names:
        gene_comp_df = calc_gene_expression(
                            mix_vec=mean_df_no_corr.loc[sample_name],
                            bulk_vec=bulk_df[sample_name],
                            sc_ref_df=sc_ref_df,
                            gene_list=gene_dict[sample_name])
        # construct export path for this sample
        sample_gene_name = """expression_{}_{}.csv""".format(
                                bulk_file_ID,
                                sample_name)
        sample_gene_path = genexpr_folder_path / sample_gene_name
        gene_comp_df.sort_index(axis=0, inplace=True)
        gene_comp_df.to_csv(sample_gene_path, header=True, index=True, sep=',')

    """ 5) Produce plot and save to folder"""
    print('\nProducing figures ...')
    # plot results, first spread of repeated annealing
    box_path = figure_folder_path / 'boxplot_{}.pdf'.format(bulk_file_ID)
    plot_repeats(parent_df, save_path=box_path)

    # next, same plots as above, for mean expression of all repeats
    pie_path = figure_folder_path / 'pies_{}.pdf'.format(bulk_file_ID)
    plot_pies_from_df(mean_df, save_path=pie_path)

    heat_path = figure_folder_path / 'heat_{}.pdf'.format(bulk_file_ID)
    plot_mix_heatmap(mean_df, rownorm=False, save_path=heat_path)

    scatter_path = figure_folder_path / 'scatter_{}.pdf'.format(bulk_file_ID)
    plot_scatter(mean_df, bulk_df, sc_ref_df, gene_dict,
                 save_path=scatter_path)

    print('\nAll done! :-)\n')
