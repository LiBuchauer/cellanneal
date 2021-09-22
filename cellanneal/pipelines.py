import pandas as pd
import time

from .general import make_gene_dictionary, deconvolve, calc_gene_expression
from .plots import (plot_pies_from_df, plot_mix_heatmap,
                    plot_repeats, plot_scatter)
from .stability import repeat_annealing

def cellanneal_pipe(
        celltype_data_path,  # path object!
        celltype_df,
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
    celltypes = celltype_df.columns.tolist()
    # produce lists of genes on which to base deconvolution
    print('2. Constructing gene sets ...')
    gene_dict = make_gene_dictionary(
                    celltype_df,
                    bulk_df,
                    disp_min=disp_min,
                    bulk_min=bulk_min,
                    bulk_max=bulk_max)

    """ 3) Run cellanneal. """
    print('\n3. Running cellanneal ...')
    all_mix_df = deconvolve(
                    celltype_df=celltype_df,
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
                            celltype_df=celltype_df,
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
    plot_scatter(all_mix_df, bulk_df, celltype_df, gene_dict,
                 save_path=scatter_path)

    print('\nAll done! :-)\n')


def repeatanneal_pipe(
        celltype_data_path,
        celltype_df,
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
    celltypes = celltype_df.columns.tolist()
    # produce lists of genes on which to base deconvolution
    print('2. Constructing base gene sets ...')
    gene_dict = make_gene_dictionary(
                    celltype_df,
                    bulk_df,
                    disp_min=disp_min,
                    bulk_min=bulk_min,
                    bulk_max=bulk_max)

    """ 3) Run cellanneal. """
    print('\n3. Running cellanneal {} times...'.format(N_repeat))
    parent_df = repeat_annealing(
                    celltype_df,
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
                            celltype_df=celltype_df,
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
    plot_scatter(mean_df, bulk_df, celltype_df, gene_dict,
                 save_path=scatter_path)

    print('\nAll done! :-)\n')
