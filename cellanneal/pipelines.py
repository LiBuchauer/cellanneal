import time

from .general import make_gene_dictionary, deconvolve, calc_gene_expression
from .plots import plot_pies, plot_mix_heatmap, plot_mix_heatmap_log, plot_scatter


def cellanneal_pipe(
    celltype_data_path,  # path object!
    celltype_df,
    bulk_data_path,  # path object!
    bulk_df,
    disp_min,
    bulk_min,
    bulk_max,
    maxiter,
    output_path,
):  # path object!
    """Serves as entrypoint into cellanneal pipeline for both gui and cli
    once all data and parameters have been collected."""

    """ 2) Identify highly variable genes and genes that pass the thresholds
    for each bulk. """
    # extract names
    bulk_names = bulk_df.columns.tolist()
    celltypes = celltype_df.columns.tolist()

    # check here for uniqueness
    if len(bulk_names) > len(set(bulk_names)):
        print(
            "Error: The names of the mixtures are not unique. Please give a unique name to each mixture."
        )
        return 0
    if len(celltypes) > len(set(celltypes)):
        print(
            "Error: The names of the cell types are not unique. Please give a unique name to each cell type."
        )
        return 0

    # produce lists of genes on which to base deconvolution
    print("\n+++ Constructing gene sets ... +++")
    gene_dict = make_gene_dictionary(
        celltype_df, bulk_df, disp_min=disp_min, bulk_min=bulk_min, bulk_max=bulk_max
    )

    """ 3) Run cellanneal. """
    print("\n+++ Running cellanneal ... +++")
    all_mix_df = deconvolve(
        celltype_df=celltype_df, bulk_df=bulk_df, maxiter=maxiter, gene_dict=gene_dict
    )

    """ 4) Write results to file."""
    print("\n+++ Writing results to file ... +++")

    # make top level folder for all results from this run
    # get timestamp for labelling
    timestamp = time.asctime().replace(" ", "_").replace(":", "-")
    bulk_file_name = bulk_data_path.name
    bulk_file_ID = bulk_file_name.split(".")[0]
    top_folder_name = "cellanneal_" + bulk_file_ID + "_" + timestamp
    top_folder_path = output_path / top_folder_name
    top_folder_path.mkdir(parents=True, exist_ok=True)

    # make subfolders for deconv, gen expr and figures
    deconv_folder_path = top_folder_path / "deconvolution_results"
    deconv_folder_path.mkdir(parents=True, exist_ok=True)
    figure_folder_path = top_folder_path / "figures"
    figure_folder_path.mkdir(parents=True, exist_ok=True)
    genexpr_folder_path = top_folder_path / "genewise_comparison"
    genexpr_folder_path.mkdir(parents=True, exist_ok=True)

    # first, write the mix matrix to csv
    deconv_name = "deconvolution_" + bulk_file_ID + ".csv"
    result_path = deconv_folder_path / deconv_name
    all_mix_df.sort_index(axis=0, inplace=True)
    all_mix_df.to_csv(result_path, header=True, index=True, sep=",")

    # next, write the actual and estimated gene expression to file,
    # this has to be done per sample as the genes are sample specific
    # a mix_df version without correlation entries is needed
    all_mix_df_no_corr = all_mix_df[celltypes]
    for sample_name in bulk_names:
        gene_comp_df = calc_gene_expression(
            mix_vec=all_mix_df_no_corr.loc[sample_name],
            bulk_vec=bulk_df[sample_name],
            celltype_df=celltype_df,
            gene_list=gene_dict[sample_name],
        )
        # construct export path for this sample
        sample_gene_name = "expression_" + bulk_file_ID + "_" + sample_name + ".csv"
        sample_gene_path = genexpr_folder_path / sample_gene_name
        gene_comp_df.sort_index(axis=0, inplace=True)
        gene_comp_df.to_csv(sample_gene_path, header=True, index=True, sep=",")

    # write a text file with all parameters
    param_file_path = top_folder_path / "parameters_{}.txt".format(timestamp)
    with open(param_file_path, "a") as file:
        file.write(
            "parameters and data used for this cellanneal run ({})\n\n".format(
                timestamp
            )
        )
        file.write("mixture data: {}\n".format(bulk_data_path))
        file.write("signature data: {}\n".format(celltype_data_path))
        file.write("minimum expression in mixture: {}\n".format(bulk_min))
        file.write("maximum expression in mixture: {}\n".format(bulk_max))
        file.write("minimum dispersion: {}\n".format(disp_min))
        file.write("maximum number of iterations: {}\n".format(maxiter))

    """ 5) Produce plots and save to folder"""
    # we only want figures if there are less than 100 samples
    if len(bulk_names) > 100:
        print(
            "\nInfo: cellanneal does not produce figures for runs with more than 100 samples. If you would like cellanneal to produce figures, consider splitting your data into several input files with less than 100 mixtures each."
        )
    # plot results
    else:
        print('\n+++ Storing figures in folder "figures" ... +++')
        try:
            pie_path = figure_folder_path / "pies_{}.pdf".format(bulk_file_ID)
            plot_pies(all_mix_df, save_path=pie_path)

            heat_path = figure_folder_path / "heat_{}.pdf".format(bulk_file_ID)
            plot_mix_heatmap(all_mix_df, rownorm=False, save_path=heat_path)

            heat_log_path = figure_folder_path / "heat_log10_{}.pdf".format(
                bulk_file_ID
            )
            plot_mix_heatmap_log(all_mix_df, rownorm=False, save_path=heat_log_path)

            scatter_path = figure_folder_path / "scatter_{}.pdf".format(bulk_file_ID)
            plot_scatter(
                all_mix_df, bulk_df, celltype_df, gene_dict, save_path=scatter_path
            )
        except:
            print("\nError: Plots could not be created.")

    print("\n+++ Finished. +++\n")


def run_cellanneal(
    celltype_df, bulk_df, disp_min, bulk_min, bulk_max, maxiter
):  # path object!
    """Combines gene set identification and deconvolution into a single
    function.

    Input:
    celltype_df  -  dataframe containing signature data
    bulk_df  -  daatframe containing mixture data
    disp_min  -  minimum scaled dispersion
    bulk_min   -  minimum expression in mixture data for genes
    bulk_max  -  maximum expression in mixture data for genes
    maxiter  -  maximum number of iterations for scipy's dual annealing

    Output:
    all_mix_df  -  a dataframe containing cell type fractions for each mixture"""

    # produce lists of genes on which to base deconvolution
    print("\n+++ Constructing gene sets ... +++")
    gene_dict = make_gene_dictionary(
        celltype_df, bulk_df, disp_min=disp_min, bulk_min=bulk_min, bulk_max=bulk_max
    )

    """ 3) Run cellanneal. """
    print("\n+++ Running cellanneal ... +++")
    all_mix_df = deconvolve(
        celltype_df=celltype_df, bulk_df=bulk_df, maxiter=maxiter, gene_dict=gene_dict
    )

    return all_mix_df
