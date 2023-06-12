import numpy as np
from pandas import DataFrame, cut, Series

# functional requirements
from scipy.spatial.distance import correlation

# personalized dual_annealing function
from .dual_annealing import dual_annealing

# we choose to ignore warnings at this stage because console output is
# part of the user experience - make sure to enable when developing
import warnings

warnings.filterwarnings("ignore")


def make_gene_dictionary(
    celltype_df,
    bulk_df,
    disp_min=0.5,
    bulk_min=1e-5,
    bulk_max=0.01,
    remove_mito=True,
):
    """Finds highly variable genes across cell types and checks for expression
    thresholds within each bulk separately, returns a dictionary of lists were
    each list contains the identified gene for the bulk given as as they appear
    __alphabetically sorted__ by column names.
    If n_high_var_genes is given, this number of highly variable genes is
    returned. If it is None, the default parameters for flavor='seurat' are
    used and the length of the resulting gene list depends on availability."""
    # order bulk columns alphabetically
    bulk_df = bulk_df.sort_index(axis=1)
    # first, find the most variable genes across cell types
    high_var_genes = find_high_var_genes(celltype_df, disp_min=disp_min)

    print(
        """{} highly variable genes identified in cell type
        reference.""".format(
            len(high_var_genes)
        )
    )

    # now, for each bulk, we find the genes which comply with our expression
    # thresholds, and then keep only those highly variable genes which do
    # to ensure usage of correct gene list later one, store in dict
    gene_dict = {}
    for bulk in bulk_df.columns:
        min_max_genes = find_thr_genes(
            bulk_df[bulk], min_thr=bulk_min, max_thr=bulk_max, remove_mito=remove_mito
        )
        thr_highvar_genes = [x for x in min_max_genes if x in high_var_genes]
        gene_dict[bulk] = thr_highvar_genes
        print(
            "\t{} of these are within thresholds for sample {}".format(
                len(thr_highvar_genes), bulk
            )
        )

    return gene_dict


def calc_gene_expression(mix_vec, bulk_vec, celltype_df, gene_list):
    """Given a mixing vector, a bulk expression vector and the list of genes
    for this sample (i.e. all data only for one of the bulk samples),
    calculates the mixed expression and the fold change compared to the
    experimentally observed bulk expression and returns a dataframe with
    these two values as well as the original bulk measurement."""

    # subset celltype_df and bulk_vec to the genes supplied in gene_list
    bulk_vec_sub = bulk_vec.loc[gene_list].values
    bulk_vec_sub_comp = bulk_vec_sub / bulk_vec_sub.sum()
    celltype_df_sub = celltype_df.loc[gene_list].values

    # calculate the mixed gene expression vector
    mixed_expression = np.dot(mix_vec, celltype_df_sub.T).T
    mixed_expression_comp = mixed_expression / mixed_expression.sum()

    # calculate fold change as experimental over mixed
    exp_over_mixed = bulk_vec_sub_comp / mixed_expression_comp
    log_exp_over_mixed = np.log10(exp_over_mixed)

    # combine all into dataframe
    data = np.vstack(
        (bulk_vec_sub_comp, mixed_expression_comp, exp_over_mixed, log_exp_over_mixed)
    ).T
    gene_comp_df = DataFrame(
        data=data,
        columns=[
            "experimental bulk",
            "cellanneal mixed bulk",
            "fold change, exp/mixed",
            "log10 fold change, exp/mixed",
        ],
        index=gene_list,
    )
    return gene_comp_df


def find_high_var_genes(
    celltype_df,
    disp_min=0.5,
):
    """Finds highly variable genes across cell types in celltype_df.
    The implementation follows scanpy's highly_variable_genes procedure for
    flavor 'Seurat'."""
    # normalize counts within each celltype to sum 1
    sc_ref_norm = celltype_df.div(celltype_df.sum(axis=0), axis=1)

    # calculate mean and variance of each gene across types
    mean = np.mean(sc_ref_norm, axis=1)
    mean_of_sq = np.multiply(sc_ref_norm, sc_ref_norm).mean(axis=1)
    var = mean_of_sq - mean**2
    # enforce R convention (unbiased estimator) for variance
    var *= np.shape(sc_ref_norm)[1] / (np.shape(sc_ref_norm)[1] - 1)
    # set entries equal to zero to small value to avoid div by 0 value
    mean[mean == 0] = 1e-12
    # caculate dispersion from var and mean
    dispersion = var / mean

    # log versions of mean and dispersion are needed
    dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)
    mean = np.log1p(mean)

    # collect in a dataframe
    df = DataFrame(index=sc_ref_norm.index)
    df["means"] = mean
    df["dispersions"] = dispersion

    # group into 20 bins
    df["mean_bin"] = cut(df["means"], bins=20)
    disp_grouped = df.groupby("mean_bin")["dispersions"]
    # mean and std of dispersion in each group
    disp_mean_bin = disp_grouped.mean()
    disp_std_bin = disp_grouped.std(ddof=1)

    # retrieve those genes that have nan std, these are the ones where
    # only a single gene fell in the bin and implicitly set them to have
    # a normalized disperion of 1
    one_gene_per_bin = disp_std_bin.isnull()
    disp_std_bin[one_gene_per_bin.values] = disp_mean_bin[
        one_gene_per_bin.values
    ].values
    disp_mean_bin[one_gene_per_bin.values] = 0

    # normalize dispersions with respect to mean and std within each
    # gene's expression bin
    df["dispersions_norm"] = (
        df["dispersions"].values - disp_mean_bin[df["mean_bin"].values].values
    ) / disp_std_bin[df["mean_bin"].values].values

    # check which genes pass dispersion and expression thresholds
    dispersion_norm = df["dispersions_norm"].values.astype("float32")
    dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
    gene_subset = dispersion_norm > disp_min

    # write to df
    df["highly_variable"] = gene_subset

    # get list of gene names
    high_var_genes = df.loc[df["highly_variable"]].index.tolist()

    return high_var_genes


def find_thr_genes(
    series,  # pandas series in which to find compliant genes
    min_thr=1e-5,  # minimum required expression
    max_thr=0.01,  # maximum allowed expression
    remove_mito=True  # if True, remove mitochondrial genes
    # (starting with "mt-" or "MT-")
):
    """Based on a pandas series containing gene expression values
    (not log!) of some kind, finds for genes which are expressed above the
    minimum relative threshold, below the maximum relative threshold, and,
    if requested, are not mitochondrial genes (start with any of "MT-", "Mt-",
    "mt-").
    """
    # check if input is a series
    if type(series) is not Series:
        raise TypeError(
            """The object you have passed to 'find_genes()' is not
                a pandas Series."""
        )

    # from original data, retain only genes which are expressed below max_thr
    colsum = series.sum()
    subset_maxclear = series < colsum * max_thr
    genes_maxclear = subset_maxclear[subset_maxclear == True].index.tolist()
    # ! above, keep '==' and dont change to 'is' - will not work

    # ... and above the min threshold
    subset_minclear = series > colsum * min_thr
    genes_minclear = subset_minclear[subset_minclear == True].index.tolist()

    # find the intersection of all lists and return these genes
    joint_genes = list(set.intersection(*map(set, [genes_minclear, genes_maxclear])))

    # if required, remove mitochondrial genes
    if remove_mito:
        import re

        mt = re.compile("mt-", re.I)  # to allow for upper and lower case
        joint_genes = [g for g in joint_genes if not bool(mt.match(g))]

    return joint_genes


# deconvolution functions
def rankdata(a):
    """
    Assign ranks to data, dealing with ties appropriately.
    Ranks begin at 1. The average of the ranks that would have been assigned to
    all the tied values is assigned to each value.

    This function is a subset of scipy.stats.rankdata, to be found here
    https://github.com/scipy/scipy/blob/v1.4.1/scipy/stats/stats.py .

    Parameters
    ----------
    a : array_like
        The array of values to be ranked. Must be 1D.

    Returns
    -------
    ranks : ndarray
         An array of length equal to the size of `a`, containing rank
         scores.

    """
    arr = np.ravel(np.asarray(a))
    sorter = np.argsort(arr, kind="quicksort")

    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)

    arr = arr[sorter]
    obs = np.r_[True, arr[1:] != arr[:-1]]
    dense = obs.cumsum()[inv]

    # cumulative counts of each unique value
    count = np.r_[np.nonzero(obs)[0], len(obs)]

    # average method
    return 0.5 * (count[dense] + count[dense - 1] + 1)


# function that takes parameters specifying distribution and returns distance
# of resulting mixed data to a single bulk composition
def calculate_distance(
    params,  # collection of independent heights of discrete dist
    comp_vec_ranked,  # the gene expression vector with which to compare the resulting mixture
    sc_data,  # single_cell data from which to mix new samples
):
    # create array of the contributions and normalise to sum 1
    conv = np.array(params)
    mixture = (conv / sum(conv)).reshape(1, len(conv))

    # using this mixture, compute the countsums of the mixture from sc data
    mixed_counts = np.dot(mixture, sc_data.T).T

    # making this data compositional (normalizing it to match the bulk) is not
    # necessary because we will only compare ranks which are not affected by this

    # calculate Spearman correlation. The bulk vector has already been ranked,
    # so we only need to rank the newly mixed vector here.
    mixed_compositional_ranked = rankdata(mixed_counts)
    # calculate Pearson correlation on ranked data
    dist = (
        1 - np.corrcoef(comp_vec_ranked, mixed_compositional_ranked, rowvar=False)[0][1]
    )

    return dist


# define a function that returns the composition given the parameters of the
# distribution
def return_mixture(params):
    # create array of the contributions and normalise to sum 1
    conv = np.array(params)
    mixture = conv / sum(conv)
    return mixture


# function to select genes according to given threshold and deconvolve the
# resulting mixture
def deconvolve(
    celltype_df,
    bulk_df,
    maxiter,
    gene_dict,
):
    # sort bulk df columns alphabetically to ensure consistency
    bulk_df = bulk_df.sort_index(axis=1)
    # for each of the bulks, subset bulk and single-cell data according to the
    # gene list, retain only raw data matrices after this (df --> np.array)
    sc_list = []
    bulk_comp_list = []  # compositional version of subset bulk data
    bulk_ranked_list = []  # ranked version of subset bulk data

    for b, bulk in enumerate(bulk_df.columns):

        # first, subset and rank bulk data
        bulk_sub = bulk_df[bulk].loc[gene_dict[bulk]].values
        bulk_comp_list.append(bulk_sub)
        bulk_ranked = rankdata(bulk_sub)
        bulk_ranked_list.append(bulk_ranked)

        # next, subset sc data
        sc_sub = celltype_df.loc[gene_dict[bulk]].values
        sc_list.append(sc_sub)

    # go through all mixtures and deconvolve them separately
    mixture_list = []
    # total number of samples for print message
    N_samples = len(bulk_df.columns)
    for i, mixt in enumerate(bulk_df.columns):
        print("Deconvolving sample {} of {} ({}) ...".format(i + 1, N_samples, mixt))
        try:
            res = dual_annealing(
                calculate_distance,
                bounds=[[0, 1] for x in range(len(celltype_df.columns))],
                maxiter=maxiter,
                args=[bulk_ranked_list[i], sc_list[i]],
                no_local_search=False,
            )
            mixture = return_mixture(res.x)
            mixture_list.append(mixture)
        except ValueError:
            print(
                "\nError: Sample {} could not be deconvolved.\nPossibly the gene set for this sample is too small.\nSee online documentation for more info.\n".format(
                    mixt
                )
            )
            mixture = np.empty(len(celltype_df.columns))
            mixture[:] = np.nan
            mixture_list.append(mixture)

    # grab the results, write them into a dataframe and return it
    spears = []
    pears = []
    for i, mixture in enumerate(mixture_list):
        # calculate final spearson correlations
        mixed_counts = np.dot(mixture, sc_list[i].T).T
        mixed_compositional = mixed_counts / mixed_counts.sum()
        mixed_ranked = rankdata(mixed_counts)
        spears.append(1 - correlation(mixed_ranked, bulk_ranked_list[i]))
        pears.append(1 - correlation(mixed_compositional, bulk_comp_list[i]))

    data_out = np.hstack((np.array(mixture_list), np.array([spears, pears]).T))
    cols_out = celltype_df.columns.tolist() + ["rho_Spearman", "rho_Pearson"]

    all_mix_df = DataFrame(data=data_out, columns=cols_out, index=bulk_df.columns)

    return all_mix_df
