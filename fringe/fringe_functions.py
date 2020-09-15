import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pylab
pylab.ion()
import pandas as pd
from anndata import AnnData
from sklearn.preprocessing import normalize
# functional requirements
from scipy.spatial.distance import euclidean, cosine, correlation, jensenshannon
from scipy.stats import spearmanr
import time

# personalized dual_annealing function
from .C_fun_dual_anneal import dual_annealing


def make_sc_ref(
        sc_anndata,
        groupby='cluster'):
    """ Given an anndata object with cluster annotations, averages the
    counts per cluster and returns a dataframe with cluster names as
    columns and genes as the index.

    Parameters
    ----------
    sc_anndata: AnnData object, contains single cell counts (not normalised, not logged),
        annotated with cell types
    groupby: string label, name of the anndata.obs column which holds the cell type
        annotations by which to group in this function
    """
    # import AnnData and scanpy if available
    try:
        import scanpy
        from anndata import AnnData
    except ModuleNotFoundError:
        return "Please install scanpy if you wish to use 'make_sc_ref()''."
    # check if the input is an anndata object
    if type(sc_anndata) is not AnnData:
        raise TypeError("The object you have passed to 'make_sc_ref()' is not an AnnData object.")

    # now, calculate the mean counts per gene and cluster
    sc_ref_df = _grouped_obs_mean(sc_anndata, groupby=groupby)

    return sc_ref_df


def make_high_var_gene_lists(
        sc_ref_df,
        bulk_df,
        n_top_genes=None,
        bulk_min_thr=1e-5,
        bulk_max_thr=0.01,
        remove_mito=True
            ):
    """ Finds highly variable genes across cell types and checks for expression
    thresholds within each bulk separately, returns a dictionary of lists were
    each list contains the identified gene for the bulk given as as they appear
    __alphabetically sorted__ by column names.
    If n_high_var_genes is given, this number of highly variable genes is returned. If
    it is None, the default parameters for flavor='seurat' are used and the length of
    the resulting gene list depends on availability."""
    # order bulk columns alphabetically
    bulk_df = bulk_df.sort_index(axis=1)
    # first, find the the most variable genes across cell types (!, not single cells)
    # using scanpy's highly variable genes function - to do this, the data needs to
    # be transferred into an AnnData object
    sc_means_ann = AnnData(X=sc_ref_df.values.T)
    sc_means_ann.obs_names = sc_ref_df.columns
    sc_means_ann.var_names = sc_ref_df.index

    # now, identify highly variable genes, but first normalize and log1p
    sc.pp.normalize_total(sc_means_ann, target_sum=1e4, inplace=True)
    sc.pp.log1p(sc_means_ann, copy=False)
    sc.pp.highly_variable_genes(sc_means_ann, flavor='seurat', n_top_genes=n_top_genes)
    high_var_genes = sc_means_ann.var_names[sc_means_ann.var['highly_variable']==True].tolist()
    print("{} highly variable genes identified".format(len(high_var_genes)))
    # now, for each bulk, we find the genes which comply with our expression thresholds,
    # and then keep only those highly variable genes which do
    # to ensure usage of correct gene list later one, store in dict
    gene_lists = {}
    for bulk in bulk_df.columns:
        min_max_genes = _find_genes(bulk_df[bulk], min_thr=bulk_min_thr, max_thr=bulk_max_thr, remove_mito=remove_mito)
        thr_highvar_genes = [x for x in min_max_genes if x in high_var_genes]
        gene_lists[bulk] = thr_highvar_genes

    return gene_lists

def _grouped_obs_mean(adata, groupby):
    """ Given an AnnData object with obs annotations in column 'groupby', groups
    cells and calculates mean expression values (raw counts) which are written
    into a dataframe with groupby categories as columns and genes as index.

    Source
    ------
    adapted from https://github.com/theislab/scanpy/issues/181, ivirshup
    """
    # group obs df
    grouped = adata.obs.groupby(groupby)

    # prepare pandas dataframe for mean gene values in groups
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    # for each group, get the data of member cells and average it,
    # then write to dataframe columns
    for group, idx in grouped.indices.items():
        X = adata[idx].X
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))

    return out


def _find_genes(
        series, # pandas series in which to find compliant genes
        min_thr=1e-5, # minimum required expression
        max_thr=0.01, # maximum allowed expression
        remove_mito=True # if True, remove mitochondrial genes (starting with "mt-" or "MT-")
               ):
    """ Based on a pandas series containing gene expression values (not log!) of some kind, finds
    for genes which are expressed above the minimum relative threshold, below the maximum relative
    threshold, and, if requested, are not mitochondrial genes (start with any of "MT-", "Mt-", "mt-").
    """
    # check if input is a series
    if type(series) is not pd.Series:
        raise TypeError("The object you have passed to 'find_genes()' is not a pandas Series.")

    # from original data, retain only genes which are expressed below the max_thr
    colsum = series.sum()
    subset_maxclear = (series < colsum*max_thr)
    genes_maxclear =  subset_maxclear[subset_maxclear == True].index.tolist()

    # ... and above the min threshold
    subset_minclear = (series > colsum*min_thr)
    genes_minclear =  subset_minclear[subset_minclear == True].index.tolist()

    # find the intersection of all lists and return these genes
    joint_genes = list(set.intersection(*map(set, [genes_minclear, genes_maxclear])))

    # if required, remove mitochondrial genes
    if remove_mito:
        import re
        print(len(joint_genes), 1)
        mt = re.compile('mt-', re.I) # to allow for upper and lower case
        joint_genes = [g for g in joint_genes if not bool(mt.match(g))]
        print(len(joint_genes), 2)


    return joint_genes

def grouped_obs_sem(adata, group_key, layer=None, gene_symbols=None):
    """ function from https://github.com/theislab/scanpy/issues/181, ivirshup
    with modifications (bugfixes)"""
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = gene_symbols
        bdata = adata[:, gene_symbols]
    else:
        new_idx = adata.var_names
        bdata = adata

    grouped = bdata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((bdata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=new_idx
    )

    for group, idx in grouped.indices.items():
        X = getX(bdata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)/np.sqrt(len(X)))
    return out


def grouped_obs_sum(adata, group_key, layer=None, gene_symbols=None):
    """ function from https://github.com/theislab/scanpy/issues/181, ivirshup
    with modifications (bugfixes)"""
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = gene_symbols
        bdata = adata[:, gene_symbols]
    else:
        new_idx = adata.var_names
        bdata = adata

    grouped = bdata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((bdata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=new_idx
    )

    for group, idx in grouped.indices.items():
        X = getX(bdata[idx])
        out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))
    return out


# deconvolution functions
# first, find a list of genes that passes expression, dynamic and specificity thresholds
def find_genes(list_of_dfs, # list of dataframes for whom to find joint genes
               min_thr=5e-6, # minimum required expression
               max_thr=0.01, # genes which make up more in at least one cell type/sample are removed
               ):
    # list to store filtered gene lists
    acceptable_gene_list = []
    # go through the data sets provided and find genes that comply with both thresholds
    for a_df in list_of_dfs:
        # from original data, retain only genes which are expressed below the max_thr (one offender enough for a no)
        a_colsum = a_df.sum(axis=0).values
        a_subset_maxclear = ((a_df > a_colsum*max_thr).sum(axis=1)) == 0 # determine whether all compositions are below the max frac
        a_genes_maxclear = [x for x in a_df.index if a_subset_maxclear.loc[x] == True] # retain only genes that are
        acceptable_gene_list.append(a_genes_maxclear)
        # ... and above the min threshold (one is enough for a yes)
        a_subset_maxclear = ((a_df > a_colsum*min_thr).sum(axis=1)) > 0 # determine whether all compositions are below the max frac
        a_genes_minclear = [x for x in a_df.index if a_subset_maxclear.loc[x] == True] # retain only genes that are
        acceptable_gene_list.append(a_genes_minclear)

    # find the intersection of all lists and return these genes
    joint_genes = list(set.intersection(*map(set, acceptable_gene_list)))

    return joint_genes


def _rankdata(a):
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
    sorter = np.argsort(arr, kind='quicksort')

    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)

    arr = arr[sorter]
    obs = np.r_[True, arr[1:] != arr[:-1]]
    dense = obs.cumsum()[inv]

    # cumulative counts of each unique value
    count = np.r_[np.nonzero(obs)[0], len(obs)]

    # average method
    return .5 * (count[dense] + count[dense - 1] + 1)


# function that takes parameters specifying distribution and returns distance
# of resulting mixed data to a single bulk composition
def dist_1mixture(params,  # collection of independent heights of discerete dist
                  comp_vec_ranked,  # the gene expression vector with which to compare the resulting mixture
                  sc_data,  # single_cell data from which to mix new samples
                  ):
    # create array of the contributions and normalise to sum 1
    conv = np.array(params)
    mixture = (conv/sum(conv)).reshape(1, len(conv))

    # using this mixture, compute the countsums of the mixture from sc data
    mixed_counts = np.dot(mixture, sc_data.T).T

    # making this data compositional (normalizing it to match the bulk) is not
    # necessary because we will only compare ranks which are not affected by this

    # calculate Spearman correlation. The bulk vector has already been ranked,
    # so we only need to rank the newly mixed vector here.
    mixed_compositional_ranked = _rankdata(mixed_counts)
    # calculate Pearson correlation on ranked data
    dist = 1 - np.corrcoef(comp_vec_ranked, mixed_compositional_ranked, rowvar=False)[0][1]

    return dist


# define a function that returns the composition given the parameters of the distribution
def return_mixture(params):
    # create array of the contributions and normalise to sum 1
    conv = np.array(params)
    mixture = conv/sum(conv)
    return mixture


# function to select genes according to given threshold and deconvolve the resulting mixture
def deconvolve(sc_ref_df,
               bulk_df,
               maxiter,
               gene_dict,
               no_local_search
               ):
    # for each of the bulks, subset bulk and single-cell data according to the
    # gene list, retain only raw data matrices after this (df --> np.array)
    sc_list = []
    bulk_comp_list = []  # compositional version of subset bulk data
    bulk_ranked_list = []  # ranked version of subset bulk data

    for b, bulk in enumerate(bulk_df.columns):

        # first, subset and rank bulk data
        bulk_sub = bulk_df[bulk].loc[gene_dict[bulk]].values
        bulk_comp_list.append(bulk_sub)
        bulk_ranked = _rankdata(bulk_sub)
        bulk_ranked_list.append(bulk_ranked)

        # next, subset sc data
        sc_sub = sc_ref_df.loc[gene_dict[bulk]].values
        sc_list.append(sc_sub)

    # go through all mixtures and deconvolve them separately
    mixt_results_list = []
    for i, mixt in enumerate(bulk_df.columns):
        res = dual_annealing(dist_1mixture, bounds=[[0,1] for x in range(len(sc_ref_df.columns))],
                             maxiter=maxiter, args=[bulk_ranked_list[i], sc_list[i]], no_local_search=no_local_search)
        mixt_results_list.append(res)

    # grab the results, write them into a dataframe and return it
    mixture_list = []
    spears = []
    pears = []
    for i, res in enumerate(mixt_results_list):
        mixture = return_mixture(res.x)
        mixture_list.append(mixture)
        # calculate final spearson correlations
        mixed_counts = np.dot(mixture, sc_list[i].T).T
        mixed_compositional = mixed_counts / mixed_counts.sum()
        mixed_ranked = _rankdata(mixed_counts)
        spears.append(1-correlation(mixed_ranked, bulk_ranked_list[i]))
        pears.append(1-correlation(mixed_compositional, bulk_comp_list[i]))

    data_out = np.hstack((np.array(mixture_list), np.array([spears, pears]).T))
    cols_out = sc_ref_df.columns.tolist() + ['rho_Spearman', 'rho_Pearson']

    all_mix_df = pd.DataFrame(data=data_out, columns=cols_out, index=bulk_df.columns)

    return all_mix_df
