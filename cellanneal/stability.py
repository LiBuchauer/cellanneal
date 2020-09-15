import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from random import sample
# functional requirements
from .general import deconvolve

"""

1)
Functions for testing reproducibility of deconvolution results vis a vis optimisation
stochasticity and gene set differences

"""

def repeat_annealing(
        sc_ref_df,
        bulk_df,
        gene_dict,
        no_local_search,
        N_repeat=10,
        sample_fraction=0.9,
        maxiter=1000
        ):
    """Function that allows to run the deconvolution process N_repeat
    times, storing and merging the results into a long-form dataframe
    for inspection regarding variability and instabilities."""
    # sort bulk df columns alphabetically to ensure consistency
    bulk_df = bulk_df.sort_index(axis=1)

    # list for collecting outputs of all runs
    df_list = []

    # repeat deconvolution n times and collect results
    for i in range(N_repeat):
        print("\nRepeat number {}".format(i+1))
        # now, draw new gene dictionary
        gene_dict_new = {}
        for key in list(gene_dict.keys()):
            number_of_genes = int(sample_fraction * len(gene_dict[key]))
            gene_dict_new[key] = sample(gene_dict[key], number_of_genes)
        
        # deconvolve with the new gene dict
        all_mix_df = all_mix_df = deconvolve(sc_ref_df=sc_ref_df,
                              bulk_df=bulk_df,
                              gene_dict=gene_dict_new,
                              maxiter=maxiter, no_local_search=no_local_search)
        try:
            b_mix_df = all_mix_df.drop(['rho_Spearman', 'rho_Pearson'], axis=1)
        except KeyError:
            b_mix_df = all_mix_df
        b_mix_df['bulk'] = b_mix_df.index
        long_df = pd.melt(b_mix_df, id_vars=['bulk'], var_name='celltype', value_name='fraction')
        long_df['run'] = i
        df_list.append(long_df)

    # finally, merge the dfs in the list, write to file and return
    master_df = pd.concat(df_list, axis=0)

    return master_df
