"""cellanneal - deconvolution of bulk RNAseq data using simulated annealing"""

__version__ = '0.1.0'
__author__ = 'Lisa Buchauer <lisa.buchauer@posteo.de>'

from .general import make_gene_dictionary, return_mixture, deconvolve
from .plots import (plot_pies_from_df, plot_mix_heatmap,
                          plot_1D_lines, plot_repeats, plot_scatter)
from .stability import repeat_annealing
