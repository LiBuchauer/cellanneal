"""cellanneal - deconvolution of bulk RNAseq data using simulated annealing"""

__version__ = "0.1.0"
__author__ = "Lisa Buchauer <lisa.buchauer@posteo.de>"

# functions which become available for import
from .general import make_gene_dictionary, return_mixture, deconvolve
from .plots import plot_pies, plot_mix_heatmap, plot_mix_heatmap_log, plot_scatter
from .pipelines import cellanneal_pipe, run_cellanneal
