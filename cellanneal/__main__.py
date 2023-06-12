"""This script defines the command line entry point for cell anneal and
implements a pipeline including identification of highly variable genes,
deconvolution of each bulk sample, ouput of deconvolution results and
production of a set of plots."""

import argparse
from pathlib import Path
from pandas import read_csv, read_excel
import time
import openpyxl  # for xlsx import
import xlrd  # for xls import

from .pipelines import cellanneal_pipe


def init_parser(parser):
    """Initialize parser arguments."""
    parser.add_argument(
        "bulk_data_path",
        type=str,
        help=(
            """Path to mixture data file; .csv, .txt, .xlsx or .xls format
        with sample names as columns and genes as rows."""
        ),
    )

    parser.add_argument(
        "celltype_data_path",
        type=str,
        help=(
            """Path to signature data file; .csv, .txt, .xlsx or .xls format
        with sample names as columns and genes as rows."""
        ),
    )

    parser.add_argument(
        "output_path",
        type=str,
        help=("""Path to folder in which to store the results. """),
    )

    parser.add_argument(
        "--bulk_min",
        type=float,
        default=1e-5,
        help=(
            """Minimum relative expression in the mixture sample for a gene
            to be considered."""
        ),
    )

    parser.add_argument(
        "--bulk_max",
        type=float,
        default=0.01,
        help=(
            """Maximum relative expression in the mixture sample for a gene
            to be considered."""
        ),
    )

    parser.add_argument(
        "--disp_min",
        type=float,
        default=0.5,
        help=(
            """Minimum scaled dispersion (var/mean) relative to other
            genes of similar expression for a gene to be considered."""
        ),
    )

    parser.add_argument(
        "--maxiter",
        type=int,
        default=1000,
        help=("""Maximum number of iterations for scipy's dual_annealing."""),
    )

    return parser


def main():
    """cellanneal. User-friendly deconvolution of RNA-Seq mixture data.

    Deconvovles mixture RNA-Seq data into cell type populations based
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
        description=("cellanneal deconvolves bulk RNA-Seq data.")
    )
    args = init_parser(my_parser).parse_args()

    # grab the individual inputs for further use
    bulk_data_path = Path(args.bulk_data_path)
    celltype_data_path = Path(args.celltype_data_path)
    output_path = Path(args.output_path)
    bulk_min = args.bulk_min
    bulk_max = args.bulk_max
    disp_min = args.disp_min
    maxiter = args.maxiter

    print("\n+++ Welcome to cellanneal! +++")
    print("{}\n".format(time.ctime()))

    """ 1) Import bulk and cell type data """
    print("\n+++ Importing mixture data ... +++ \n")
    try:
        bfile = bulk_data_path
        # depending on extension, use different import function
        if bfile.name.split(".")[-1] in ["csv", "txt"]:
            bulk_df = read_csv(bfile, index_col=0, sep=None)
        elif bfile.name.split(".")[-1] in ["xlsx"]:
            bulk_df = read_excel(bfile, index_col=0, engine="openpyxl")
        elif bfile.name.split(".")[-1] in ["xls"]:
            bulk_df = read_excel(bfile, index_col=0, engine="xlrd")
        else:
            raise ImportError
        # here, in order to make further course case insensitive,
        # change all gene names to uppercase only
        bulk_df.index = bulk_df.index.str.upper()
        # also, if there are duplicate genes, the are summed here
        bulk_df = bulk_df.groupby(bulk_df.index).sum()
        # finally, if there are nan's after import, set them to 0 to
        # avoid further issues
        bulk_df = bulk_df.fillna(0)
    except ValueError:
        print(
            """Your bulk data file could not be imported.
        Please check the documentation for format requirements
        and look at the example bulk data files.\n"""
        )
        print("+++ Aborted. +++")
        return 0

    print("\n+++ Importing signature data ... +++ \n")
    # import single cell based reference
    try:
        # depending on extension, use different import function
        cfile = celltype_data_path
        if cfile.name.split(".")[-1] in ["csv", "txt"]:
            celltype_df = read_csv(cfile, index_col=0, sep=None)
        elif cfile.name.split(".")[-1] in ["xlsx"]:
            celltype_df = read_excel(cfile, index_col=0, engine="openpyxl")
        elif cfile.name.split(".")[-1] in ["xls"]:
            celltype_df = read_excel(cfile, index_col=0, engine="xlrd")
        else:
            raise ImportError
        # here, in order to make further course case insensitive,
        # change all gene names to uppercase only
        celltype_df.index = celltype_df.index.str.upper()
        # also, if there are duplicate genes, the are summed here
        celltype_df = celltype_df.groupby(celltype_df.index).sum()
        # finally, if there are nan's after import, set them to 0 to
        # avoid further issues
        celltype_df = celltype_df.fillna(0)
    except ValueError:
        print(
            """Your celltype data file could not be imported.
        Please check the documentation for format requirements
        and look at the example celltype data files."""
        )
        print("+++ Aborted. +++")
        return 0

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
        output_path,
    )
