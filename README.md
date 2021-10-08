## Welcome to `cellanneal` - The user-friendly application for deconvolving omics data sets.

`cellanneal` is an application for deconvolving biological mixture data into constituting cell types. It comes both as a python package which includes a command line interface (CLI) and as a graphical software (graphical user interface, GUI) with the entire application bundled into a single executable. The python package with CLI can be downloaded from this repository; the graphical version is available for Microsoft Windows and MacOS and can be downloaded from the [Itzkovitz group website](http://shalevlab.weizmann.ac.il/resources/).

### Contents
1. [How does `cellanneal` work?](#1-how-does-cellanneal-work)  
2. [Installation](#2-installation)  
    a. [python package and command line interface](#2a-installing-the-python-package-and-cli)  
    b. [graphical software](#2b-installing-the-gui)  
3. [Requirements for input data files](#3-requirements-for-input-data)  
4. [Parameters](#4-parameters)  
5. [Using `cellanneal`](#5-using-cellanneal)  
    a. [python package](#5-using-the-python-package)  
    b. [command line interface](#5b-using-the-command-line-interface)  
    c. [graphical software](#5c-using-the-graphical-software)
6. [`cellanneal` output](#6-cellanneal-output)  
    a. [deconvolution results](#6a-folder-deconvolution-results)  
    b. [figures](#6b-folder-figures)  
    c. [genewise comparison](#6c-folder-genewise-comparison)  
7. [FAQs](#7-frequently-asked-questions)  


***

### 1. How does `cellanneal` work?
Given a gene-expression vector of a cellular mixture (for example derived from bulk RNA sequencing, the "mixture data") and gene-expression vectors characterising individual cell types (for example derived from clustered single-cell RNA sequencing data, the "signature data"), `cellanneal` provides an estimate of what fraction of each cell type is present in the bulk sample.  

During the deconvolution process, a computational mixture sample is constructed from a set of cell type fractions and the signature data. The resulting synthetic gene expression vector is compared to the gene expression vector of the real mixture by calculating Spearman's correlation coefficient between the two. Cell type fractions are then changed until this correlation is maximised using the optimisation algorithm simulated annealing as implemented in `scipy`'s `dual_annealing`. The cell type fractions associated with the highest Spearman correlation between the gene expression data of the experimental mixture (bulk sample gene expression) and the computational mixture are the `cellanneal` estimate for the mixture composition in terms of the cell types supplied in the signature file.

***

### 2. Installation
The python package comes with a set of functions which can be included in python workflows, scripts and notebooks as well as with a command-line entry point to `cellanneal`. The required code can be downloaded from this repository. The graphical software is available for Microsoft Windows and MacOS operating systems and can be downloaded from the [Itzkovitz group website](http://shalevlab.weizmann.ac.il/resources/).

#### 2a. Installing the python package and CLI
Clone this code repository or download the zipped version and unpack it into a location fo choice.  

Installing `cellanneal` into a virtual environment, for example via `anaconda`, is recommended. `cellanneal` has been tested with `python 3.7` and `python 3.8`.

It is recommended to install `cellanneal`'s dependencies first; if using `conda`:

    conda install numpy scipy matplotlib pandas seaborn xlrd openpyxl

If using pip:

    pip install numpy scipy matplotlib pandas seaborn xlrd openpyxl

Next, navigate into the `cellanneal` directory (`cellanneal-master`, the directory containing the file `setup.py`) on the command line. There, execute the command

    pip install .

That's it. Now you should be able to use `cellanneal` in your python projects via

```python
import cellanneal
```

and via the command line as

    cellanneal mixture_data.csv signature_data.csv output_folder

For more details on how to use it, see [Using `cellanneal`](#5-using-cellanneal).

#### 2b. Installing the GUI

Installing the graphical software is as simple as downloading it from the [Itzkovitz group website](http://shalevlab.weizmann.ac.il/resources/). **Please note that the GUI has an initial start-up time of up to one minute**. For more information on how to use the GUI, see [Using the graphical software](#5c-using-the-graphical-software).

***

### 3. Requirements for input data
`cellanneal` accepts text files (\*.csv and \*.txt) as well as excel files (\*.xlsx and \*.xls) as inputs for both mixture and signature data provided that they are formatted correctly. Specifically, gene names need to appear in the first column for both mixture and signature data files, and sample names (for mixture data file) or cell type names (for signature data file) need to appear in the first row. **ADD SCREENSHOTS OF CORRECT DATA HERE**. Example data files can be found in this repository in the directory [examples](https://github.com/LiBuchauer/cellanneal/tree/master/examples).
***

### 4. Parameters
`cellanneal` allows the user to set four parameters. The first three govern the set of genes underlying the deconvolution process for each sample; the fourth parameter (iteration number) specifies for how long to run the optimisation process. Each parameter is discussed below.

* **Minimum expression in mixture** (`bulk_min`) - minimum required expression level in the mixture sample (where total expression is normalised to sum up to 1) for a gene to be considered, `default=1e-5`. Allowed values are in the range `[0, 1)` but must be smaller than the maximum allowed expression. This parameter allows to exclude lowly expressed and potentially noisy genes.

* **Maximum expression in mixture** (`bulk_max`) - maximum allowed expression level in the mixture sample (where total expression is normalised to sum up to 1) for a gene to be considered, `default=0.01`. Allowed values are in the range `(0, 1]` but must be larger than the minimum allowed expression. This parameter allows to exclude very highly expressed, potential contaminant, genes.

* **Minimum scaled dispersion** (`min_disp`) - minimum scaled dispersion (variance/mean) over cell types for a gene to be considered, `default=0.5`. The value indicates the number of standard deviations which the dispersion of a specific gene lies above or below the mean when compared to genes of similar expression. All numerical values are allowed, but reasonable values for most cases lie between 0 and 3 as this parameter is used to select genes which vary strongly across cell types in the signature file.

* **Maximum number of iterations** (`maxiter`) - the maximum number of iterations through the logical chain of the underlying optimisation algorithm, [scipy’s dual annealing](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.dual_annealing.html). `default=1000`, after which typical problems have converged. Problems with a very high number of celltypes may require a higher number of iterations.

***

### 5. Using `cellanneal`
`cellanneal` can be used as part of a `python` workflow or individually via the command line or the graphical software. All three use cases are explained below.

#### 5a. Using the python package
The python package provides functions for the three main steps of a deconvolution anlysis with `cellanneal`: identification of a gene set for deconvolution, deconvolution using simulated annealing, and plotting the results. A [quick start workflow](https://github.com/LiBuchauer/cellanneal/blob/master/examples/cellanneal_quickstart.ipynb) is available in the examples folder.  

In order to use `cellanneal` in your `python` workflow, you need to import it:

```python
import cellanneal
```

As a first step, a gene set on which to base deconvolution has to be identified for each mixture sample. This step uses the parameters `bulk_min`, `bulk_max` and `disp_min` which are explained in the section [Parameters]((#4-parameters). The function `make_gene_dictionary` takes these inputs and produces a `dictionary` holding a gene list for each mixture sample:

```python
gene_dict = cellanneal.make_gene_dictionary(
                    signature_df,
                    mixture_df,
                    disp_min=0.5,
                    bulk_min=1e-5,
                    bulk_max=0.01)
```

Next, deconvolution is run and a `pandas.DataFrame` holding the results is returned:
```python
all_mix_df = cellanneal.deconvolve(
                signature_df,
                mixture_df,
                maxiter=1000,
                gene_dict=gene_dict)
```
Finally, four plotting options for deconvolution results are provided with `cellanneal` - pie charts, two heatmaps, and a scatter plot showing correlations between computational and real mixture samples.

```python
cellanneal.plot_pies(all_mix_df)
cellanneal.plot_mix_heatmap(all_mix_df)
cellanneal.plot_mix_heatmap_log(all_mix_df)
cellanneal.plot_scatter(all_mix_df, mixture_df, signature_df, gene_dict)
```

#### 5b. Using the command line interface
After installing the python package, a single command line command, `cellanneal`,
becomes available. Note that if you are using `conda` environments, this command will only be available inside the environment into which you installed it and you need to activate this environment via `conda activate my_env` before you can make calls to `cellanneal`.  

`cellanneal` requires three arguments,
* the path to the mixture data file (\*.csv, \*.txt, \*.xlsx, or \*.xls)  
* the path to the signature data file (\*.csv, \*.txt, \*.xlsx, or \*.xls)  
* the path to the folder in which the results are to be stored  

and allows the user to set four parameters,  
* `bulk_min`, the minimum required gene expression in the mixture  
* `bulk_max`, the maximum allowd gene expression in the mixture  
* `min_disp`, the minimum required scaled dispersion across cell types  
* `maxiter`, the maximum iteration number for `scipy`'s `dual_annealing`

resulting in the following call signature:
```
cellanneal [-h] [--bulk_min BULK_MIN] [--bulk_max BULK_MAX]
                [--disp_min DISP_MIN] [--maxiter MAXITER]
                bulk_data_path celltype_data_path output_path
```
Further information about each parameter can be found in section [Parameters](#4-parameters).


#### 5c. Using the graphical software

***

### 6. `cellanneal` output
Independent of whether it was started from the command line or the graphical user interface, each `cellanneal` run creates a timestamped directory containing  three folders with tabular results and figures into the user-specifed output folder. Their contents are discussed below. Additionally, a text file containing the names of mixture and signature files and the parameters of the run is stored at the top level of the results folder.

#### 6a. Folder "deconvolution results"

#### 6b. Folder "figures"

#### 6c. Folder "genewise comparison"

***

### 6. Frequently Asked Questions
Notes - include antivirus issue for windows  
mention how many genes one should aim for
mention what lower cut-off (min expression) makes sense based on the number of counts
