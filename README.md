## Welcome to `cellanneal` - The user-friendly application for deconvolving omics data sets.

<p align="center" width="100%">
    <img width="17%" src="https://github.com/LiBuchauer/cellanneal/blob/master/img/logo.png">
</p>

`cellanneal` is an application for deconvolving biological mixture data into constituting cell types. It comes both as a python package which includes a command line interface (CLI) and as a graphical software (graphical user interface, GUI) with the entire application bundled into a single executable. The python package with CLI can be downloaded from this repository; the graphical version is available for Microsoft Windows and macOS and can be downloaded from [zenodo](https://doi.org/10.5281/zenodo.5559545). If `cellanneal` is useful for your research, please cite [the preprint](https://arxiv.org/abs/2110.08209).



| [Download `cellanneal` graphical software for  Windows](https://zenodo.org/record/5559545/files/cellanneal_windows.zip?download=1) |
| --- |

| [Download `cellanneal` graphical software for  macOS](https://zenodo.org/record/5559545/files/cellanneal_macos.zip?download=1) |
| --- |

| IMPORTANT: The graphical software has a startup time of up to one minute. |
| --- |

| IMPORTANT: On macOS, when opening the graphical software for the first time,  you must do so via right-click --> "Open" and then choose "Open" in the emerging dialogue, [see below for more info](#2b-installing-the-gui). |
| --- |

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
The python package comes with a set of functions which can be included in python workflows, scripts and notebooks as well as with a command-line entry point to `cellanneal`. The required code can be downloaded from this repository. The graphical software is available for Microsoft Windows and macOS operating systems and can be downloaded from [zenodo](https://doi.org/10.5281/zenodo.5559545).

#### 2a. Installing the python package and CLI
Clone this code repository or download the zipped version and unpack it into a location of choice.  

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

Installing the graphical software is as simple as downloading the correct version for your operating system [zenodo](https://doi.org/10.5281/zenodo.5559545) and unzipping the content. The archive contains an executable file (the `cellanneal` application), an example mixture data file and an example signature data file.  **Please note that the GUI has an initial start-up time of up to one minute**.

* **macOS**: Under most security settings, macOS does not allow to open software from unidentified developers via double-click. To circumvent this, right-click (secondary click) onto the cellanneal executable and choose "Open" at the top of the emerging context menu. This primes a dialogue in which you can then press "Open". Subsequently, the software will also be accessible via double-click.

* **Windows**: Antiviral software may inhibit the launching of the software - it may be necessary to set an exception or click "Allow" when asked whether to procede.

For more information on how to use the GUI, see [Using the graphical software](#5c-using-the-graphical-software).

***

### 3. Requirements for input data
`cellanneal` accepts text files (\*.csv and \*.txt) as well as excel files (\*.xlsx and \*.xls) as inputs for both mixture and signature data provided that they are formatted correctly. Specifically, gene names need to appear in the first column for both mixture and signature data files, and sample names (for mixture data file) or cell type names (for signature data file) need to appear in the first row. Example data files can be found in this repository in the [example directory](https://github.com/LiBuchauer/cellanneal/tree/master/examples). The top of an exemplary mixture.csv file may look like this ![mixture csv file example](/img/data_example_csv.png) and the top of an exemplary signature.xlsx file looks like this ![signature xlsx file example](/img/data_example_excel.png)   

Further important points regarding the input data:

* It is not required that mixture and signature data sets contain exactly the same genes, or that these genes appear in the same order (or in alphabetical order).
* Please do not logarithmise the input data before passing it into `cellanneal`.  
* Normalisation of mixture data: it is not required that the individual sample columns are normalised to a specific sum value; the normalisation will not affect the outcome.
* Normalisation of signature data: normalising the individual cell type columns to the same count sum or not will lead to different results and whether you wish to normalise or not may depend on your biological question and available data. Specifically, if you do normalise all cell types to the same count sum, the output of `cellanneal` will tell you which fraction of the overall RNA was contributed by each cell type. This will not take into account size differences between cell types. In a toy example, if you analyse a mixture of one cell of type A and one cell of type B, where cell A at the base had ten times more RNA than cell B, after normalisation you will obtain the result that 10/11=91% of the RNA stem from type A. If all your reference data stems from the same data set, and you think that the average count sum of cells of a given type is a good proxy for their size, you may instead not normalise the values (in this case, the sum of all counts for cell type A would be 10 times higher than for B). Then, `cellanneal`’s output can be interpreted as cell fractions, i.e. the above example would return 50% type A and 50% type B. Following this concept, you can also think about normalising your cell types to different sum values based on known or estimated size factors.

***

### 4. Parameters
`cellanneal` allows the user to set four parameters. The first three govern the set of genes underlying the deconvolution process for each sample; the fourth parameter (iteration number) specifies for how long to run the optimisation process. Each parameter is discussed below.

* **Minimum expression in mixture** (`bulk_min`) - minimum required expression level in the mixture sample (where total expression is normalised to sum up to 1) for a gene to be considered, `default=1e-5`. Allowed values are in the range `[0, 1)` but must be smaller than the maximum allowed expression. This parameter allows to exclude lowly expressed and potentially noisy genes.

* **Maximum expression in mixture** (`bulk_max`) - maximum allowed expression level in the mixture sample (where total expression is normalised to sum up to 1) for a gene to be considered, `default=0.01`. Allowed values are in the range `(0, 1]` but must be larger than the minimum allowed expression. This parameter allows to exclude very highly expressed, potential contaminant, genes.

* **Minimum scaled dispersion** (`min_disp`) - minimum scaled dispersion (variance/mean) over cell types for a gene to be considered, `default=0.5`. The value indicates the number of standard deviations which the dispersion of a specific gene lies above or below the mean when compared to genes of similar expression. All numerical values are allowed, but reasonable values for most cases lie between 0 and 1 as this parameter is used to select for genes which vary across cell types in the signature file while still keeping a broad gene base for robust deconvolution.

* **Maximum number of iterations** (`maxiter`) - the maximum number of iterations through the logical chain of the underlying optimisation algorithm, [scipy’s dual annealing](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.dual_annealing.html). `default=1000`, after which typical problems have converged. Problems with a very high number of celltypes may require a higher number of iterations.

***

### 5. Using `cellanneal`
`cellanneal` can be used as part of a `python` workflow or individually via the command line or the graphical software. All three use cases are explained below.

#### 5a. Using the python package
The python package provides functions for the three main steps of a deconvolution analysis with `cellanneal`: identification of a gene set for deconvolution, deconvolution using simulated annealing, and plotting the results. A [quick start workflow](https://github.com/LiBuchauer/cellanneal/blob/master/examples/cellanneal_quickstart.ipynb) is available in the examples folder.  

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

After download, the graphical user interface can be opened by double-clicking the executable. A console directly opens up; the graphical user interface follows with a delay of up to one minute as the bundled python packages have to be unpacked into a temporary directory. Please be patient and do not close the console. Once started, the interface looks like this: ![cellanneal GUI](/img/cellanneal_gui.png)  

The user can now select mixture data, signature data and an output folder from the file system using the three `Browse file system` buttons in the upper half of the interface. Optionally, the four parameters (see also section [Parameters](#4-parameters)) can be changed via the `Change parameters` button which opens a separate window for entering parameter values. Parameter value defaults can be restored by clicking on the `Reset to default values` button.  

Finally, a deconvolution run is started by pressing the button `run cellanneal` at the bottom of the interface. Once running, the interface becomes unresponsive until the process finishes. While  `cellanneal` is running, detailed progress updates are  printed into  the accompanying console. When the run has finished, all results can be found in a directory labelled with the name of the mixture file and a timestamp inside the user-defined output folder. For further information on the output created by `cellanneal`, see section [Output](#6-cellanneal-output).

In order to shut down the application, the console window needs to be closed.

***

### 6. `cellanneal` output
`cellanneal` runs which were started from either the command line or the graphical user interface create a timestamped directory containing  three folders with tabular results and figures into the user-specifed output folder. Their contents are discussed below. Additionally, a text file containing the names of mixture and signature files and the parameters of the run is stored at the top level of the results folder.

#### 6a. Folder "deconvolution results"
This folder contains a CSV file with the main result of `cellanneal`: the fractional composition of each mixture in terms of cell types. Cell type names are shown in the first row; mixture sample names in the first column. Each numerical value in the table indicates the fraction the corresponding cell type occupies in the corresponding sample.

#### 6b. Folder "figures"
If the input mixture file contains up to 100 samples, a standard `cellanneal` run produces four figures:
* A figure with one pie chart per sample with each part of the pie representing the size of a cell type fraction.
* A heatmap in which mixture samples run across the horizontal axis and cell types along the vertical one, each coloured square indicating the corresponding cell type fraction.
* A second heatmap, similar to the first one, but showing log10(cell type fractions) instead in order to display small cell populations more clearly.
* A figure with one scatter plot per sample. In each scatter plot, each dot represents one gene and a dot's location is determined by its expression in the real mixture (x-axis) and its expression in the optimal computational mixture (i.e. the `cellanneal`result, y-axis). This figure helps judge how well `cellanneal` was able to approximate the real mixture sample by producing a computational mixture of the supplied cell types.

#### 6c. Folder "genewise comparison"
This folder contains one CSV file per mixture sample in the input data. Based on the deconvolution gene set for each sample , the file shows the normalised gene-wise expression in the experimental mixture (user input) in the first column and the corresponding expression in the optimal computational mixture in the second. The third column gives the ratio between the two (experimental/computational); the fourth the  logarithm of this fold change. The purpose of this file is to allow to search for genes with particularly high discrepancies between experimental and computational mixtures. Such genes may be of biological or medical interest: as an example, if the signature data stemmed from healthy people, but the mixture file from a pathology,  genes with high fold change between experiment and deconvolution result may have implications in the disease.

***

### 7. Frequently Asked Questions
* How small/large should my gene sets be?  
`cellanneal` draws its robustness and ability to identify small cell populations from its permissive gene selection strategy. Ideally, we want to include every gene which we believe to carry information (i.e. which is expressed above the noise level and has at least some meaningful variability between cell types). Very large gene sets (~5000 genes and above) may lead to long runtimes, but if the data quality allows, the user should aim to have gene set sizes upwards of 2000 or even 3000 genes.

* Why can my data not be imported?  
Please make sure that your data is formatted as described in the [Input data section](#3-requirements-for-input-data). Common pitfalls include:
    * excel files downloaded from publications contain a title row (e. g. "Supplementary Table 2")
    * excel may have converted some of your gene names to dates (e.g. "MAR1", "SEPT9"...)
    * CSV files have an unequal number of columns in the first row (the row with the sample or cell type names) compared to subsequent rows because the first row looks like this `sample1, sample2, sample3` instead of `gene_name, sample1, sample2, sample3` as it should be (in the subsequent data rows, the first column contains the gene name).

* What happens to mitochondrial genes? I noticed they are not part of my output.  
Mitochondrial genes (gene names starting with "MT-", "Mt-", "mt-") are removed from the gene list on which deconvolution is based in the `cellanneal` workflow. This happens after genes are selected based on minimum and maximum expression thresholds.

* My run fails with error `"Error: Sample XXX could not be deconvolved.Possibly the gene set for this sample is too small. See online documentation for more info."` when using the python module even though the selected gene set is large.  
This can be due to gene name duplications in either the signature file or the mixture file. Try to run `mixture_df = mixture_df.groupby(mixture_df.index).sum()` and/or  `signature_df = signature_df.groupby(signature_df.index).sum()` prior to running `cellanneal.deconvolve()`

