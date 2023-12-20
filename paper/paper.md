---
title: '`cellanneal`: A user-friendly deconvolution software for transcriptomics data'
tags:
  - bioinformatics
  - computational biology
  - mixture deconvolution
  - bulk deconvolution
  - transcriptomics
  - RNA sequencing
  - omics methods
  
authors:
  - name: Lisa Buchauer
    corresponding: true # (This is how to denote the corresponding author)
    orcid: 0000-0002-4722-8390
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Shalev Itzkovitz
    orcid: 0000-0003-0685-2522
    affiliation: 1
affiliations:
 - name: Department of Molecular Cell Biology, Weizmann Institute of Science, Rehovot, Israel
   index: 1
 - name: Department of Infectious Diseases and Respiratory Medicine, Charité-Universitätsmedizin Berlin, Berlin, Germany
   index: 2
date: 12 June 2023
bibliography: paper.bib
---


# Summary
Single-cell sequencing methods enable precise characterization of gene expression patterns in individual cells. However, they may provide inaccurate information about the cell type composition of samples, as required preprocessing procedures such as tissue dissociation or cell sorting affect viability of different cell types to varying extent [@erdmann2021likelihood]. Further, especially in the clinical context, single-cell sequencing of patient samples is currently not routinely applied because of high cost and required expertise, while bulk sequencing is more prevalent. 

For these reasons, computational deconvolution methods are gaining popularity in basic and clinical research. Computational deconvolution approaches infer the cell type proportions constituting a given bulk RNA sample based on separately obtained cell type reference data. Several computational deconvolution methods have been developed in the last decade and have contributed to our understanding of tissue composition [@cobos2020benchmarking; @sturm2019benchmarking]. Generally, during deconvolution, the computational mixture is constructed from a set of cell type fractions and reference gene expression vectors for each of the participating cell types, most commonly derived from single-cell data. The cell type fractions are then iteratively changed until agreement between the _in silico_ gene expression vector and the observed bulk sample gene expression vector is optimal by a measure of choice. Here, published methods rely almost exclusively on minimizing the sum of squared residuals between bulk and computationally mixed vectors. Algorithms for such optimization problems are readily available and include variants of least squares regression (e.g. weighted least squares regression [@racle2017simultaneous], non-negative least squares regression [@wang2019bulk; @jew2020accurate] or least trimmed squares [@hao2019fast]) and support vector regression [@newman2015robust; @newman2019determining]. 

However, least squares-based optimization is faced with a particular challenge in bulk RNAseq deconvolution because of the highly skewed nature of mRNA copy number distributions, ranging from less than 1 to more than 10,000 average mRNA copies per cell [@li2016comprehensive; @schwanhausser2011global]. In such settings, optimization results may be strongly influenced by few highly expressed genes and are thus not robust to noise or platform effects influencing the readout of these genes. Support vector regression based models like CIBERSORT [@newman2015robust] perform gene feature selection out of a user-defined signature gene list, the contents of which can strongly affect the cell proportion estimates. Overall, identifying the right genes for deconvolution becomes a task in itself [@aliee2021autogenes]. As a result, deconvolution methods may yield inferred mixed gene expression vectors that do not correlate well with measured bulk gene expression.

Here, we introduce `cellanneal`, a python-based software for deconvolving bulk RNA sequencing data. `cellanneal` relies on the optimization of Spearman's rank correlation coefficient between experimental and computational mixture gene expression vectors using simulated annealing. Transforming gene expression values into ranks prior to optimization allows genes of different expression magnitudes to contribute similarly to deconvolution; further, `cellanneal` employs a permissive gene selection procedure that includes as many informative genes as possible. Together, these approaches limit the influence of highly expressed genes on the one hand and reduce dependency on specific gene list choices.  `cellanneal` can be used as a python package or via a command line interface, but importantly also provides a simple graphical user interface which is distributed as a single executable file for user convenience.

# Statement of need

Making sense of bulk RNA sequencing datasets often requires analysis of the cell type composition of the samples. This is particularly relevant in clinical samples that analyze the transcriptome of tissues or tumors which consist of epithelial, stromal and immune cell types. In parallel, publicly available single-cell data sets enable precise characterization of the expression signature of multiple individual cell types. However, software tools for computational bulk deconvolution are often slow, non-robust and not easy to use. Some existing methods address the aspect of user-friendliness by providing graphical web interfaces, but submitting sensitive medical data to an external web server is not always compatible with privacy legislation.

To address these challenges, we have developed `cellanneal`, a deconvolution approach that uses Spearman’s rank correlation coefficient between synthetic and bulk gene expression vectors as the optimization procedure’s objective function. Because this correlation measure is calculated from ranks rather than absolute data values, each gene influences the optimization result to a similar extent. Users are encouraged to include as many informative genes as possible in the analysis. `cellanneal` optimizes cell type fractions by simulated annealing, a flexible, rapid and robust algorithm for global optimization [@kirkpatrick1983optimization; @2020SciPy-NMeth]. `cellanneal` can be used as a python package, via its command line interface or via a user-friendly graphical software which runs locally. Its typical processing time for one mixture sample is below one minute on a desktop machine (MacBook Pro 2020, 2.3 GHz Quad-Core Intel Core i7, 16 GB RAM).

# Availability and Features

The python package and command line interface are available at [https://github.com/LiBuchauer/cellanneal](https://github.com/LiBuchauer/cellanneal) and can be installed using `pip`. The graphical software for Microsoft Windows and MacOS can be downloaded at [http://shalevlab.weizmann.ac.il/resources](http://shalevlab.weizmann.ac.il/resources) and does not require installation. Instructions for installation and use as well as general documentation is available at [https://github.com/LiBuchauer/cellanneal](https://github.com/LiBuchauer/cellanneal).

The python package provides functions for the three main steps of a deconvolution analysis with `cellanneal`: identification of a gene set for deconvolution, deconvolution using simulated annealing, and plotting the results. A quick start workflow is available as part of the documentation. For the command line interface and the graphical user interface, these three steps are combined into one call (click).

`cellanneal` runs which were started from either the command line or the graphical user interface produce a collection of result files including tabular deconvolution results (cell type fractions for each sample) and figures illustrating these cell type distributions. Further, `cellanneal` computes and stores the gene-wise fold change between the observed bulk expression and the estimated expression based on the inferred cell type composition. This enables identifying genes for which expression may be specifically induced or inhibited in the bulk sample compared to the single cell reference. Such genes may be of biological or medical interest.  
Figures produced by `cellanneal` include a heatmap showing sample compositions (\autoref{fig:heatmap}), pie charts showing sample compositions (\autoref{fig:pie}), and scatter plots showing correlation between experimental bulk gene expression values and their `cellanneal`-derived counterparts from the best identified computational mixture (\autoref{fig:scatter}). The examples presented in this manuscript use data from [@massalha2020liver].


![A heatmap produced by `cellanneal`. Constituting cell types are on the y-axis, deconvolved bulk sample names on the x-axis. The colour scale shows the fractional presence of cell type in each bulk.\label{fig:heatmap}](heatmap_example.png){ width=50% }

![Pie charts produced by `cellanneal`. Each pie corresponds to one deconvolved bulk sample from the input data.\label{fig:pie}](piechart_example.png)

![Gene correlation scatter plots produced by `cellanneal` Each panel corresponds to one deconvolved bulk sample from the input data. Each dot represents a gene used during deconvolution. The x-axis shows the experimentally measured expression of each gene after normalizing so that the total count sum is 1. The y-axis shows the normalized expression of each gene in the best identified synthetic bulk mixed from cell type signature data .\label{fig:scatter}](scatter_example.png)


`cellanneal` relies on the python packages `scipy` [@2020SciPy-NMeth], `numpy` [@harris2020array], `pandas` [@reback2020pandas], `seaborn` [@Waskom2021] and `matplotlib` [@Hunter:2007].

# Citations
Examples of published research projects using cellanneal include [@egozi2023single] and [@berkova2022terminal].

# Acknowledgements

We thank all members of the Itzkovitz lab for valuable feedback and for testing the software.

L.B. was supported by the European Molecular Biology Organization under EMBO Long-Term Fellowship ALTF 724-2019. S.I. is supported by the Wolfson Family Charitable Trust, the Edmond de Rothschild Foundations, the Fannie Sherr Fund, the Helen and Martin Kimmel Institute for Stem Cell Research grant, the Minerva Stiftung grant, the Israel Science Foundation grant No. 1486/16, the Broad Institute‐Israel Science Foundation grant No. 2615/18, the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation program grant No. 768956, the Chan Zuckerberg Initiative grant No. CZF2019‐002434, the Bert L. and N. Kuggie Vallee Foundation and the Howard Hughes Medical Institute (HHMI) international research scholar award. 

# References
