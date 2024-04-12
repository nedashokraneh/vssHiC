# vssHiC

## Overview

```vssHiC``` is an extension of the [VSS](https://github.com/faezeh-bayat/VSS), a signal transformation approach used for variance stabilization of epigenomic signals to support 
Hi-C modality. Genome-wide chromatin conformation capture assay, Hi-C, allows studying chromatin 3D structural features
and their functional implications, but the variance instability of raw read counts hinders the visualization and downstream analysis of this data modality. ```vssHiC``` stabilizes 
the variance of signals across a dynamic range, makes heatmap visualization of contact maps more appealing, and improves the performance of subcompartment callers 
relying on Gaussian observed variables. ```vssHiC``` is compatible with a well-developed R class ```InteractionSet``` and file format ```.cool``` for Hi-C data. There are also 
wrapper functions implemented for three TAD callers to support ```InteractionSet``` object and output TAD domains as a ```GRanges``` object which can be easily visualized on top of the 
heatmap of Hi-C contact maps by the ```HiContacts``` package.

## Installation

You can install the development version of vssHiC from Github:

``` r
library(devtools)
devtools::install_github("nedashokraneh/vssHiC")
```

## Vignette

A vignette of the vssHiC manual is provided at [https://rpubs.com/nshokran/vssHiC](https://rpubs.com/nshokran/vssHiC).

## Citation

Neda Shokraneh Kenari, Faezeh Bayat, Maxwell Libbrecht. [VSS-Hi-C: Variance-stabilized signals for chromatin contacts](https://doi.org/10.1101/2021.10.19.465027). bioRxiv 2021.10.19.
