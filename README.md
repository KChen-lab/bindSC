
<image src="./doc/images/logo.png" width="400">
  
**bindSC** (**B**i-order  **IN**tegration of multi-omics **D**ata from **S**ingle **C**ell sequencing technologies) is an R package for single cell multi-omic integration analysis, developed and maintained by [Ken chen's lab](https://sites.google.com/view/kchenlab/Home) in MDACC. `bindSC` is developed to address the challenge of single-cell multi-omic data integration that consists of unpaired cells measured with unmatched features across modalities. Previous methods such as [Seurat](https://satijalab.org/seurat/), [Liger](https://github.com/MacoskoLab/liger), [Harmony](https://github.com/immunogenomics/harmony) did not work on this case unless match feature empricallcy. For example, integration of scRNA-seq and scATAC-seq data requires to calculate the gene/promoter activity by counting peaks in gene body, which always loses information. This strategy also did not work on integrating scRNA-seq and cytof data becasue gene pression and protein abundance level is not always correlated due to sparsity of scRNA-seq data or post translational modification. 

The core algorithm implemented in `bindSC` package is `BiCCA` (**Bi**-order **C**anonical **C**orrelation **A**nalysis), which utilizes a transition matrix **Z** (*M* features by *L* samples) to bridge observed **X** (*M* features by *K* cells) with **Y** (*N* features by *L* cells). Initialized from prior knowldge, the matrix **Z** is solved iteratively by maximizing correlation of pair (**X**, **Z**) and correlation of pair (**Y**, **Z**) simultaneously. Under estimated matrix **Z**, the cell/feature correspondence across modalities can be obtained by implementing standard `CCA` on pair (**X**, **Z**) and pair (**Y**, **Z**) respectively. 


Once multiple datasets are integrated, `bindSC` provides functionality for further data exploration, analysis, and visualization. User can: 


* Jointly defining cell types from multi-omic datasets
* Identifying comprehensive molecular multi-view of biological processes in cell type level.

Improvements and new features will be added on a regular basis, please contact jinzhuangdou198706@gmail.com or kchen3@mdanderson.org with any question.

## Version History 

### v1.1.1 [11/14/2020]

* Add integraion of scRNA-seq and cytof data demo from CITE-seq technology 
* Add `fast-SVD` module to integrate large scale datasets (cell size > 1M)  

### v1.1 [9/9/2020]

* Update parameter optimization module
* Provide joint profiles of gene expression, chromatin accessibility, and TF activity on pseudocell level. 

### v1.0.0 [7/7/2020]
* Release `bindSC`.

## System Requirements

### Hardware requirements
The `bindSC` package requires only a standard computer with enough RAM to support the in-memory operations. For minimal performance, please make sure that the computer has at least about `10 GB` of RAM. For optimal performance, we recommend a computer with the following specs:

* RAM: 10+ GB
* CPU: 4+ cores, 2.3 GHz/core

Before setting up the `bindSC` package, users should have `R` version 3.6.0 or higher, and several packages set up from CRAN and other repositories. The user can check the dependencies in `DESCRIPTION`.

## Installation 

`bindSC` is written in `R` and can be installed by following `R` commands:

``` bash
$ R
> install.packages('devtools')
> library(devtools)
> install_github('KChen-lab/bindSC')
```

Users can also install `bindSC` from source code: 
``` bash 
$ git clone https://github.com/KChen-lab/bindSC.git
$ R CMD INSTALL bindSC
```
## Usage 

For usage examples and guided walkthroughs, check the `vignettes` directory of the repo.

*  [Testing unmatched feature integration using simulated data](https://htmlpreview.github.io/?https://github.com/KChen-lab/bindSC/blob/master/vignettes/sim/sim.html)

* [Jointly Defining Cell Types from scRNA-seq and scATAC-seq](https://htmlpreview.github.io/?https://github.com/KChen-lab/bindSC/blob/master/vignettes/A549/A549.html)


* [Integraion of scRNA-seq and cytof](https://htmlpreview.github.io/?https://github.com/KChen-lab/bindSC/blob/master/vignettes/CITE-seq/main.html)


## Bug report

## License
This project is covered under the **GNU General Public License 3.0**.

## Citation
Preprint: [Unbiased integration of single cell multi-omics data](https://www.biorxiv.org/content/10.1101/2020.12.11.422014v1)
  

