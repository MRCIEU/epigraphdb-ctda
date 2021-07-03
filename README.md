# EpiGraphDB CTDA MR project repository 
This is the GitHub repository to provide all information, news and updates of the EpiGraphDB COVID-19 target-disease atlas project, which is a phenome-wide Mendelian randomization study of 353 targets on COVID-19 and other diseases. 

<a href="http://epigraphdb.org/pqtl/"><img src="https://epigraphdb.org/img/epigraphdb-logo.ed38e02a.svg" alt="" height="60" style="padding:10px"/></a> <span class="pull-right">&nbsp;&nbsp;&nbsp; <a href="http://www.bris.ac.uk"><img src="https://epigraphdb.org/img/uob.16744ca9.svg" alt="" height="60" style="padding:10px"/></a>&nbsp;&nbsp;&nbsp; <a href="http://www.bris.ac.uk/ieu"><img src="http://www.bristol.ac.uk/media-library/sites/integrative-epidemiology/images/mrc-ieu-logo.png" alt="" height="60" style="padding:10px"/></a> </span>

<!-- badges: start -->

<!--
[![CRAN status](https://www.r-pkg.org/badges/version/epigraphdb)](https://cran.r-project.org/package=epigraphdb)
[![Travis build status](https://travis-ci.org/MRCIEU/epigraphdb-r.svg?branch=master)](https://travis-ci.org/MRCIEU/epigraphdb-ctda)
-->

<!-- badges: end -->

[`epigraphdb-ctda`](https://github.com/MRCIEU/epigraphdb-ctda/) is a GitHub repo to provide the following information: 
1. news and updates of the EpiGraphDB CTAD project, and issue reporting functionality for users. 
2. easy to use R scripts to run MR and colocalization analysis for the pQTL MR project. 

The MR analyses results from the EpiGraphDB CTDA project have been pre-calculated and stored in the [EpiGraphDB CTDA browser](https://epigraphdb.org/ctda/). 

## Installation of related R packages

The scripts in this repository require the following dependencies to be installed:

[`devtools`](https://devtools.r-lib.org/)
is required to install from github:

```r
###install the Two sample MR package (just need do this once) 
source("https://bioconductor.org/biocLite.R")
install.packages("devtools")

##to install/update the R package (once there is a )
library(devtools)
install_github("MRCIEU/TwoSampleMR")

#example of use the older version of the package
devtools::install_github("MRCIEU/TwoSampleMR@0.3.2")
```

## Run pQTL MR analysis

`epigraphdb-pqtl` provides a simple and intuitive way to run pQTL MR analysis using the "TwoSampleMR" R package

For more information on how to run the MR analysis, please check out R script `MR-script-covid19.R`

For more information on how to run the multi-trait colocalization analysis, please check out R script `MR-script-covid19.R`


## EpiGraphDB pQTL PheWAS browser 

The [EpiGraphDB Covid-19 target-disease atlas](https://epigraphdb.org/ctda/) currently contains the Mendelian randomization and sensitivity analyses results for 353 drug targets effects on COvid-19 and 622 other phenotypes, i.e. diseases and risk factors. To start using this browser, simply type a target or phenotype name into the "search" field, for example, [JAK2](https://epigraphdb.org/ctda/JAK2) or [Lung cancer](https://epigraphdb.org/ctda/Lung%20cancer). The full list of targets can be found by following the link on top of the "search" field.

## Citation

Please cite the pQTL MR analysis as

> Multi-omics study revealing putative drug targets of COVID-19 severity and other viral infection diseases

> Jie Zheng, Yuemiao Zhang, Yi Liu, Denis Baird, Mohd Anisul Karim, Maya Ghoussaini, Jeremy Schwartzentruber, Ian Dunham, Benjamin Elsworth, Katherine Roberts, Hannah Compton, Felix Miller-Molloy, Xingzi Liu, Lin Wang, Hong Zhang, George Davey Smith, Tom R Gaunt. Multi-omics study revealing putative drug targets of COVID-19 severity and other viral infection diseases. MedRxiv. 

```
@article {Covid-19 target-disease atlas paper,
  author = {Jie Zheng, Yuemiao Zhang, Yi Liu, Denis Baird, Mohd Anisul Karim, Maya Ghoussaini, Jeremy Schwartzentruber, Ian Dunham, Benjamin Elsworth, Katherine Roberts, Hannah Compton, Felix Miller-Molloy, Xingzi Liu, Lin Wang, Hong Zhang, George Davey Smith, Tom R Gaunt},
  title = {Multi-omics study revealing putative drug targets of COVID-19 severity and other viral infection diseases},
  url = {https://www.medrxiv.org/content/10.1101/2020.05.07.20093286v1}
}
```

## Contact

Please get in touch with us for issues, comments, suggestions, etc. via the following methods:

- [The issue tracker on the `epigrapdb-ctda` repo](https://github.com/MRCIEU/epigraphdb-ctda/issues)
- [The support email](mailto:feedback@epigraphdb.org)
- [The EpiGraphDB twitter](https://twitter.com/epigraphdb)
