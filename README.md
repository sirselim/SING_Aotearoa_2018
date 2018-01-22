# SING Aotearoa Bioinformatics Workshop 2018

## Bioinformatics – A Genetic and Epigenetic Taster

Welcome to the bioinformatics 'taster' workshop for the SING Aotearoa 2018 program. This is aimed to be a basic introduction where we will interrogate and analyse a combination of genetic and epigenetic data sets, explore variation at both the single marker and multi-marker level and discuss what this means, how it is interpreted and what it may be used for.

There are two aspects of what you've learnt so far in the program that we will explore with some basic bioinformatic tools:

  - Populations Genetics
  - Epigenetics
  
### 1. Population Genetics

We will investigate genetic variation (SNPs) from several populations which have had individuals sequenced as part of the 1000 Genomes project.

***Note:*** *The majority of code and data for this section is borrowed/modified from the course **GENE360** maintained by [Associate Professor Mik Black](http://www.otago.ac.nz/genetics/staff/members/black.html). You can find the full resource here:* https://github.com/mikblack/GENE360-PopDiv

### 2. Epigenetics

We will explore some DNA methylation data and demonstrate how specific markers may allow the discrimination of both tissues and cells.

***Note:*** *This section borrows from previous epigenetic workshops, see below for the full resources:*

  - GitHub repository: https://github.com/sirselim/methylation_EWAS_workshop  
  - Full workshop: https://sirselim.github.io/methylation_EWAS_workshop/
  
## Getting started and following along

This workshop requires you to have installed both R and RStudio, basic knowledge in either/both is nice but not essential. They can be obtained from here:

  - R: https://cran.rstudio.com/
  - RStudio: https://www.rstudio.com/product/rstudio/download/  

An RMarkdown document has been created and is part of this repository, you can either download the repository and open the `SING_workshop.Rmd` file in RStudio and follow along or feel free to work through the online version.

The online version of this workshop can be viewed here (**best viewed in Google Chrome**): https://rawgit.com/sirselim/SING_Aotearoa_2018/master/SING_workshop.html

### Required Packages

This workshop requires a few specific R packages to be installed. The below contains all packages required to follow along with the RMarkdown file, and how to install them.

```R
# install required packages from CRAN
install.packages(printr)
install.packages(scales)
install.packages(FactoMineR)
install.packages(scatterplot3d)
install.packages(plotly)
install.packages(knitr)

# install package from Github
# need to first install devtools
install.packages(devtools)
# then grab and install package from Github
devtools::install_github("gabraham/flashpca/flashpcaR")

# install Bioconductor (https://www.bioconductor.org/install/)
source("https://bioconductor.org/biocLite.R")
# get basic packages and minfi
biocLite()
biocLite('minfi')
```