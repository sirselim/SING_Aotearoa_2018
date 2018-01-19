# SING Aotearoa Bioinformatics Workshop 2018

## Bioinformatics – A Genetic and Epigenetic Taster

There are two aspects of what you've learnt so far that we will explore with some basic bioinformatic tools:

  - Populations Genetics
  - Epigenetics
  
## 1. Population Genetics

***Note:*** *The majority of code and data for this section is borrowed from the course GENE360 maintained by Associate Professor Mik Black. You can find the full resource here:* https://github.com/mikblack/GENE360-PopDiv

## 2. Epigenetics

***Note:*** *This section borrows from previous epigenetic workshops, see below for the full resources:*

  - GitHub repository: https://github.com/sirselim/methylation_EWAS_workshop  
  - Full workshop: https://sirselim.github.io/methylation_EWAS_workshop/
  

The online version of this workshop can be viewed here (**best viewed in Google Chrome**): https://rawgit.com/sirselim/SING_Aotearoa_2018/master/SING_workshop.html

**Note:** if the above is too slow loading or doesn't work at all here is a link to the html document without the fancy interactive plots: https://rawgit.com/sirselim/SING_Aotearoa_2018/master/SING_workshop_no3d.html

## Required Packages

The below contains all packages required to follow along with the RMarkdown file, and how to install them.

```R
# install required packages from CRAN
install.packages(printr)
install.packages(scales)
install.packages(formatR)
install.packages(FactoMineR)
install.packages(scatterplot3d)
install.packages(plotly)

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