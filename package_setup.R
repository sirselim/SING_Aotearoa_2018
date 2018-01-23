####
## Package setup for SING 2018 bioinformatics workshop
## date: 190118

# install required packages from CRAN
install.packages(printr)
install.packages(scales)
install.packages(FactoMineR)
install.packages(scatterplot3d)
install.packages(plotly)
install.packages(knitr)
install.packages(prettydoc)

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
##/END