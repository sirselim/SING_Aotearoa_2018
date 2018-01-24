####
## workshop script only

# source all required packages
# install and load required packages from CRAN
# install.packages(printr)
require(printr)
# install.packages(scales)
require(scales)
# install.packages(FactoMineR)
require(FactoMineR)
# install.packages(scatterplot3d)
require(scatterplot3d)
# install.packages(plotly)
require(plotly)
# install.packages(knitr)
require(knitr)
# install.packages(prettydoc)
require(prettydoc)

# install package from Github
# need to first install devtools
# install.packages(devtools)
# then grab and install package from Github
# devtools::install_github("gabraham/flashpca/flashpcaR")
require(flashpcaR)

# install Bioconductor (https://www.bioconductor.org/install/)
# source("https://bioconductor.org/biocLite.R")
# get basic packages and minfi
# biocLite()
# biocLite('minfi')
require(minfi)


# Read in the SNP data
snpData = read.table('data/GENE360snpData.txt',sep='\t',header=T)


# Make a genotype frequency table for the first SNP
table(snpData$rs3826656)


# Calculate proportions
prop.table(table(snpData$rs3826656))


# Create contingeny table - genotypes across populations
table(snpData$Population, snpData$rs3826656)

# And calculate proportions (rounded)
round(prop.table(table(snpData$Population, snpData$rs3826656),1),2)


snpFreqs = t(prop.table(table(snpData$Population, snpData$rs3826656), 1))
barplot(snpFreqs, beside=TRUE, legend.text=TRUE, ylim=c(0,1))


# Load ancestry data
snpAns = read.table('data/GENE360snpAncestry.txt', header=T, sep='\t')
# Check dimensionality of data
dim(snpAns)


# Create object containing only the SNP data - remove the first three columns
snpAnsDat = snpAns[,-c(1,2,3)]


# Load a custom function to convert to allele counts
source('alleleCounts.R')

# Apply the function to the genotype data, one column (SNP) at a time
# This will take 30 seconds or so...
snpAnsCount = apply(snpAnsDat, 2, alleleCounts)


# use flashpcaR package for fast PCA analysis
# perform PCA
r <- flashpca(snpAnsCount, do_loadings=TRUE, verbose=TRUE, stand="binom", ndim=3)


# plot component 1 vs 2
plot(r$vectors[,1], r$vectors[,2], col = alpha(as.numeric(snpAns$Population), 0.5), pch = 16, 
     main = "PCA of SNP data showing dim 1 vs dim 2")
legend('bottomleft', legend = c('AFR', 'SAS', 'EAS', 'AMR', 'EUR'), col = c('black', 'green', 'cyan', 'red', 'blue'), pch = 16)


# plot component 1 vs 3
plot(r$vectors[,1], r$vectors[,3], col = alpha(as.numeric(snpAns$Population), 0.5), pch = 16, 
     main = "PCA of SNP data showing dim 1 vs dim 3")
legend('topleft', legend = c('AFR', 'SAS', 'EAS', 'AMR', 'EUR'), col = c('black', 'green', 'cyan', 'red', 'blue'), pch = 16)


# plot component 2 vs 3
plot(r$vectors[,2], r$vectors[,3], col = alpha(as.numeric(snpAns$Population), 0.5), pch = 16, 
     main = "PCA of SNP data showing dim 2 vs dim 3")
legend('topleft', legend = c('AFR', 'SAS', 'EAS', 'AMR', 'EUR'), col = c('black', 'green', 'cyan', 'red', 'blue'), pch = 16)


# 3D plotting of PCA results
scatterplot3d(r$vectors[,1], r$vectors[,2], r$vectors[,3], color=alpha(as.numeric(snpAns$Population), 0.5), pch=16,
              cex.symbols=0.5, xlab="PC1", ylab="PC2", zlab="PC3")
legend(-3.5, 2.5, legend = c('AFR', 'SAS', 'EAS', 'AMR', 'EUR'), col = c('black', 'green', 'cyan', 'red', 'blue'), pch = 16, bty = 'n')


# use plotly for interactive graph to explore PCA results
# set margins
m <- list(
  l = 0,
  r = 0,
  b = 40,
  t = 40,
  pad = 1
)
# create plot
p <- plot_ly(x = r$vectors[,1], y = r$vectors[,2], z = r$vectors[,3], width = 900, height = 500) %>%
  add_markers(sizes=c(100,100), opacity = 0.5, color = snpAns$Population, size = r$vectors[,1]) %>%
  layout(autosize = F, margin = m,
         legend = list(x = 0.06, y = 0.04),
         scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

# embed plot
p


# Load custom function for calculating major homozygote frequency:
source('calcMajorFreq.R')

# Apply this function to the SNP count data
# Takes a few seconds
popFreqs = calcMajorFreq(snpAnsCount, snpAns$Population)


# Make a cluster tree based on these frequencies
plot(hclust(dist(popFreqs)),hang=-1)


# create some tissue methylation data
tissue_data <- data.frame(CG_1 = c(100, 50, 0), CG_2 = c(50, 50, 100), CG_3 = c(0, 100, 50))
rownames(tissue_data) <- c("blood", "heart", "brain")
tissue_data <- as.matrix(tissue_data)

barplot(tissue_data, beside=TRUE, legend.text=TRUE, ylim=c(0,100), col = c('cadetblue', 'darkred', 'grey65'), ylab = '% methylation',
        args.legend = list(bty = 'n'))



# load beta/methylation data for tissues
betas <- readRDS('data/blood_buccal.RDS')


# define samples
tissues <- c(rep(0, 5), rep(1, 5))
tissues <- as.factor(tissues)
levels(tissues) <- c('Blood', 'Buccal')


# Hierarchical Clustering
d <- t(betas)
rownames(d) <- tissues
d <- dist(d, method = "minkowski", p = 2)
hc <- hclust(d)
plot(hc, main="Cluster on all probes", labels = tissues, xlab = '', sub="")


# minfi dmpfinder
results <- dmpFinder(betas, pheno = tissues, type = "categorical", qCutoff = 1)
results <- na.omit(results)


# filter results
top.results <- results[results$qval <= 0.00005,]
head(top.results)


# create boxplots for 4 results
par(mfrow = c(2,2))
boxplot(betas['cg26297005',]*100 ~ tissues, ylim = c(50,100), main = 'cg26297005', ylab = '% methylation')
boxplot(betas['cg16780847',]*100 ~ tissues, ylim = c(50,100), main = 'cg16780847', ylab = '% methylation')
boxplot(betas['cg05350268',]*100 ~ tissues, ylim = c(0,60), main = 'cg05350268', ylab = '% methylation')
boxplot(betas['cg07759042',]*100 ~ tissues, ylim = c(50,100), main = 'cg07759042', ylab = '% methylation')


# subset the beta matrix by the top markers and cluster
sig.betas <- betas[rownames(betas) %in% rownames(top.results),]


# Hierarchical Clustering 
d <- t(sig.betas)
rownames(d) <- tissues
d <- dist(d, method = "minkowski", p = 2)
hc <- hclust(d)
plot(hc, main="Cluster on selected top tissue discrimination CpG sites", labels = tissues, xlab = '', sub="")


include_graphics('images/Bauer_Fig4.jpg')


# read in the methylation data for 10 different cell populations
blood_data <- read.csv('data/publicbeta_data.csv', head = T, row.names = 1, as.is = F)
colnames(blood_data) <- gsub('\\.', '', colnames(blood_data))


# Hierarchical Clustering 
d <- t(blood_data)
rownames(d) <- colnames(blood_data)
d <- dist(d, method = "minkowski", p = 2)
hc <- hclust(d)
plot(hc, main="Cluster on selected top cell type discrimination CpG sites", labels = colnames(blood_data), xlab = '', sub="")


# create a data frame
blood_data2 <- data.frame(cells = as.factor(colnames(blood_data)), t(blood_data))
# clean cell type groups
blood_data2$cells <- gsub('_.*', '', blood_data2$cells)
blood_data2$cells <- as.factor(blood_data2$cells)
# perfom the PCA
cell_pca <- PCA(blood_data2, quali.sup = 1, graph = F)
# create the plots# Create object containing only the SNP data - remove the first three columns
snpAnsDat = snpAns[,-c(1,2,3)]
plot.PCA(cell_pca, habillage = 'cells', label = 'quali', axes = c(1,2))
plot.PCA(cell_pca, habillage = 'cells', label = 'quali', axes = c(1,3))
plot.PCA(cell_pca, habillage = 'cells', label = 'quali', axes = c(2,3))
plot.PCA(cell_pca, habillage = 'cells', label = 'quali', axes = c(3,4))


# use plotly for interactive graph to explore PCA results
# set margins
m <- list(
  l = 0,
  r = 0,
  b = 40,
  t = 40,
  pad = 1
)
# create plot
p <- plot_ly(blood_data2, x = cell_pca$ind$coord[,1], y = cell_pca$ind$coord[,2], z = cell_pca$ind$coord[,3], width = 900, height = 550) %>%
  add_markers(sizes=c(10,100), opacity = 0.5, color = ~cells) %>%
  layout(autosize = F, margin = m,
         legend = list(x = 0.06, y = 0.04),
         scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

# embed plot
p


# load data (ms vs hc)
ms_data <- readRDS('data/disease_data.RDS')
# set variable with disease status
disease <- colnames(ms_data)


# Hierarchical Clustering 
d <- t(ms_data)
rownames(d) <- colnames(disease)
d <- dist(d, method = "minkowski", p = 2)
hc <- hclust(d)
plot(hc, main="Cluster on all CpG sites", labels = disease, xlab = '', sub="")


# minfi dmpfinder
results <- dmpFinder(as.matrix(ms_data), pheno = disease, type = "categorical", qCutoff = 1)
results <- na.omit(results)
head(results)


# create boxplots of 4 CpG sites from results
par(mfrow = c(2,2))
boxplot(as.numeric(ms_data['cg25822783',])*100 ~ disease, col = alpha(c('cadetblue', 'darkred'), 0.5), main = 'cg25822783', ylab = '% meth')
boxplot(as.numeric(ms_data['cg23365293',])*100 ~ disease, col = alpha(c('cadetblue', 'darkred'), 0.5), main = 'cg23365293', ylab = '% meth')
boxplot(as.numeric(ms_data['cg24274665',])*100 ~ disease, col = alpha(c('cadetblue', 'darkred'), 0.5), main = 'cg24274665', ylab = '% meth')
boxplot(as.numeric(ms_data['cg27343976',])*100 ~ disease, col = alpha(c('cadetblue', 'darkred'), 0.5), main = 'cg27343976', ylab = '% meth')


# create boxplots of 4 CpG sites from results
par(mfrow = c(2,2))
boxplot(as.numeric(ms_data['cg25822783',])*100 ~ disease, ylim = c(0,100),
        col = alpha(c('cadetblue', 'darkred'), 0.5), main = 'cg25822783', ylab = '% meth')
boxplot(as.numeric(ms_data['cg23365293',])*100 ~ disease, ylim = c(0,100),
        col = alpha(c('cadetblue', 'darkred'), 0.5), main = 'cg23365293', ylab = '% meth')
boxplot(as.numeric(ms_data['cg24274665',])*100 ~ disease, ylim = c(0,100),
        col = alpha(c('cadetblue', 'darkred'), 0.5), main = 'cg24274665', ylab = '% meth')
boxplot(as.numeric(ms_data['cg27343976',])*100 ~ disease, ylim = c(0,100),
        col = alpha(c('cadetblue', 'darkred'), 0.5), main = 'cg27343976', ylab = '% meth')


# get a list of cpgs passing set threshold
cpgs <- rownames(results[results$pval < 0.005,])
# filter the ms data for these sites
data_filtered <- ms_data[rownames(ms_data) %in% cpgs,]


# Hierarchical Clustering 
d <- t(data_filtered)
rownames(d) <- colnames(disease)
d <- dist(d, method = "canberra", p = 2)
hc <- hclust(d)
plot(hc, main="Cluster on selected top discriminatory CpG sites", labels = disease, xlab = '', sub="")


# load tissue data
load('/data/PostDoc/tissue_typing/data/tissue_typing_betas.RData')


