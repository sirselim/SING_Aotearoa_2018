## Function for creating 0,1,2 genotype data, with major allele homozygote as 0
alleleCounts <- function(x){
  gt <- names(table(x))
  alleles <- unique(unlist(strsplit(gt,'')))
  oo <- order(sapply(alleles,function(z) sum(grep(z,x))),decreasing=T)
  alleles <- alleles[oo]
  geno<-c(paste(alleles[1],alleles[1],collapse='',sep=''),
          paste(sort(alleles),collapse='',sep=''),  
          paste(alleles[2],alleles[2],collapse='',sep=''))
  return(sapply(x,match,geno) - 1)
}