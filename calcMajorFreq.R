calcMajorFreq<-function(x,y){
## x is count-based genotypes
## y is population data  
  pops<-names(table(y))
  gtFreq<-matrix(0,length(pops),ncol(x))
  for(i in seq(pops))
    gtFreq[i,] <- apply(x[as.vector(y)==pops[i],],2,function(z) 1-prop.table(table(z))[1])
  gtFreq[is.na(gtFreq)]<-0
  rownames(gtFreq)<-pops
  colnames(gtFreq)<-colnames(x)
  return(gtFreq)
}

