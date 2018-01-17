pcaGenotypes <- function(x){
  sm<-colMeans(x/2)
  gg <- t(t(x/2) - rowMeans(t(x/2))) / sqrt(sm*(1-sm))
  hh <- 1/ncol(x) * gg%*%t(gg)
  return(eigen(hh)$vectors[,1:3])
}