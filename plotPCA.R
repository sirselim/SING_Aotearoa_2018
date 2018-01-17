plotPCA <- function(x,y){
  par(mfrow=c(2,2))
  plot(0,0,col='white',xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab='',ylab='')
#   legend(0.1,0.9,fill=rainbow(length(table(y))),names(table(names(y))))
  colTab = table(popCol,names(popCol))
  colTab = colTab[apply(colTab,2,which.max),]
  pCol = rownames(colTab)
  names(pCol) = colnames(colTab)
  legend(0.1,0.9,fill=pCol,names(pCol))
  plot(x[,1],x[,2],pch=20,col=popCol,xlab="PC 1", ylab="PC 2")
  plot(x[,1],x[,3],pch=20,col=popCol,xlab="PC 1", ylab="PC 3")
  plot(x[,2],x[,3],pch=20,col=popCol,xlab="PC 2", ylab="PC 3")
}
