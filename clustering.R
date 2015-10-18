spectral = function(data, centers=c(4:10)) {
  simmatrix = 1 - as.matrix(dist(data, method="euclidean"))
  D=diag(rowSums(simmatrix))
  Dinv=diag(1/rowSums(simmatrix))
  L=diag(rep(1,nrow(simmatrix)))-Dinv %*% simmatrix
  U=D-simmatrix
  
  evL=eigen(L,symmetric=TRUE)
  evU=eigen(U,symmetric=TRUE)
  
  for (i in centers) {
    kmL=kmeans(evL$vectors[,(ncol(simmatrix)-1):(ncol(simmatrix)-0)], centers=i, nstart=5)
   # segmatL=matrix(kmL$cluster-1,dimmat[1],dimmat[2],byrow=T)
    
    #if(max(segmatL) & sum(segmatL==1)<sum(segmatL==0)) {segmatL=abs(segmatL-1)}
  
  kmU=kmeans(evU$vectors[,(ncol(simmatrix)-1):(ncol(simmatrix)-0)],centers=i,nstart=5)
  #segmatU=matrix(kmU$cluster-1,newdim[1],newdim[2],byrow=T)
  #if(max(segmatU) &sum(segmatU==1)<sum(segmatU==0)) {segmatU=abs(segmatU-1)}
  
}}

kmeans_range = function(data, range=c(2:10)) {
  d = lapply(range, function(x) {
    kmeans(data, centers = x)$cluster
  })
  d = do.call(cbind, d)
  colnames(d) = paste("k", range, sep="-")
  d
}