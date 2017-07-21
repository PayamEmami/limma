normalizeCyclicLoessBlock<-function (x,block, weights = NULL, span = 0.7, iterations = 3, method = "fast",cln=1) 
{
  require(limma)
  require(lsmeans)
  if(cln>1)
  require(parallel)
  x <- as.matrix(x)
  method <- match.arg(method, c("normal", "block"))
  n <- ncol(x)
  if (method == "normal") {
    for (k in 1:iterations) {
      a <- rowMeans(x, na.rm = TRUE)
      for (i in 1:n) {
        m <- x[, i] - a
        f <- loessFit(m, a, weights = weights, span = span)$fitted
        x[, i] <- x[, i] - f
      }
    }
  }
  if (method == "block") {
    totalgroup<-1
    block<-as.factor(block)
    for (k in 1:iterations) {
a<-NA
      if(cln>1)
{
 cl <- makeCluster(cln)
 clusterEvalQ(cl, {library(lsmeans)})
 #clusterExport(cl=cl, varlist=c("x", "block"))
 a<-parApply(cl,x,1,function(row){as.numeric(summary(lsmeans(lm(data ~ global+block,data = data.frame(block=block,global=1,data=row)),"global"))["lsmean"])})
stopCluster(cl)
}else{
      a <- apply(x,1,function(row){as.numeric(summary(lsmeans(lm(data ~ global+block,data = data.frame(block=block,global=1,data=row)),"global"))["lsmean"])})
}
      for (i in 1:n) {
        m <- x[, i] - a
        f <- loessFit(m, a, weights = weights, span = span)$fitted
        x[, i] <- x[, i] - f
      }
    }
  }
  x
}
