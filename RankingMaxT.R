####### RankingMaxT #######

RankingMaxT<-function(x, y, type="unpaired", B=1000, gene.names=NULL, ...){
  require(GeneSelector, quietly=TRUE)
  require(multtest, quietly=TRUE)
  gc(reset=TRUE)
  mode(x) <- "numeric"
  if(length(y) != ncol(x)){
    stop("Length of y is not equal to the number of columns of the expression matrix \n.")
  } 
  type <- match.arg(type)
  if( !is.element(type, eval(formals(RankingMaxT)$type))){ 
    stop("Permutation test is only possible for the two class
         unpaired setting \n")
  }
  ll <- eval(substitute(list(...)))
  taby <- table(y)
  if(length(taby) != 2)
  {stop("y has not exactly two levels ! \n")}
  if(!hasArg(test)) ll$test <- "wilcoxon"
  if(B>1000){
    B <- 1000
    warning("number of permutations > 1000; reset to 1000 \n")
  }
  ll$B <- B
  ll$classlabel <- ifelse(y==names(taby[2]),0,1)
  ll$X=x
  res <- do.call(mt.maxT, ll)
  idx= order(res$index)
  statistic <- res$teststat[idx]
  pvals <- res$adjp[idx]
  
  ranking <- match(1:nrow(x), res$index)
  
  if(!is.null(gene.names))
    names(pvals) <- names(statistic) <- gene.names
  else{
    if(!is.null(rownames(x)))
      names(pvals) <- names(statistic) <- rownames(x)
  }
  new("GeneRanking", x=x, y=as.factor(y), statistic=statistic,
      ranking=ranking, pval=pvals, type=type, method="maxT")
}