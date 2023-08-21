#nflods-cross validation of lambda1 and lambda2
Sig.cv<- function(y, nlam=5, flmin=0.1, lam1list=NULL, lam2list=NULL,folds=NULL,nfolds=5){
  n=nrow(y)
  S=cov(y)*(n-1)/n
  if(is.null(folds)){
    folds=makefolds(n, nfolds)
    nfolds=length(folds)
  }
  if (is.null(lam1list) &&  (is.null(lam2list)) ){
    lambda1max=2*q*sqrt(log(m*q)/n)
    lambda2max=2*sqrt(log(q)/n)
    lam1list=sort(lambda1max*seq(flmin, 5*flmin, length= nlam))
    lam2list=sort(lambda2max*seq(flmin, 5*flmin, length= nlam))
  }
   P=array(0,c(m*q,m*q,nlam*nlam))
   errs=matrix(0, nlam*nlam, nfolds) 
    y.tr = y[-folds[[1]], ]
    S.tr = crossprod(y.tr)/(dim(y.tr)[1])
    for(i in seq(nlam)){
      for(j in seq(nlam)){
        lambda1 = lam1list[i]
        lambda2 = lam2list[j]
        P[ , ,(i-1)*nlam+j]=Siginv_C(S, q, lambda1, lambda2)
      }
    }
    for (k in seq(nfolds)) {
      
      y.te=y[folds[[k]], ]
      S.te=crossprod(y.te)/(dim(y.te)[1])
      for(j in seq(nlam*nlam)){
        errs[j,k]=AICS(S.te, P[, , j])
      }
      
    }
   me=rowMeans(errs)
  ibest=which.min(me)
  for(i in seq(nlam)){
    for(j in seq(nlam)){
      if(ibest == (i-1)*nlam+j){
        ibest1=i
        ibest2=j
      }
    }
  }
  lam1best = lam1list[ibest1]
  lam2best = lam2list[ibest2]
  S.cv = Siginv_C(S, q, lambda1, lambda2)
  list(errs=errs, me=me, lam1best = lam1list[ibest1], lam2best = lam2list[ibest2],
       ibest = c(ibest1,ibest2), S.cv=S.cv)

}

