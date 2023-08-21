#AR(p)
ar1=function(a,n){
  M=matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n) 
      M[i,j]=outer(a,abs(i-j),"^")
  }
  return(M)
}
#counting nonzero elements
nonzerocount<-function(M){
  n1=nrow(M)
  n2=ncol(M)
  count1=0
  eps=1e-4
  for(i in 1:n1){
    for(j in 1:n2){
      if(M[i,j]>eps){
        count1=count1+1
      } 
    }
  }
  return(count1)
}
#computing TPR,TNR
S.tpn<-function(S, S.est, eps=1e-4){
  p=nrow(S)
  count0=0
  count0.est=0
  count1=0
  count1.est=0
  for(i in 1:p){
    for(j in 1:p){
      if(abs(S[i,j]) < eps){
        count0=count0+1
        if(abs(S.est[i,j]) < eps){
          count0.est=count0.est+1
        }
      }
      if(abs(S[i,j]) > eps){
        count1=count1+1
        if(abs(S.est[i,j]) > eps){
          count1.est=count1.est+1
        }
      }
    }
  }
  tnr=count0.est/count0
  tpr=count1.est/count1
  return(list(tnr=tnr, tpr=tpr)) 
}
#operator norm loss
opnloss<-function(S1,S2){
  temp=t(S1-S2)%*%(S1-S2)
  res=sqrt(max(eigen(temp)$values))
  return(res)
}
#scaled Fnorm loss
Fnloss <- function(S1, S2) {
  p = ncol(S1)
  res = as.numeric(sqrt(crossprod(as.vector(S1-S2))))
  return(res/p)
}

#scaled Kullback-Leibler loss
KLloss<-function(S1,S2){
  p = ncol(S1)
  res = sum(diag(S1%*%S2))-determinant(S1%*%S2,logarithm = TRUE)$modulus[1]-p
  return(res/p)
}
#generate positive defined sparse and not diagnoal Omega
SparseOmega<-function(q){
  Omega0 = as.matrix(rsparsematrix(q,q,nnz=q/2,symmetric = TRUE, rand.x = runif))
  diag(Omega0)=2
  return(Omega0)
}

BS<-function(A,B,rho,C){
  m=nrow(C)
  n=ncol(C)
  tol=1e-8
  if(norm(C,"2")<tol){
    X=matrix(0,m,n)
    return(X)
  }
  Y=matrix(0,m,n)
  Q1=qz(A)$Q
  Q2=qz(B)$Q
  D=t(Q1)%*%C%*%Q2 
  H=qz(A)$T
  S=qz(B)$T
  I=diag(m)
  DD=rep(0, m)
  Y[, n]=backsolve(S[n, n]*H+rho*I, D[, n])
  for(k in (n-1):1){
    for(j in (k+1):n){
      temp=S[k,j]*Y[,j]
      DD=DD+temp
    }
    Y[, k]=backsolve(S[k, k]*H+rho*I, D[, k]-H%*%DD)
  }
  X=Q1%*%Y%*%t(Q2)
  return(X)
}

makefolds<-function(n,nfolds){
  nn=round(n/nfolds)
  sizes=rep(nn,nfolds)
  sizes[nfolds]=sizes[nfolds]+n-nn*nfolds
  b=c(0,cumsum(sizes))
  ii=sample(n)
  folds=list()
  for(i in seq(nfolds)){
    folds[[i]]=ii[seq(b[i]+1,b[i+1])]
  }
  return(folds)
}

AICS<-function(S,Sigma.inv){
  res = sum(Sigma.inv%*%S) - determinant(Sigma.inv,logarithm = TRUE)$modulus[1]+
    nonzerocount(Sigma.inv)
  return(res)
}

BICS<-function(S,Sigma.inv){
  res = sum(Sigma.inv%*%S) - determinant(Sigma.inv,logarithm = TRUE)$modulus[1]+
      log(n)*nonzerocount(Sigma.inv)
  return(res)
}

nlikelihood<-function(S,Sigma.inv){
  res=sum(Sigma.inv%*%S)-determinant(Sigma.inv,logarithm = TRUE)$modulus[1]
  return(res)
}

lam1max<-function(S){
  
  p = ncol(S)
  temp = rep(0,p-1)
  for(j in 2:p){
    temp[j-1] = max(abs(S[1:(j-1), j]))/sqrt(S[j, j])
  }
  res = 2*max(temp)
  return(res)
  
}

pathGen<-function(S,nlam,flmin){
  #lambda1max = c(lam1max(S))
  lambda1max=2.5*q*sqrt(log(m*q)/n)
  lambda2max = 2*sqrt(log(p)/n)
  lam1list = sort(lambda1max*exp(seq(0, log(flmin), length= nlam)))
  lam2list = sort(lambda2max*exp(seq(0, log(flmin), length= nlam)))
  list(lam1list=lam1list,lam2list=lam2list)
}


######----------For parLapply
sumQ <- function(j, paras){
  library(glasso)
  library(Rcpp)
  sourceCpp("ABRfun.cpp")
  Qnj_C(paras$S, j, paras$q, paras$lambda1, paras$lambda2)
}

#####image matrix
matimage <- function(Mat, main = NULL){
  
  
  tmppar <- par(pty = "s")
  
  
  image(sign(t(apply(Mat, 2, rev))), axes = FALSE, col = c("gray50","white","black"), main = main)
  
  
  par(tmppar)
  
  
}

####AROLS-------------------------------------

AROLS <- function(q, Y){
  
  SS = (n-1)*cov(Y)
  
  p = ncol(SS)
  
  m = p/q
  
  bT = diag(p)
  
  bOmega = matrix(0, p, p)
  
  bD = matrix(0, p, p)
  
  bOmega[1:q, 1:q] = solve(SS[1:q, 1:q]/(n - 1))
  
  for (j in 2:m) {
    
    PPhij = diag(j*q)
    
    SSj = ginv(SS[1:((j-1)*q),1:((j-1)*q)])
    
    Phij = SS[((j-1)*q+1):(j*q),1:((j-1)*q)] %*% SSj
    
    PPhij[1:((j-1)*q),1:((j-1)*q)] = t(Phij)%*%Phij
    
    PPhij[((j-1)*q+1):(j*q),1:((j-1)*q)] = (-1) * Phij
    
    PPhij[1:((j-1)*q), ((j-1)*q+1):(j*q)] = (-1) * t(Phij)
    
    Dj = sum(diag(PPhij%*%SS[1:(j*q), 1:(j*q)]))/(n - 1)
    
    bT[((j-1)*q+1):(j*q),1:((j-1)*q)] = (-1) * Phij
    
    bD[((j-1)*q+1):(j*q), ((j-1)*q+1):(j*q)] = Dj
    
    bOmega[((j-1)*q+1):(j*q), ((j-1)*q+1):(j*q)] = ginv(Dj)
    
  }
  
  bS = t(bT)%*%bOmega%*%bT
  
  list(bT = bT, bOmega = bOmega, bD = bD, bS = bS)
}

dmax <- function(q, Y){
  
  m = ncol(Y)/q
  
  bD = AROLS(q, Y)$bD
  
  d = rep(0, m)
  
  for (j in 1:m) {
    
    d[j] = norm(bD[((j-1)*q+1):(j*q), ((j-1)*q+1):(j*q)])
    
  }
  
  dm = max(d)
  
  return(dm)
}














