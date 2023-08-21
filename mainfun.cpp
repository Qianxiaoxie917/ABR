// [[Rcpp::depends(RcppArmadillo)]]
#include<cmath>
#include <RcppArmadillo.h> 
#include<Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]] 
arma::mat lowertri(const int l) {
  arma::mat L=zeros<mat>(l,l);
  for (int i = 0; i < l; i++) {
    for (int j = 0; j <= i; j++) {
      L(i, j) = 1;
    }
  }
  return L;
}

// [[Rcpp::export]] 
arma::mat blockupdate(const int l, const int q, const double tau, arma::mat Rl) {
  arma::mat el = lowertri(l);
  arma::mat Q  = ones<mat>(q,q);
  arma:: mat indl = arma::kron(el, Q);
  arma::mat Rindl = Rl;
  for (int i = 0; i < l*q; i++) {
    for (int j = 0; j < l*q; j++) {
      if (indl(i, j) == 0) {
        Rindl(i, j) = 0;
      }
    }
  }
  double tmpnorm = arma::norm(Rindl, "fro");
  double c = std::max(1 - tau / tmpnorm, 0.00);
  for (int i = 0; i < l*q; i++) {
    for (int j = 0; j < l*q; j++) {
      if (indl(i, j) != 0) {
        Rl(i, j) = c*Rl(i, j);
      }
    }
  }
  return Rl;
}

// [[Rcpp::export]] 
arma::mat BS_C(arma::mat A, arma::mat B, const double rho, arma::mat C) {
  int m = C.n_rows;
  int n = C.n_cols;
  double tol = 1e-8;
  arma::mat X = zeros<mat>(m, n);
  arma::mat Y = zeros<mat>(m, n);
  arma::mat Q1 = zeros<mat>(m, m);
  arma::mat H = zeros<mat>(m, m);
  arma::mat Q2 = zeros<mat>(n, n);
  arma::mat S = zeros<mat>(n, n);
  if (arma::norm(C, 2) < tol) {
    return X;
  }
  arma::schur(Q1, H, A);
  arma::schur(Q2, S, B);	
  arma::mat D = Q1.t()*C*Q2;
  arma::mat I = eye<mat>(m,m);
  Y.col(n-1) = solve(trimatu(S(n-1,n-1)*H+rho*I), D.col(n-1));
  arma::mat DD = zeros<mat>(m, 1);
  for(int k = n-2; k >= 0; k--){
    for(int j = k+1; j < n; j++){
      DD +=S(k,j)*Y.col(j);
    }
    Y.col(k)=solve(trimatu(S(k,k)*H+rho*I), D.col(k)-H*DD);
  }
  X=Q1*Y*Q2.t();
  return X;
}

// [[Rcpp::export]] 
arma::mat glasso_C(arma::mat S, double lambda2){
  Environment gla("package:glasso");
  Function gl = gla["glasso"];
  double thr = 1e-4; 
  double maxit = 1e4;
  bool approx = 0; 
  bool diag = 1;
  List bc(7);
  bc = gl(S, lambda2, R_NilValue, thr, maxit, approx, diag);
  return bc["wi"];
}

// [[Rcpp::export]] 
arma::mat elliproj_C(const arma::mat y, const double tau, 
                const int m, const int q) {
  arma::mat R = y;
  for (int l = 1; l < m; l++) {
    arma:: mat Rl = R.submat((m - l)*q, 0, m*q - 1, l*q - 1);
    arma::mat Rlnew = blockupdate(l, q, tau, Rl);
    for (int i = 0; i < l*q; i++) {
      for (int j = 0; j < l*q; j++) {
        R((m - l)*q + i, j) = Rlnew(i, j);
      }
    }
  }
  return R;
}

// [[Rcpp::export]] 
arma::mat bTupdate_C( const arma::mat S, const arma::mat hatOmega, 
                const int m, const int q, const double rho, 
                const arma::mat u, const arma::mat bgamma ) {
  int p = m*q;
  arma::mat res =  eye<mat>(p, p);
  for (int j = 2; j <= m; j++) {
    arma::mat A = 2* hatOmega.submat((j - 1)*q, (j - 1)*q, j*q - 1, j*q - 1);
    arma::mat B = S.submat(0, 0, (j - 1)*q - 1, (j - 1)*q - 1);
    arma::mat C = 2 * hatOmega.submat((j - 1)*q, (j - 1)*q, j*q - 1, j*q - 1)*
                  S.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1) + 
                  rho*bgamma.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1)-
                  u.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1);
    res.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1) = BS_C(A, B, rho, C);
  }
  return res;
}

// [[Rcpp::export]] 
arma::mat bTadmm_C(const arma::mat S, const arma::mat hatOmega, 
             const arma::mat init_bT, const int m, const int q, 
             const double lambda1, double tol = 1e-4, const int itermax = 1e+4) {
  int p = m*q;
  double tolabs = tol;
  double tolrel = tol;
  double rho = 2.0;
  double mu = 10.0;
  double inc = 2.0;
  double dec = 2.0;
  
  double pres = 0.0;
  double dres = 0.0;
  double peps = 0.0;
  double deps = 0.0;
  
  arma::mat bT = init_bT;
  arma::mat bgamma = init_bT;
  arma::mat bT_new = eye(p, p);
  arma::mat bgamma_new = eye(p, p);
  arma::mat u = zeros<mat>(p, p);
  for (int i = 0; i < itermax; i++) {
    bT_new = bTupdate_C(S, hatOmega, m, q, rho, u, bgamma);
    bgamma_new = elliproj_C(bT_new + u / rho, lambda1 / rho, m, q);
    u = u + rho*(bT_new - bgamma_new);
    pres = arma::norm(bT_new - bgamma_new, "fro");
    dres = rho*arma::norm(bgamma_new - bgamma, "fro");
    peps = tolabs*p + tolrel*std::max(arma::norm(bT_new, "fro"), arma::norm(bgamma_new, "fro"));
    deps = tolabs*p + tolrel*arma::norm(u, "fro");
    if (pres <= peps && dres <= deps)
      return bgamma_new;
    else {
      bT = bT_new;
      bgamma = bgamma_new;
      if (pres>mu*dres) {
        rho *= inc;
      }
      else if (dres>mu*pres) {
        rho /= dec;
      }
    }
  }
  Rcpp::Rcout << "ADMM fails to converge" << std::endl;
  return bgamma_new;
  exit(0);
}

// [[Rcpp::export]] 
arma::mat hatOmega_update_C(arma::mat S, arma::mat bT, const int m, const int q, const double lambda2) {
  
  int p = m*q;
  arma::mat I = eye<mat>(p, p);
  arma::mat res = zeros<mat>(p, p);
  bT = I - bT;
  res.submat(0, 0, q-1, q-1) = glasso_C(S.submat(0, 0, q-1, q-1), lambda2);
  for(int j = 2; j <= m; j++){
    arma::mat bT1 = bT.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1);
    arma::mat S0 = S.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1)*bT1.t();
    arma::mat S1 = S.submat((j - 1)*q, (j - 1)*q, j*q - 1, j*q - 1)-S0-S0.t()+
                   bT1*S.submat(0, 0, (j - 1)*q - 1, (j - 1)*q - 1)*bT1.t();
    res.submat((j - 1)*q, (j - 1)*q, j*q - 1, j*q - 1) = glasso_C(S1, lambda2);
  }
  return res;
  
}

// [[Rcpp::export]] 
arma::mat Sigmainvband_CC(arma::mat S, const int m, const int q,
                          const double lambda1, const double lambda2){
  int p = m*q;
  double eps = 1e-4;
  int itermax = 1e2;
  arma::mat I = eye<mat>(p, p);
  arma::mat T0 =  eye<mat>(p, p);
  arma::mat Omega0 = hatOmega_update_C(S, T0, m, q, lambda2);
  arma::mat T1 = 2*I - bTadmm_C(S, Omega0, T0, m, q, lambda1);
  arma::mat Omega1 =  hatOmega_update_C(S, T1, m, q, lambda2);
  int itercount = 0;
  double teps = accu(abs(T1 - T0))/accu(abs(T0));
  double oeps = accu(abs(Omega1 - Omega0))/accu(abs(Omega0));
  while( teps > eps || oeps > eps ){
    T0 = T1;
    Omega0 = Omega1;
    Omega1 = hatOmega_update_C(S, T0, m, q, lambda2);
    T1 = 2*I - bTadmm_C(S, Omega1, T0, m, q, lambda1);
    teps = accu(abs(T1 - T0))/accu(abs(T0));
    oeps = accu(abs(Omega1 - Omega0))/accu(abs(Omega0));
    itercount += 1;
    if (itercount > itermax){
      Rcpp::Rcout << "ACS fails to converge" << std::endl;
      break;
    }
  }
  arma::mat Sigmainv = T1.t()*Omega1*T1;
  return Sigmainv;
}  



  






