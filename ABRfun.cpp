// [[Rcpp::depends(RcppArmadillo)]]
#include<cmath>
#include <RcppArmadillo.h> 
#include<Rcpp.h>
using namespace Rcpp;
using namespace arma;


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
arma::mat elliproj_C(const arma::mat y, const int q, const double tau) {
  
  int j = y.n_cols/q + 1;
  arma::mat Rj = y;
  for (int l = 1; l < j; l++) {
    arma:: mat Rjl = Rj.cols(0, l*q - 1);
    double Rnorm = arma::norm(Rjl, "fro");
    double c = std::max(1 - tau/Rnorm, 0.00);
    Rj.cols(0, l*q-1) = c*Rj.cols(0, l*q-1);
  }
  return Rj;
}

// [[Rcpp::export]] 
arma::mat Phijupdate_C( const arma::mat S, const arma::mat hatOmegaj, const int q, 
                        const double rho, const arma::mat Uj, const arma::mat Psij) {
  
  int j = Uj.n_cols/q +1;
  arma::mat A = 2* hatOmegaj;
  arma::mat B = S.submat(0, 0, (j - 1)*q - 1, (j - 1)*q - 1);
  arma::mat C = 2 * hatOmegaj*S.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1)+ rho*Psij-Uj;
  arma::mat res = BS_C(A, B, rho, C);
  return res;
}

// [[Rcpp::export]] 
arma::mat Phijadmm_C(const arma::mat S, const arma::mat hatOmegaj, 
                   const arma::mat init_Phij, const int q, const double lambda1, 
                   double tol = 1e-4, const int itermax = 1e+4) {
 
  int r = init_Phij.n_cols/q;
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
  
  arma::mat Phij = init_Phij;
  arma::mat Psij = init_Phij;
  arma::mat Phij_new = zeros<mat>(q, r*q);
  arma::mat Psij_new = zeros<mat>(q, r*q);
  arma::mat Uj  =     zeros<mat>(q, r*q);
  for (int i = 0; i < itermax; i++) {
    arma::mat Phij_new = Phijupdate_C(S, hatOmegaj, q, rho, Uj, Psij);
    arma::mat Psij_new = elliproj_C(Phij_new + Uj / rho, q, lambda1 / rho);
    Uj= Uj + rho*(Phij_new  - Psij_new );
    pres = arma::norm(Phij_new  - Psij_new , "fro");
    dres = rho*arma::norm(Psij_new  - Psij, "fro");
    peps = tolabs*sqrt(r*q*q) + tolrel*std::max(arma::norm(Phij_new , "fro"), arma::norm(Psij_new , "fro"));
    deps = tolabs*sqrt(r*q*q) + tolrel*arma::norm(Uj, "fro");
    if (pres <= peps && dres <= deps)
      return Psij_new ;
    else {
      Phij = Phij_new;
      Psij= Psij_new;
      if (pres>mu*dres) {
        rho *= inc;
      }
      else if (dres>mu*pres) {
        rho /= dec;
      }
    }
  }
  Rcpp::Rcout << "ADMM fails to converge" << std::endl;
  return Psij_new;
}

// [[Rcpp::export]] 
arma::mat hatOmegaj_update_C(arma::mat S, arma::mat Phij, const int q, const double lambda2) {
  
  int j = Phij.n_cols/q + 1;
  arma::mat S0 = S.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1)*Phij.t();
  arma::mat S1 = S.submat((j - 1)*q, (j - 1)*q, j*q - 1, j*q - 1)-S0-S0.t()+
  Phij*S.submat(0, 0, (j - 1)*q - 1, (j - 1)*q - 1)*Phij.t();
  arma::mat res = glasso_C(S1, lambda2);
  return res;
  
}

// [[Rcpp::export]] 
arma::mat Siginv_C(arma::mat S, const int q, const double lambda1,
                   const double lambda2){
  
  int p = S.n_cols;
  int m = p/q;
  double eps = 1e-4;
  int itermax = 1e2;
  arma::mat bT0 = zeros<mat>(p,p);
  arma::mat bOmega0 = zeros<mat>(p, p);
  arma::mat bT1 = zeros<mat>(p,p);
  arma::mat bOmega1 = zeros<mat>(p, p);
  arma::mat I = eye<mat>(p, p);
  bOmega1.submat(0, 0, q - 1, q - 1) = glasso_C(
    S.submat(0, 0, q - 1, q - 1), lambda2);
  for (int j = 2; j <= m; j++){
    arma::mat Phij0 = bT0.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1);
    arma::mat hatOmegaj0 = bOmega0.submat((j - 1)*q, (j - 1)*q, j*q - 1, j*q - 1);
    arma::mat Phij1 =  Phijadmm_C(S, hatOmegaj0, Phij0,  q, lambda1);
    bT1.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1) = Phij1;
    arma::mat hatOmegaj1 = hatOmegaj_update_C(S, Phij1, q, lambda2);
    bOmega1.submat((j - 1)*q, (j - 1)*q, j*q - 1, j*q - 1) = hatOmegaj1;
  }
  int itercount = 0;
  double teps = accu(abs(bT1 - bT0))/accu(abs(bT0));
  double oeps = accu(abs(bOmega1 - bOmega0))/accu(abs(bOmega0));
  while( teps > eps || oeps > eps ){
    bT0 = bT1;
    bOmega0 = bOmega1;
    for (int j = 2; j <= m; j++){
      arma::mat Phij0 =   bT0.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1);
      arma::mat hatOmegaj1 = hatOmegaj_update_C(S, Phij0, q, lambda2);
      bOmega1.submat((j - 1)*q, (j - 1)*q, j*q - 1, j*q - 1) = hatOmegaj1;
      arma::mat Phij1 =  Phijadmm_C(S, hatOmegaj1, Phij0,  q, lambda1);
      bT1.submat((j - 1)*q, 0, j*q - 1, (j - 1)*q - 1) = Phij1;
    }
    teps = accu(abs(bT1 - bT0))/accu(abs(bT0));
    oeps = accu(abs(bOmega1 - bOmega0))/accu(abs(bOmega0));
    itercount += 1;
    if (itercount > itermax){
      Rcpp::Rcout << "ACS fails to converge" << std::endl;
      break;
    }
  } 
  bT1 = I - bT1;
  arma::mat Siginv = bT1.t()*bOmega1*bT1;
  return Siginv ;
}



  
  
  
  
  
   
  
  
  
  
  
  
  
  




