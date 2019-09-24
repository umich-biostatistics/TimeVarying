//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int sign_num(double x){
  if(x>0) return 1;
  if(x<0) return -1;
  return 0;
}




// This function gives statrified likelihood
// [[Rcpp::export]]
double dloglik_likelihood_stratify(int knot,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat &theta)
{
int n     = delta.n_rows;
int N     = z.n_rows;

double partial_likelihood = 0;
arma::rowvec S0  = arma::zeros<arma::rowvec>(n);
arma::colvec unique_facility = unique(facility);
int number_facility    = unique_facility.n_elem;
for(int i = 0; i < number_facility; i++){
  arma::uvec loc   = find(facility == i+1);
  arma::uvec loc1  = find(facility == i+1 && delta == 1);
  arma::mat z_temp = z.rows(loc);
  int size_temp     = z_temp.n_rows;
  arma::mat b_spline_temp = b_spline.rows(loc);
  int size          = loc1.n_elem;
  int q             = 0;
  arma::mat aa;
  for(int j = 0; j < size_temp; j++)
  {
    aa         = exp(z_temp.rows(j,size_temp-1)*theta*b_spline_temp.row(j).t());
    S0(loc(j)) = accu(aa);
    if(q < size&&loc(j) == loc1(q))
    {
      partial_likelihood = partial_likelihood + (accu(z_temp.row(j)*theta*(b_spline.row(loc(j)).t()))-log(S0(loc(j))));
      q++;}
  }
}
partial_likelihood = partial_likelihood/(N*1.0);

return partial_likelihood;
}


// [[Rcpp::export]]
List ddloglik_gradient(int knot,arma::colvec &facility,arma::mat &z,arma::colvec &delta,arma::mat &b_spline,arma::mat &theta,int number_facility=1){
  int n         = delta.n_rows;
  int p         = theta.n_rows;
  arma::uvec delta_1  = find(delta==1);
  arma::mat L         = arma::zeros<arma::mat>(n,p);
  arma::rowvec L1     = arma::zeros<arma::rowvec>(p*knot);
  double partial_likelihood = 0;
  arma::rowvec S0     = arma::zeros<arma::rowvec>(n);
  arma::rowvec S1;
  for(int i = 0; i < number_facility; i++){
    arma::uvec loc   = find(facility==i+1);
    arma::uvec loc1  = find(facility==i+1&&delta==1);
    arma::mat z_temp = z.rows(loc);
    int size_temp     = z_temp.n_rows;
    arma::mat b_spline_temp = b_spline.rows(loc);
    int size          = loc1.n_elem;
    int q             = 0;
    arma::mat aa;
    arma::mat bb;
    for(int j = 0 ; j < size_temp; j++)
    { bb    = z_temp.rows(j,size_temp-1);
      aa    = exp(bb * theta * b_spline_temp.row(j).t());
      S0(loc(j)) = accu(aa);
      aa       = bb%repmat(aa,1,p);
      S1         = sum(aa,0);
      if(q<size&&loc(j)==loc1(q))
      {
        L.row(loc(j)) = z_temp.row(j)-S1/S0(loc(j));
        L1            = L1-kron(L.row(loc(j)),b_spline.row(loc(j)));
        partial_likelihood = partial_likelihood+(accu(z_temp.row(j) * theta * (b_spline.row(loc(j)).t()))-log(S0(loc(j))));
        q++;}
    }
  }
  List result;
  result["L"]  = L;
  result["L1"] = L1;
  result["partial_likelihood"] = partial_likelihood;
  return result;
}

// [[Rcpp::export]]
List GDboost_gradient(int knot,double rate,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat theta){
  int n = delta.n_rows;
  int p = z.n_cols;
  arma::colvec unique_facility = unique(facility);
  int number_facility    = unique_facility.n_elem;
  List result;
  arma::rowvec dist            = arma::zeros<arma::rowvec>(knot);
  arma::mat beta_tuta          = arma::zeros<arma::mat>(p,knot);
  arma::mat beta_tuta2         = arma::zeros<arma::mat>(p,knot);
  beta_tuta2.fill(1/sqrt(knot));
  arma::mat beta_tuta_normalize= arma::zeros<arma::mat>(p,knot);
  arma::mat beta_temp2         = arma::zeros<arma::mat>(n,p);
  arma::rowvec GD              = arma::zeros<arma::rowvec>(p);
  arma::rowvec GD1             = arma::zeros<arma::rowvec>(p);
  arma::rowvec GD2             = arma::zeros<arma::rowvec>(p);
  arma::rowvec test            = arma::zeros<arma::rowvec>(p);
  arma::rowvec rate_check      = arma::zeros<arma::rowvec>(p);
  result = ddloglik_gradient(knot,facility,z,delta,b_spline,theta,number_facility);
  arma::mat diff    = result["L"];
  arma::rowvec L1   = result["L1"];
  for(int i = 0; i < p; i++){
    arma::rowvec gd = arma::zeros<arma::rowvec>(2);
    beta_tuta_normalize.row(i) = L1.cols(i*knot,(i+1)*knot-1)/sqrt(accu(L1.cols(i*knot,(i+1)*knot-1)%L1.cols(i*knot,(i+1)*knot-1)));
    GD1(i)    = accu(L1.cols(i*knot,(i+1)*knot-1)%L1.cols(i*knot,(i+1)*knot-1));
    gd(0)     = GD1(i);
    
    beta_temp2.col(i) = diff.col(i);
    GD2(i)            = fabs(accu(L1.cols(i*knot,(i+1)*knot-1)%L1.cols(i*knot,(i+1)*knot-1)));
    gd(1)             = GD2(i);
    
    GD(i)             = gd.max();
    test(i)           = 1-gd.index_max();
  }
  int j_star          = GD.index_max();
  arma::rowvec temp;
  temp                = rate * beta_tuta_normalize.row(j_star);
  theta.row(j_star)   = theta.row(j_star)-temp;
  double partial_likelihood = result["partial_likelihood"];
  
  List Result;
  Result["theta"]   = theta;
  Result["j_star"]  = j_star;
  Result["GD"]      = GD;
  Result["test"]    = test;
  Result["L"]       = result["L"];
  Result["L1"]      = result["L1"];
  Result["partial_likelihood"] = partial_likelihood;
  Result["temp"]    = temp;

  return Result;
}

// [[Rcpp::export]]
double dloglik_likelihood_gradient(int knot,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat &theta)
{
  int N     = z.n_rows;
  int n     = delta.n_rows;
  double partial_likelihood = 0;
  arma::rowvec S0 = arma::zeros<arma::rowvec>(n);
  arma::colvec unique_facility = unique(facility);
  int number_facility = unique_facility.n_elem;

for(int i = 0; i < number_facility; i++){
  arma::uvec loc   = find(facility==i+1);
  arma::uvec loc1  = find(facility==i+1 && delta==1);
  arma::mat z_temp = z.rows(loc);
  int size_temp     = z_temp.n_rows;
  arma::mat b_spline_temp = b_spline.rows(loc);
  int size          = loc1.n_elem;
  int q             = 0;
  arma::mat aa;
  for(int j = 0; j < size_temp; j++)
  {
    aa       = exp(z_temp.rows(j,size_temp-1)*theta*b_spline_temp.row(j).t());
    S0(loc(j)) = accu(aa);
    if(q < size && loc(j) == loc1(q))
    {
      partial_likelihood = partial_likelihood+(accu(z_temp.row(j) * theta * (b_spline.row(loc(j)).t()))-log(S0(loc(j))));
      q++;
    }
  }
}
partial_likelihood = partial_likelihood/(N*1.0);
return partial_likelihood;
}



// [[Rcpp::export]]
List ddloglik(int knot,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat &theta,int number_facility){
  int n          = delta.n_rows;
  int p          = theta.n_rows;
  arma::uvec delta_1   = find(delta == 1);              // record the events
  arma::mat schoenfeld = arma::zeros<arma::mat>(n,p);               
  //double partial_likelihood = 0; 
  arma::mat GVG        = arma::zeros<arma::mat>(p*knot,p*knot);     
  arma::rowvec GR_test = arma::zeros<arma::rowvec>(p*knot);
  arma::rowvec GR      = arma::zeros<arma::rowvec>(p*knot);
  arma::rowvec S0      = arma::zeros<arma::rowvec>(n);
  arma::mat S1         = arma::zeros<arma::mat>(n,p);
  arma::mat S2         = arma::zeros<arma::mat>(p,p);
  for(int i = 0; i < number_facility; i++){
    arma::uvec loc     = find(facility==i+1);
    arma::uvec loc1    = find(facility==i+1&&delta==1);
    arma::mat z_temp   = z.rows(loc);
    int size_temp= z_temp.n_rows;
    arma::mat b_spline_temp = b_spline.rows(loc);
    int size     = loc1.n_elem;
    int q        = 0;
    arma::mat aa;
    arma::mat bb;
    for(int j = 0; j < size_temp; j++)
    { bb         = z_temp.rows(j,size_temp-1);
      aa       = exp(bb * theta * b_spline_temp.row(j).t());
      S0(loc(j)) = accu(aa);
      aa       = bb%repmat(aa,1,p);
      S1.row(loc(j)) = sum(aa,0);
      if(q<size&&loc(j) == loc1(q))
      {for(int k = 0; k < p; k++)
      {S2.col(k) = sum(aa % repmat(bb.col(k),1,p),0).t();}
      GVG        = GVG + kron(S2/S0(loc(j))-(S1.row(loc(j)).t()*S1.row(loc(j))/(pow(S0(loc(j)),2))),
                   b_spline_temp.row(j).t() * b_spline_temp.row(j));
      q++;
      }
    }
  }

  schoenfeld     = repmat(delta,1,p)%(z-S1/repmat(S0.t(),1,p));
  int number_1   = delta_1.n_elem;
  arma::mat gr_test    = arma::zeros<arma::mat>(p*knot,number_1);
  for(int i = 0; i < number_1; i++){
    gr_test.col(i)     = kron(schoenfeld.row(delta_1(i)),b_spline.row(delta_1(i))).t();
    // partial_likelihood = partial_likelihood+
    //   (accu(z.row(delta_1(i))*theta*(b_spline.row(delta_1(i)).t()))-log(S0(delta_1(i))));
  } 
  GR_test              = sum(gr_test,1).t();
  arma::mat dist             = solve(GVG,GR_test.t());
  dist.reshape(knot,p);
  dist=dist.t();
  List result;
  result["GVG"]        = GVG;
  result["GR_test"]    = GR_test;
  //result["partial_likelihood"] = partial_likelihood;
  result["dist"]       = dist;
  return result;
}


// different version for NR
// [[Rcpp::export]]
List NRloop(int knot,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat &theta, int M_stop, int number_facility, double rate, double tol){
  
  arma::colvec likelihood_NR_all(M_stop+1);
  likelihood_NR_all(0) = dloglik_likelihood_stratify(knot, facility, delta, z, b_spline, theta);
  

    List temp;
    int i = 1;
    while(i <= M_stop)
  {
    temp = ddloglik(knot,facility,delta,z,b_spline,theta,number_facility);
    arma::mat dist = temp["dist"];
    arma::mat GVG = temp["GVG"];
    arma::rowvec GR_test = temp["GR_test"];
    double gamma = 1;
    arma::mat theta_temp = theta + gamma * dist;

    double likelihood_stratify = dloglik_likelihood_stratify(knot,facility, delta, z, b_spline, theta_temp);

    while(likelihood_stratify < likelihood_NR_all(i-1) + gamma * rate){
      gamma = gamma / 2.0;
      theta_temp = theta + gamma * dist;
      likelihood_stratify = dloglik_likelihood_gradient(knot, facility, delta, z, b_spline, theta_temp);
    }
    theta = theta + gamma*dist;

    likelihood_stratify = dloglik_likelihood_stratify(knot, facility, delta, z, b_spline, theta);
    likelihood_NR_all(i) = likelihood_stratify;

    if(i >= 2){
      double llk_diff  = likelihood_NR_all(i) - likelihood_NR_all(i-1);
      double llk_diff2 = likelihood_NR_all(i) -likelihood_NR_all(1);
      arma::mat GVG_Q, GVG_R;
      qr(GVG_Q,GVG_R,GVG);
      arma::mat incrementm = GR_test * (solve(GVG_R,(GVG_Q.t())*GR_test.t()));
      double increment = incrementm(0,0);
      if(abs(llk_diff/llk_diff2) < tol || increment < tol){
        break;
      }
    }
    i++;
  }


  List result;
  result["theta"]  = theta;
  result["likelihood_NR_all"] = likelihood_NR_all;
  result["GVG"]    = temp["GVG"];

  return result;

}




















