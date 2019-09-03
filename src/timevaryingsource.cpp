//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
int sign_num(double x){
  if(x>0) return 1;
  if(x<0) return -1;
  return 0;
}


// [[Rcpp::export]]
List ddloglik_stratify(int knot,colvec &facility,mat &z,colvec &delta,mat &b_spline,mat &theta,int number_facility = 1){
  int n       = delta.n_rows;
  int p       = theta.n_rows;
  uvec delta_1= find(delta==1);
  mat L       = zeros<mat>(n,p);
  rowvec L1   = zeros<rowvec>(p*knot);
  mat L2      = zeros<mat>(p*knot,knot);
  double partial_likelihood = 0;
  rowvec S0   = zeros<rowvec>(n);
  rowvec S1;
  rowvec S2;
  for(int i = 0; i < number_facility ; i++){
    uvec loc       = find(facility == i+1);
    uvec loc1      = find(facility == i+1 && delta == 1);
    mat z_temp     = z.rows(loc);
    int size_temp  = z_temp.n_rows;
    mat b_spline_temp = b_spline.rows(loc);
    int size       = loc1.n_elem;
    int q          = 0;
    field<mat> a(size_temp);
    field<mat> b(size_temp);
    for(int j=0; j < size_temp; j++)
    { b(j)         = z_temp.rows(j,size_temp-1);
      a(j)         = exp(b(j)*theta*b_spline_temp.row(j).t());
      S0(loc(j))   = accu(a(j));
      a(j)         = b(j)%repmat(a(j),1,p);
      S1           = sum(a(j),0);
      S2           = sum(b(j)%a(j),0);
      if(q < size && loc(j)==loc1(q))
      {
        L.row(loc(j))      = z_temp.row(j) - S1/S0(loc(j));
        L1                 = L1 - kron(L.row(loc(j)),b_spline.row(loc(j)));
        partial_likelihood = partial_likelihood + (accu(z_temp.row(j) * theta * (b_spline.row(loc(j)).t()))-log(S0(loc(j))));
        L2                 = L2 + kron((S2/S0(loc(j))-(S1%S1)/(S0(loc(j)) * S0(loc(j)))).t(),(b_spline.row(loc(j)).t()) * b_spline.row(loc(j)));
        q++;
      }
    }
  }
  List result;
  result["L"]  = L;
  result["L1"] = L1;
  result["L2"] = L2;
  result["partial_likelihood"] = partial_likelihood;
  return result;
}

// [[Rcpp::export]]
List GDboost_stratify(int knot,double rate,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat theta){
  int n                   = delta.n_rows;
  int p                   = z.n_cols;
  colvec unique_facility  = unique(facility);
  int number_facility     = unique_facility.n_elem;
  List result;
  rowvec dist             = zeros<rowvec>(knot);
  mat beta_tuta           = zeros<mat>(p,knot);
  mat beta_tuta2          = zeros<mat>(p,knot);
  beta_tuta2.fill(1/sqrt(knot));
  mat beta_tuta_normalize = zeros<mat>(p,knot);
  mat beta_temp2          = zeros<mat>(n,p);
  rowvec GD               = zeros<rowvec>(p);
  rowvec GD1              = zeros<rowvec>(p);
  rowvec GD2              = zeros<rowvec>(p);
  rowvec test             = zeros<rowvec>(p);
  rowvec rate_check       = zeros<rowvec>(p);
  result      = ddloglik_stratify(knot,facility,z,delta,b_spline,theta,number_facility);
  mat diff    = result["L"];
  rowvec L1   = result["L1"];
  mat L2      = result["L2"];
  for(int i = 0; i < p; i++){
    rowvec gd = zeros<rowvec>(2);
    dist      = solve(L2.rows(i*knot,(i+1)*knot-1),L1.cols(i*knot,(i+1)*knot-1).t()).t();
    rate_check(i)    = (abs(dist)).max();
    beta_tuta.row(i) = beta_tuta.row(i) - dist;
    beta_tuta_normalize.row(i) = beta_tuta.row(i)/sqrt(accu(beta_tuta.row(i) % beta_tuta.row(i)));
    GD1(i)           = accu(diff.col(i) % (b_spline*(beta_tuta_normalize.row(i).t())));
    gd(0)            = GD1(i);
    
    beta_temp2.col(i) = b_spline * (beta_tuta2.row(i).t());
    GD2(i)           = fabs(accu(diff.col(i)%beta_temp2.col(i)));
    gd(1)            = GD2(i);
    
    GD(i)            = gd.max();
    test(i)          = 1-gd.index_max();
  }
  int j_star         = GD.index_max();
  rowvec temp;
  if(GD1(j_star) >= GD2(j_star)){
    temp = rate * beta_tuta_normalize.row(j_star);
  }
  else{
    temp = rate*beta_tuta2.row(j_star) * sign_num(accu(diff.col(j_star) % beta_temp2.col(j_star)));
  }

  theta.row(j_star) = theta.row(j_star) + temp;
  double partial_likelihood = result["partial_likelihood"];
  List Result;
  Result["theta"]  = theta;
  Result["j_star"] = j_star;
  Result["GD"]     = GD;
  Result["GD1"]    = GD1;
  Result["GD2"]    = GD2;
  Result["test"]   = test;
  Result["rate_check"] = rate_check;
  Result["L"]      = result["L"];
  Result["L1"]     = result["L1"];
  Result["L2"]     = result["L2"];
  Result["partial_likelihood"] = partial_likelihood;
  return Result;
}


// [[Rcpp::export]]
List GDboost_stratify_no_GD2(int knot,double rate,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat theta){
  int n   = delta.n_rows;
  int p   = z.n_cols;
  colvec unique_facility = unique(facility);
  int number_facility    = unique_facility.n_elem;
  List result;
  rowvec dist   = zeros<rowvec>(knot);
  mat beta_tuta = zeros<mat>(p,knot);
  mat beta_tuta_normalize = zeros<mat>(p,knot);
  mat beta_temp2= zeros<mat>(n,p);
  rowvec GD     = zeros<rowvec>(p);
  rowvec test   = zeros<rowvec>(p);
  rowvec rate_check = zeros<rowvec>(p);
  result            = ddloglik_stratify(knot,facility,z,delta,b_spline,theta,number_facility);
  mat diff          = result["L"];
  rowvec L1         = result["L1"];
  mat L2            = result["L2"];
  for(int i = 0; i < p; i++){
    rowvec gd       = zeros<rowvec>(2);
    dist            = solve(L2.rows(i * knot,(i+1) * knot-1),L1.cols(i * knot,(i+1)*knot-1).t()).t();
    rate_check(i)   = (abs(dist)).max();
    beta_tuta.row(i)= beta_tuta.row(i)-dist;
    beta_tuta_normalize.row(i) = beta_tuta.row(i)/sqrt(accu(beta_tuta.row(i)%beta_tuta.row(i)));
    GD(i)           = accu(diff.col(i)%(b_spline*(beta_tuta_normalize.row(i).t())));
  }
  int j_star        = GD.index_max();
  rowvec temp;
  temp              = rate*beta_tuta_normalize.row(j_star);
  theta.row(j_star) = theta.row(j_star) + temp;
  double partial_likelihood = result["partial_likelihood"];
  List Result;
  Result["theta"]   = theta;
  Result["j_star"]  = j_star;
  Result["GD"]      = GD;
  Result["rate_check"] = rate_check;
  Result["L"]       = result["L"];
  Result["L1"]      = result["L1"];
  Result["L2"]      = result["L2"];
  Result["partial_likelihood"] = partial_likelihood;
  return Result;
}

// This function gives statrified likelihood
// [[Rcpp::export]]
double dloglik_likelihood_stratify(int knot,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat &theta)
{int n     = delta.n_rows;
double partial_likelihood = 0;
rowvec S0  = zeros<rowvec>(n);
colvec unique_facility = unique(facility);
int number_facility    = unique_facility.n_elem;
for(int i = 0; i < number_facility; i++){
  uvec loc   = find(facility == i+1);
  uvec loc1  = find(facility == i+1 && delta == 1);
  mat z_temp = z.rows(loc);
  int size_temp     = z_temp.n_rows;
  mat b_spline_temp = b_spline.rows(loc);
  int size          = loc1.n_elem;
  int q             = 0;
  field<mat> a(size_temp);
  for(int j = 0; j < size_temp; j++)
  {
    a(j)       = exp(z_temp.rows(j,size_temp-1)*theta*b_spline_temp.row(j).t());
    S0(loc(j)) = accu(a(j));
    if(q < size&&loc(j) == loc1(q))
    {
      partial_likelihood = partial_likelihood + (accu(z_temp.row(j)*theta*(b_spline.row(loc(j)).t()))-log(S0(loc(j))));
      q++;}
  }
}
return partial_likelihood;
}

// [[Rcpp::export]]
List ddloglik_gradient(int knot,colvec &facility,mat &z,colvec &delta,mat &b_spline,mat &theta,int number_facility=1){
  int n         = delta.n_rows;
  int p         = theta.n_rows;
  uvec delta_1  = find(delta==1);
  mat L         = zeros<mat>(n,p);
  rowvec L1     = zeros<rowvec>(p*knot);
  double partial_likelihood = 0;
  rowvec S0     = zeros<rowvec>(n);
  rowvec S1;
  for(int i = 0; i < number_facility; i++){
    uvec loc   = find(facility==i+1);
    uvec loc1  = find(facility==i+1&&delta==1);
    mat z_temp = z.rows(loc);
    int size_temp     = z_temp.n_rows;
    mat b_spline_temp = b_spline.rows(loc);
    int size          = loc1.n_elem;
    int q             = 0;
    field<mat> a(size_temp);
    field<mat> b(size_temp);
    for(int j = 0 ; j < size_temp; j++)
    { b(j)    = z_temp.rows(j,size_temp-1);
      a(j)    = exp(b(j) * theta * b_spline_temp.row(j).t());
      S0(loc(j)) = accu(a(j));
      a(j)       = b(j)%repmat(a(j),1,p);
      S1         = sum(a(j),0);
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
List GDboost_gradient(int knot,double rate,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat theta){
  int n = delta.n_rows;
  int p = z.n_cols;
  colvec unique_facility = unique(facility);
  int number_facility    = unique_facility.n_elem;
  List result;
  rowvec dist            = zeros<rowvec>(knot);
  mat beta_tuta          = zeros<mat>(p,knot);
  mat beta_tuta2         = zeros<mat>(p,knot);
  beta_tuta2.fill(1/sqrt(knot));
  mat beta_tuta_normalize= zeros<mat>(p,knot);
  mat beta_temp2         = zeros<mat>(n,p);
  rowvec GD              = zeros<rowvec>(p);
  rowvec GD1             = zeros<rowvec>(p);
  rowvec GD2             = zeros<rowvec>(p);
  rowvec test            = zeros<rowvec>(p);
  rowvec rate_check      = zeros<rowvec>(p);
  result = ddloglik_gradient(knot,facility,z,delta,b_spline,theta,number_facility);
  mat diff    = result["L"];
  rowvec L1   = result["L1"];
  for(int i = 0; i < p; i++){
    rowvec gd = zeros<rowvec>(2);
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
  rowvec temp;
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
double dloglik_likelihood_gradient(int knot,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat &theta)
{
  int n     = delta.n_rows;
  double partial_likelihood = 0;
  rowvec S0 = zeros<rowvec>(n);
  colvec unique_facility = unique(facility);
  int number_facility = unique_facility.n_elem;

for(int i = 0; i < number_facility; i++){
  uvec loc   = find(facility==i+1);
  uvec loc1  = find(facility==i+1 && delta==1);
  mat z_temp = z.rows(loc);
  int size_temp     = z_temp.n_rows;
  mat b_spline_temp = b_spline.rows(loc);
  int size          = loc1.n_elem;
  int q             = 0;
  field<mat> a(size_temp);
  for(int j = 0; j < size_temp; j++)
  {
    a(j)       = exp(z_temp.rows(j,size_temp-1)*theta*b_spline_temp.row(j).t());
    S0(loc(j)) = accu(a(j));
    if(q < size && loc(j) == loc1(q))
    {
      partial_likelihood = partial_likelihood+(accu(z_temp.row(j) * theta * (b_spline.row(loc(j)).t()))-log(S0(loc(j))));
      q++;
    }
  }
}
return partial_likelihood;
}



// [[Rcpp::export]]
double dloglik_likelihood(colvec &par,int knot,colvec &facility,colvec &delta,mat &z,mat &b_spline)
{ mat theta   = par;
  int p       = z.n_cols;
  int n       = delta.n_rows;
  theta.reshape(p,knot);

  double partial_likelihood  = 0;
  rowvec S0                  = zeros<rowvec>(n);
  colvec unique_facility     = unique(facility);
  int number_facility        = unique_facility.n_elem;
  for(int i = 0; i < number_facility; i++){
    uvec loc          = find(facility==i+1);
    uvec loc1         = find(facility==i+1 && delta==1);
    mat z_temp        = z.rows(loc);
    int size_temp     = z_temp.n_rows;
    mat b_spline_temp = b_spline.rows(loc);
    int size          = loc1.n_elem;
    int q             = 0;
    field<mat> a(size_temp);
    for(int j = 0; j < size_temp; j++)
    {
      a(j)       = exp(z_temp.rows(j,size_temp-1)*theta*b_spline_temp.row(j).t());
      S0(loc(j)) = accu(a(j));
      if(q < size && loc(j) == loc1(q))
      {
        partial_likelihood = partial_likelihood + (accu(z_temp.row(j) * theta * (b_spline.row(loc(j)).t())) - log(S0(loc(j))));
        q++;
      }
    }
  }
  return -partial_likelihood;
}


// [[Rcpp::export]]
List ddloglik(int knot,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat &theta,int number_facility=1){
  int n          = delta.n_rows;
  int p          = theta.n_rows;
  uvec delta_1   = find(delta == 1);              // record the events
  mat schoenfeld = zeros<mat>(n,p);               
  double partial_likelihood = 0; 
  mat GVG        = zeros<mat>(p*knot,p*knot);     
  rowvec GR_test = zeros<rowvec>(p*knot);
  rowvec GR      = zeros<rowvec>(p*knot);
  rowvec S0      = zeros<rowvec>(n);
  mat S1         = zeros<mat>(n,p);
  mat S2         = zeros<mat>(p,p);
  for(int i = 0; i < number_facility; i++){
    uvec loc     = find(facility==i+1);
    uvec loc1    = find(facility==i+1&&delta==1);
    mat z_temp   = z.rows(loc);
    int size_temp= z_temp.n_rows;
    mat b_spline_temp = b_spline.rows(loc);
    int size     = loc1.n_elem;
    int q        = 0;
    field<mat> a(size_temp);
    field<mat> b(size_temp);
    for(int j = 0; j < size_temp; j++)
    { b(j)       = z_temp.rows(j,size_temp-1);
      a(j)       = exp(b(j) * theta * b_spline_temp.row(j).t());
      S0(loc(j)) = accu(a(j));
      a(j)       = b(j)%repmat(a(j),1,p);
      S1.row(loc(j)) = sum(a(j),0);
      if(q<size&&loc(j) == loc1(q))
      {for(int k = 0; k < p; k++)
      {S2.col(k) = sum(a(j) % repmat(b(j).col(k),1,p),0).t();}
      GVG        = GVG + kron(S2/S0(loc(j))-(S1.row(loc(j)).t()*S1.row(loc(j))/(pow(S0(loc(j)),2))),
                   b_spline_temp.row(j).t() * b_spline_temp.row(j));

      q++;
      }
    }
  }
  schoenfeld     = repmat(delta,1,p)%(z-S1/repmat(S0.t(),1,p));
  int number_1   = delta_1.n_elem;
  mat gr_test    = zeros<mat>(p*knot,number_1);
  for(int i = 0; i < number_1; i++){
    gr_test.col(i)     = kron(schoenfeld.row(delta_1(i)),b_spline.row(delta_1(i))).t();
    partial_likelihood = partial_likelihood+
      (accu(z.row(delta_1(i))*theta*(b_spline.row(delta_1(i)).t()))-log(S0(delta_1(i))));
  } 
  GR_test              = sum(gr_test,1).t();
  mat dist             = solve(GVG,GR_test.t());
  dist.reshape(knot,p);
  dist=dist.t();
  List result;
  result["GVG"]        = GVG;
  result["GR_test"]    = GR_test;
  result["partial_likelihood"] = partial_likelihood;
  result["dist"]       = dist;
  return result;
}

// [[Rcpp::export]]
List stratify(int knot,double tol,colvec &facility,colvec &delta,mat &z,mat &b_spline,int max_iteration=5000)
{ int p                = z.n_cols;
  mat theta            = zeros<mat>(p,knot);
  colvec unique_facility = unique(facility);
  int number_facility    = unique_facility.n_elem;
  List result;
  double difference      = 1;
  mat l2;
  mat L2;
  rowvec theta_NR_i;
  
  colvec partial_likelihood_all(max_iteration);
  
  int iteration        = 1;
  int key              = 0;
  while(difference >= tol && iteration <= max_iteration)
  { 
    result             = ddloglik(knot,facility,delta,z,b_spline,theta,number_facility);
    iteration          = iteration+1;
    mat GVG            = result["GVG"];
    rowvec GR_test     = result["GR_test"];
    mat dist           = result["dist"];
    double risk        = result["partial_likelihood"];
    partial_likelihood_all[key] = risk;
    key                = key+1;
    theta              = theta + dist;
    difference         = (abs(dist)).max();
  }
  
  
  
  rowvec test_constrast = zeros<rowvec>(p);
  mat constrain         = -diff(diagmat(ones<rowvec>(knot)));
  mat GVG               = result["GVG"];
  double partial_likelihood = result["partial_likelihood"];
  for(int i = 0; i < p; i++){
    rowvec theta_NR_i   = theta.row(i);
    mat l2              = inv(GVG);
    mat L2              = l2.rows(i*knot,(i+1)*knot-1).cols(i*knot,(i+1)*knot-1);
    test_constrast(i)   = accu((constrain*theta_NR_i.t()).t()*inv(constrain*L2*constrain.t())
                             *(constrain*theta_NR_i.t()));
  }
  List Result;
  Result["test_constrast"] = test_constrast;
  Result["H"]              = GVG;
  Result["theta"]          = theta;
  Result["knot"]           = knot;
  Result["partial_likelihood"] = partial_likelihood;
  Result["partial_likelihood_all"] = partial_likelihood_all;
  Result["GR_test"]        = result["GR_test"];
  Result["Iteration"]      = key;
  return Result;
}






