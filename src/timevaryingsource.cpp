//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
//using namespace arma;

// [[Rcpp::export]]
int sign_num(double x){
  if(x>0) return 1;
  if(x<0) return -1;
  return 0;
}


// [[Rcpp::export]]
List ddloglik_stratify(int knot,arma::colvec &facility,arma::mat &z,arma::colvec &delta,arma::mat &b_spline,arma::mat &theta,int number_facility = 1){
  int n       = delta.n_rows;
  int p       = theta.n_rows;
  arma::uvec delta_1= find(delta==1);
  arma::mat L       = zeros<arma::mat>(n,p);
  arma::rowvec L1   = zeros<arma::rowvec>(p*knot);
  arma::mat L2      = zeros<arma::mat>(p*knot,knot);
  double partial_likelihood = 0;
  arma::rowvec S0   = zeros<arma::rowvec>(n);
  arma::rowvec S1;
  arma::rowvec S2;
  for(int i = 0; i < number_facility ; i++){
    arma::uvec loc       = find(facility == i+1);
    arma::uvec loc1      = find(facility == i+1 && delta == 1);
    arma::mat z_temp     = z.rows(loc);
    int size_temp  = z_temp.n_rows;
    arma::mat b_spline_temp = b_spline.rows(loc);
    int size       = loc1.n_elem;
    int q          = 0;
    field<arma::mat> a(size_temp);
    field<arma::mat> b(size_temp);
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
List GDboost_stratify(int knot,double rate,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat theta){
  int n                   = delta.n_rows;
  int p                   = z.n_cols;
  arma::colvec unique_facility  = unique(facility);
  int number_facility     = unique_facility.n_elem;
  List result;
  arma::rowvec dist             = zeros<arma::rowvec>(knot);
  arma::mat beta_tuta           = zeros<arma::mat>(p,knot);
  arma::mat beta_tuta2          = zeros<arma::mat>(p,knot);
  beta_tuta2.fill(1/sqrt(knot));
  arma::mat beta_tuta_normalize = zeros<arma::mat>(p,knot);
  arma::mat beta_temp2          = zeros<arma::mat>(n,p);
  arma::rowvec GD               = zeros<arma::rowvec>(p);
  arma::rowvec GD1              = zeros<arma::rowvec>(p);
  arma::rowvec GD2              = zeros<arma::rowvec>(p);
  arma::rowvec test             = zeros<arma::rowvec>(p);
  arma::rowvec rate_check       = zeros<arma::rowvec>(p);
  result      = ddloglik_stratify(knot,facility,z,delta,b_spline,theta,number_facility);
  arma::mat diff    = result["L"];
  arma::rowvec L1   = result["L1"];
  arma::mat L2      = result["L2"];
  for(int i = 0; i < p; i++){
    arma::rowvec gd = zeros<arma::rowvec>(2);
    dist      = arma::solve(L2.rows(i*knot,(i+1)*knot-1),L1.cols(i*knot,(i+1)*knot-1).t()).t();
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
  arma::rowvec temp;
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
List GDboost_stratify_no_GD2(int knot,double rate,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat theta){
  int n   = delta.n_rows;
  int p   = z.n_cols;
  arma::colvec unique_facility = unique(facility);
  int number_facility    = unique_facility.n_elem;
  List result;
  arma::rowvec dist   = zeros<arma::rowvec>(knot);
  arma::mat beta_tuta = zeros<arma::mat>(p,knot);
  arma::mat beta_tuta_normalize = zeros<arma::mat>(p,knot);
  arma::mat beta_temp2= zeros<arma::mat>(n,p);
  arma::rowvec GD     = zeros<arma::rowvec>(p);
  arma::rowvec test   = zeros<arma::rowvec>(p);
  arma::rowvec rate_check = zeros<arma::rowvec>(p);
  result            = ddloglik_stratify(knot,facility,z,delta,b_spline,theta,number_facility);
  arma::mat diff          = result["L"];
  arma::rowvec L1         = result["L1"];
  arma::mat L2            = result["L2"];
  for(int i = 0; i < p; i++){
    arma::rowvec gd       = zeros<arma::rowvec>(2);
    dist            = arma::solve(L2.rows(i * knot,(i+1) * knot-1),L1.cols(i * knot,(i+1)*knot-1).t()).t();
    rate_check(i)   = (abs(dist)).max();
    beta_tuta.row(i)= beta_tuta.row(i)-dist;
    beta_tuta_normalize.row(i) = beta_tuta.row(i)/sqrt(accu(beta_tuta.row(i)%beta_tuta.row(i)));
    GD(i)           = accu(diff.col(i)%(b_spline*(beta_tuta_normalize.row(i).t())));
  }
  int j_star        = GD.index_max();
  arma::rowvec temp;
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
double dloglik_likelihood_stratify(int knot,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat &theta)
{int n     = delta.n_rows;
double partial_likelihood = 0;
arma::rowvec S0  = zeros<arma::rowvec>(n);
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
  field<arma::mat> a(size_temp);
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
List ddloglik_gradient(int knot,arma::colvec &facility,arma::mat &z,arma::colvec &delta,arma::mat &b_spline,arma::mat &theta,int number_facility=1){
  int n         = delta.n_rows;
  int p         = theta.n_rows;
  arma::uvec delta_1  = find(delta==1);
  arma::mat L         = zeros<arma::mat>(n,p);
  arma::rowvec L1     = zeros<arma::rowvec>(p*knot);
  double partial_likelihood = 0;
  arma::rowvec S0     = zeros<arma::rowvec>(n);
  arma::rowvec S1;
  for(int i = 0; i < number_facility; i++){
    arma::uvec loc   = find(facility==i+1);
    arma::uvec loc1  = find(facility==i+1&&delta==1);
    arma::mat z_temp = z.rows(loc);
    int size_temp     = z_temp.n_rows;
    arma::mat b_spline_temp = b_spline.rows(loc);
    int size          = loc1.n_elem;
    int q             = 0;
    field<arma::mat> a(size_temp);
    field<arma::mat> b(size_temp);
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
List GDboost_gradient(int knot,double rate,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat theta){
  int n = delta.n_rows;
  int p = z.n_cols;
  arma::colvec unique_facility = unique(facility);
  int number_facility    = unique_facility.n_elem;
  List result;
  arma::rowvec dist            = zeros<arma::rowvec>(knot);
  arma::mat beta_tuta          = zeros<arma::mat>(p,knot);
  arma::mat beta_tuta2         = zeros<arma::mat>(p,knot);
  beta_tuta2.fill(1/sqrt(knot));
  arma::mat beta_tuta_normalize= zeros<arma::mat>(p,knot);
  arma::mat beta_temp2         = zeros<arma::mat>(n,p);
  arma::rowvec GD              = zeros<arma::rowvec>(p);
  arma::rowvec GD1             = zeros<arma::rowvec>(p);
  arma::rowvec GD2             = zeros<arma::rowvec>(p);
  arma::rowvec test            = zeros<arma::rowvec>(p);
  arma::rowvec rate_check      = zeros<arma::rowvec>(p);
  result = ddloglik_gradient(knot,facility,z,delta,b_spline,theta,number_facility);
  arma::mat diff    = result["L"];
  arma::rowvec L1   = result["L1"];
  for(int i = 0; i < p; i++){
    arma::rowvec gd = zeros<arma::rowvec>(2);
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
  int n     = delta.n_rows;
  double partial_likelihood = 0;
  arma::rowvec S0 = zeros<arma::rowvec>(n);
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
  field<arma::mat> a(size_temp);
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
double dloglik_likelihood(arma::colvec &par,int knot,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline)
{ arma::mat theta   = par;
  int p       = z.n_cols;
  int n       = delta.n_rows;
  theta.reshape(p,knot);

  double partial_likelihood  = 0;
  arma::rowvec S0                  = zeros<arma::rowvec>(n);
  arma::colvec unique_facility     = unique(facility);
  int number_facility        = unique_facility.n_elem;
  for(int i = 0; i < number_facility; i++){
    arma::uvec loc          = find(facility==i+1);
    arma::uvec loc1         = find(facility==i+1 && delta==1);
    arma::mat z_temp        = z.rows(loc);
    int size_temp     = z_temp.n_rows;
    arma::mat b_spline_temp = b_spline.rows(loc);
    int size          = loc1.n_elem;
    int q             = 0;
    field<arma::mat> a(size_temp);
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
List ddloglik(int knot,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,arma::mat &theta,int number_facility=1){
  int n          = delta.n_rows;
  int p          = theta.n_rows;
  arma::uvec delta_1   = find(delta == 1);              // record the events
  arma::mat schoenfeld = zeros<arma::mat>(n,p);               
  double partial_likelihood = 0; 
  arma::mat GVG        = zeros<arma::mat>(p*knot,p*knot);     
  arma::rowvec GR_test = zeros<arma::rowvec>(p*knot);
  arma::rowvec GR      = zeros<arma::rowvec>(p*knot);
  arma::rowvec S0      = zeros<arma::rowvec>(n);
  arma::mat S1         = zeros<arma::mat>(n,p);
  arma::mat S2         = zeros<arma::mat>(p,p);
  for(int i = 0; i < number_facility; i++){
    arma::uvec loc     = find(facility==i+1);
    arma::uvec loc1    = find(facility==i+1&&delta==1);
    arma::mat z_temp   = z.rows(loc);
    int size_temp= z_temp.n_rows;
    arma::mat b_spline_temp = b_spline.rows(loc);
    int size     = loc1.n_elem;
    int q        = 0;
    field<arma::mat> a(size_temp);
    field<arma::mat> b(size_temp);
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
  schoenfeld     = arma::repmat(delta,1,p)%(z-S1/arma::repmat(S0.t(),1,p));
  int number_1   = delta_1.n_elem;
  arma::mat gr_test    = zeros<arma::mat>(p*knot,number_1);
  for(int i = 0; i < number_1; i++){
    gr_test.col(i)     = arma::kron(schoenfeld.row(delta_1(i)),b_spline.row(delta_1(i))).t();
    partial_likelihood = partial_likelihood+
      (accu(z.row(delta_1(i))*theta*(b_spline.row(delta_1(i)).t()))-log(S0(delta_1(i))));
  } 
  GR_test              = sum(gr_test,1).t();
  arma::mat dist             = arma::solve(GVG,GR_test.t());
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
List stratify(int knot,double tol,arma::colvec &facility,arma::colvec &delta,arma::mat &z,arma::mat &b_spline,int max_iteration=5000)
{ int p                = z.n_cols;
  arma::mat theta            = zeros<arma::mat>(p,knot);
  arma::colvec unique_facility = unique(facility);
  int number_facility    = unique_facility.n_elem;
  List result;
  double difference      = 1;
  arma::mat l2;
  arma::mat L2;
  arma::rowvec theta_NR_i;
  
  arma::colvec partial_likelihood_all(max_iteration);
  
  int iteration        = 1;
  int key              = 0;
  while(difference >= tol && iteration <= max_iteration)
  { 
    result             = ddloglik(knot,facility,delta,z,b_spline,theta,number_facility);
    iteration          = iteration+1;
    arma::mat GVG            = result["GVG"];
    arma::rowvec GR_test     = result["GR_test"];
    arma::mat dist           = result["dist"];
    double risk        = result["partial_likelihood"];
    partial_likelihood_all[key] = risk;
    key                = key+1;
    theta              = theta + dist;
    difference         = (abs(dist)).max();
  }
  
  
  
  arma::rowvec test_constrast = zeros<arma::rowvec>(p);
  arma::mat constrain         = -diff(diagmat(ones<arma::rowvec>(knot)));
  arma::mat GVG               = result["GVG"];
  double partial_likelihood = result["partial_likelihood"];
  for(int i = 0; i < p; i++){
    arma::rowvec theta_NR_i   = theta.row(i);
    arma::mat l2              = arma::inv(GVG);
    arma::mat L2              = l2.rows(i*knot,(i+1)*knot-1).cols(i*knot,(i+1)*knot-1);
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






