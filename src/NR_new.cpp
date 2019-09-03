//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;



// [[Rcpp::export]]
double dloglik_likelihood_stratify(int knot,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat &theta)
{
int n     = delta.n_rows;
int N     = z.n_rows;

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
  mat aa;
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
List ddloglik(int knot,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat &theta,int number_facility){
  int n          = delta.n_rows;
  int p          = theta.n_rows;
  uvec delta_1   = find(delta == 1);              // record the events
  mat schoenfeld = zeros<mat>(n,p);               
  //double partial_likelihood = 0; 
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
    mat aa;
    mat bb;
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
  mat gr_test    = zeros<mat>(p*knot,number_1);
  for(int i = 0; i < number_1; i++){
    gr_test.col(i)     = kron(schoenfeld.row(delta_1(i)),b_spline.row(delta_1(i))).t();
    // partial_likelihood = partial_likelihood+
    //   (accu(z.row(delta_1(i))*theta*(b_spline.row(delta_1(i)).t()))-log(S0(delta_1(i))));
  } 
  GR_test              = sum(gr_test,1).t();
  mat dist             = solve(GVG,GR_test.t());
  dist.reshape(knot,p);
  dist=dist.t();
  List result;
  result["GVG"]        = GVG;
  result["GR_test"]    = GR_test;
  //result["partial_likelihood"] = partial_likelihood;
  result["dist"]       = dist;
  return result;
}




// [[Rcpp::export]]
double dloglik_likelihood_gradient(int knot,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat &theta)
{
  int N     = z.n_rows;
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
  mat aa;
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



// different version for NR
// [[Rcpp::export]]
List NRloop(int knot,colvec &facility,colvec &delta,mat &z,mat &b_spline,mat &theta, int M_stop, int number_facility, double rate, double tol){
	
	colvec likelihood_NR_all(M_stop+1);
	likelihood_NR_all(0) = dloglik_likelihood_stratify(knot, facility, delta, z, b_spline, theta);
	

    List temp;
    int i = 1;
    while(i <= M_stop)
	{
		temp = ddloglik(knot,facility,delta,z,b_spline,theta,number_facility);
		mat dist = temp["dist"];
		mat GVG = temp["GVG"];
		rowvec GR_test = temp["GR_test"];
		double gamma = 1;
		mat theta_temp = theta + gamma * dist;

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
			mat GVG_Q, GVG_R;
			qr(GVG_Q,GVG_R,GVG);
			mat incrementm = GR_test * (solve(GVG_R,(GVG_Q.t())*GR_test.t()));
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














