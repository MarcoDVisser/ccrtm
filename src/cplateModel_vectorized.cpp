#include <Rcpp.h>

using namespace Rcpp;
// plate model that accepts a matrix of transmittance values
//
// @param(theta) Angle (in radians!)
// @param(ref) refractive index
// 
// @references Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
// Properties Model Separating Photosynthetic Pigments, Remote Sensing of
// Environment
// @references Stern F. (1964), Transmission of isotropic radiation across an
// interface between two dielectrics, Appl. Opt., 3(1):111-113.
// @references Allen W.A. (1973), Transmission of isotropic light across a
// dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
// 63(6):664-666.
// 
// @export
// [[Rcpp::export]]
List cplateModel_vectorized(NumericVector r12, NumericVector t12,
		 NumericVector r21, NumericVector t21,
		 NumericVector x, NumericVector y,
     NumericMatrix trans,
     NumericVector N){
		 

  NumericMatrix rho(trans.rows(),r12.size());
  NumericMatrix tau(trans.rows(),r12.size());

  double trans2,ra,ta,r90,t90,r902,t902,delta,beta,va,vb;
  double vbNN,vbNNinv,vainv,s1,s2,s3;

  for(int j =0; j < trans.rows(); j++){
    
    for(int i =0; i < r12.size(); i++){
      
      // reflectance and transmittance of the elementary layer N = 1
      trans2 = trans(j,i)*trans(j,i);
      ra = r12[i]+(t12[i]*t21[i]*r21[i]*trans2)/(1.0-(r21[i]*r21[i])*trans2);
      ta = (t12[i]*t21[i]*trans(j,i))/(1.0-(r21[i]*r21[i])*trans2);
      r90 = (ra-y[i])/x[i];
      t90 = ta/x[i];
      r902 = r90*r90;
      t902 = t90*t90;
      
      // reflectance and transmittance of N layers
      delta = sqrt(pow(t902-r902-1,2)-4.0*r902);
      beta = (1.0+r902-t902-delta)/(2.0*r90);
      va = (1.0+r902-t902+delta)/(2.0*r90);
      vb = sqrt(beta*(va-r90)/(va*(beta-r90)));
      
      vbNN = pow(vb,N[j]-1.0);
      vbNNinv = 1.0/vbNN;
      vainv = 1.0/va;
      s1 = ta*t90*(vbNN-vbNNinv);
      s2 = ta*(va-vainv);
      s3 = va*vbNN-vainv*vbNNinv-r90*(vbNN-vbNNinv);
      
      rho(j,i) = ra+s1/s3; // leaf rho
      tau(j,i) = s2/s3; // leaf tau
    }
  }
  
  return Rcpp::List::create(rho,tau);
} // cplateModel_vectorized
