#include <Rcpp.h>

using namespace Rcpp;
// tav function
// 
// Stern F. (1964), Transmission of isotropic radiation across an
// interface between two dielectrics, Appl. Opt., 3(1):111-113.
// Allen W.A. (1973), Transmission of isotropic light across a
// dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
// 63(6):664-666.
// @references Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
// Properties Model Separating Photosynthetic Pigments, Remote Sensing of
// Environment
// 
// @param(theta) Angle (in degrees!)
// @param(ref) refractive index
// @export
// [[Rcpp::export]]
NumericVector ctav(double theta, NumericVector n){

  static const double pi = 3.14159265;
  double n2,np,nm,a,k,ds,ds2,k2,nm2,b1,b2,b,ts,tp1,tp2;
  double tp3,n22,np3,tp4,tp5,tp;
  
  NumericVector f(n.size());

    if(theta==0){
      
      NumericVector  f0 = 4.0*n/((n+1.0)*(n+1.0));
      return f0;    

    } else {
      
      for(int i =0; i < n.size(); i++){
	n2 = n[i]*n[i];
	np = n2+1.0;
	nm = n2-1.0;
	a = ((n[i]+1.0)*(n[i]+1.0))/2.0;
	k = -((n2-1)*(n2-1))/4.0;
	ds = sin(theta*pi/180.0);
	ds2 = ds*ds;
	k2 = k*k;
	nm2 = nm*nm;

	if(theta!=90){
	  b1 = sqrt(pow(ds2-np/2.0,2)+k);
	} else {

	  b1 = 0.0;
	}

	//	Rcpp::Rcout << "value is " << std::endl << b1 << std::endl;

	b2 = ds2-np/2.0;
	b = b1-b2;
	ts = (k2/(6.0*b*b*b)+k/b-b/2.0)-(k2/(6.0*a*a*a)+k/a-a/2.0);
	tp1 = -2.0*n2*(b-a)/(np*np);
	tp2 = -2.0*n2*np*log(b/a)/nm2;
	tp3 = n2*((1.0/b)-(1.0/a))/2.0;
	n22 = n2*n2;
	np3 = np*np*np;
	tp4 = 16.0*n22*(n22+1.0)*log((2.0*np*b-nm2)/(2.0*np*a-nm2))/(np3*nm2);
	tp5 = 16.0*(n2*n2*n2)*(pow(2.0*np*b-nm2,-1)-pow(2.0*np*a-nm2,-1))/np3;
	tp = tp1+tp2+tp3+tp4+tp5;
	f[i] = (ts+tp)/(2.0*ds2);
	
      }
    }
    
    return f;    
} // ctav 

