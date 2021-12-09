#include <Rcpp.h>

using namespace Rcpp;

//' Leaf inclination distribution function
//' cummulative lagden function from Wout Verhoef's dissertation
//' Extended here for any angle  
//'
//' @param a parameter
//' @param b parameter
//' @param theta angle in degrees
//' @return angle fraction value
//'
//' @export
// [[Rcpp::export]]
NumericVector cdcum(double a, double b, NumericVector theta) {

  NumericVector rad(theta.size());
  NumericVector dcum(theta.size());
  double eps = 1e-6;
  double x,dx,y;
  static const double pi = 3.14159265;

  if(a>1) {
    // Rcpp suger for the special case 
    rad = theta * (pi / 180.0); // theta to radians
    dcum=1-cos(rad);
  } else {
  
    for(int i =0; i < theta.size(); i++){
    
      rad[i] = theta[i] * (pi / 180.0);
      x = 2.0 * rad[i]; // to radians
      dx = 1.0;
      while (dx > eps) {
	y = a * sin(x) + 0.5 * b * sin(2.0 * x);
	dx = 0.5 * (y - x + (2.0 * rad[i]));
	x = x + dx;
	dx = fabs(dx);
      }
    
      dcum[i] = 2.0 * (y + rad[i]) / pi;
    }
  }

  return dcum;
}  
  

//##############################################################################
//                          Campbell                            
//     
//    Computation of the leaf angle distribution function value (freq) 
//    Ellipsoidal distribution function caracterised by the average leaf 
//    inclination angle in degree (ala)                                     
//    Campbell 1986                                                      
//                                                                              
//##############################################################################
//' Leaf inclination distribution function
//' Ellipsoidal distribution function
//'
//' @param ala average leaf angle parameter
//' @param tx1 angle in degree
//' @param tx2 angle in degree
//' @return angle fraction value
//'
//' @export
// [[Rcpp::export]]
NumericVector cambell(double ala, NumericVector tx1, NumericVector tx2){ 

  static const double pi = 3.14159265;
  double excent = exp(-1.6184e-5*(ala*ala*ala)+2.1145e-3*(ala*ala)-1.2390e-1*ala+3.2491);
  double excent2 = excent*excent;
  double alpha = excent/sqrt(fabs(1-excent2));
  double alpha2 = alpha * alpha;
  double tl1,tl2,x1,x2,x12,x22,alpx1,alpx2,almx1,almx2,dum1,dum2;


  NumericVector freq(tx1.size());

    for(int i =0; i < tx1.size(); i++){
        
      tl1 = tx1[i]*pi/180.0;
      tl2 = tx2[i]*pi/180.0;
      x1 =  excent/(sqrt(1+excent2*(tan(tl1)*tan(tl1))));
      x2 =  excent/(sqrt(1+excent2*(tan(tl2)*tan(tl2))));
      
      if(excent==1) {
	
        freq[i]  =  fabs(cos(tl1)-cos(tl2));
	
	  } else {
	
       x12  =  x1*x1;
       x22  =  x2*x2;
       
       if(excent>1){
	 
	 alpx1  =  sqrt(alpha2+x12);
	 alpx2  =  sqrt(alpha2+x22);
	 dum1   =  x1*alpx1+alpha2*log(x1+alpx1);
	 dum2   =  x2*alpx2+alpha2*log(x2+alpx2);
	 freq[i]  =  fabs(dum1-dum2);
	 //	 Rcpp::Rcout << "value is " << std::endl << dum1-dum2 << std::endl;     
	 
       } else {
	 
	 almx1  =  sqrt(alpha2-x12);
	 almx2  =  sqrt(alpha2-x22);
	 dum1   = x1*almx1+alpha2*asin(x1/alpha);
	 dum2   = x2*almx2+alpha2*asin(x2/alpha);
	 freq[i] =  fabs(dum1-dum2);
       }

      }
    }
    
    return freq;
}


    
