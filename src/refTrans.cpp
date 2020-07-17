#include <Rcpp.h>

using namespace Rcpp;


// J1 function with avoidance of singularity problem
// [[Rcpp::export]]
double cJfunc1(double k,
	      double l,
	      double t){
  
  double Jout, del;

  del = (k-l)*t;
    
  if(fabs(del)>1e3) {
      
    Jout = ((exp(-l*t)-exp(-k*t))/(k-l)); // for all del>1e-3
  } else {  
      
    Jout = 0.5*t*(exp(-k*t)+exp(-l*t))*(1.0-del*del/12.0);
	
  }
    
    
  return Jout;
  
} //Jfunc1


// J2 function
// some version of the foursail code have an alternative j3 function 
// [[Rcpp::export]]
double cJfunc2or3(double k,
		 double l,
		 double t) {

  double  Jout;
  Jout = (1-exp(-(k+l)*t))/(k+l);
  return Jout;

} // Jfunc2 or 3



// J4 function for treating (near) conservative scattering
// used in foursail2
// [[Rcpp::export]]
double cJfunc4(double m, double t){
  
  double del =m*t;
  double e, Jout;
  
  if(del>1e-3){

      // sub1
    e = exp(-del);
    Jout = (1.0-e)/(m*(1.0+e));
    
  } else {

    // sub2
    e = m*0;
    Jout = 0.5*t*(1.0-del*del/12.0);
      
      }
  
  
  return Jout;
  
} // Jfunc4


// Single layer reflection and transmission function 
// 4sail
// \code{cReflTrans} calculates
// @param rho numeric vector of particle reflectance
// @param tau numeric vector of particle transmission
// @param lai leaf area index
// @param att to be documented
// @param m to be documented
// @param sigb to be documented
// @param ks to be documented
// @param ko to be documented
// @param sf to be documented
// @param sb to be documented
// @param vf to be documented
// @param vb to be documented
// 
// @author Marco D. Visser
// @return Returns a list of reflectances
// @export
// [[Rcpp::export]]
List cReflTrans(NumericVector rho,
			 NumericVector tau,
			 double lai,
			 NumericVector att,
			 NumericVector m,
			 NumericVector sigb,
			 double ks,
			 double ko,
			 NumericVector sf,
			 NumericVector sb,
			 NumericVector vf,
			 NumericVector vb){

  // result vectors
  NumericVector rdd(rho.size());
  NumericVector tdd(rho.size());
  NumericVector tsd(rho.size());
  NumericVector rsd(rho.size());
  NumericVector tdo(rho.size());
  NumericVector rdo(rho.size());
  NumericVector rsod(rho.size());

  // result doubles
  double tss   = exp(-ks*lai);
  double too   = exp(-ko*lai);
  double z     = cJfunc2or3(ks,ko,lai);

  
  // return list
  List ret;

  // helper doubles
  double e1,e2,rinf,rinf2,re,denom;
  double J1ks,J2ks,J1ko,J2ko;
  double Ps,Qs,Pv,Qv;
  double g1,g2;
  double Tv1,Tv2,T1,T2,T3;
    
  //  NumericVector xx = exp( x );

  for(int i =0; i < rho.size(); i++){

    e1    = exp(-m[i]*lai);
    e2    = e1*e1;
    rinf  = (att[i]-m[i])/sigb[i];
    rinf2 = rinf*rinf;
    re    = rinf*e1;
    denom = 1-rinf2*e2;
    	   
    // Jfunctions sun and observer
    J1ks = cJfunc1(ks,m[i],lai);
    J2ks = cJfunc2or3(ks,m[i],lai);
    J1ko = cJfunc1(ko,m[i],lai);
    J2ko = cJfunc2or3(ko,m[i],lai);

    Ps    = (sf[i]+sb[i]*rinf)*J1ks;
    Qs    = (sf[i]*rinf+sb[i])*J2ks;
    Pv    = (vf[i]+vb[i]*rinf)*J1ko;
    Qv    = (vf[i]*rinf+vb[i])*J2ko;

    rdd[i]   = rinf*(1-e2)/denom;
    tdd[i]   = (1-rinf2)*e1/denom;
    tsd[i]   = (Ps-re*Qs)/denom;
    rsd[i]   = (Qs-re*Ps)/denom;
    tdo[i]   = (Pv-re*Qv)/denom;
    rdo[i]   = (Qv-re*Pv)/denom;

    g1    = (z-J1ks*too)/(ko+m[i]);
    g2    = (z-J1ko*tss)/(ks+m[i]);

    Tv1   = (vf[i]*rinf+vb[i])*g1;
    Tv2   = (vf[i]+vb[i]*rinf)*g2;
    T1    = Tv1*(sf[i]+sb[i]*rinf);
    T2    = Tv2*(sf[i]*rinf+sb[i]);
    T3    = (rdo[i]*Qs+tdo[i]*Ps)*rinf;

    // Multiple scattering contribution to bidirectional canopy reflectance
    rsod[i] = (T1+T2-T3)/(1-rinf2);
      }


  return Rcpp::List::create(rdd,tdd,tsd,rsd,tdo,rdo,tss,too,rsod);

} // ReflTrans






// Single layer for two layer prospect
// 4sail2
// [[Rcpp::export]]
List cReflTransSingleLayer(NumericVector rho, NumericVector tau,
			   double lai,
			   double ks, double ko,
			   double sdf,double sdb, double dof, double dob,
			   double sob,double sof,double  ddb,double  ddf){

  //output
  NumericVector rdd(rho.size());
  NumericVector rsd(rho.size());
  NumericVector rdo(rho.size());
  NumericVector rsod(rho.size());
  NumericVector tdd(rho.size());
  NumericVector tsd(rho.size());
  NumericVector tdo(rho.size());
  NumericVector w(rho.size());

  // constants
  double tss = exp(-ks*lai);
  double too = exp(-ko*lai);

  //helper
  double sb,sf,vb,vf,sigb,sigf,att,m2,m;  
  double J1ks,J2ks,J1ko,J2ko,z,J4,e1,e2,rinf,rinf2,re,denom,Ps,Qs;
  double Pv,Qv,g1,g2,Tv1,Tv2,T1,T2,T3,amsig,apsig,rtp,rtm,dns,dno;
  double cks,cko,dks,dko,ho;
  
  for(int i =0; i < rho.size(); i++){
    
    //Here the reflectance and transmission come in
    //with the suits coefficients

    //RTgeomRes function in R
    sb = sdb*rho[i]+sdf*tau[i];
    sf = sdf*rho[i]+sdb*tau[i];
    vb = dob*rho[i]+dof*tau[i];
    vf = dof*rho[i]+dob*tau[i];
    w[i] = sob*rho[i]+sof*tau[i];
    
    sigb = ddb*rho[i]+ddf*tau[i];
    sigf = ddf*rho[i]+ddb*tau[i];
    att = 1.0-sigf;

    m2 = (att+sigb)*(att-sigb);
    m;

    if(m2<0){
      m = 0;
    } else {
      m = sqrt(m2);

    }
       
    J1ks = cJfunc1(ks,m,lai);
    J2ks = cJfunc2or3(ks,m,lai);
    J1ko = cJfunc1(ko,m,lai);
    J2ko = cJfunc2or3(ko,m,lai);
    z  = cJfunc2or3(ks,ko,lai);
    J4 = cJfunc4(m,lai);

     // condition 1
    if(m>0.01) {

      // Normal case
      e1 = exp(-m*lai);
      e2 = e1*e1;
      rinf = ((att-m)/sigb);
      rinf2 = rinf*rinf;
      re = rinf*e1;
      denom = 1.0-rinf2*e2;

      Ps  =  (sf+sb*rinf)*J1ks;
      Qs = (sf*rinf+sb)*J2ks;
      Pv = (vf+vb*rinf)*J1ko;
      Qv = (vf*rinf+vb)*J2ko;
    
      rdd[i] = rinf*(1.0-e2)/denom;
      tdd[i] = (1.0-rinf2)*e1/denom;
    
      tsd[i] = (Ps-re*Qs)/denom;
      rsd[i] = (Qs-re*Ps)/denom;
      tdo[i] = (Pv-re*Qv)/denom;
      rdo[i] = (Qv-re*Pv)/denom;
    
      g1 = (z-J1ks*too)/(ko+m);
      g2 = (z-J1ko*tss)/(ks+m);

      Tv1 = (vf*rinf+vb)*g1;
      Tv2 = (vf+vb*rinf)*g2;
      T1 = Tv1*(sf+sb*rinf);
      T2 = Tv2*(sf*rinf+sb);
      T3 = (rdo[i]*Qs+tdo[i]*Ps)*rinf;

      // Multiple scattering contribution to
      //bidirectional canopy reflectance
      
      rsod[i] = (T1+T2-T3)/(1.0-rinf2);
      
    } else {

      // else m<0.01
      // Near or complete conservative scattering
      amsig = att-sigb;
      apsig = att+sigb;
	
      rtp = (1.0-amsig*J4)/(1+amsig*J4);
      rtm = (-1.0+apsig*J4)/(1+apsig*J4);
      rdd[i] = 5.0*(rtp+rtm);
      tdd[i] = 5.0*(rtp-rtm);
	
      dns = ks*ks-m*m;
      dno = ko*ko-m*m;
      cks = (sb*(ks-att)-sf*sigb)/dns;
      cko = (vb*(ko-att)-vf*sigb)/dno;
      dks = (-sf*(ks+att)-sb*sigb)/dns;
      dko = (-vf*(ko+att)-vb*sigb)/dno;
      ho  = (sf*cko+sb*dko)/(ko+ks);
	
      rsd[i] = cks*(1.0-tss*tdd[i])-dks*rdd[i];
      rdo[i] = cko*(1.0-too*tdd[i])-dko*rdd[i];
      tsd[i] = dks*(tss-tdd[i])-cks*tss *rdd[i];
      tdo[i] = dko*(too -tdd[i])-cko*too*rdd[i];
	
      //  Multiple scattering contribution to
      // bidirectional canopy reflectance
      rsod[i] = ho*(1.0-tss*too)-cko*tsd[i]*too-dko*rsd[i];
    }
  
    //    Rcpp::Rcout << "value is " << std::endl << rsod << std::endl;

    
  }
  
  return Rcpp::List::create(rdd,rsd,rdo,rsod,tdd,tsd,tdo,too,tss,w);
  
} //  ReflTransSingleLayer


