 ################################################################################
## Translation of Prospect_5B into R
## Marco D. Visser
################################################################################
## Note from original code follow below
## _______________________________________________________________________
##
## prospect_5B.m (including carotenoids and brown pigments)
## version 5b (october, 20th 2009)
## _______________________________________________________________________
## for any question or request, please contact: 
##
## Jean-Baptiste FERET
## UMR-TETIS, IRSTEA Montpellier
## Maison de la Teledetection
## 500 rue Jean-Fracois Breton
## 34093 Montpellier cedex 5
## E-mail: jb.feret[at]teledetection.fr
##
## Stephane JACQUEMOUD
## Universite Paris Diderot / Institut de Physique du Globe de Paris
## 35 rue Helene Brion
## 75013 Paris, France
## E-mail: jacquemoud[at]ipgp.fr
##
## http://teledetection.ipgp.fr/prosail/
## Plant leaf reflectance and transmittance are calculated from 400 nm to
## 2500 nm (1 nm step) with the following parameters:
##
##       - N     = leaf structure parameter
##       - Cab   = chlorophyll a+b content 
##       - Car   = carotenoids content 
##       - Cbrown= brown pigments content in arbitrary units
##       - Cw    = equivalent water thickness
##       - Cm    = dry matter content 
##
## Here are some examples observed during the LOPEX'93 experiment on
## fresh (F) and dry (D) leaves :
##
## ---------------------------------------------
##                N     Cab     Cw        Cm    
## ---------------------------------------------
## min          1.000    0.0  0.004000  0.001900
## max          3.000  100.0  0.040000  0.016500
## corn (F)     1.518   58.0  0.013100  0.003662
## rice (F)     2.275   23.7  0.007500  0.005811
## clover (F)   1.875   46.7  0.010000  0.003014
## laurel (F)   2.660   74.1  0.019900  0.013520
## ---------------------------------------------
## min          1.500    0.0  0.000063  0.0019
## max          3.600  100.0  0.000900  0.0165
## bamboo (D)   2.698   70.8  0.000117  0.009327
## lettuce (D)  2.107   35.2  0.000244  0.002250
## walnut (D)   2.656   62.8  0.000263  0.006573
## chestnut (D) 1.826   47.7  0.000307  0.004305
## ---------------------------------------------
## _______________________________________________________________________


#' The PROSPECT model of leaf optical properties (version 5B)
#'
#' Note from the original code; 
#' Jacquemoud S., Baret F. (1990), PROSPECT: a model of leaf optical
#' properties spectra, Remote Sens. Environ., 34:75-91.
#' Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
#' Properties Model Separating Photosynthetic Pigments, Remote Sensing of
#' Environment, 112:3030-3043
#' The specific absorption coefficient corresponding to brown pigment is
#' provided by Frederic Baret (EMMAH, INRA Avignon, baret[at]avignon.inra.fr)
#' and used with his autorization.
#' ***********************************************************************
#' @param(N) leaf structure parameter
#' @param(Cab) chlorophyll a+b content in ug/cm2
#' @param(Car) carotenoids content in ug/cm2
#' @param(Cbrown) brown pigments content in arbitrary units
#' @param(Cw) equivalent water thickness in g/cm2 or cm
#' @param(Cm) dry matter content in g/cm2
#' @import pracma
#' @export
prospect_5BR <- function(N=1.5,Cab=60,Car=10,Cbrown=0.01,Cw=0.013,Cm=0.003){

    CoefMat <- data_prospect5 ## get reflective and absorbtion coefficient data
    l <- CoefMat[,1] # wavelength (nm)
    n <- CoefMat[,2] # refractive index
    ## specific absorption coefficient for each element at each wl (k= total absorbtion coefficient)
    k <- (Cab*CoefMat[,3]+Car*CoefMat[,4]+Cbrown*CoefMat[,5]+Cw*CoefMat[,6]+Cm*CoefMat[,7])/N;
#    k[which(k==0)] <- .Machine$double.eps ## precision of zero in R
#    trans <- (1-k)*exp(-k)+k^2*pracma::expint(k) ## transmittance
    trans <- (1-k)*exp(-k)+k^2*expint::expint_E1(k) ## transmittance

    ## ***********************************************************************
    ## reflectance and transmittance of one layer
    ## ***********************************************************************
    ## Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
    ## Interaction of isotropic ligth with a compact plant leaf, J. Opt.
    ## Soc. Am., 59(10):1376-1379.
    ## ***********************************************************************
    ## reflectivity and transmissivity at the interface
    ##-------------------------------------------------
    alpha <- 40
    t12 <- tav(alpha,n) ## Transmission of isotropic light from 40degrees
    tav90n <- tav(90,n)
    t21 <- tav90n/n^2 ## 90degrees
    r12 <- 1-t12
    r21 <- 1-t21
    x <- t12/tav90n
    y <- x*(tav90n-1)+1-t12

    ## reflectance and transmittance of the elementary layer N = 1
    #########################################################################

    ra <- r12+(t12*t21*r21*trans^2)/(1-r21^2*trans^2)
    ta <- (t12*t21*trans)/(1-r21^2*trans^2)
    r90 <- (ra-y)/x
    t90 <- ta/x

        
    ## ***********************************************************************
    ## reflectance and transmittance of N layers
    ## ***********************************************************************
    ## Stokes G.G. (1862), On the intensity of the light reflected from
    ## or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
    ## 11:545-556.
    ## ***********************************************************************
    delta <- sqrt((t90^2-r90^2-1)^2-4*r90^2)
    beta <- (1+r90^2-t90^2-delta)/(2*r90)
    va <- (1+r90^2-t90^2+delta)/(2*r90)
    vb <- sqrt(beta*(va-r90)/(va*(beta-r90)))

    vbNN <- vb^(N-1)
    vbNNinv <- 1/vbNN
    vainv <- 1/va
    s1 <- ta*t90*(vbNN-vbNNinv)
    s2 <- ta*(va-vainv)
    s3 <- va*vbNN-vainv*vbNNinv-r90*(vbNN-vbNNinv)

    RN <- ra+s1/s3
    TN <- s2/s3


        
    LTR <- matrix(c(l,RN,TN),ncol=3)

    colnames(LTR) <- c("wl","ref","trans")

    return(LTR)
    
}



#' tav function
#' 
#' Stern F. (1964), Transmission of isotropic radiation across an
#' interface between two dielectrics, Appl. Opt., 3(1):111-113.
#' Allen W.A. (1973), Transmission of isotropic light across a
#' dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
#' 63(6):664-666.
#' @section references:
#' Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
#' Properties Model Separating Photosynthetic Pigments, Remote Sensing of
#' Environment
#' 
#' @param(teta) Angle in degrees
#' @param(ref) refractive index
#' @export
tav <- function(teta,ref) {

    s <- length(ref)
    teta <- teta*pi/180 ## degrees to radians
    r2 <- ref^2
    rp <- r2+1
    rm <- r2-1
    a <- (ref+1)^2/2
    k <- -(r2-1)^2/4
    ds <- sin(teta)

    k2 <- k^2
    rm2 <- rm^2

    if(teta==0){
        f <- 4*ref/(ref+1)^2
        return(f)
    } else if (teta==pi/2) {
        b1 <- numeric(s)
    }  else{
        b1 <- sqrt((ds^2-rp/2)^2+k)
    }

    b2 <- ds^2-rp/2
    b <- b1-b2
    ts <- (k2/(6*b^3)+k/b-b/2)-(k2/(6*a^3)+k/a-a/2)
    tp1 <- -2*r2*(b-a)/(rp^2)
    tp2 <- -2*r2*rp*log(b/a)/rm2
    tp3 <- r2*(b^(-1)-a^(-1))/2
    tp4 <- 16*r2^(2)*(r2^2+1)*log((2*rp*b-rm2)/(2*rp*a-rm2))/(rp^(3)*rm2)
    tp5 <- 16*r2^(3)*((2*rp*b-rm2)^(-1)-(2*rp*a-rm2)^(-1))/rp^3
    tp <- tp1+tp2+tp3+tp4+tp5
    f <- (ts+tp)/(2*ds^2)

    return(f)

}


