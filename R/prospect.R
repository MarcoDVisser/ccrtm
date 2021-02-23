#' PROSPECT model version 5 and 5B
#'
#' The PROSPECT5(b) leaf reflectance model. The model was implemented based on
#' Jacquemoud and Ustin (2019), and is further described in detail in Feret et al (2008).
#' PROSPECT models use the plate models developed in 
#' Allen (1969) and Stokes (1862). Set Cbrown to 0 for prospect version 5.
#' 
#' @param param A named vector of PROSPECT parameters (note: program ignores case):
#' 
#'  \itemize{
#' \item [1] = leaf structure parameter (N)
#' \item [2] = chlorophyll a+b content in ug/cm2 (Cab)
#' \item [3] = carotenoids content in ug/cm2 (Car)
#' \item [4] = brown pigments content in arbitrary units (Cbrown)
#' \item [5] = equivalent water thickness in g/cm2  (Cw)
#' \item [6] = leaf dry matter content in g/cm2 - lma - (Cm)
#' }
#' 
#' @return spectra matrix with leaf reflectance and transmission 
#' for wavelengths 400 to 2500nm:
#'  \itemize{
#' \item [1] = leaf reflectance (rho)
#' \item [2] = leaf transmission (tau)
#' }
#'
#' @importFrom expint expint_E1
#' 
#' @export
#' @references Jacquemoud, S., and Ustin, S. (2019). Leaf optical properties. 
#'   Cambridge University Press.
#' @references Feret, J.B., Francois, C., Asner, G.P., Gitelson, A.A., 
#'   Martin, R.E., Bidel, L.P.R., Ustin, S.L., le Maire, G., Jacquemoud, S. (2008),
#'   PROSPECT-4 and 5: Advances in the leaf optical properties model separating photosynthetic 
#'   pigments. Remote Sens. Environ. 112, 3030-3043. 
#' @references Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R.
#'   (1969), Interaction of isotropic ligth with a compact plant leaf,
#'   Journal of the Optical Society of American, 59:1376-1379.
#' @references Stokes G.G. (1862), On the intensity of the light
#'   reflected from or transmitted through a pile of plates,
#'   Proceedings of the Royal Society of London, 11:545-556.
#' @useDynLib ccrtm    
prospect5 <- function(param){

    ## force case to lower
    names(param) <- tolower(names(param))
    
    ## input paramters
    N <-    as.numeric(param["n"])
    Cab <-    as.numeric(param["cab"])
    Car <- as.numeric(param["car"])
    Cbrown <- as.numeric(param["cbrown"])
    Cw <- as.numeric(param["cw"])
    Cm <- as.numeric(param["cm"])
    rhoTau <- matrix(0, 2101, 2)

    CoefMat <- data_prospect5 ## get reflective and absorbtion coefficient data
    
    l <- CoefMat[,1] # wavelength (nm)
    n <- CoefMat[,2] # refractive index
    ## specific absorption coefficient for each element at each wl (k= total absorbtion coefficient)

    alpha <- c(Cab,Car,Cbrown,Cw,Cm)/N 
    ## k <- (Cab*CoefMat[,3]+Car*CoefMat[,4]+Cbrown*CoefMat[,5]+Cw*CoefMat[,6]+Cm*CoefMat[,7])/N;

    k <- CoefMat[,-c(1:2)]%*%alpha ## 1% faster
    ##    k[which(k==0)] <- .Machine$double.eps ## precision of zero in R

    ## using expint::expint_E1 instead of pracma::expint due to performance
    ## test
    trans <- (1-k)*exp(-k)+k^2*expint::expint_E1(k) ## transmittance

    ## the plate model 

    ## reflectance and transmittance though interface, one layer and N layers
    alpha <- 40
    t12 <- ctav(alpha,n) ## Transmission of isotropic light from 40degrees
    tav90n <- ctav(90,n)
    t21 <- tav90n/n^2 ## 90degrees
    r12 <- 1-t12
    r21 <- 1-t21
    x <- t12/tav90n
    y <- x*(tav90n-1)+1-t12

    pm <-  cplateModel(r12,t12,r21,t21, x, y, trans, N)

    rhoTau[,1] <- pm[[1]]
    rhoTau[,2] <- pm[[2]]
    colnames(rhoTau) <- c("rho", "tau")
    return(rhoTau)
}

#' PROSPECT model version D
#'
#' The PROSPECTD leaf reflectance model. The model was implemented based on
#' Jacquemoud and Ustin (2019), and is further described in detail in Feret et al (2017).
#' PROSPECT models use the plate models developed in 
#' Allen (1969) and Stokes (1862).
#'
#' @param param A named vector of PROSPECT parameters (note: program ignores case):
#'  \itemize{
#' \item [1] = leaf structure parameter (N)
#' \item [2] = chlorophyll a+b content in ug/cm2 (Cab)
#' \item [3] = carotenoids content in ug/cm2 (Car)
#' \item [4] = Leaf anthocyanin content (ug/cm2) (Canth)
#' \item [5] = brown pigments content in arbitrary units (Cbrown)
#' \item [6] = equivalent water thickness in g/cm2  (Cw)
#' \item [7] = leaf dry matter content in g/cm2 - lma - (Cm)
#' }
#' 
#' @return spectra matrix with leaf reflectance and transmission 
#' for wavelengths 400 to 2500nm:
#'  \itemize{
#' \item [1] = leaf reflectance (rho)
#' \item [2] = leaf transmission (tau)
#' }
#' @references Jacquemoud, S., and Ustin, S. (2019). Leaf optical properties. 
#'   Cambridge University Press.
#' @references Feret, J.B., Gitelson, A.A., Noble, S.D., Jacquemoud, S. (2017).
#'   PROSPECT-D: Towards modeling leaf optical properties through a complete lifecycle. 
#'   Remote Sens. Environ. 193, 204-215.  
#' @references Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R.
#'   (1969), Interaction of isotropic ligth with a compact plant leaf,
#'   Journal of the Optical Society of American, 59:1376-1379.
#' @references Stokes G.G. (1862), On the intensity of the light
#'   reflected from or transmitted through a pile of plates,
#'   Proceedings of the Royal Society of London, 11:545-556.
#' @importFrom expint expint_E1
#' @export
#' @useDynLib ccrtm    
prospectd <- function(param){
 
    ## force case to lower
    names(param) <- tolower(names(param))
    
    ## input paramters
    N <-    as.numeric(param["n"])
    Cab <-    as.numeric(param["cab"])
    Car <- as.numeric(param["car"])
    Canth <- as.numeric(param["canth"])
    Cbrown <- as.numeric(param["cbrown"])
    Cw <- as.numeric(param["cw"])
    Cm <- as.numeric(param["cm"])
    rhoTau <- matrix(0, 2101, 2)

#   CoefMat <- ccrtm:::data_prospectd ## get reflective and absorbtion coefficient data
    CoefMat <- data_prospectd ## get reflective and absorbtion coefficient data
     
    l <- CoefMat[,1] # wavelength (nm)
    n <- CoefMat[,2] # refractive index
    ## specific absorption coefficient for each element at each wl (k= total absorbtion coefficient)
    ## should be able to switch this to CoefMat%*%Input
    ## the formula is over N as concentrations are assumed uniform over the 
    ## leaf
    k <- (Cab*CoefMat[,3]+Car*CoefMat[,4]+Canth*CoefMat[,5]+Cbrown*CoefMat[,6]+
          Cw*CoefMat[,7]+Cm*CoefMat[,8])/N;
    
    ##    k[which(k==0)] <- .Machine$double.eps ## precision of zero in R

    ## using expint::expint_E1 instead of pracma::expint due to performance
    ## test
    trans <- (1-k)*exp(-k)+k^2*expint::expint_E1(k) ## transmittance

    ## the plate model 

    ## reflectance and transmittance though interface, one layer and N layers
    alpha <- 40
    t12 <- ctav(alpha,n) ## Transmission of isotropic light from 40degrees
    tav90n <- ctav(90,n)
    t21 <- tav90n/n^2 ## 90degrees
    r12 <- 1-t12
    r21 <- 1-t21
    x <- t12/tav90n
    y <- x*(tav90n-1)+1-t12

    pm <-  cplateModel(r12,t12,r21,t21, x, y, trans, N)

    rhoTau[,1] <- pm[[1]]
    rhoTau[,2] <- pm[[2]]

    colnames(rhoTau) <- c("rho", "tau")
    return(rhoTau)
}




# Low level prospect D implementation
.prospect5 <- function(param,
                       expected=c("n","cab","car","cbrown","cw","cm")){
 
    ## force case to lower
    names(param) <- tolower(names(param))
    
    ## input paramters
    N <-    as.numeric(param[expected[1]])
    Cab <-    as.numeric(param[expected[2]])
    Car <- as.numeric(param[expected[3]])
    Cbrown <- as.numeric(param[expected[4]])
    Cw <- as.numeric(param[expected[5]])
    Cm <- as.numeric(param[expected[6]])
    rhoTau <- matrix(0, 2101, 2)

    CoefMat <- data_prospect5 ## get reflective and absorbtion coefficient data
    n <- CoefMat[,2] # refractive index

    alpha <- c(Cab,Car,Cbrown,Cw,Cm)/N 
    ## k <- (Cab*CoefMat[,3]+Car*CoefMat[,4]+Cbrown*CoefMat[,5]+Cw*CoefMat[,6]+Cm*CoefMat[,7])/N;

    k <- CoefMat[,-c(1:2)]%*%alpha ## 1% faster

    
    ##    k[which(k==0)] <- .Machine$double.eps ## precision of zero in R

    ## using expint::expint_E1 instead of pracma::expint due to performance
    ## test
    trans <- (1-k)*exp(-k)+k^2*expint::expint_E1(k) ## transmittance

    ## the plate model 

    ## reflectance and transmittance though interface, one layer and N layers
    alpha <- 40
    t12 <- ctav(alpha,n) ## Transmission of isotropic light from 40degrees
    tav90n <- ctav(90,n)
    t21 <- tav90n/n^2 ## 90degrees
    r12 <- 1-t12
    r21 <- 1-t21
    x <- t12/tav90n
    y <- x*(tav90n-1)+1-t12

    pm <-  cplateModel(r12,t12,r21,t21, x, y, trans, N)

    rhoTau[,1] <- pm[[1]]
    rhoTau[,2] <- pm[[2]]

    colnames(rhoTau) <- c("rho", "tau")
    return(rhoTau)
}


# Low level prospect D implementation
.prospectd <- function(param,
                       expected=c("n","cab","car",
                                  "canth","cbrown","cw","cm")){
 
    ## force case to lower
    names(param) <- tolower(names(param))
    
    ## input paramters
    N <-    as.numeric(param[expected[1]])
    Cab <-    as.numeric(param[expected[2]])
    Car <- as.numeric(param[expected[3]])
    Canth <- as.numeric(param[expected[4]])
    Cbrown <- as.numeric(param[expected[5]])
    Cw <- as.numeric(param[expected[6]])
    Cm <- as.numeric(param[expected[7]])
    rhoTau <- matrix(0, 2101, 2)

    CoefMat <- data_prospectd ## get reflective and absorbtion coefficient data
    n <- CoefMat[,2] # refractive index

    ## specific absorption coefficient for each element at each wl (k= total absorbtion coefficient)
    alpha <- c(Cab,Car,Canth,Cbrown,Cw,Cm)/N 

    k <- CoefMat[,-c(1:2)]%*%alpha ## 1% faster
    
    ## using expint::expint_E1 instead of pracma::expint due to performance
    ## test
    trans <- (1-k)*exp(-k)+k^2*expint::expint_E1(k) ## transmittance

    ## the plate model 

    ## reflectance and transmittance though interface, one layer and N layers
    alpha <- 40
    t12 <- ctav(alpha,n) ## Transmission of isotropic light from 40degrees
    tav90n <- ctav(90,n)
    t21 <- tav90n/n^2 ## 90degrees
    r12 <- 1-t12
    r21 <- 1-t21
    x <- t12/tav90n
    y <- x*(tav90n-1)+1-t12

    pm <-  cplateModel(r12,t12,r21,t21, x, y, trans, N)

    rhoTau[,1] <- pm[[1]]
    rhoTau[,2] <- pm[[2]]

    colnames(rhoTau) <- c("rho", "tau")
    return(rhoTau)
}




