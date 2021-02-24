#' Leaf inclination distribution models
#' s3 method for calling leaf models
#' @param pars a parameter vector with a class 
#' lidf.[modelnumber]. Models include:
#'  \itemize{
#' \item [1] = Dlagden distribution (1, lidf.1)
#' \item [2] = Ellipsoid (Campebll) distribution (2, lidf.2)
#' \item [3] = Beta distribution (3, lidf.3)
#' \item [4] = One parameter beta distribution (4, lidf.4)
#' }
#'  Models 1 and 2 are the standard models from the SAIL model
#' @return a vector of proportions for each leaf angle calculated 
#' from each leaf inclination model  
#' @export
lidf <- function(pars) {
    
    lidf <- UseMethod("lidf",pars)
    return(lidf)

} # LIDF_fun

################################################################################
##                          Campbell                            
##     
##    Computation of the leaf angle distribution function value (freq) 
##    Ellipsoidal distribution function caracterised by the average leaf 
##    inclination angle in degree (ala)                                     
##    Campbell 1986                                                      
##                                                                              
################################################################################
##  old pure R version kept here for testing 
r_lidf.lidf.2 <- function(param){

    na <- param[1]
    ala <- param[2]

    ## 
    tx2 <- as.double(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 82, 84, 86, 88))
    tx1 <- as.double(c(10, 20, 30, 40, 50, 60, 70, 80, 82, 84, 86, 88, 90))

    tl1 <- tx1*pi/180
    tl2 <- tx2*pi/180

    excent <- exp(-1.6184e-5*ala^3+2.1145e-3*ala^2-1.2390e-1*ala+3.2491)
   
    x1   <-  excent/(sqrt(1+excent^2*tan(tl1)^2))
    x2   <-  excent/(sqrt(1+excent^2*tan(tl2)^2))
    
    if(excent==1) {
        freq  <-  abs(cos(tl1)-cos(tl2))
    } else {
        alpha   <-  excent/sqrt(abs(1-excent^2))
        alpha2  <-  alpha^2
        x12  <-  x1^2
        x22  <-  x2^2
        
    }
    
    if(excent>1){
        alpx1  <-  sqrt(alpha2+x12)
        alpx2  <-  sqrt(alpha2+x22)
        dum   <-  x1*alpx1+alpha2*log(x1+alpx1)
        freq  <-  abs(dum-(x2*alpx2+alpha2*log(x2+alpx2)))
        
    } else {
        almx1  <-  sqrt(alpha2-x12)
        almx2  <-  sqrt(alpha2-x22)
        dum    <- x1*almx1+alpha2*asin(x1/alpha)
        freq <-  abs(dum-(x2*almx2+alpha2*asin(x2/alpha)))
       }
    
   
    finalfreq <- freq/sum(freq)
    
    return(finalfreq)
} #campbell

################################################################################
##                          Campbell                            
##     
##    Computation of the leaf angle distribution function value (freq) 
##    Ellipsoidal distribution function caracterised by the average leaf 
##    inclination angle in degree (ala)                                     
##    Campbell 1986                                                      
##                                                                              
################################################################################
lidf.lidf.2 <- function(param){

    na <- param[1]
    ala <- param[2]

    ## using c++ cambell function
    tx2 <- as.double(c(0, 10, 20, 30, 40, 50, 60, 70, 80, 82, 84, 86, 88))
    tx1 <- as.double(c(10, 20, 30, 40, 50, 60, 70, 80, 82, 84, 86, 88, 90))
    freq <- cambell(ala,tx1,tx2) 
   
    finalfreq <- freq/sum(freq)
    
    return(finalfreq)
} #campbell


## dlagden function
##  old pure R version kept here for testing 
lidf.lidf.1 <- function(param){

    na <- param[1]
    a <- param[2]
    b <- param[3]

    freq <- numeric(na)
    
    t1 <- c(1:8)*10
    t2 <- 80+(c(9:12)-8)*2
    
    freq[1:12] <- cdcum(a,b,c(t1,t2))
    freq[13] <- 1
    freq[13:2] <- freq[13:2]-freq[c(13:2)-1]
    return(freq)
} # dladgen

## dlagden function
##  old pure R version kept here for testing 
r_lidf.lidf.1 <- function(param){

    na <- param[1]
    a <- param[2]
    b <- param[3]

    freq <- numeric(na)
    
    t <- c(1:8)*10
    freq[1:8] <- sapply(t,function(X) dcum(a,b,X))
  

    t <- 80+(c(9:12)-8)*2
    freq[9:12] <- sapply(t,function(X) dcum(a,b,X))

    freq[13] <- 1

    freq[13:2] <- freq[13:2]-freq[c(13:2)-1]

    return(freq)
} # dladgen

## cummulative lagden function
##  old pure R version kept here for testing 
dcum <- function(a,b,t){

    rd <- pi/180
    
    if(a>1){
        
        dcum <- 1-cos(rd*t)
    } else{
        
        eps <- 1e-8
        delx <- 1
        
        x <- 2*rd*t
        p <- x
        
        while(delx>eps){
            y  <-  a*sin(x)+.5*b*sin(2*x)
            dx <- .5*(y-x+p) 
            x <- x+dx
            delx <- abs(dx)
        }
        
        dcum <- (2*y+p)/pi
        
    }
    
    return(dcum)
} #dcum

## the beta leaf inclination function
## @param a number of success parameter
## @param b number of failures parameter
lidf.lidf.3 <- function(param){

    na <- param[1]
    a <- param[2] * param[3]
    b <- (1-param[2]) * param[3]
    
    freq <- numeric(na)

    t1 <- c(0:7)*10
    t2 <- c(1:8)*10
    
    freq[1:8] <- pbeta(t2/90,a,b)-pbeta(t1/90,a,b)
  
    t1 <- 80+(c(8:12)-8)*2
    t2 <- 80+(c(9:13)-8)*2
    freq[9:13] <- pbeta(t2/90,a,b)-pbeta(t1/90,a,b)

    return(freq)
} # dbeta


## one parameter beta
lidf.lidf.4<-function(param){

    na <- param[1]
    theta <- param[2]

    freq <- numeric(na)

    t1 <- c(0:7)*10
    t2 <- c(1:8)*10

    ## density function = theta*x^(theta-1)
    in2 <- (t2/90)^(theta)
    in1 <- (t1/90)^(theta)
    freq[1:8] <- in2-in1
  
    t1 <- 80+(c(8:12)-8)*2
    t2 <- 80+(c(9:13)-8)*2
    in2 <- (t2/90)^(theta)
    in1 <- (t1/90)^(theta)
   
    freq[9:13] <- in2-in1

 
 return(freq)
} # dbeta1
    


