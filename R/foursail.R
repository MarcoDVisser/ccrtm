#' Optimized R implementation of foursail (4SAIL) 
#'
#' The foursail (or 4SAIL) radiative transfer model is commonly used to simulate bidirectional 
#' reflectance distribution functions within vegetation canopies. Foursail (4SAIL) refers 
#' to "Scattering by Arbitrary Inclined Leaves" in a 4-stream model. The four-streams represents 
#' the scattering and absorption of upward, downward and two directional radiative fluxes with 
#' four linear differential equations in a 1-D canopy. The model was initially developed by
#' Verhoef (1984), who extended work by Suits (1971) 4-steam model.  
#' 
#' @param rho input leaf reflectance from 400-2500nm (can be measured or modeled)
#' @param tau input leaf transmittance from 400-2500nm (can be measured or modeled)
#' @param bgr background reflectance. Usual input is
#' soil reflectance spectra from 400-2500nm (can be measured or modeled)
#' @param param A named vector of SAIL parameter values (note: program ignores case):
#'  \itemize{
#' \item [1] = Leaf angle distribution function parameter a (LIDFa)
#' \item [2] = Leaf angle distribution function parameter b (LIDFb)
#' \item [3] = Leaf angle distribution function type (see ?lidfFun)
#' \item [4] = Leaf area index (LAI)
#' \item [5] = Hot spot effect parameter (hspot)
#' \item [6] = Solar zenith angle (tts)
#' \item [7] = Observer zenith angle (tto)
#' \item [8] = Sun-sensor azimuth angle (psi)
#' }
#' 
#' @return spectra matrixwith 4 reflectance factors and canopy transmission 
#' for wavelengths 400 to 2500nm:
#'  \itemize{
#' \item [1] = bi-hemispherical reflectance (rddt). White-sky albedo: the reflectance of the canopy
#' under diffuse illumination. The BRDF integrated over all viewing and illumination directions.  
#' \item [2] = hemispherical directional reflectance (rsdt). Black-sky albedo: reflectance of a surface
#' under direct light without a diffuse component. It is the integral of the BRDF over all viewing
#' directions.
#' \item [3] = directional hemispherical reflectance (rdot). Diffuse reflectance in the vieweing
#' direction. 
#' \item [4] = bi-directional reflectance (rsot). The ratio of reflected radiance in the viewing direction
#' to the incoming radiant flux in the solar direction. 
#' \item [5] = Canopy transmission of diffuse light through the canopy (taud).
#' \item [6] = transmission of direct light through the canopy (taus).
#' }
#'
#' @examples
#' ## lower-level implementation example
#' ## see ?fRTM for the typical mode of simulation
#' ## e.g. fRTM(rho~prospectd+foursail) 
#'
#' ## 1) get parameters
#' params<-getDefaults(rho~prospectd+foursail) 
#' ## getDefaults("foursail") will also work
#' bestpars<-params$foursail$best
#' ## ensure the vector is named
#' names(bestpars) <- rownames(params$foursail)
#' 
#' ## 2) get leaf reflectance and transmission 
#' rt<-fRTM(rho+tau~prospectd)
#'
#' ## 3) get soil reflectance to model background reflectance
#' data(soil)
#' 
#' ## a linear mixture soil model 
#' bgRef<- bestpars["psoil"]*soil[,"drySoil"] + (1-bestpars["psoil"])*soil[,"wetSoil"]
#' 
#' ## 4) run 4SAIL
#' foursail(rt[,"rho"],rt[,"tau"],bgRef,bestpars)
#' 
#' @references Suits, G.H., 1971. The calculation of the directional reflectance of a 
#'  vegetative canopy. Remote Sens. Environ. 2, 117-125.
#' @references Verhoef, W. (1984). Light scattering by leaf layers with application to 
#'   canopy reflectance modeling: The SAIL model. Remote Sens. Environ. 16, 125-141.
#' @export
foursail <- function(rho, tau, bgr,param){

    names(param) <- tolower(names(param))
    
    ## input paramters
    LIDFa <-    as.numeric(param["lidfa"])
    LIDFb <-    as.numeric(param["lidfb"])
    TypeLIDF <- as.numeric(param["typelidf"])
    lai <- as.numeric(param["lai"])
    q <- as.numeric(param["hspot"])
    tts <- as.numeric(param["tts"])
    tto <- as.numeric(param["tto"])
    psi <- as.numeric(param["psi"])
    
    
    ## leaf angles
    litab <- c(5,15,25,35,45,55,65,75,81,83,85,87,89)
    ## number of angles
    na <- length(litab)

    ## degrees to radians
    rd <-  pi / 180
    
    ##  Sensor geometry: Compute geometric quantities
    ## output = list(cts,cto,ctscto,dso)
    sensGeom <- sail_sensgeom(tts,tto,psi,rd)
    cts <- sensGeom[[1]]
    cto <- sensGeom[[2]]
    ctscto <- sensGeom[[3]]
    dso <- sensGeom[[4]]
    
    ##  Leaf angle distribution function
    lidfpar <- c(na,LIDFa, LIDFb)
    names(lidfpar) <- c("na","a","b")
    class(lidfpar) <- paste0("lidf.",TypeLIDF)
    lidFun <- lidf(lidfpar)

    ## Angular distance - shadow length compensation
    ## output = list(ks,ko,sob,sof,sdb,sdf,dob,dof,ddb,ddf)
    suitRes <- SUITS(na,litab,lidFun,tts,tto,cts,cto,psi,ctscto)

    ks  <- suitRes[[1]]
    ko  <- suitRes[[2]] 
    sob <- suitRes[[3]] 
    sof <- suitRes[[4]] 
    sdb <- suitRes[[5]]  
    sdf <- suitRes[[6]]  
    dob <- suitRes[[7]]  
    dof <- suitRes[[8]]  
    ddb <- suitRes[[9]]  
    ddf <- suitRes[[10]] 
    
    ## if 3 or leaf refl-- Geometric leaf reflectance
    ## Input rho and tau from PROSPECT/ other particle model
    ## Correct using SUITS factors
    RTgeomRes <- RTgeom(rho,tau,ddb,ddf,sdb,sdf,dob,dof,sob,sof)

    sigb <- RTgeomRes[[1]]
    att  <- RTgeomRes[[2]]
    m    <- RTgeomRes[[3]]
    sb   <- RTgeomRes[[4]]
    sf   <- RTgeomRes[[5]]
    vb   <- RTgeomRes[[6]]
    vf   <- RTgeomRes[[7]]
    w    <- RTgeomRes[[8]]

    ##  Canopy reflectance and transmittance
    ## Input LAI and leaf rho and tau
    ## Also give geometric factors
    refTransRes <- cReflTrans(rho,tau,lai,att,m,sigb,ks,ko,sf,sb,vf,vb)

    rdd  <- refTransRes[[1]]
    tdd  <- refTransRes[[2]]
    tsd  <- refTransRes[[3]]
    rsd  <- refTransRes[[4]]
    tdo  <- refTransRes[[5]]
    rdo  <- refTransRes[[6]]
    tss  <- refTransRes[[7]]
    too  <- refTransRes[[8]]
    rsod <- refTransRes[[9]]


    ##  Hot spot effect
    hspotRes <- HotSpot(lai,q,tss,ks,ko,dso)
    tsstoo <- hspotRes[[1]]
    sumint <- hspotRes[[2]] 

     ## Bidirectional reflectance
     finalOut <- sail_BDRF(w,lai,sumint,tsstoo,bgr,
                           rdd,tdd,tsd,rsd,tdo,rdo,tss,too,rsod)

    
  ## rddt    Bi-hemispherical reflectance
  ## rsdt    Directional-hemispherical reflectance for solar incident flux
  ## rdot    Hemispherical-directional reflectance in viewing direction
  ## rsot    Bi-directional reflectance factor

    reflMat <- as.matrix(do.call(cbind,finalOut))
    colnames(reflMat) <- c('rddt','rsdt','rdot','rsot','tau')
   
    return(reflMat)

} # foursail



#' R implementation of foursail (pure R)
#'
#' The pure R version of foursail is included in the package as an easy
#' way to review the code, and to check more optimized versions against.
#' Model originally developed by Wout Verhoef.
#' 
#' @param rho input leaf reflectance from 400-2500nm (can be measured or modeled)
#' @param tau input leaf transmittance from 400-2500nm (can be measured or modeled)
#' @param bgr background reflectance. Usual input is
#' soil reflectance spectra from 400-2500nm (can be measured or modeled)
#' @param param A named vector of SAIL parameter values (note: program ignores case):
#'  \itemize{
#' \item [1] = Leaf angle distribution function parameter a (LIDFa)
#' \item [2] = Leaf angle distribution function parameter b (LIDFb)
#' \item [3] = Leaf angle distribution function type (see ?lidfFun)
#' \item [4] = Leaf area index (LAI)
#' \item [5] = Hot spot effect parameter (hspot) - The foliage hot spot size
#'             parameter is equal to the ratio of the correlation length of leaf
#'             projections in the horizontal plane and the canopy height (Verhoef & Bach 2007).
#' \item [6] = Solar zenith angle (tts)
#' \item [7] = Observer zenith angle (tto)
#' \item [8] = Sun-sensor azimuth angle (psi)
#' }
#' 
#' @return spectra matrixwith 4 reflectance factors and canopy transmission 
#' for wavelengths 400 to 2500nm:
#'  \itemize{
#' \item [1] = bi-hemispherical reflectance (rddt). White-sky albedo: the reflectance of the canopy
#' under diffuse illumination. The BRDF integrated over all viewing and illumination directions.  
#' \item [2] = hemispherical directional reflectance (rsdt). Black-sky albedo: reflectance of a surface
#' under direct light without a diffuse component. It is the integral of the BRDF over all viewing
#' directions.
#' \item [3] = directional hemispherical reflectance (rdot). Diffuse reflectance in the vieweing
#' direction. 
#' \item [4] = bi-directional reflectance (rsot). The ratio of reflected radiance in the viewing direction
#' to the incoming radiant flux in the solar direction. 
#' \item [5] = Canopy transmission of diffuse light through the canopy (taud).
#' \item [6] = transmission of direct light through the canopy (taus).
#' }
#'
#' @author   Marco D. Visser (R implementation) 
#'
r_foursail <- function(rho, tau,bgr,param){

    names(param) <- tolower(names(param))
    
    ## input paramters
    LIDFa <-    as.numeric(param["lidfa"])
    LIDFb <-    as.numeric(param["lidfb"])
    TypeLIDF <- as.numeric(param["typelidf"])
    lai <- as.numeric(param["lai"])
    q <- as.numeric(param["hspot"])
    tts <- as.numeric(param["tts"])
    tto <- as.numeric(param["tto"])
    psi <- as.numeric(param["psi"])
    
    
    ## leaf angles
    litab <- c(5,15,25,35,45,55,65,75,81,83,85,87,89)
    ## number of angles
    na <- length(litab)

    ## degrees to radians
    rd <-  pi / 180
    
    ##  Sensor geometry: Compute geometric quantities
    ## output = list(cts,cto,ctscto,dso)
    sensGeom <- sail_sensgeom(tts,tto,psi,rd)
    cts <- sensGeom[[1]]
    cto <- sensGeom[[2]]
    ctscto <- sensGeom[[3]]
    dso <- sensGeom[[4]]
    
    ##  Leaf angle distribution function
    lidfpar <- c(na,LIDFa, LIDFb)
    names(lidfpar) <- c("na","a","b")
    class(lidfpar) <- paste0("lidf.",TypeLIDF)
    lidFun <- lidf(lidfpar)

    ## Angular distance - shadow length compensation
    ## output = list(ks,ko,sob,sof,sdb,sdf,dob,dof,ddb,ddf)
    suitRes <- SUITS(na,litab,lidFun,tts,tto,cts,cto,psi,ctscto)

    ks  <- suitRes[[1]]
    ko  <- suitRes[[2]] 
    sob <- suitRes[[3]] 
    sof <- suitRes[[4]] 
    sdb <- suitRes[[5]]  
    sdf <- suitRes[[6]]  
    dob <- suitRes[[7]]  
    dof <- suitRes[[8]]  
    ddb <- suitRes[[9]]  
    ddf <- suitRes[[10]] 
    
    ## if 3 or leaf refl-- Geometric leaf reflectance
    ## Input rho and tau from PROSPECT/ other particle model
    ## Correct using SUITS factors
    RTgeomRes <- RTgeom(rho,tau,ddb,ddf,sdb,sdf,dob,dof,sob,sof)

    sigb <- RTgeomRes[[1]]
    att  <- RTgeomRes[[2]]
    m    <- RTgeomRes[[3]]
    sb   <- RTgeomRes[[4]]
    sf   <- RTgeomRes[[5]]
    vb   <- RTgeomRes[[6]]
    vf   <- RTgeomRes[[7]]
    w    <- RTgeomRes[[8]]

    ##  Canopy reflectance and transmittance
    ## Input LAI and leaf rho and tau
    ## Also give geometric factors
    refTransRes <- cReflTrans(rho,tau,lai,att,m,sigb,ks,ko,sf,sb,vf,vb)

    rdd  <- refTransRes[[1]]
    tdd  <- refTransRes[[2]]
    tsd  <- refTransRes[[3]]
    rsd  <- refTransRes[[4]]
    tdo  <- refTransRes[[5]]
    rdo  <- refTransRes[[6]]
    tss  <- refTransRes[[7]]
    too  <- refTransRes[[8]]
    rsod <- refTransRes[[9]]


    ##  Hot spot effect
    hspotRes <- HotSpot(lai,q,tss,ks,ko,dso)
    tsstoo <- hspotRes[[1]]
    sumint <- hspotRes[[2]] 

     ## Bidirectional reflectance
     finalOut <- sail_BDRF(w,lai,sumint,tsstoo,bgr,
                           rdd,tdd,tsd,rsd,tdo,rdo,tss,too,rsod)

    
  ## rddt    Bi-hemispherical reflectance
  ## rsdt    Directional-hemispherical reflectance for solar incident flux
  ## rdot    Hemispherical-directional reflectance in viewing direction
  ## rsot    Bi-directional reflectance factor

    reflMat <- as.matrix(do.call(cbind,finalOut))
    colnames(reflMat) <- c('rddt','rsdt','rdot','rsot','tau')
   
    return(reflMat)

} # foursail



## Sensor geometry function
sail_sensgeom <- function(tts,tto,psi,rd){

    cts     <-  cos(rd*tts)
    cto     <-  cos(rd*tto)
    ctscto  <-  cts*cto
    tants   <-  tan(rd*tts)
    tanto   <-  tan(rd*tto)
    cospsi  <-  cos(rd*psi)

    ## angular distance measure also used in the hotspot effect
    dso     <-  sqrt(tants*tants+tanto*tanto-2*tants*tanto*cospsi)
        
    return(list(cts,cto,ctscto,dso))

} # sail_sensgeom



## SUITS function to calculate SUITS coefficients
SUITS <- function(na,litab,lidv,tts,tto,cts,cto,psi,ctscto){

    ## Calculate geometric factors associated with extinction and scattering 
    ##	Initialise sums

    ## R FLAG 1 : THIS LIKELY DOESNT NEED TO BE A LOOP!!
    rd <- pi/180
    sof <- sob <- bf <- ko <- ks <- 0

    ##	Weighted sums over LIDF
    for(i in 1:na) { ## loop over leaf inclination discrete values
        
        ttl <- litab[i]     
        ctl = cos(rd*ttl)
        ## SAIL volume scattering phase function gives interception and portions to be 
        ##  multiplied by rho and tau

        volscatRes <- volscatt(tts,tto,psi,ttl,chi_s,chi_o,frho,ftau)

        chi_s <- volscatRes[[1]]
        chi_o <- volscatRes[[2]]
        frho  <- volscatRes[[3]]
        ftau  <- volscatRes[[4]]
        
    ############################################################################
    ##                   SUITS SYSTEM COEFFICIENTS 
    ##
    ##	ks  : Extinction coefficient for direct solar flux
    ##	ko  : Extinction coefficient for direct observed flux
    ##	att : Attenuation coefficient for diffuse flux
    ##	sigb : Backscattering coefficient of the diffuse downward flux
    ##	sigf : Forwardscattering coefficient of the diffuse upward flux
    ##	sf  : Scattering coefficient of the direct solar flux for downward
    ##  diffuse flux
    ##	sb  : Scattering coefficient of the direct solar flux for upward diffuse
    ##  flux
    ##	vf   : Scattering coefficient of upward diffuse flux in the observed
    ##  direction
    ##	vb   : Scattering coefficient of downward diffuse flux in the observed
    ##  direction
    ##	w   : Bidirectional scattering coefficient
    ############################################################################
    
        ## Extinction coefficients
        ksli <- chi_s/cts
        koli <- chi_o/cto

        ## Area scattering coefficient fractions
        sobli <- frho*pi/ctscto
        sofli <- ftau*pi/ctscto
        bfli <- ctl*ctl
        ks    <- ks+ksli*lidv[i]
        ko    <- ko+koli*lidv[i]
        bf    <- bf+bfli*lidv[i]
        sob   <- sob+sobli*lidv[i]
        sof   <- sof+sofli*lidv[i]

    }

    ## Geometric factors to be used later with rho and tau
    sdb <- 0.5*(ks+bf)
    sdf <- 0.5*(ks-bf)
    dob <- 0.5*(ko+bf)
    dof <- 0.5*(ko-bf)
    ddb <- 0.5*(1+bf)
    ddf <- 0.5*(1-bf)

    return(list(ks,ko,sob,sof,sdb,sdf,dob,dof,ddb,ddf))

} # SUITS


## RTgeom: correct rho and tau using SUITS factors
RTgeom <- function(rho,tau,ddb,ddf,sdb,sdf,dob,dof,sob,sof){

    
    sigb <- ddb*rho+ddf*tau
    sigf <- ddf*rho+ddb*tau
    att <- 1-sigf
    m2 <- (att+sigb)*(att-sigb)
    ##  Old WHERE (m2 < 0) statement 
    m2[m2<0] <- 0
    m <- sqrt(m2)
    sb <- sdb*rho+sdf*tau
    sf <- sdf*rho+sdb*tau
    vb <- dob*rho+dof*tau
    vf <- dof*rho+dob*tau
    w  <- sob*rho+sof*tau

    return(list(sigb,att,m,sb,sf,vb,vf,w))
} #RTgeom


## volume scattering function
volscatt <- function(tts,tto,psi,ttl,chi_s,chi_o,frho,ftau){

    ###########################################################################
    ## tts= solar zenith
    ## tto= viewing zenith
    ## psi= azimuth
    ## ttl= leaf inclination angle
    ## chi_s= interception functions
    ## chi_o= interception functions
    ## frho= function to be multiplied by leaf reflectance rho
    ## ftau= functions to be multiplied by leaf transmittance tau
    ##

    ## Compute volume scattering functions and interception coefficients
    ## for given solar zenith, viewing zenith, azimuth and
    ## leaf inclination angle.

    ## chi_s and chi_o are the interception functions.
    ## frho and ftau are the functions to be multiplied by leaf reflectance rho and
    ## leaf transmittance tau, respectively, in order to obtain the volume scattering
    ## function.

    ## Wout Verhoef, april 2001, for CROMA
    rd <- pi/180

    costs <- cos(rd*tts)
    costo <- cos(rd*tto)
    sints <- sin(rd*tts)
    sinto <- sin(rd*tto)
    cospsi <- cos(rd*psi)

    psir <- rd*psi

    costl <- cos(rd*ttl)
    sintl <- sin(rd*ttl)
    cs <- costl*costs
    co <- costl*costo
    ss <- sintl*sints
    so <- sintl*sinto

    ##.........................................................................
    ##    betas -bts- and betao -bto- computation
    ##    Transition angles (beta) for solar (betas) and view (betao) directions
    ##    if thetav+thetal>pi/2, bottom side of the leaves is observed for leaf
    ##    azimut interval betao+phi<leaf azimut<2pi-betao+phi.
    ##    if thetav+thetal<pi/2, top side of the leaves is always observed,
    ##    betao=pi  same consideration for solar direction to compute betas
    ##  .......................................................................


    cosbts <- 5
    if(abs(ss)>1e-6){
        cosbts <- -cs/ss
    }

    cosbto <- 5
    if(abs(so)>1e-6){
        cosbto <- -co/so
    }

    if(abs(cosbts)<1){ ## no need for .d0 as R is double by default
        bts <- acos(cosbts)
        ds <- ss
    } else {
        
        bts <- pi
        ds <- cs
    }

    ## sun interception function
    chi_s <- 2/pi*((bts-pi*.5)*cs+sin(bts)*ss)

    if(abs(cosbto)<1){  ## no need for .d0 as R is double by default
        bto <- acos(cosbto)
        doo <- so
    } else if(tto<90){
        bto <- pi
        doo <- co
    } else {
        bto <- 0
        doo <- -co
    }

    ## observed interception function
    chi_o <- 2/pi*((bto-pi*.5)*co+sin(bto)*so)

    ##......................................................................
    ##    Computation of auxiliary azimut angles bt1, bt2, bt3 used          
    ##  for the computation of the bidirectional scattering coefficient w   
    ##  .....................................................................


    btran1 <- abs(bts-bto)
    btran2 <- pi-abs(bts+bto-pi)

    if(psir<=btran1){
        bt1 <- psir
        bt2 <- btran1
        bt3 <- btran2
    } else {
        bt1 <- btran1
        if(psir<=btran2){
            bt2 <- psir
            bt3 <- btran2
        } else {
            bt2 <- btran2
            bt3 <- psir
        }
    }

    t1 <- 2*cs*co+ss*so*cospsi
    t2 <- 0
    
    if(bt2>0){
        t2 <- sin(bt2)*(2*ds*doo+ss*so*cos(bt1)*cos(bt3)) #
    }

    denom <- 2*pi*pi
    frho <- ((pi-bt2)*t1+t2)/denom
    ftau <- (-bt2*t1+t2)/denom

    if(frho<0){
        frho <- 0
    }

    if(ftau<0){
        ftau <- 0
    }

    return(list(chi_s,chi_o,frho,ftau))
} # volscatt



ReflTrans <- function(rho,tau,lai,att,m,sigb,ks,ko,sf,sb,vf,vb){

    e1    <- exp(-m*lai)
    e2    <- e1*e1
    rinf  <- (att-m)/sigb
    rinf2 <- rinf*rinf
    re    <- rinf*e1
    denom <- 1-rinf2*e2

    ## Jfunctions sun and observer
    J1ks <- Jfunc1(ks,m,lai)
    J2ks <- Jfunc2(ks,m,lai)
    J1ko <- Jfunc1(ko,m,lai)
    J2ko <- Jfunc2(ko,m,lai)

    Ps    <- (sf+sb*rinf)*J1ks
    Qs    <- (sf*rinf+sb)*J2ks
    Pv    <- (vf+vb*rinf)*J1ko
    Qv    <- (vf*rinf+vb)*J2ko

    rdd   <- rinf*(1-e2)/denom
    tdd   <- (1-rinf2)*e1/denom
    tsd   <- (Ps-re*Qs)/denom
    rsd   <- (Qs-re*Ps)/denom
    tdo   <- (Pv-re*Qv)/denom
    rdo   <- (Qv-re*Pv)/denom

    tss   <- exp(-ks*lai)
    too   <- exp(-ko*lai)
    z     <- Jfunc3(ks,ko,lai)
    g1    <- (z-J1ks*too)/(ko+m)
    g2    <- (z-J1ko*tss)/(ks+m)

    Tv1   <- (vf*rinf+vb)*g1
    Tv2   <- (vf+vb*rinf)*g2
    T1    <- Tv1*(sf+sb*rinf)
    T2    <- Tv2*(sf*rinf+sb)
    T3    <- (rdo*Qs+tdo*Ps)*rinf

    ## Multiple scattering contribution to bidirectional canopy reflectance
    rsod <- (T1+T2-T3)/(1-rinf2)

    return(list(rdd,tdd,tsd,rsd,tdo,rdo,tss,too,rsod))

} # ReflTrans


## J1 function with avoidance of singularity problem
Jfunc1 <- function(k,l,t){
    del <- (k-l)*t
    Jout <- ((exp(-l*t)-exp(-k*t))/(k-l)) # for all del>1e-3
    inc <- abs(del)<=1e-3
    Jout[inc] <- (0.5*t*(exp(-k*t)+exp(-l*t))*(1-del*del/12))[inc]
    return(Jout)
} # Jfunc1



## J2 function
Jfunc2 <- function(k,l,t) {

    (1-exp(-(k+l)*t))/(k+l)

} # Jfunc2



Jfunc3 <- function(k,l,t){
    
    (1-exp(-(k+l)*t))/(k+l)
 
} # Jfunc3


## J4 function for treating (near) conservative scattering
## used in foursail2
Jfunc4<-function(m,t){

    del<-m*t
    sub1<-which(del>1e-3)
    sub2<-which(del<1e-3)
    Jout  <- numeric(length(del))
    
    e<-m*0
    e[sub1]<-exp(-del[sub1])
    Jout[sub1]<-(1-e[sub1])/(m[sub1]*(1+e[sub1]))
    Jout[sub2] <- 0.5*t[sub2]*(1.-del[sub2]*del[sub2]/12)
    return(Jout)
    
} # Jfunc4


## Hotspot function
HotSpot <- function(lai,q,tss,ks,ko,dso){
    
    ## Treatment of the hotspot-effect
    
    alf <- 1e6
    
    ## Apply correction 2/(K+k) suggested by F.M. Breon
    
    if(q>0) alf <- (dso/q)*2/(ks+ko)

    if(alf>200) alf <- 200 ## inserted H. Bach 1/3/04 

    if(alf==0){
        ## The pure hotspot - no shadow
        tsstoo <- tss
        sumint <- (1-tss)/(ks*lai)
    } else {
        ## Outside the hotspot
        fhot <- lai*sqrt(ko*ks)

        ## Integrate by exponential Simpson method in 20 steps
        ## the steps are arranged according to equal partitioning
        ## of the slope of the joint probability function

        x1     <- 0
        y1     <- 0
        f1     <- 1
        ca <- exp(-1*alf)
        fint   <- (1-ca)*.05
        sumint <- 0

        for(i in 1:20){
            
        if(i<20) {
            x2 <- -log(1-i*fint)/alf
        } else  {
            x2 <- 1
        }
        
        y2     <- -(ko+ks)*lai*x2+fhot*(1-exp(-alf*x2))/alf
        f2     <- exp(y2)
        sumint <- sumint+(f2-f1)*(x2-x1)/(y2-y1)
        x1     <- x2
        y1     <- y2
        f1     <- f2
        }
        
        tsstoo <- f1
    }

    return(list(tsstoo,sumint))
} # HotSpot


#' The SAIL BDRF function
#' 
#' @param w goemeric reflectance parameter
#' @param lai leaf area index
#' @param sumint exp int
#' @param tsstoo Bi-directional gap fraction
#' @param rsoil background reflectance
#' @param rdd Bi-hemispherical reflectance over all in & outgoing directions
#' (white-sky albedo).
#' @param tdd Bi-hemispherical transmittance (diffuse transmittance in all directions)
#' @param tsd Directional hemispherical transmittance for solar flux
#' @param rsd Directional hemispherical reflectance for solar flux (diffuse solar angle)
#' @param tdo Directional hemispherical transmittance (diffuse in viewing direction)
#' @param rdo Directional hemispherical reflectance in viewing direction
#' @param tss Direct transmittance of solar flux
#' @param too Direct transmittance in viewing direction
#' @param rsod Multi scattering contribution 
#' @return spectra matrixwith 4 reflectance factors and canopy transmission 
#' for wavelengths 400 to 2500nm:
#'  \itemize{
#' \item [1] = bi-hemispherical reflectance (rddt). White-sky albedo: the reflectance of the canopy
#' under diffuse illumination. The BRDF integrated over all viewing and illumination directions.  
#' \item [2] = hemispherical directional reflectance (rsdt). Black-sky albedo: reflectance of a surface
#' under direct light without a diffuse component. It is the integral of the BRDF over all viewing
#' directions.
#' \item [3] = directional hemispherical reflectance (rdot). Diffuse reflectance in the vieweing
#' direction. 
#' \item [4] = bi-directional reflectance (rsot). The ratio of reflected radiance in the viewing direction
#' to the incoming radiant flux in the solar direction. 
#' \item [5] = Canopy transmission of diffuse light through the canopy (taud).
#' \item [6] = transmission of direct light through the canopy (taus).
#' }
sail_BDRF <- function(w,lai,sumint,tsstoo,rsoil,
        rdd,tdd,tsd,rsd,tdo,rdo,tss,too,rsod) {

    rsos  <- w*lai*sumint        # Single scattering contribution
    rso   <- rsos+rsod           # Total canopy contribution
    dn    <- 1-rsoil*rdd         # Interaction with the soil

    ## rddt: bi-hemispherical reflectance factor
    rddt <- rdd+tdd*rsoil*tdd/dn 

    ## rsdt: directional-hemispherical reflectance factor
    ## for solar incident flux
    rsdt <- rsd+(tsd+tss)*rsoil*tdd/dn

    ## rdot: hemispherical-directional reflectance factor in viewing direction
    rdot <- rdo+tdd*rsoil*(tdo+too)/dn

    ## rsot: bi-directional reflectance factor
    rsodt <- rsod+((tss+tsd)*tdo+(tsd+tss*rsoil*rdd)*too)*rsoil/dn 
    rsost <- rsos+tsstoo*rsoil
    rsot <- rsost+rsodt # total

    ## transmission of light through the canopy
    ## taus = tranmission of diffuse light from solar direction
    ## plus gap fraction in solar direction
    taus <- (tsd+tss) ## direct and diffuse transmission of light

    return(list(rddt,rsdt,rdot,rsot,taus))

} ## sail_BDRF


   
