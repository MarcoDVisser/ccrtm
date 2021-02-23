#' R implementation of the foursail2 model with 2 canopy layers.
#'
#' The foursail2 model is a two layer implementation of the 
#' foursail model described in Verhoef and Bach (2007).
#' Layers are assumed identical in particle inclination and
#' hotspot, but may differ in the relative density and types of
#' particles (see foursail2b for a layer specific inclination angle).
#' In comparison to foursail, the background (soil),
#' can now be non-Lambertain, having it own 4-stream
#' BDRF (not implemented here but may be input by the user).
#' There are two types of particles, generalized
#' to primary and secondary (originally termed "green"
#' and "brown" particles). The  realtive abundance of
#' the secondary particle in the top canopy is regulated by
#' the dissociation paramerter.The model 4SAIL2 combines with
#' prospect, libery or procosine for the reflectance
#' and transmittance of the particles, and with the the foursail
#' or Hapke elements for the background reflectance.
#' If run alone, these require direct inputs which could be 
#' measured leaf reflectance. 
#' 
#' @param rhoA primary particle reflectance from 400-2500nm (can be measured or modeled)
#' @param tauA primary particle transmittance from 400-2500nm (can be measured or modeled)
#' @param rhoB secondary particle reflectance from 400-2500nm (can be measured or modeled)
#' @param tauB secondary particle reflectance from 400-2500nm (can be measured or modeled)
#' @param rsobgr : background bidirectional reflectance (rso)
#' @param rdobgr : background directional hemispherical reflectance (rdo)
#' @param rsdbgr : background hemispherical directional reflectance (rsd)
#' @param rddbgr : background bi-hemispherical diffuse reflectance (rdd)

#' @param bgr background reflectance. Usual input is
#' soil reflectance spectra from 400-2500nm (can be measured or modeled)
#' @param param A named vector of 4SAIL2 parameter values (note:
#' program ignores case):
#'  \itemize{
#' \item [1] = Leaf angle distribution function parameter a (LIDFa)
#' \item [2] = Leaf angle distribution function parameter b (LIDFb)
#' \item [3] = Leaf angle distribution function type (TypeLidf, see ?lidfFun)
#' \item [4] = Total Leaf Area Index (LAI), including primary and secondary 
#' particles (brown and green leafs).
#' \item [5] = fraction secondary particles ("brown leaf fraction", fb)
#' \item [6] = Canopy dissociation factor for secondary particles ("diss")
#' \item [7] = Hot spot effect parameter (hspot). Often defined as the
#' ratio of mean leaf width and canopy height.
#' \item [7] = vertical crown coverage fraction (Cv), models clumping in combination
#' with parameter zeta. 
#' \item [7] = tree shape factor (zeta), defined as the ratio of crown diameter and height.
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
#'  Diffuse reflectance for diffuse incidence.
#' \item [2] = hemispherical directional reflectance (rsdt). Black-sky albedpo: reflectance of a surface
#' under direct light without a diffuse component. It is the integral of the BRDF over all viewing
#' directions. Diffuse reflectance for direct solar incidence.
#' \item [3] = directional hemispherical reflectance (rdot). Diffuse reflectance in the vieweing
#' direction. 
#' \item [4] = bi-directional reflectance (rsot). The ratio of reflected radiance in the viewing direction
#' to the incoming radiant flux in the solar direction. 
#' }
#'
#' @examples
#' ## see ?foursail for lower-level implementations
#' fRTM(rho~prospect5+foursail2)
#'
#' @references Verhoef, W., Bach, H. (2007). Coupled soil-leaf-canopy and atmosphere radiative transfer 
#'   modeling to simulate hyperspectral multi-angular surface reflectance and TOA radiance data. 
#'   Remote Sens. Environ. 109, 166-182. 
#' @export
foursail2 <- function(rhoA,tauA,rhoB=NULL,tauB=NULL,bgr,
                      rsobgr=NULL,rdobgr=NULL,rsdbgr=NULL,rddbgr=NULL,
                      param){

    if(!is.null(rsobgr)) stop("SOIL BDRF not implemented currently")

    if(is.null(rhoB)) { ## standard if B is empty

        rhoB <- rhoA
        tauB <- tauA
   
    }
    
    names(param) <- tolower(names(param))
    
    ## input paramters
    LIDFa <-    as.numeric(param["lidfa"])
    LIDFb <-    as.numeric(param["lidfb"])
    TypeLIDF <- as.numeric(param["typelidf"])
    lai <- as.numeric(param["lai"])
    q <- as.numeric(param["hspot"])
    Cv <- as.numeric(param["cv"])
    Zeta0 <- as.numeric(param["zeta"])
    fb <- as.numeric(param["fb"])
    diss <- as.numeric(param["diss"])
    psoil <- as.numeric(param["psoil"]) ## change to background
    #skyl <- as.numeric(param["skyl"])
    tts <- as.numeric(param["tts"])
    tto <- as.numeric(param["tto"])
    psi <- as.numeric(param["psi"])

    ##  Leaf angle distribution function
    ## TO DO: build generic N angle LIA models
    na  <-  18 # number of angle
    rd <- pi/180 
    lna <- seq(0,17,length.out=18)
    thetal <- 2.5+5*lna
    ## reventing back to 4sail standard
    ## leaf angles
    lna <- c(5,15,25,35,45,55,65,75,81,83,85,87,89)
    ## number of angles
    na <- length(lna)


    lidfpar <- c(na,LIDFa, LIDFb)
    names(lidfpar) <- c("na","a","b")
    class(lidfpar) <- paste0("lidf.",TypeLIDF)
    lidFun <- lidf(lidfpar)
    
    ## special case when lai = 0
    if(lai<=0){
        ## gap fractions (direct transmittance)
	tss <-  1
        too <- 1
        tsstoo	 <- 1 # (bi-directional direct transmittance)
        
        ## reflectance and transmission
        ## hemispherical diffuse
        rdd <- 0
        tdd <- 1
        ## directional hemispheric for direct solar (s) and in viewing angle/observer (o)
        rsd <- 0
	tsd <- 0
        rdo <- 0
        tdo <- 0
        ## bi-directional in viewing angle
	rso <-  0
        rsos <- 0
        rsod <- 0
	rddt <-  bgr
        rsdt <-  bgr
        rdot <- bgr
	rsodt <-  0
        rsost <- bgr
        rsot <-  bgr
        
    }   else {

        ##  Sensor geometry: Compute geometric quantities
        ## output = list(cts,cto,ctscto,dso)
        sensGeom <- sail_sensgeom(tts,tto,psi,rd)
        cts <- sensGeom[[1]]
        cto <- sensGeom[[2]]
        ctscto <- sensGeom[[3]]
        dso <- sensGeom[[4]]
        
	## Crown and vegatation clumping effects
        ## similar to flim implementation
	Cs <- 1 # initialize 
        Co <- 1 # initialize 

	if(Cv<=1){
            Cs <- 1-(1-Cv)^(1/cts)
            Co <- 1-(1-Cv)^(1/cto)
          }

	Overlap <- 0
        
	if(Zeta0>0){
            Overlap <- min(Cs*(1-Co),Co*(1-Cs))*exp(-dso/Zeta0)
	}
        
	Fcd <- Cs*Co+Overlap
	Fcs <- (1-Cs)*Co-Overlap
	Fod <- Cs*(1-Co)-Overlap
	Fos <- (1-Cs)*(1-Co)+Overlap
	Fcdc <- 1-(1-Fcd)^(.5/cts+.5/cto)

       
	##	Part depending on diss, fb, and leaf optical properties
	##	First save the input fb as the old fb, as the following change is only artificial
	##	Better define an fb that is actually used: fbu, so that the input is not modified!

	fbu <- fb
        
	if(fb==0){ # if only green leaves
		fbu <- 0.5
		rhoB <- rhoA
		tauB <- tauA
	}
	if(fb==1){ # only brown
		fbu <- 0.5
		rhoA <- rhoB
		tauA <- tauB
	}
        
        s <- (1-diss)*fbu*(1-fbu)
        ## mixing of primary and secondary particles
	## rho1 & tau1 : green foliage
	## rho2 & tau2 : brown foliage (bottom layer)
      rho1 <- ((1-fbu-s)*rhoA+s*rhoB)/(1-fbu)
      tau1 <- ((1-fbu-s)*tauA+s*tauB)/(1-fbu)
      rho2 <- (s*rhoA+(fbu-s)*rhoB)/fbu
      tau2 <- (s*tauA+(fbu-s)*tauB)/fbu


        ## Angular distance - shadow length compensation
        ## output = list(ks,ko,sob,sof,sdb,sdf,dob,dof,ddb,ddf)
        suitRes <- SUITS(na,lna,lidFun,tts,tto,cts,cto,psi,ctscto)
        
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

        ## devide the leaf area index in two layers
	lai1 <- (1-fbu)*lai
	lai2 <- fbu*lai
        
        ##  2 layer hot spot effect
        hspotRes <- HotSpot2(lai,lai1,fbu,q,tss,ks,ko,dso)
        tsstoo <- hspotRes[[1]]
        sumint1 <- hspotRes[[2]]
        sumint2 <- hspotRes[[3]]
        
        
        ##     Calculate reflectances and transmittances for 2 layers
        ## Input LAI and leaf rho and tau
        ## and geometric factors
        
        refTransRes <- ReflTrans2(rho1,rho2,tau1,tau2,ks,ko,
                                  lai1,lai2,sdf,sdb,dof,dob,
                                  sob,sof,ddb,ddf,
                                  sumint1,sumint2)
        
        dn <- refTransRes[[1]]
        tup <- refTransRes[[2]]
        tdn <- refTransRes[[3]]
        rsdt <- refTransRes[[4]]
        rdot <- refTransRes[[5]]
        rsodt <- refTransRes[[6]]
        rsost <- refTransRes[[7]]
        rsot <- refTransRes[[8]]
        rddt_t <- refTransRes[[9]]
        rddt_b <- refTransRes[[10]]
        tsst <- refTransRes[[11]]
        toot <- refTransRes[[12]]
        tsdt <- refTransRes[[13]]
        tdot <- refTransRes[[14]]
        tddt <- refTransRes[[15]]
        
        
        ## Apply clumping effects to vegetation layer
	rddcb <- Cv*rddt_b
	rddct <- Cv*rddt_t
	tddc <- 1-Cv+Cv*tddt
	rsdc <- Cs*rsdt
	tsdc <- Cs*tsdt
	rdoc <- Co*rdot
	tdoc <- Co*tdot
	tssc <- 1-Cs+Cs*tsst
	tooc <- 1-Co+Co*toot


        ## New weight function Fcdc for crown contribution
        ## (W. Verhoef, 22-05-08)
	rsoc <- Fcdc*rsot
	tssooc <- Fcd*tsstoo+Fcs*toot+Fod*tsst+Fos
        
        ##	Add the soil background 
	dn <- 1-rddcb*bgr
	tup <- (tssc*bgr+tsdc*bgr)/dn
	tdn <- (tsdc+tssc*bgr*rddcb)/dn
        
	rddt <- rddct+tddc*bgr*tddc/dn
	rsdt <- rsdc+tup*tddc
	rdot <- rdoc+tddc*(bgr*tdoc+bgr*tooc)/dn
	rsot <- rsoc+tssooc*bgr+tdn*bgr*tooc+tup*tdoc
    }

    taus <- (tsdc+tssc) ## direct and diffuse transmission of light

    ## rddt    Bi-hemispherical reflectance
    ## rsdt    Directional-hemispherical reflectance for solar incident flux
    ## rdot    Hemispherical-directional reflectance in viewing direction
    ## rsot    Bi-directional reflectance factor
    finalOut  <- list(rddt,rsdt,rdot,rsot,taus)
    reflMat <- as.matrix(do.call(cbind,finalOut))
    colnames(reflMat) <- c('rddt','rsdt','rdot','rsot','tau')
    
    return(reflMat)
  
}



## 2 layer Hotspot function
HotSpot2 <- function(lai,lai1,fbu,q,tss,ks,ko,dso){
    
    ## Treatment of the hotspot-effect for 2 layers
    tss <- exp(-ks*lai) ## full canopy
    ck <- exp(-ks*lai1) ## only top
    
    alf <- 1e6
    
    ## Apply correction 2/(K+k) suggested by F.M. Breon
    
    if(q>0) alf <- (dso/q)*2/(ks+ko)

    if(alf>200) alf <- 200 ## inserted H. Bach 1/3/04 

    if(alf==0){
        ## The pure hotspot - no shadow
        tsstoo <- tss
        sumint1 <- (1-ck)/(ks*lai)
        sumint2 <- (ck-tss)/(ks*lai)

    } else {
        ## Outside the hotspot
        fhot <- lai*sqrt(ko*ks)

        ## Integrate by exponential Simpson method in 20 steps
        ## the steps are arranged according to equal partitioning
        ## of the slope of the joint probability function

        x1     <- 0
        y1     <- 0
        f1     <- 1
        ca <- exp(alf*(fbu-1))
        fint   <- (1-ca)*.05
        sumint1 <- 0

        for(i in 1:20){
            
        if(i<20) {
            x2 <- -log(1-i*fint)/alf
        } else  {
            x2 <- 1-fbu
        }
        
        y2     <- -(ko+ks)*lai*x2+fhot*(1-exp(-alf*x2))/alf
        f2     <- exp(y2)
        sumint1 <- sumint1+(f2-f1)*(x2-x1)/(y2-y1)
        x1     <- x2
        y1     <- y2
        f1     <- f2
        }

	fint <- (ca-exp(-alf))*.05
	sumint2 <- 0
        for(i in 1:20){

            if(i<20) {
                x2 <- -log(ca-i*fint)/alf
            } else  {
                x2 <- 1
            }

            y2     <- -(ko+ks)*lai*x2+fhot*(1-exp(-alf*x2))/alf
            f2     <- exp(y2)

            sumint2 <- sumint2+(f2-f1)*(x2-x1)/(y2-y1)
            x1     <- x2
            y1     <- y2
            f1     <- f2
           
        }

        
        tsstoo <- f1
    }

    return(list(tsstoo,sumint1,sumint2))
} # HotSpot2
     
    
ReflTrans2 <- function(rho1,rho2,tau1,tau2,ks,ko,
                       lai1,lai2,sdf,sdb,dof,dob,
                       sob,sof,ddb,ddf,
                       sumint1,sumint2){
    
    ##	Bottom layer (rho2,tau2,lai2)
    bottom <- cReflTransSingleLayer(rho2,tau2,lai2,
                                    ks,ko,
                                    sdf,sdb,dof,dob,
                                    sob,sof,ddb,ddf)
    
    rddb <- bottom[[1]]
    rsdb <- bottom[[2]]
    rdob <- bottom[[3]]
    rsodb <- bottom[[4]]
    tddb <- bottom[[5]]
    tsdb <- bottom[[6]]
    tdob <- bottom[[7]]
    toob <- bottom[[8]]
    tssb <- bottom[[9]]
    w2 <- bottom[[10]]

    
    ##	Top layer (rho1,tau1,lai1)
    
    top <- cReflTransSingleLayer(rho1,tau1,lai1,
                                ks,ko,
                                sdf,sdb,dof,dob,
                                sob,sof,ddb,ddf)
    
    rdd <- top[[1]]
    rsd <- top[[2]]
    rdo <- top[[3]]
    rsod <- top[[4]]
    tdd <- top[[5]]
    tsd <- top[[6]]
    tdo <- top[[7]]
    too <- top[[8]]
    tss <- top[[9]]
    w1 <-  top[[10]]

    ##  Combine with bottom layer reflectances and
    ##  transmittances (adding method)
    lai <- lai1+lai2 ## total LAI
    
    dn <- 1-rdd*rddb
    tup<-(tss*rsdb+tsd*rddb)/dn
    tdn<-(tsd+tss*rsdb*rdd)/dn
    rsdt<-rsd+tup*tdd
    rdot<-rdo+tdd*(rddb*tdo+rdob*too)/dn
    rsodt<-rsod+(tss*rsodb+tdn*rdob)*too+tup*tdo

    rsost<-(w1*sumint1+w2*sumint2)*lai
    rsot<-rsost+rsodt

    ## Diffuse reflectances at the top and the bottom are now different 
    rddt_t<-rdd+tdd*rddb*tdd/dn
    rddt_b<-rddb+tddb*rdd*tddb/dn

    ## Transmittances of the combined canopy layers
    tsst<-tss*tssb
    toot<-too*toob
    tsdt<-tss*tsdb+tdn*tddb
    tdot<-tdob*too+tddb*(tdo+rdd*rdob*too)/dn
    tddt<-tdd*tddb/dn

    return(list(dn,tup,tdn,rsdt,rdot,rsodt,rsost,rsot,rddt_t,
                rddt_b,tsst,toot,tsdt,tdot,tddt))
}







## calculate transmittance and reflectance of a single layer
ReflTransSingleLayer <- function(rho,tau,lai,
                                 ks,ko,
                                 sdf,sdb,dof,dob,
                                 sob,sof,ddb,ddf){


    ## if 3 or leaf refl-- Geometric leaf reflectance
    ## Input rho and tau from PROSPECT/ other particle model
    ## Correct using SUITS factors
    ##	Here the reflectance and transmission come in
    ## with the suits coefficients
    RTgeomRes <- RTgeom(rho,tau,ddb,ddf,sdb,sdf,dob,dof,sob,sof)
    
    sigb <- RTgeomRes[[1]]
    att  <- RTgeomRes[[2]]
    m    <- RTgeomRes[[3]]
    sb   <- RTgeomRes[[4]]
    sf   <- RTgeomRes[[5]]
    vb   <- RTgeomRes[[6]]
    vf   <- RTgeomRes[[7]]
    w    <- RTgeomRes[[8]]

    tss <- exp(-ks*lai);
    too <- exp(-ko*lai);

    J1ks <- Jfunc1(ks,m,lai)
    J2ks <- Jfunc2(ks,m,lai)
    J1ko <- Jfunc1(ko,m,lai)
    J2ko <- Jfunc2(ko,m,lai)
    z  <- Jfunc3(ks,ko,lai)
    J4 <- Jfunc4(m,lai)
   
    m_val1 <- which(m>0.01)
    m_val2 <- which(m<0.01)

    ## Normal case
    rinf<-rinf2<-e1<-e2<-denom<-re<-rdd<-tdd <- 0*m
    e1 <- exp(-m[m_val1]*lai)
    e2[m_val1] <- e1[m_val1]*e1[m_val1]
    rinf[m_val1] <- ((att[m_val1]-m[m_val1])/sigb[m_val1])
    rinf2[m_val1] <- rinf[m_val1]*rinf[m_val1]
    re[m_val1] <- rinf[m_val1]*e1[m_val1]
    denom[m_val1] <- 1-rinf2[m_val1]*e2[m_val1]


    Ps<-Qs<-Pv<-Qv<-rdo<-tdo<-rds<-tds<-tsd<-rsd<- 0*m
    Ps[m_val1]  <-  (sf[m_val1]+sb[m_val1]*rinf[m_val1])*J1ks[m_val1]
    Qs[m_val1] <- (sf[m_val1]*rinf[m_val1]+sb[m_val1])*J2ks[m_val1]
    Pv[m_val1] <- (vf[m_val1]+vb[m_val1]*rinf[m_val1])*J1ko[m_val1]
    Qv[m_val1] <- (vf[m_val1]*rinf[m_val1]+vb[m_val1])*J2ko[m_val1]
    
    rdd[m_val1]	<- rinf[m_val1]*(1-e2[m_val1])/denom[m_val1]
    tdd[m_val1]	<- (1-rinf2[m_val1])*e1[m_val1]/denom[m_val1]
    
    tsd[m_val1]	<- (Ps[m_val1]-re[m_val1]*Qs[m_val1])/denom[m_val1]
    rsd[m_val1]	<- (Qs[m_val1]-re[m_val1]*Ps[m_val1])/denom[m_val1]
    tdo[m_val1]	<- (Pv[m_val1]-re[m_val1]*Qv[m_val1])/denom[m_val1]
    rdo[m_val1]	<- (Qv[m_val1]-re[m_val1]*Pv[m_val1])/denom[m_val1]
    
    g1 <- g2 <- 0*m
    g1[m_val1] <- (z-J1ks[m_val1]*too)/(ko+m[m_val1])
    g2[m_val1] <- (z-J1ko[m_val1]*tss)/(ks+m[m_val1])

    Tv1 <- Tv2 <- T1 <- T2 <- T3 <- 0*m
    Tv1[m_val1] <- (vf[m_val1]*rinf[m_val1]+vb[m_val1])*g1[m_val1]
    Tv2[m_val1] <- (vf[m_val1]+vb[m_val1]*rinf[m_val1])*g2[m_val1]
    T1[m_val1] <- Tv1[m_val1]*(sf[m_val1]+sb[m_val1]*rinf[m_val1])
    T2[m_val1] <- Tv2[m_val1]*(sf[m_val1]*rinf[m_val1]+sb[m_val1])
    T3[m_val1] <- (rdo[m_val1]*Qs[m_val1]+
                   tdo[m_val1]*Ps[m_val1])*rinf[m_val1]
    ## Multiple scattering contribution to
    ## bidirectional canopy reflectance

    rsod <- 0*m
    rsod[m_val1] <- (T1[m_val1]+T2[m_val1]-T3[m_val1])/(1-rinf2[m_val1])

    ## Near or complete conservative scattering
    amsig <- apsig <- rtp <- rtm <- 0*m
    amsig[m_val2] <- att[m_val2]-sigb[m_val2]
    apsig[m_val2] <- att[m_val2]+sigb[m_val2]
    
    rtp[m_val2] <- (1-amsig[m_val2]*J4[m_val2])/
        (1+amsig[m_val2]*J4[m_val2])
    rtm[m_val2]<-(-1+apsig[m_val2]*J4[m_val2])/
        (1+apsig[m_val2]*J4[m_val2])
    rdd[m_val2]<-5*(rtp[m_val2]+rtm[m_val2])
    tdd[m_val2]<-5*(rtp[m_val2]-rtm[m_val2])
    
    dns<-dno<-cks<-cko<-dks<-dko<-ho<-0*m
    dns[m_val2]<-ks*ks-m[m_val2]*m[m_val2]
    dno[m_val2]<-ko*ko-m[m_val2]*m[m_val2]
    cks[m_val2]<-(sb[m_val2]*(ks-att[m_val2])-
                  sf[m_val2]*sigb[m_val2])/dns[m_val2]
    cko[m_val2]<-(vb[m_val2]*(ko-att[m_val2])-
                  vf[m_val2]*sigb[m_val2])/dno[m_val2]
    dks[m_val2]<-(-sf[m_val2]*(ks+att[m_val2])-
                  sb[m_val2]*sigb[m_val2])/dns[m_val2]
    dko[m_val2]<-(-vf[m_val2]*(ko+att[m_val2])-
                  vb[m_val2]*sigb[m_val2])/dno[m_val2]
    ho[m_val2]<-(sf[m_val2]*cko[m_val2]+sb[m_val2]*dko[m_val2])/(ko+ks)
    
    rsd[m_val2]<-cks[m_val2]*(1-tss*tdd[m_val2])-dks[m_val2]*rdd[m_val2]
    rdo[m_val2]<-cko[m_val2]*(1-too*tdd[m_val2])-dko[m_val2]*rdd[m_val2]
    tsd[m_val2]<-dks[m_val2]*(tss-tdd[m_val2])-
        cks[m_val2]*tss *rdd[m_val2]
    tdo[m_val2]<-dko[m_val2]*(too -tdd[m_val2])-
        cko[m_val2]*too*rdd[m_val2]
    
    ##  Multiple scattering contribution to
    ## bidirectional canopy reflectance
    rsod[m_val2]<-ho[m_val2]*(1-tss*too)-
        cko[m_val2]*tsd[m_val2]*too-dko[m_val2]*rsd[m_val2]
    
    return(list(rdd,rsd,rdo,rsod,tdd,tsd,tdo,too,tss,w))
} #  ReflTransSingleLayer


