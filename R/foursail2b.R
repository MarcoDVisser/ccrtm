#' R implementation of the foursail2 model with 2 canopy layers.
#'
#' The foursail2b model is a two layer implementation of the 
#' foursail model described in Zhang et al (2005).
#' Layers are assumed identical in hotspot,
#' but may differ in the relative density, inclination and types of
#' particles. In comparison to foursail, the background (soil),
#' can now be non-Lambertain, having it own 4-stream
#' BDRF. There are two types of particles, generalized
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
#' fRTM(rho~prospectd+foursail2b)
#'
#' @references Zhang, Q., Xiao, X., Braswell, B., Linder, E., Baret, F., Moore, B. (2005).
#'   Estimating light absorption by chlorophyll, leaf and canopy in a deciduous broadleaf forest
#'   using MODIS data and a radiative transfer model. Remote Sens. Environ. 99, 357-371. 
#' @export
foursail2b <- function(rhoA,tauA,rhoB=NULL,tauB=NULL,bgr,
                      rsobgr=NULL,rdobgr=NULL,rsdbgr=NULL,rddbgr=NULL,
                      param){

    if(!is.null(rsobgr)) stop("SOIL BDRF not implemented currently")

    if(is.null(rhoB)) { ## standard if B is empty

        rhoB <- rhoA
        tauB <- tauA
   
    }
    
    names(param) <- tolower(names(param))
    
    
    ## input paramters
    LIDFat <-    as.numeric(param["lidfa"])
    LIDFbt <-    0 #as.numeric(param["lidfa"])
    LIDFab <-    as.numeric(param["lidfb"])
    LIDFbb <-    0 #as.numeric(param["lidfb"])

    ## TypeLIDF <- as.numeric(param["typelidf"])
    TypeLIDF <- 2  ## hard coded Cambell!
    
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


    ## top leaf angles
    lidfpart <- c(na,LIDFat, LIDFbt)
    names(lidfpart) <- c("na","a","b")
    class(lidfpart) <- paste0("lidf.",TypeLIDF)
    lidFunt <- lidf(lidfpart)

    ## bottom leaf angles
    lidfparb <- c(na,LIDFab, LIDFbb)
    names(lidfparb) <- c("na","a","b")
    class(lidfparb) <- paste0("lidf.",TypeLIDF)
    lidFunb <- lidf(lidfparb)
    
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


        ## top canopy 
        ## Angular distance - shadow length compensation
        ## output = list(ks,ko,sob,sof,sdb,sdf,dob,dof,ddb,ddf)
        suitRes <- SUITS(na,lna,lidFunt,tts,tto,cts,cto,psi,ctscto)
        
        kst  <- suitRes[[1]]
        kot  <- suitRes[[2]] 
        sobt <- suitRes[[3]] 
        soft <- suitRes[[4]] 
        sdbt <- suitRes[[5]]  
        sdft <- suitRes[[6]]  
        dobt <- suitRes[[7]]  
        doft <- suitRes[[8]]  
        ddbt <- suitRes[[9]]  
        ddft <- suitRes[[10]]

        ## bottom canopy
        ## Angular distance - shadow length compensation
        ## output = list(ks,ko,sob,sof,sdb,sdf,dob,dof,ddb,ddf)
        suitRes <- SUITS(na,lna,lidFunb,tts,tto,cts,cto,psi,ctscto)
        
        ksb  <- suitRes[[1]]
        kob  <- suitRes[[2]] 
        sobb <- suitRes[[3]] 
        sofb <- suitRes[[4]] 
        sdbb <- suitRes[[5]]  
        sdfb <- suitRes[[6]]  
        dobb <- suitRes[[7]]  
        dofb <- suitRes[[8]]  
        ddbb <- suitRes[[9]]  
        ddfb <- suitRes[[10]]

        ## devide the leaf area index in two layers
	lai1 <- (1-fbu)*lai
	lai2 <- fbu*lai
        
        ##  2 layer hot spot effect
        hspotRes <- HotSpot2b(lai,fbu,q,tss,kst,kot,ksb,kob,dso)
        tsstoo <- hspotRes[[1]]
        sumint1 <- hspotRes[[2]]
        sumint2 <- hspotRes[[3]]
        
        
        ## Calculate reflectances and transmittances for 2 layers
        ## Input LAI and leaf rho and tau
        ## and geometric factors
        
        refTransRes <- ReflTrans2b(rho1,rho2,tau1,tau2,
                                  lai1,lai2,
                                  kst,kot,
                                  sdft,sdbt,doft,dobt,
                                  sobt,soft,ddbt,ddft,
                                  ksb,kob,
                                  sdfb,sdbb,dofb,dobb,
                                  sobb,sofb,ddbb,ddfb,
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
HotSpot2b <- function(lai,fbu,q,tss,kst,kot,ksb,kob,dso){
    
    ## Treatment of the hotspot-effect for 2 layers

    ks <- (ksb*fbu + (1-fbu)*kst) ## full canopy ks
    ko <- (kob*fbu + (1-fbu)*kot) ## full canopy ko
    
    tss <- exp(-ks*lai) ## full canopy
    ck <- exp(-kst*(1-fbu)*lai) ## only top!!!!!!
    
    alf <- 1e6
    
    ## Apply correction 2/(K+k) suggested by F.M. Breon
    ## discussed in Verhoef and Bach 2007

    ## alf = q in Verhoef and Bach 2007 (there sl = q)
    if(q>0) alf <- (dso/q)*2/(ks+ko)

    if(alf>200) alf <- 200 ## inserted H. Bach 1/3/04 

    if(alf==0){
        ## The pure hotspot - no shadow
        tsstoo <- tss
        sumint1 <- (1-ck)/(ks*lai) ## bottom
        sumint2 <- (ck-tss)/(ks*lai) ## top

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
     
    
ReflTrans2b <- function(rho1,rho2,tau1,tau2,
                        lai1,lai2,
                        kst,kot,
                        sdft,sdbt,doft,dobt,
                        sobt,soft,ddbt,ddft,
                        ksb,kob,
                        sdfb,sdbb,dofb,dobb,
                        sobb,sofb,ddbb,ddfb,
                        sumint1,sumint2){
    
    ##	Bottom layer (rho2,tau2,lai2)
    bottom <- cReflTransSingleLayer(rho2,tau2,lai2,
                                    ksb,kob,
                                    sdfb,sdbb,dofb,dobb,
                                    sobb,sofb,ddbb,ddfb)
    
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
                                kst,kot,
                                sdft,sdbt,doft,dobt,
                                sobt,soft,ddbt,ddft)
    
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



