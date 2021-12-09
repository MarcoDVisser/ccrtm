#' run a requested RTM (internal function)
#'
#' @param modReq  model request object built in fRTM
#' @param \dots additional plot arguments
#' @return prediction from the requested model
runRTM <- function(modReq){

    result <- UseMethod("rtm",pars)
    return(result)

}

rtm.prospect5 <- function(pars){
    
    prospect5(pars)

}

rtm.prospectd <- function(pars){
    
    prospectd(pars)
    
}




rtm.prosail5 <- function(pars){
    
    LS <- prospect5(pars[["prosail5"]])
    foursail(LS[,"rho"],LS[,"tau"],pars[["foursail"]])

}


rtm.prosaild <- function(pars){
    
    LS <- prospectd(pars[["prosaild"]])
    foursail(LS[,"rho"],LS[,"tau"],pars[["foursail"]])
        

}



rtm.prosail2_55 <- function(pars){
    
    LS <- prospect5(pars[["prosail5"]])
    LS <- prospect5(pars[["prosail5"]])

    foursail(LS[,"rho"],LS[,"tau"],pars[["foursail"]])

}


rtm.prosail2_dd <- function(pars){
    
    LS <- prospectd(pars[["prosaild"]])
    foursail(LS[,"rho"],LS[,"tau"],pars[["foursail"]])
        

}


rtm.prosail2_d5 <- function(pars){
    
    LS <- prospect5(pars[["prosail5"]])
    foursail(LS[,"rho"],LS[,"tau"],pars[["foursail"]])

}


rtm.prosail2_5d <- function(pars){
    
    LS <- prospectd(pars[["prosaild"]])
    foursail(LS[,"rho"],LS[,"tau"],pars[["foursail"]])
        

}


## example
## require(ccrtm)
## informpars <-  ccrtm:::defaults.inform5()
## R <- ccrtm:::rtm.inform5(informpars)
rtm.inform5 <- function(pars){

    ## setup inform parameters
    parsc <- pars[["foursail"]][["canopy"]]
    parsu <- pars[["foursail"]][["understorey"]]

    skyl <- pars[["skyl"]]

    leafpar <- pars[["prospect5"]]
    lRTc <- prospect5(leafpar[["canopy"]]) ## canopy particles 
    lRTu <- prospect5(leafpar[["understorey"]]) ## unstorey particles

    ## calculate background (soil) reflectance
    bgRef<- parsc["psoil"]*soil[,"drySoil"] +
        (1-parsc["psoil"])*soil[,"wetSoil"]

    tmpar <- parsc
    tmpar["LAI"] <- 15
    RTinf <- foursail(lRTc[,"rho"], lRTc[,"tau"], bgRef,tmpar) ## canopy reflectance infinate depth 
    RTu <- foursail(lRTu[,"rho"], lRTu[,"tau"], bgRef,parsu) ## understorey reflectance depth

    ## !expected implementation following available code!
    rhoinf <- skyl(RTinf[,"rddt"], RTinf[,"rsdt"], RTinf[,"rdot"], RTinf[,"rsot"],
                   skyl,Es=1,Ed=1)$directional
    rhou <- skyl(RTu[,"rddt"], RTu[,"rsdt"], RTu[,"rdot"], RTu[,"rsot"],
                 skyl,Es=1,Ed=1)$directional
    
    ## Crown transmittance in sun and observer direction (taus and tauo)
    RTs <- foursail(lRTc[,"rho"], lRTc[,"tau"],rhou,parsc)

    tss <- RTs[,"tss"]
    tsd <- RTs[,"tsd"]
    rdd <- RTs[,"rdd"]
    tdd <- RTs[,"tdd"]
   
    dn    <- 1-bgRef*rdd         # Interaction with the soil/background
    tsd0 <- tss+(tsd+tss*bgRef*rdd)/dn
    tdd0 <- tdd/dn
 
    ## Transmission in observed direction
    ## inform implementation of tauo - not robust in ALL CASES!
    ## when tts=\= tto &  psi --> 0 bias is maximized
    ## only use when you are 100% certain you want inform
    ## otherwise use SAIL2 as this is done correctly/unbiased here
    tmpar <- parsc
    tmpar["tts"] <- tmpar["tto"] 
    RTo <- foursail(lRTc[,"rho"], lRTc[,"tau"],rhou,parsc)
    
    tss <- RTo[,"tss"]
    tsd <- RTo[,"tsd"]
    rdd <- RTo[,"rdd"]
    tdd <- RTo[,"rdd"]

    dn    <- 1-bgRef*rdd         # Interaction with the soil/background
    tsd1 <- tss+(tsd+tss*bgRef*rdd)/dn
    tdd1 <- tdd/dn

    ## diffuse light part with sky light model
    ## in original code: trans_hemi=(TSD*ES+TDD*ED)./(ES+ED);
    trans  <- skyl(tsd0,tdd0,tsd1,tdd1,skyl,Es=1,Ed=1) ## Es,Ed =1 should be replaced with solar rediance!
    taus <- trans[[1]]
    tauo <- trans[[2]]
   
    ## apply flim
    R <- flim(rhoinf,rhou,tauo,taus,params=pars[["flim"]])$rho
    return(R)
    
}



rtm.informd <- function(pars){

## setup inform parameters
    parsc <- pars[["foursail"]][["canopy"]]
    parsu <- pars[["foursail"]][["understorey"]]

    skyl <- pars[["skyl"]]

    leafpar <- pars[["prospectd"]]
    lRTc <- prospectd(leafpar[["canopy"]]) ## canopy particles 
    lRTu <- prospectd(leafpar[["understorey"]]) ## unstorey particles

    ## calculate background (soil) reflectance
    bgRef<- parsc["psoil"]*soil[,"drySoil"] +
        (1-parsc["psoil"])*soil[,"wetSoil"]

    tmpar <- parsc
    tmpar["LAI"] <- 15
    RTinf <- foursail(lRTc[,"rho"], lRTc[,"tau"], bgRef,tmpar) ## canopy reflectance infinate depth 
    RTu <- foursail(lRTu[,"rho"], lRTu[,"tau"], bgRef,parsu) ## understorey reflectance depth

    ## !expected implementation following available code!
    rhoinf <- skyl(RTinf[,"rddt"], RTinf[,"rsdt"], RTinf[,"rdot"], RTinf[,"rsot"],
                   skyl,Es=1,Ed=1)$directional
    rhou <- skyl(RTu[,"rddt"], RTu[,"rsdt"], RTu[,"rdot"], RTu[,"rsot"],
                 skyl,Es=1,Ed=1)$directional
    
    ## Crown transmittance in sun and observer direction (taus and tauo)
    RTs <- foursail(lRTc[,"rho"], lRTc[,"tau"],rhou,parsc)

    tss <- RTs[,"tss"]
    tsd <- RTs[,"tsd"]
    rdd <- RTs[,"rdd"]
    tdd <- RTs[,"tdd"]
   
    dn    <- 1-bgRef*rdd         # Interaction with the soil/background
    tsd0 <- tss+(tsd+tss*bgRef*rdd)/dn
    tdd0 <- tdd/dn
 
    ## Transmission in observed direction
    ## inform implementation of tauo - not robust in ALL CASES!
    ## when tts=\= tto &  psi --> 0 bias is maximized
    ## only use when you are 100% certain you want inform
    ## otherwise use SAIL2 as this is done correctly/unbiased here
    tmpar <- parsc
    tmpar["tts"] <- tmpar["tto"] 
    RTo <- foursail(lRTc[,"rho"], lRTc[,"tau"],rhou,parsc)
    
    tss <- RTo[,"tss"]
    tsd <- RTo[,"tsd"]
    rdd <- RTo[,"rdd"]
    tdd <- RTo[,"rdd"]

    dn    <- 1-bgRef*rdd         # Interaction with the soil/background
    tsd1 <- tss+(tsd+tss*bgRef*rdd)/dn
    tdd1 <- tdd/dn

    ## diffuse light part with sky light model
    ## in original code: trans_hemi=(TSD*ES+TDD*ED)./(ES+ED);
    trans  <- skyl(tsd0,tdd0,tsd1,tdd1,skyl,Es=1,Ed=1) ## Es,Ed =1 should be replaced with solar rediance!
    taus <- trans[[1]]
    tauo <- trans[[2]]
   
    ## apply flim
    R <- flim(rhoinf,rhou,tauo,taus,params=pars[["flim"]])$rho
    return(R)
}





