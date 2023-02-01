#' run a requested RTM (internal function)
#' 
#' List of aliases: prospect5, prospectd, prosail5,
#' prosaild, prosail2_55,prosail2_dd, prosail2_5d,
#' prosail2_d5, rtm.inform5, rtm.informd
#' @param modReq  model request object built in fRTM
#' @param pars the required parameters (vector or list)
#' @return prediction from the requested model
runRTM <- function(pars){
    result <- UseMethod("rtm",pars)
    return(result)
}

## leaf models

rtm.prospect5 <- function(pars){
  prospect5(pars[["prospect5"]])

}

rtm.prospectd <- function(pars){
  prospectd(pars[["prospectd"]])
}

## 4SAIL2

rtm.prosail5 <- function(pars){

  SS <- pars[["foursail"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LS <- prospect5(pars[["prospect5"]]) ## leaf spectra

  foursail(LS[,"rho"],LS[,"tau"],SS,pars[["foursail"]])
}

rtm.prosaild <- function(pars){

  SS <- pars[["foursail"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LS <- prospectd(pars[["prosaild"]])

  foursail(LS[,"rho"],LS[,"tau"],SS,pars[["foursail"]])
}

## 4SAIL2

rtm.prosail2_55 <- function(pars){

  SS <- pars[["foursail2"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail2"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LSa <- prospect5(pars[["prospect5.a"]]) ## leaf spectra
  LSb <- prospect5(pars[["prospect5.b"]]) ## leaf spectra

  foursail2(LSa[,"rho"],LSa[,"tau"],LSb[,"rho"],LSb[,"tau"],SS,param=pars[["foursail2"]])

}

rtm.prosail2_dd <- function(pars){
  SS <- pars[["foursail2"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail2"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LSa <- prospect5(pars[["prospectd.a"]]) ## leaf spectra
  LSb <- prospect5(pars[["prospectd.b"]]) ## leaf spectra

  foursail2(LSa[,"rho"],LSa[,"tau"],LSb[,"rho"],LSb[,"tau"],SS,param=pars[["foursail2"]])

}

rtm.prosail2_d5 <- function(pars){

  SS <- pars[["foursail2"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail2"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LSa <- prospect5(pars[["prospectd"]]) ## leaf spectra
  LSb <- prospect5(pars[["prospect5"]]) ## leaf spectra

  foursail2(LSa[,"rho"],LSa[,"tau"],LSb[,"rho"],LSb[,"tau"],SS,param=pars[["foursail2"]])

}

rtm.prosail2_5d <- function(pars){
  SS <- pars[["foursail2"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail2"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LSa <- prospect5(pars[["prospect5"]]) ## leaf spectra
  LSb <- prospect5(pars[["prospectd"]]) ## leaf spectra

  foursail2(LSa[,"rho"],LSa[,"tau"],LSb[,"rho"],LSb[,"tau"],SS,param=pars[["foursail2"]])
}

## 4SAIL2b

rtm.prosail2b_55 <- function(pars){

  SS <- pars[["foursail2b"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail2b"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LSa <- prospect5(pars[["prospect5.a"]]) ## leaf spectra
  LSb <- prospect5(pars[["prospect5.b"]]) ## leaf spectra

  foursail2b(LSa[,"rho"],LSa[,"tau"],LSb[,"rho"],LSb[,"tau"],SS,param=pars[["foursail2b"]])

}

rtm.prosail2b_dd <- function(pars){
  SS <- pars[["foursail2b"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail2b"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LSa <- prospect5(pars[["prospectd.a"]]) ## leaf spectra
  LSb <- prospect5(pars[["prospectd.b"]]) ## leaf spectra

  foursail2b(LSa[,"rho"],LSa[,"tau"],LSb[,"rho"],LSb[,"tau"],SS,param=pars[["foursail2b"]])

}

rtm.prosail2b_d5 <- function(pars){

  SS <- pars[["foursail2b"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail2b"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LSa <- prospect5(pars[["prospectd"]]) ## leaf spectra
  LSb <- prospect5(pars[["prospect5"]]) ## leaf spectra

  foursail2b(LSa[,"rho"],LSa[,"tau"],LSb[,"rho"],LSb[,"tau"],SS,param=pars[["foursail2b"]])

}

rtm.prosail2b_5d <- function(pars){
  SS <- pars[["foursail2b"]]["psoil"]*soil[,"drySoil"] +
    (1-pars[["foursail2b"]]["psoil"])*soil[,"wetSoil"] ## soil spectra

  LSa <- prospect5(pars[["prospect5"]]) ## leaf spectra
  LSb <- prospect5(pars[["prospectd"]]) ## leaf spectra

  foursail2b(LSa[,"rho"],LSa[,"tau"],LSb[,"rho"],LSb[,"tau"],SS,param=pars[["foursail2b"]])
}


## INFORM

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
