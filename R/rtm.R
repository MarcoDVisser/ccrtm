### TO DO: 
## 1) internal parameter and model order and existence check function
## 2) internal prediction function like fit.lm so bRTM can call fRTM 
## without the need of double parameter and model checks - this  
## could be done via and orderedRTM prediction function.

## Species actions
## rho + tau ~ sp(foursail, rho=ref,tau=trans,bgr=rsoil) + geo

## All models should have a named/ numbered list as input
## so all models can just accept a parameter vector
## and the output of the previous model directly
## RTM2(..., param=pars[[2]],chain=predictions1)

#' Forward implementation of coupled Radiative Transfer Models.
#'
#' @param fm A formula specifying which rtm to run
#' @param pars a list of _named_ parameter vectors for all models.
#' The parameter list for a model call as rho ~ prospect + foursail
#' therefore contains two named vectors the first with parameters for
#' prospect and the second with parameters for foursail 
#' if left empty default parameters are generated
#' @param wl wavelengths (in nm) add only if
#' certain wavelengths are required as output.
#' Input is expected to integers between 400 and 2500, or will 
#' be forced to be an integer. Integers outside the 400:2500 range
#' will not be returned.
#'
#' @return spectra matrix with reflectance (and transmission,
#' depending on the formula inputs). 
#' See seperate model helpfiles for details. 
#'
#' @examples
#' ## setup graphics for plots 
#' oldpar<-par()
#' par(mfrow=c(3,2))
#' 
#' ## get reflectance for a leaf 
#' ref <- fRTM(rho~prospect5)
#' plot(ref,main="Prospect 5")
#'      
#' ## get reflectance and transmission for a leaf 
#' reftrans <- fRTM(rho+tau~prospect5)
#' plot(reftrans,main="Prospect 5")
#'      
#' ## get reflectance for a single layered canopy 
#' ref <- fRTM(rho~prospect5+foursail)
#' plot(ref,main="Prospect 5 + 4SAIL")
#' 
#' ## get reflectance for a 2 layered canopy with two leaf types 
#' ref <- fRTM(rho~prospectd+prospect5+foursail2)
#' plot(ref,main="Prospect D + Prospect 5  + 4SAIL2")
#' 
#' ## edit the parameters: sparse vegatation LAI 
#' parlist<- list(prospect5=NULL,prospectd=NULL,foursail2=c(LAI=0.05))
#' 
#' ## update reflectance
#' ref <- fRTM(rho~prospect5+prospectd+foursail2,parlist)
#' plot(ref,main="LAI=0.05")
#' 
#' ## change leaf area index to dense vegetation
#' parlist$foursail2["LAI"]<-8.5
#' 
#' ## update reflectance
#' ref <- fRTM(rho~prospect5+prospectd+foursail2,parlist)
#' plot(ref,main="LAI=8.5")
#'
#' par(oldpar)   
#' @export
fRTM <- function(fm= rho + tau ~ prospect5 + foursail , pars=NULL,
                 wl=400:2500){

    rtmModels <- getModels()

    wl <- as.integer(wl) # force integer
    
    rowInc <- 400:2500%in%wl

    if(!any(rowInc)) stop("check wavelengths: wl")
    
    tmp <- checkForm(fm)
    reqMods <- tmp[[1]]
    ordN <- tmp[[2]]

    ## get parameters or check parameters
    pars <- checkPars(pars,reqMods,ordN)

    ## get desired predictions 
    tmp<-getPredictions(fm, ordN)
    
    returnCol<-tmp[[1]]
    pred<-tmp[[2]]

    model<-orderedRTM(reqMods,pars)
    
    final<-model(pars)[rowInc,returnCol]

    ## finalize output
    attr(final,"main")<-paste0(pred, " ")
    attr(final,"wavelength")<-wl
    attr(final,"models") <- reqMods
    class(final) <- c(class(final),"rtm.spectra")
    return(final)
  
}

 


#' Generates an invertable model for backward implementation
#' of Radiative Transfer Models
#' @inheritParams fRTM
#' @param fixed a list of parameters to fix
#' @param data ignored as of yet
#' 
# @export
# plot(bRTM(rho ~ prospect,parvec=c(Car=20,Cw=.002,Cm=0.001,Cbrown=0.001),fixed=c(N=2,Cab=9)))
bRTM <- function(fm = rho ~ prospect5, data=NULL, pars=NULL, fixed=NULL,
                 wl=400:2500){

    
    wl <- as.integer(wl) # force integer
    
    rowInc <- 400:2500%in%wl

    if(!any(rowInc)) stop("check wavelengths: wl")
    
    
    tmp <- checkForm(fm)
    reqMods <- tmp[[1]]
    ordN <- tmp[[2]]
    
    ## get parameters or check parameters
    pars <- checkPars(pars,reqMods,ordN)


    ## remove duplication id
    dub <- duplicated(reqMods)
    if(any(dub)){
        for(i in 1:length(fixed)){
            tmp <- names(fixed[[i]])
            if(!is.null(tmp)) names(fixed[[i]]) <- gsub("[0-9._//-]","",tmp)
        }
    }
    

    ## replace pars with fixed pars
 
    if(!is.null(fixed)) {

        if(!class(fixed)=="list") fixed <- list(fixed)

        if(length(fixed)!=length(pars)){
            
            stop("fixed parameters is not a list of expected length")
            
        }
        
        for(i in 1:length(pars)){

            if(!is.null(fixed[[i]])){
            names(fixed[[i]]) <- tolower(names(fixed[[i]])) 	
            stopifnot(all(names(fixed[[i]])%in%names(pars[[i]])))
            pars[[i]][names(fixed[[i]])] <- fixed[[i]]
            }
        }
        

        
    }
    
    
    ## get desired predictions 
    tmp<-getPredictions(fm, ordN) 	
    returnCol<-tmp[[1]]
    pred<-tmp[[2]]
    
    
    model<-orderedRTM(reqMods,pars)
    final<-model(pars)[rowInc,returnCol]
    
    ## finalize output
    attr(final,"main")<-paste0(pred, " ")
    attr(final,"wavelength")<-wl
    attr(final,"models") <- reqMods
    class(final) <- c(class(final),"rtm.spectra")
    return(final)
    
#ref ~ prospect + sail + geo | N = 2, LAI=4 .,,
#ref ~ (prospect | N = 2) + (sail | LAI =3) + (geo | X= 2)
#predict(prospect + sail + geo)
}


## internal prediction building function
orderedRTM<-function(reqMods,prms){

ordN<-length(reqMods)

if(ordN==1){
       
       model1<-function(prms,RTM=get(reqMods[1])){ 
                
        ord1pred <- RTM(prms[[1]])

        return(ord1pred)   
       }
    return(model1)

    }
    
    if(ordN==2){

	
       
        bgRef<- prms[[2]]["psoil"]*soil[,"drySoil"] +
            (1-prms[[2]]["psoil"])*soil[,"wetSoil"]

	model2<-function(prms,RTM1=get(reqMods[1]),
			RTM2=get(reqMods[2]),
			Rsoil=bgRef){	
        
        
        ## get first order prediction: particles
        ord1pred <- RTM1(param=prms[[1]])
        
        ## get second order prediction: medium 
        ord2pred <- RTM2(ord1pred[,"rho"],
                         ord1pred[,"tau"],
                         bgr=Rsoil,
                         param=prms[[2]])
            
        return(ord2pred)

        }
        
        return(model2)

    }
    
    
    if(ordN==3){
        
        bgRef<- prms[[3]]["psoil"]*soil[,"drySoil"] +
            (1-prms[[3]]["psoil"])*soil[,"wetSoil"]
        
        model3<-function(prms,
                        RTM1=get(reqMods[1]),
                        RTM2=get(reqMods[2]),
                        RTM3=get(reqMods[3]),
                        Rsoil=bgRef){	
            
            ## get first order prediction
            ord1pred <- RTM1(param=prms[[1]])

            ## get second order prediction
            ord2pred <- RTM2(param=prms[[2]])
            
            ## get third order prediction
            ord3pred <- RTM3(ord1pred[,"rho"],
                             ord1pred[,"tau"],
                             ord2pred[,"rho"],
                             ord2pred[,"tau"],
                             bgr=Rsoil,
                             param=prms[[3]])
            
            return(ord3pred)
            
            
        }
        
        return(model3)
        
    }
    
 
}

## get prediction columns
getPredictions<-function(fm,ordN){
    
    pred <- as.character(fm[[2]])
    pred <- pred[!grepl("\\+",as.character(pred))]
    
    if(!any(c("rho","tau")%in%pred)) stop("Unknown predictions requested")
    
    if(length(pred)>1){

        if(ordN>1){
        returnCol <- 1:5 ## return reflectance and transmittance
        } else{
        returnCol <- 1:2 ## return reflectance and transmittance
        }
        
    }  else if(any("rho"%in%pred)){

        if(ordN>1){
            returnCol <- 1:4 ## return reflectance
        } else {
            returnCol <- 1 ## return reflectance
        }
        
    } else {
        
        if(ordN>1){
            returnCol <- 5 ## return transmittance
        } else {
            returnCol <- 2 ## return transmittance
        }
        
    }

return(list(returnCol, pred))

}


## Get list of implemented models
getModels <- function(){

    mods <- c("prospect5",
              "prospectd",
              "foursail",
              "foursail2",
              "foursail2b")
    
    data.frame(model=mods,
               level=c(1,1,2,3,3),
                row.names = mods)
    
}




## internal function to check rtm formulas 
checkForm <- function(fm=NULL){

    if(!class(fm)=="formula") stop("Invalid formula")
    
    ## get order of RTM
    reqMods <- attr(terms(fm),"term.labels")

    ## assert that repeated terms are not lost!
    test <-unlist(strsplit(as.character(fm)[[3]]," \\+ "))
    if(length(test)>length(reqMods)) reqMods <- test

    
    ordN <- length(reqMods)

    ## get model information
    rtmModels <- getModels()

    ## check models
    if(any(!reqMods%in%rtmModels$model)){
        
        stop(reqMods[!reqMods%in%rtmModels$model], " not implemented")
    }

    suppOrder <- rtmModels[reqMods,]$level
    expOrder <- sort(suppOrder) ## must increase in order
    
    if(!all.equal(suppOrder,expOrder)){
        
        stop("Incorrect ordering of models")

    }

    
    if(sum(suppOrder==1)>2|sum(suppOrder==2)>1|sum(suppOrder==3)>1) {
        
        warning("models appear repeated check formula")

    }

    ## check if models are compatible
    if(max(suppOrder)<ordN) stop("Non compatible models: check formula")

    return(list(reqMods,ordN,suppOrder))
    
}


## function to check and return parameters
checkPars <- function(pars,reqMods,ordN){

    ## generate expected parameters
    parlist <- lapply(reqMods, getDefaults)
    expparnames <- lapply(parlist, rownames)

    ## add duplicated parameters
    dub <- duplicated(reqMods)
    if(any(dub)){

        for(i in 1:sum(dub)){
            tmp <- expparnames[[which(dub)[i]]]
            expparnames[[which(dub)[i]]] <- paste0(tmp,i)
        }

    }
    
    expparlist <- lapply(parlist, function(X) X$best)

    for(i in 1:length(expparnames)){
        names(expparlist[[i]]) <- tolower(expparnames[[i]])
    }
    
    ## ADD fixed parameter setting here! ( | N=1)
    
    if(is.null(pars)) {

        if(class(expparlist)!="list"){
            
            pars <- list(expparlist)
            
        } else {

            pars <- expparlist
        }
        
    } else {
        
        if(!class(pars)=="list") {
            
            if(ordN==1){
                
                pars <- list(pars)
                
            } else{
                
                stop("parameters must be supplied as a list")
            }
        }
        
        if(length(pars)!=ordN){
            stop("Looks like not each model was supplied parameter vector")
        }
        
        parnames <- lapply(pars,function(X) names(X))
        
        ## enforce lowercase
        for(i in 1:length(parnames)){
            
            if(!is.null(pars[[i]])) names(pars[[i]]) <- tolower(parnames[[i]]) 
        }
        
        test <- tolower(unlist(parnames))%in%tolower(unlist(expparnames))
        
        if(!all(test)) {
            errorvec<-unlist(parnames)[!test]
            msg<-paste0(paste(errorvec,collapse=", "), " not recognized \n")
            stop(msg)
        }

        test <- sapply(expparlist,length)!=sapply(pars,length)
        
        if(any(test)){

            ## fill in defaults where parameters are missing
            repl <- which(test)
            for(i in 1:length(repl)){
                
                inc <- names(pars[[repl[i]]])
                expparlist[[repl[i]]][inc] <- pars[[repl[i]]]
                expparlist[[repl[i]]] -> pars[[repl[i]]]
            }
        }
        
    }

    ## parmeters should be good now - remove any duplication id
    if(any(dub)){
        for(i in 1:length(pars)){
            tmp <- names(pars[[i]])
            names(pars[[i]]) <- gsub("[0-9._//-]","",tmp)
        }
    }
    return(pars)
    
}

## function to convert a parameter list to a vector or viceversa
## depending on a formula
## for use in the orderedRTM function from bRTM so it
## is compatible with MCMC samplers
vec2list <- function(fm=NULL, pars=NULL){


    if(class(fm)!="formula"&!is.null(fm)){
        stop("first argument must be a formula (fm)")
    }
    
    if(class(pars)=="ccrtm.priors"){

        mu <- lapply(pars,"[[","mu")
        sigma <- lapply(pars,"[[","sigma")
        prior <- lapply(pars,"[[","prior")
        names(mu) <- names(prior) <- names(sigma) <- NULL
        mu <- unlist(mu)
        sigma <- unlist(sigma)
        prior <- unlist(prior)

        ## fix duplicated names
        nms <- names(mu)
        nms[duplicated(nms)] <- paste0(nms[duplicated(nms)],".",1)
        
        names(mu) <-names(sigma) <- names(prior) <- nms
        prepPriors <- list(mu,sigma,prior)
        names(prepPriors) <- c("mu","sigma","prior")
        
        return(prepPriors)
    }


    if(class(pars)=="list"){

        parvec <- unlist(pars)
        nms <- unlist(lapply(pars,names))

        tmp <- 1

        while(anyDuplicated(nms)>1){
            
            nms[duplicated(nms)] <- paste0(nms[duplicated(nms)],".",tmp)
            tmp <- tmp+1

        }

        names(parvec) <- nms

        return(parvec)
      

    } else{

        if(is.null(fm)|class(fm)!="formula"){
        stop("when pars is not a list, a valid formula must be supplied")
        }

        reqMods <- checkForm(fm)[[1]]

        if(length(reqMods)==1) {
            parlist <- list(getDefaults(reqMods))
        } else {
            parlist <- lapply(reqMods, getDefaults)
        }

        names(pars) <- tolower(names(pars))
        
        expnms <- lapply(1:length(parlist),
                         function(X)
                             data.frame(n=tolower(rownames(parlist[[X]]))
                                       ,l=X))
        expnms <- do.call(rbind,expnms)

        newparlist <- vector("list",length(parlist))

        for(i in 1:length(pars)){

            ## take fist hit 
            hit <- match(gsub("[0-9._//-]","",names(pars[i])),expnms$n)[1]

            if(length(hit)<1 | is.na(hit)) {
                stop(paste("unexpected parameter:",names(pars)[i]))
            }

            listn <- expnms[hit,]$l
            newparlist[[listn]] <- c(newparlist[[listn]],pars[i])
            expnms <- expnms[-hit,]
        }

        names(newparlist) <- reqMods
        return(newparlist)
    }
    

}
