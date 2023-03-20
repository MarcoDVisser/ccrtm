### TO DO: 
## 1) internal prediction function like fit.lm so bRTM can call fRTM 
## without the need of double parameter and model checks - this  
## could be done via and orderedRTM prediction function.

## Species actions
## rho + tau ~ sp(foursail, rho=ref,tau=trans,bgr=rsoil) + geo

#' Forward implementation of coupled Radiative Transfer Models.
#'
#' @param fm A formula specifying which rtm to run (see details).
#' @param pars a _named_ list of _named_ parameter vectors for all models.
#' The parameter list for a model call as rho ~ prospect + foursail,
#' therefore, contains two vectors: the first with parameters for
#' prospect and the second with parameters for foursail.
#' See ?getDefaults for an example of a parameter list.
#' If left empty default parameters are generated.
#' @param wl wavelengths (in nm) add only if
#' certain wavelengths are required as output.
#' Input is expected to integers between 400 and 2500, or will 
#' be forced to be an integer. Integers outside the 400:2500 range
#' will not be returned.
#' 
#' @details
#'
#' In general the form of the formula specifies the
#' both the expected output (left hand) and the different models
#' you would like to couple to generate the output (right hand).
#'
#' At current the following radiative transfer models are
#' implemented
#'
#' | Example of a formula                     | Model     |
#' | :--------------------------------------  | :-------: |
#' | rho ~ prospect5                          | prospect5 |
#' | rho ~ prospectd                          | prospectd |
#' | rho ~ prospect5 + foursail               | PROSAIL   |
#' | rho ~ prospect5 + foursail               | PROSAIL   |
#' | rho ~ prospectd + foursail2              | PROSAIL   |
#' | rho ~ prospectd + prospect5 + foursail2  | PROSAIL2* |
#' | rho ~ prospectd + foursail2b             | PROSAIL2b*|
#' | rho ~ prospectd + foursail + flim + sky  | INFORM**  |
#' | rho ~ prospectd + foursail + flim        | INFORM**  |
#'
#' * Note that the formula "rho ~ prospectd + foursail2" is the
#' same as "rho ~ prospectd + prospectd + foursail2" and
#' both expect a names list of 3 parameter vectors 
#' (leaf type 1, leaf type 2, and the canopy parameters).
#'
#' ** INFORM is currently restricted to a lower-level function
#' only. See the cctrm github readme page on how to use it.
#'
#' In the above examples with additive components, prospectd can
#' be replaced with prospect5 - or vice versa - if so desired.
#' See the help files for details on each right hand component.
#' For instance, ?foursail provides more elaboration on the 4SAIL model
#' and gives an example for lower-level implementations of each component
#' model.
#'
#' Tranmission can also be returned if specified in
#' the left-hand component of the formula:
#'
#' | Formula                             | Model     |
#' | :---------------------------------  | :-------: |
#' | rho + tau ~ prospect5               | prospect5 |
#' | rho + tau ~ prospectd + foursail    | 4SAIL     |
#'
#' The examples above indicate that the users wishes to predict
#' transmission next to reflectance. More specifically, The first
#' returns leaf reflectance and transmission while
#' the second returns 4 components of canopy reflectance and
#' canopy transmission in the solar and viewing direction.
#'
#' More details are given in ?cctrm.
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
#' ## edit the parameters: sparse vegetation LAI
#' parlist <- getDefaults(rho~prospectd+prospect5+foursail2)
#' parlist$foursail2["LAI"] <- 0.05
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
#' @md
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
  pars <- checkPars(pars,fm,ordN)

  ## get desired predictions
  tmp<-getPredictions(fm, ordN)

  returnCol<-tmp[[1]]
  pred<-tmp[[2]]

  alias<-getAlias(fm) ## get model alias

  class(pars) <- alias ## set alias

  final <- runRTM(pars)[rowInc,returnCol]

  ## finalize output
  attr(final,"main")<-paste0(pred, " ")
  attr(final,"wavelength")<-wl
  attr(final,"models") <- unique(reqMods)
  class(final) <- c(class(final),"rtm.spectra")
  return(final)

  }



## get prediction columns
getPredictions<-function(fm,ordN){

  pred <- as.character(fm[[2]])
  pred <- pred[!grepl("\\+",as.character(pred))]

  if(!any(c("rho","tau")%in%pred)) stop("Unknown predictions requested")

  if(length(pred)>1){

    if(ordN>1){
      returnCol <- 1:6 ## return reflectance and transmittance
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
      returnCol <- 5:6 ## return transmittance
    } else {
      returnCol <- 2 ## return transmittance
    }

  }

  return(list(returnCol, pred))

}


## Get list of implemented models
getModels <- function(){

  mods <- c("prospect5","prospectd",
            "foursail","prosail5","prosaild",
            "foursail2","prosail2_55","prosail2_dd","prosail2_d5","prosail2_5d",
            "foursail2b","prosail2b_55","prosail2b_dd","prosail2b_d5","prosail2b_5d",
            "flim",
            "inform5",
            "informd",
            "skyl"
             )

   data.frame(model=mods,
              level=c(1,1,
                      2,2,2,
                      3,3,3,3,3,
                      3,3,3,3,3,
                      3,
                      3,
                      3,
                      4),
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

## new alias method for model prediction
getAlias <- function(fm){

  request <- as.character(fm)[[3]]

  models  <-  c("prospect5","prospectd",
                "prospect5 + foursail","prospectd + foursail",
                "prospect5 + foursail2","prospectd + foursail2",
                "prospect5 + prospect5 + foursail2","prospectd + prospectd + foursail2",
                "prospectd + prospect5 + foursail2","prospect5 + prospectd + foursail2",
                "prospect5 + foursail2b","prospectd + foursail2b",
                "prospect5 + prospect5 + foursail2b","prospectd + prospectd + foursail2b",
                "prospectd + prospect5 + foursail2b","prospect5 + prospectd + foursail2b",
                "prospect5 + foursail + flim","prospectd + foursail + flim"
                )

  alias <- c("prospect5","prospectd",
             "prosail5","prosaild",
             "prosail2_55","prosail2_dd","prosail2_55","prosail2_dd","prosail2_d5","prosail2_5d",
             "prosail2b_55","prosail2b_dd","prosail2b_55","prosail2b_dd","prosail2b_d5","prosail2b_5d",
             "inform5","informd"
             )

  if(!any(models%in%request)){
    stop("Requested model not implemented. Did you maybe mean \"",
         models[which.min(adist("propect + forisail2",models))], "\"?")
  }

  ind <- which(models%in%request)

  if((alias[ind]=="inform5")|(alias[ind]=="informd")){
    stop("INFORM only available as lower-level model, see ccrtm github page")
  }

  alias[ind]
}



#' Function to check and return parameters
#' @importFrom utils adist
checkPars <- function(pars,fm,ordN){

  ## generate expected parameters
  expparlist <- getDefaults(fm)
  expparnames <- lapply(expparlist, function(X) names(X))
  expmodnames <- names(expparlist)

  for(i in 1:length(expparnames)){
    names(expparlist[[i]]) <- tolower(expparnames[[i]])
  }

  if(is.null(pars)) {

    if(class(expparlist)!="list"){
      pars <- list(expparlist)

    } else {

      pars <- expparlist
    }

  } else {

    if(!class(pars)=="list") {

      stop("parameters must be supplied as a named list")

    }

    if(any(!expmodnames%in%names(pars))) {
      msg<-paste0(paste(expmodnames,collapse=", "), " expected models in parameter list but not all found. See ?getDefaults. \n")
      stop(msg)
    }

    parnames <- lapply(pars,function(X) names(X))

    ## enforce lowercase
    for(i in 1:length(parnames)){

      if(!is.null(pars[[i]])) names(pars[[i]]) <- tolower(parnames[[i]])
    }

    ## exception 1,2,3: too few, too many or unknown.
    if(length(unlist(parnames))<length(unlist(expparnames))){

      test <-  tolower(unlist(expparnames))%in%tolower(unlist(parnames))

      if(!all(test)) {
        errorvec<-unlist(expparnames)[!test]
        msg<-paste0(paste(errorvec,collapse=", "), " expected parameters not found. See ?getDefaults. \n")
        stop(msg)
      }
    } else if(length(unlist(parnames))>length(unlist(expparnames))){

      test <-tolower(unlist(parnames))%in%tolower(unlist(expparnames))

      if(!all(test)) {
        errorvec<-unlist(parnames)[!test]
        msg<-paste0(paste(errorvec,collapse=", "), " parameters supplied but not expected. See ?getDefaults. \n")
        stop(msg)
      }
    } else {

      test <- tolower(unlist(expparnames))%in%tolower(unlist(parnames))

      if(!all(test)) {
        errorvec<-unlist(parnames)[!test]
        msg<-paste0(paste(errorvec,collapse=", "), " parameters not recognized. See ?getDefaults. \n")
        stop(msg)
      }


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
