#' S3- methods for Generate defaults settings and parameters
#' for all supported models. See ?ccrtm for details.
#' 
#' @param model a ccrtm formula (rho ~ prospectd) or character vector of modelnames
#' (e.g. "prospect5")
#' @param \\dots not used.
#' 
#' @return a data.frame with default model parameters
#' @export
getDefaults<-function(model=NULL, ...){

    rtmModels <- getModels()

  if(class(model)=="formula"){

    reqMods <- getAlias(model)

        ## check models
        if(any(!reqMods%in%rtmModels$model)){

            stop(reqMods[!reqMods%in%rtmModels$model],
                 " not implemented")
        }
    } else if(class(model)=="character"){

        test <- any(sapply(model,function(X) !X%in%rtmModels$model))
        if(test) stop("model not recognized")

        reqMods <- model
        }

    def <- vector("list",length(reqMods))

    for(i in seq_len(length(reqMods))){

        mod <- reqMods[i]
        class(mod) <- reqMods[i]
        def[[i]] <- defaults(mod)
    }

    names(def) <- reqMods
    if(length(def)==1) return(def[[1]])

    return(def)
}

## defualts function
defaults <- function(x, ...){
    UseMethod("defaults",x)
}

## get some default values for the prospect model 
##      - N     = leaf structure parameter
##      - Cab   = chlorophyll a+b content 
##      - Car   = carotenoids content 
##      - Cbrown= brown pigments content in arbitrary units
##      - Cw    = equivalent water thickness
##      - Cm    = dry matter content 
defaults.prospect5<-function(x,simple=TRUE){
   
    ## typical values are the rice values
    typical = c(2.698,70.8,20, 0.0117,0.009327,0.001)
    ## N, Cab,Car, Cw, Cm, Cbrown, resid
    omega <- c(5,250,80,.4,0.3,10) 
    alpha <- c(1,10,1,1e-5,1e-5,0) 
    
    names(typical)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown")
    names(omega)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown")
    names(alpha)<-c("N", "Cab","Car", "Cw", "Cm", "Cbrown")

    if(simple){
         
      return(list("prospect5"=typical))
   
    } else {
        
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega
    
    return(list("prospect5"=def))
    }
}

## get some default values for the prospect model 
##       - N     = leaf structure parameter
##       - Cab   = chlorophyll a+b content 
##       - Car   = carotenoids content 
##       - Canth = Leaf anthocyanin content
##       - Cbrown= brown pigments content in arbitrary units
##       - Cw    = equivalent water thickness
##       - Cm    = dry matter content 
##
defaults.prospectd<-function(x,simple=TRUE){
    
    
    ## typical values are the rice values
    typical = c(2.698,70.8,20, 0.0117,0.009327,5,0.001)
    ## N, Cab,Car, Cw, Cm, Cbrown, resid
    omega <- c(5,250,80,.4,0.3,50,10) 
    alpha <- c(1,10,1,1e-5,1e-5,0,0) 
    
    names(typical)<-c("N", "Cab","Car", "Cw", "Cm","Canth", "Cbrown")
    names(omega)<-c("N", "Cab","Car", "Cw", "Cm", "Canth","Cbrown")
    names(alpha)<-c("N", "Cab","Car", "Cw", "Cm","Canth","Cbrown")

    if(simple){
        
      return(list("prospectd"=typical))
    
    } else {

        def = data.frame(best = typical)
        def$lower = alpha
        def$upper = omega

      return(list("prospectd"=def))
    }
}




## Generate default values for the foursail model
##       - psoil = linear mixing weight for 2 soil
##                 reflectance types
##       - LAI   =  Leaf area index
##       - TypeLidf= Type of leaf angle distribution. 
##       - lidfa = Leaf angle distribution parameter a. 
##       - lidfb = Leaf angle distribution parameter b.
##       - hspot = Hotspot parameter
##       - tts   = Solar zenith angle (degrees)
##       - tto   = Observer zenith angle (degrees)
##       - psi   = Relative azimuth angle (degrees)
defaults.foursail<-function(x,simple=TRUE){

    nms <- c("psoil","LAI","TypeLidf","lidfa","lidfb",
             "hspot","tts","tto","psi")
#             "sigma")
    
    typical = c(0,4.77,1,-0.35,-0.15,
                0.01,30,10,0)
#                0.1)
    
    omega <- c(1,50,1,1,1,
               1,90,90,180)
#               4) 
    alpha <- c(0,0.1,1,-1,-1,
               1e-6,0,0,0)
#               1e-8) 
    
    names(typical)<-nms
    names(omega)<-nms
    names(alpha)<-nms

    if(simple){
        
      return(list("foursail"=typical))
    
    } else {
    
    def = data.frame(best = typical)
    def$lower = alpha
    def$upper = omega

  return(list("foursail"=def))
      
    }
}   

## Generate default values for the foursail2 model 
##       - psoil = psoil: Dry/Wet soil factor
##       - LAI   =  Leaf area index
##       - TypeLidf= Type of leaf angle distribution. 
##       - lidfa = Leaf angle distribution parameter a. 
##       - lidfb = Leaf angle distribution parameter b.
##       - Cv    = vertical crown coverage fraction.
##      - Zeta0 = tree shape factor (D/H).
##      - fb    = fraction secondary particles (fraction "brown").
##       - diss  = Canopy dissociation factor for secondary particles.
##       - skyl  = linear mixing weight for diffuse & direct radiation.
##                 (ratio of diffuse to total incident radiation)
##       - hspot = Hotspot parameter.
##       - tts   = Solar zenith angle (degrees).
##       - tto   = Observer zenith angle (degrees).
##       - psi   = Relative azimuth angle (degrees).
defaults.foursail2<-function(x,simple=TRUE){

    nms <- c("psoil","LAI","TypeLidf","lidfa","lidfb",
             "hspot",
             "Cv","Zeta","fb","diss",
             "tts","tto","psi")
#             "sigma")
    
    typical  <-  c(0,8.77,1,-0.35,-0.15,
                   0.01,                   
                   0.95,1/3,0.5,0.5,
                   30,10,0)
#                0.1)
    
    omega <- c(1,50,1,1,1,
               1,
               1,1,1,1,
               90,90,180)
#               4) 
    alpha <- c(0,0.1,1,-1,-1,
               1e-6,
               0,0,0,0,
               0,0,0)
#               1e-8) 
    
    names(typical)<-nms
    names(omega)<-nms
    names(alpha)<-nms

    if(simple){
        
      return(list("foursail2"=typical))

    } else {
    
        def = data.frame(best = typical)
        def$lower = alpha
        def$upper = omega
        
        
      return(list("foursail2"=def))
    }
    
}


## Generate default values for the foursail2b model 
##       - psoil = psoil: Dry/Wet soil factor
##       - LAI   =  Leaf area index
##       - TypeLidf= Type of leaf angle distribution. 
##       - lidfa = Leaf angle distribution parameter a. 
##       - lidfb = Leaf angle distribution parameter b.
##      - Cv    = vertical crown coverage fraction.
##      - Zeta0 = tree shape factor (D/H).
##      - fb    = fraction secondary particles (fraction "brown").
##       - diss  = Canopy dissociation factor for secondary particles.
##       - skyl  = linear mixing weight for diffuse & direct radiation.
##                 (ratio of diffuse to total incident radiation)
##       - hspot = Hotspot parameter.
##       - tts   = Solar zenith angle (degrees).
##       - tto   = Observer zenith angle (degrees).
##       - psi   = Relative azimuth angle (degrees).
defaults.foursail2b<-function(x,simple=TRUE){

    nms <- c("psoil","LAI","TypeLidf","lidfa","lidfb",
             "hspot",
             "Cv","Zeta","fb","diss",
             "tts","tto","psi")
#             "sigma")
    
    typical  <-  c(0,8.77,2,35,45,
                   0.01,                   
                   0.95,1/3,0.5,0.5,
                   30,10,0)
#                0.1)
    
    omega <- c(1,50,2,0,0,
               1,
               1,1,1,1,
               90,90,180)
#               4) 
    alpha <- c(0,0.1,2,90,90,
               1e-6,
               0,0,0,0,
               0,0,0)
#               1e-8) 
    
    names(typical)<-nms
    names(omega)<-nms
    names(alpha)<-nms

    if(simple){
        
      return(list("foursail2b"=typical))
        
    } else {
        
        def = data.frame(best = typical)
        def$lower = alpha
        def$upper = omega
        
        return("foursail2b"=def)
    }
}

## get some default values for the prospect model 
##       - skyl     = diffuse pattern
defaults.skyl<-function(x,simple=TRUE){
   
    typical = c(0.1)
    omega <- c(1) 
    alpha <- c(0) 
    
    names(alpha)<-names(omega)<-names(typical)<-c("skyl")

    if(simple){
        
      return(list("skyl"=typical))
       
    } else {
        
        def = data.frame(best = typical)
        def$lower = alpha
        def$upper = omega

    return(list("skyl"=def))

    }
}


## Generate default values for the flim model
#' d = stand density (d)  
#' cd = crown diameter (cd) 
#' h = mean crown height (h)
#' lai = leaf area index (lai)  
#' alpha = light extinction coefficient (alpha)  
#' tts = Solar zenith angle (tts)
#' tto = Observer zenith angle (tto)
#' psi = Sun-sensor azimuth angle (psi)
defaults.flim<-function(x,simple=TRUE){

    nms <- c("cd", "h", "d", "area","tto","tts","psi")
    
    typical = c(10,20,30,10000,0,15,45)
    omega <- c(Inf,Inf,Inf,Inf,90,90,180) 
    alpha <- c(0,0,0,0,0,0,0) 
    
    names(alpha)<-names(omega)<-names(typical)<-nms

    if(simple){
        
      return(list("flim"=typical))
        
    } else {
    
        def = data.frame(best = typical)
        def$lower = alpha
        def$upper = omega
        
      return(list("flim"=def))
        
    }
}


## Generate default values for the inform-5 model
defaults.inform5<-function(x,simple=TRUE){

    ## make parameter list
    leafpars <- getDefaults("prospect5",simple)[[1]]
    canopypars <- getDefaults("foursail",simple)[[1]]
    flimpars <- getDefaults("flim",simple)[[1]]
    skylpars <- getDefaults("skyl",simple)[[1]]

    def <- list("prospect5"=list("canopy"=leafpars,
                                 "understorey"=leafpars),
                "foursail"=list("canopy" = canopypars,
                                 "understorey" = canopypars),
                "flim"=flimpars,
                "skyl"=skylpars
                )
                
    return(def)
}


## Generate default values for the inform-d model
defaults.informd<-function(x,simple=TRUE){

    ## make parameter list
    leafpars <- getDefaults("prospectd",simple)[[1]]
    canopypars <- getDefaults("foursail",simple)[[1]]
    flimpars <- getDefaults("flim",simple)[[1]]
    skylpars <- getDefaults("skyl",simple)[[1]]

    def <- list("prospectd"=list("canopy"=leafpars,
                                 "understorey"=leafpars),
                "foursail"=list("canopy" = canopypars,
                                 "understorey" = canopypars),
                "flim"=flimpars,
                "skyl"=skylpars
                )
    return(def)
}


## Generate default values for the prosail5 model
defaults.prosail5<-function(x,simple=TRUE){

  ## make parameter list
  leafpars <- getDefaults("prospect5",simple)[[1]]
  canopypars <- getDefaults("foursail",simple)[[1]]

  def <- list("prospect5"=leafpars,
              "foursail"=canopypars
              )

  return(def)
}

## Generate default values for the prosaild model
defaults.prosaild<-function(x,simple=TRUE){

  ## make parameter list
  leafpars <- getDefaults("prospectd",simple)[[1]]
  canopypars <- getDefaults("foursail",simple)[[1]]

  def <- list("prospectd"=leafpars,
              "foursail"=canopypars
              )

  return(def)
}

## Generate default values for the prosail5 model
defaults.prosail2_55<-function(x,simple=TRUE){

  ## make parameter list
  leafparsa <- getDefaults("prospect5",simple)[[1]]
  leafparsb <- getDefaults("prospect5",simple)[[1]]
  canopypars <- getDefaults("foursail2",simple)[[1]]

  def <- list("prospect5.a"=leafparsa,
              "prospect5.b"=leafparsb,
              "foursail2"=canopypars
              )

  return(def)
}

## Generate default values for the prosail5 model
defaults.prosail2_dd<-function(x,simple=TRUE){

  ## make parameter list
  leafparsa <- getDefaults("prospectd",simple)[[1]]
  leafparsb <- getDefaults("prospectd",simple)[[1]]
  canopypars <- getDefaults("foursail2",simple)[[1]]

  def <- list("prospectd.a"=leafparsa,
              "prospectd.b"=leafparsb,
              "foursail2"=canopypars
              )

  return(def)
}


## Generate default values for the prosail5 model
defaults.prosail2_5d<-function(x,simple=TRUE){

  ## make parameter list
  leafparsa <- getDefaults("prospect5",simple)[[1]]
  leafparsb <- getDefaults("prospectd",simple)[[1]]
  canopypars <- getDefaults("foursail2",simple)[[1]]

  def <- list("prospect5"=leafparsa,
              "prospectd"=leafparsb,
              "foursail2"=canopypars
              )

  return(def)
}

## Generate default values for the prosail5 model
defaults.prosail2_d5<-function(x,simple=TRUE){

  ## make parameter list
  leafparsa <- getDefaults("prospectd",simple)[[1]]
  leafparsb <- getDefaults("prospect5",simple)[[1]]
  canopypars <- getDefaults("foursail2",simple)[[1]]

  def <- list("prospectd"=leafparsa,
              "prospect5"=leafparsb,
              "foursail2"=canopypars
              )

  return(def)
}

## Generate default values for the prosail5 model
defaults.prosail2b_55<-function(x,simple=TRUE){

  ## make parameter list
  leafparsa <- getDefaults("prospect5",simple)[[1]]
  leafparsb <- getDefaults("prospect5",simple)[[1]]
  canopypars <- getDefaults("foursail2b",simple)[[1]]

  def <- list("prospect5.a"=leafparsa,
              "prospect5.b"=leafparsb,
              "foursail2b"=canopypars
              )

  return(def)
}

## Generate default values for the prosail5 model
defaults.prosail2b_dd<-function(x,simple=TRUE){

  ## make parameter list
  leafparsa <- getDefaults("prospectd",simple)[[1]]
  leafparsb <- getDefaults("prospectd",simple)[[1]]
  canopypars <- getDefaults("foursail2b",simple)[[1]]

  def <- list("prospectd.a"=leafparsa,
              "prospectd.b"=leafparsb,
              "foursail2b"=canopypars
              )

  return(def)
}


## Generate default values for the prosail5 model
defaults.prosail2b_5d<-function(x,simple=TRUE){

  ## make parameter list
  leafparsa <- getDefaults("prospect5",simple)[[1]]
  leafparsb <- getDefaults("prospectd",simple)[[1]]
  canopypars <- getDefaults("foursail2b",simple)[[1]]

  def <- list("prospect5"=leafparsa,
              "prospectd"=leafparsb,
              "foursail2b"=canopypars
              )

  return(def)
}

## Generate default values for the prosail5 model
defaults.prosail2b_d5<-function(x,simple=TRUE){

  ## make parameter list
  leafparsa <- getDefaults("prospectd",simple)[[1]]
  leafparsb <- getDefaults("prospect5",simple)[[1]]
  canopypars <- getDefaults("foursail2b",simple)[[1]]

  def <- list("prospectd"=leafparsa,
              "prospect5"=leafparsb,
              "foursail2b"=canopypars
              )

  return(def)
}


