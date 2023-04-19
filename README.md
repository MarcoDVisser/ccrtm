## ccrtm: Coupled Chain Radiative Transfer Models (0.4.1)

<!-- badges: start -->
[![cran version](http://www.r-pkg.org/badges/version/ccrtm)](http://cran.rstudio.com/web/packages/ccrtm)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/ccrtm?color=E664A4)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/ccrtm?color=333FFF)](https://github.com/metacran/cranlogs.app)
[![R-CMD-check](https://github.com/ropensci/allodb/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/ccrtm/actions)
<!-- badges: end -->

A set of radiative transfer models to quantitatively describe the absorption, reflectance and transmission of solar energy in vegatation,
and model remotely sensed spectral signatures of vegetation at distinct spatial scales. The main principle behind ccrtm is that many 
radiative transfer models can form a coupled chain, basically models that feed into each other in a linked chain (from leaf, to canopy, to stand, to atmosphere). Included models are prospect family, and 2 stream and 4 stream models for canopies, with planned inclusion atmospheric models.

Starting from version 0.4.0 inversion routines are included (backward prediction: spectra to traits). Currently, prospect5 and prospectd can be run in backward model (inversion results shown below). 

The package will slowly be extended as more models are added, and tested, and optimized. Please send requests and bug reports.

Currently the following models are implemented from scratch (and refactored in c++):
- PROSPECT5 (forward and backward mode)
- PROSPECT5B (forward and backward mode) 
- PROSPECTD (forward and backward mode)
- 4SAIL 
- 4SAIL2 (leaf angles fixed for each layer; sensu Verhoef and Bach 2007)
- 4SAIL2B (leaf angles differ for each layer; sensu Zhang et al 2005)
- FLIM
- PROSAIL (5/5B/D)
- PROSAIL2(b) (5/5B/D)
- INFORM (5/5B/D) - only lower level implementation as of yet (see below)

Implementation planned in the near future:
- SMAC
- PROCOSINE
- LIBERTY
- SOILSPECT


### Quicklinks

-   [Quick start and tutorials](#quick-start-and-tutorials)
-   [Installation](#installation)
	-   [Dependencies](#dependencies)
-   [Examples](#examples)
    -   [Examples of output](#examples-of-output)
    -	[Model inversion/Backward mode](#inversion)
-   [Inform](#inform) 
-   [Thanks](#thanks)
  

## Quick start and tutorials

A "10 minute" quickstart guide is will be implemented in due time (somewhere close to the version 1.0 release) - for now see the examples below which give an overview.

## Installation
The pacakge is on [CRAN](https://cran.rstudio.com/web/packages/ccrtm/). 
Newer developmental versions will first be available from github and can be downloaded as [zip](https://github.com/MarcoDVisser/ccrtm/zipball/master) 
or [tar ball](https://github.com/MarcoDVisser/ccrtm/tarball/master).
To install decompress these and run R CMD INSTALL on the contents of the archives, or use the **devtools** package to install the current development version from R.


```r
## devtools is required
require(devtools)
install_github("MarcoDVisser/ccrtm")
```
### Dependencies

ccrtm depends on the Rcpp, pracma and expint packages (for now).  

## Examples (versions > 0.4)

The basic functionality of ccrtm, and the coupled chain nature of the forward modelling components is shown below. Backward modelling has now been tested and included (since April 2023) for leaf models.

```r
require(ccrtm)

## setup graphics for plots 
par(mfrow=c(3,2))

## get reflectance for a leaf 
ref <- fRTM(rho~prospect5)
plot(ref,main="Prospect 5")
     
## get reflectance and transmission for a leaf 
reftrans <- fRTM(rho+tau~prospect5)
plot(reftrans,main="Prospect 5")
     
## get reflectance for a single layered canopy 
ref <- fRTM(rho~prospect5+foursail)
plot(ref,main="Prospect 5 + 4SAIL")

## get reflectance for a 2 layered canopy with two leaf types 
ref <- fRTM(rho~prospectd+prospect5+foursail2)
plot(ref,main="Prospect D + Prospect 5  + 4SAIL2")

## edit the parameters: sparse vegatation LAI 
parlist <- getDefaults(rho~prospectd+prospect5+foursail2)
parlist$foursail2["LAI"] <- 0.05

## update reflectance
ref <- fRTM(rho~prospect5+prospectd+foursail2,parlist)
plot(ref,main="LAI=0.05")

## change leaf area index to dense vegetation
parlist$foursail2["LAI"]<-8.5

## update reflectance
ref <- fRTM(rho~prospect5+prospectd+foursail2,parlist)
plot(ref,main="LAI=8.5")


```
	 
### Examples of output
The package provdes some standard plots, the output from the above example code for instance is:

![](https://i.imgur.com/alouDkJ.png)


Printing any ccrtm object will return basic information:
```r
ref <- fRTM(rho~prospect5+prospectd+foursail2)
ref
```

```
RTM predicted spectra:  reflectance 
Generating model(s):  prospect5, prospectd, foursail2 
Wavelength range  400-2500 (nm) 
```

### Inversion

Backward mode (inversion) now works for prospect5 and prospectd given leaf reflectance spectra from 400 to 2400 nm (1 nm step).

#### Inversion of a single leaf spectra

```r
require(ccrtm)

## Normally you would use your own measured leaf reflectance.
## In this example we get reflectance for a leaf with forward mode
## step 1: make a parameter list
parameters<-list(prospectd=c(N=3,Cab=40,Car=15,Cw=0.01,Cm=0.025,Canth=26,Cbrown=4))

## get reflectance for a single leaf on simulated spectra
## at the inversion requirements with wavelengths betwen 400 and 2400
ref <- fRTM(rho~prospectd,pars=parameters,wl=400:2400)   
plot(ref,main="Prospect D")

## reorder with replicates measurements over rows, and make into matrix
## this is the typical format for packages such as hsdar
refdata<-t(as.matrix(ref))

fit<-bRTM(rho~prospectd,data=refdata)
summary(fit)
     
## compare fit with simulation on log-scale so all parameter are visible
plot(parameters$prospectd,fit$mu,xlab="expected",ylab="inverted",pch=16,log="xy")
  
## add uncertainty
segments(parameters$prospectd,fit$lower.ci,
parameters$prospectd,fit$upper.ci,lwd=2)
  
## 1 to 1 line
abline(0,1)
```

```
------------------------------------------------------
Inverted  model(s):  prospectd 
------------------------------------------------------
Inverted (mean) parameters and inversion uncertainty: 
 Estimate N : 3.03 (CI: 2.7-3.4)
 Estimate Cab : 39.8 (CI: 35-45)
 Estimate Car : 12.6 (CI: 2.2-23)
 Estimate Cw : 0.00973 (CI: 0.0091-0.01)
 Estimate Cm : 0.0242 (CI: 0.022-0.026)
 Estimate Canth : 22.3 (CI: 15-29)
 Estimate Cbrown : 3.61 (CI: 2.6-4.6)
 

Head of inverted values matrix: 
             N     Cab      Car          Cw         Cm    Canth  Cbrown
[1,] 3.032301 39.7715 12.57556 0.009731867 0.02422989 22.31715 3.61153

 A total of 1 RTM inversions performed. 

```

![](https://i.imgur.com/R3msl0K.png)

#### Inversion of many leaf spectra

Users will more likely be interested in the inversion of multiple replicate spectra. Hundreds of spectra can be inverted 
within a few seconds: 

```r   
## Inversion for multiple leaf spectra
     
## As above we simulate spectra, here using a lower-level vectorized prospect function
## that rapidly simulates spectra over a matrix of input parameters
set.seed(1234)

## we simulate all spectra at once
nsim<-300 ## number of leaves
     
## Generate random leaf parameters
parmat<-cbind(
     N=runif(nsim,1,6),
     Cab=runif(nsim,5,80),
     Car=runif(nsim,1,40),
     Cw=runif(nsim,0.001,.02),
     Cm=runif(nsim,0.002,0.03)+0.01,
     Canth=runif(nsim,0,6),
     Cbrown=runif(nsim,0,4)
)
     
## simulate with the lower level prospect for rapid simulation
## of many leaves
ref<-ccrtm:::.prospectdv(parmat)[[1]][,1:2001] ## subset to 400:2400 wl
     
## invert the simulations
fit<-bRTM(rho~prospectd,data=ref)
summary(fit)
     
## check inversion performace for LMA
plot(parmat[,"Cm"],fit$mu[,"Cm"],xlab="expected",ylab="inverted",pch=16)
     
## add uncertainty
segments(parmat[,"Cm"],fit$lower.ci[,"Cm"],
          parmat[,"Cm"],fit$upper.ci[,"Cm"],lwd=2)
abline(0,1)

```

```
------------------------------------------------------
Inverted  model(s):  prospectd 
------------------------------------------------------
Inverted (mean) parameters and inversion uncertainty: 
 Estimate N : 3.5 (CI: 3.2-3.8)
 Estimate Cab : 44.8 (CI: 40-50)
 Estimate Car : 20.6 (CI: 12-30)
 Estimate Cw : 0.0103 (CI: 0.0096-0.011)
 Estimate Cm : 0.026 (CI: 0.024-0.028)
 Estimate Canth : 6.33 (CI: 2.1-13)
 Estimate Cbrown : 2.63 (CI: 1.7-3.6)
 

Head of inverted values matrix: 
             N      Cab        Car          Cw         Cm     Canth    Cbrown
[1,] 1.454031 36.75596 22.3264856 0.013613748 0.02991352  9.807317 3.6700023
[2,] 4.379834 83.06349  6.2814088 0.008684775 0.02227380  7.177197 3.1339196
[3,] 4.067006 43.62785  0.9717786 0.018761527 0.03370315  3.727199 1.6434374
[4,] 4.350049 35.88337 28.9690067 0.009501127 0.01806905  0.000000 0.8650389
[5,] 5.437518 69.62519 18.8491896 0.018808470 0.01538833 22.009023 0.0000000
[6,] 4.355071 33.46050 17.5216464 0.001314402 0.03929835  2.919029 1.1305926

 A total of 300 RTM inversions performed. 

```
![](https://i.imgur.com/iqx4gWi.png)

## Inversion methods and performance
Inversion procedure uses a multivariate neural network (MANN), a partial least squares regression (PLSR) and a Bayesian weighting model to invert from spectra to traits. The MANN and PLSR are both used because their performace differs, with the one outperforming the other depending on the parameter.   The Bayesian weighting model combines the predictions from the MANN and PLSR, and ensures that inversion uncertainty can be included (see examples above). Simulation results (below over 10000 simulations) show that the inversion routine works well. 

![](https://i.imgur.com/xsQzaew.png)


## INFORM
The INFORM functions are not exported yet as this code remain untested "formally" - although it appears to behave well. The best "formal test" is to test againts original published (often fortran) code. For INFORM this isn't possible because no code is publically available to test against (i.e. made available by the authors of INFORM) - at  least as far as I am aware of. Therefore, inform remains lower-level for now (not exported by defualt).  You can however use INFORM if you wish simply by doing:   

```r
require(ccrtm)

## to use INFORM with prospect5
informpars <-  ccrtm:::defaults.inform5()
R <- ccrtm:::rtm.inform5(informpars)

## to use INFORM with prospectd
informpars <-  ccrtm:::defaults.informd()
R <- ccrtm:::rtm.informd(informpars)


```

## to do list
- Add liberty, and procosine
- Add 2 stream model
- Add inform model to fRTM
- Add GEOSAIL model

## Thanks
Special thanks to Matteo Detto for all the suggestions while I was developing this package. 
Special thanks to (in no particular order) Zavud Baghirov and Martina Schmirl. 
Their comments helped improve documentation and remove bugs.


## Version history
- 0.01 Initial package with everything in R 
- 0.1 Basic tests complete, basic optimization, and refactoring in c++. 
- 0.1.1 Documentation update, more code refactored in c++
- 0.1.1 Documentation update, more code refactored in c++
- 0.1.4 refactored code now running via Rcpp, tested against fortran
- 0.1.5 Documentation updated, code reviewed and CRAN tests passed.
- 0.1.6 CRAN new package review issues fixed.
- 0.1.7 fixed bugs identified by Zavud, and improved documentation for leaf angle model
- 0.2.0 Start of major overhual of fRTM using model aliases, and inclusion of lower-level implementation of INFORM
- 0.3.2 Ongoing of major overhual now including a 2x faster full vectorized prospect (takes matrix of parameters)
- 0.3.3 Working aliases, inform implemented at lower-level and fully vectorized prospect 
- 0.3.4 Expected: functions accept matrix of parameters for vectorized prospect 
- 0.4.0 Leaf models prospect5 and prospectd can now be inverted, and inversion uncertainty is returned.
- 0.4.1 Documentation updated for backward modelling
 

