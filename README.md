## ccrtm: Coupled Chain Radiative Transfer Models 

A set of radiative transfer models to quantitatively describe the absorption, reflectance and transmission of solar energy in vegatation,
and model remotely sensed spectral signatures of vegetation at distinct spatial scales. The main principle behind ccrtm is that many 
radiative transfer models can form a coupled chain, basically models that feed into each other in a linked chain (from leaf, to canopy, to stand, to atmosphere). Included models are prospect family, and 2 stream and 4 stream models for canopies, with planned inclusion atmospheric models.

The package will slowly be extended as more models are added, and tested, and optimized. Please send requests and bug reports.

[![cran version](http://www.r-pkg.org/badges/version/ccrtm)](http://cran.rstudio.com/web/packages/ccrtm)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/ccrtm?color=E664A4)](https://github.com/metacran/cranlogs.app)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/ccrtm?color=333FFF)](https://github.com/metacran/cranlogs.app)

### Quicklinks

-   [Quick start and tutorials](#quick-start-and-tutorials)
-   [Installation](#the-online-code-files-from-s1-text)
	-   [Dependencies](#dependencies)
-   [Examples](#examples)
    -   [Examples of output](#examples-of-output)
-   [Thanks](#thanks)
  

## Quick start and tutorials

A "10 minute" quickstart guide is will be implemented in due time. 

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

ccrtm depends on the pracma and expint packages (for now). 

## Examples

The basic functionality of ccrtm, and the coupled chain nature of the forward modelling components is shown below.
Backward modelling will be included in later stages.

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
parlist<- list(prospect5=NULL,prospectd=NULL,foursail2=c(LAI=0.05))

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

## to do list
- Add liberty, and procosine
- Add 2 stream model
- Add inform model
- Add GEOSAIL model

## Thanks
Special thanks to Matteo Detto for all the suggestions while I was developing this package. 
Special thanks to (in no particular order) Zavud Baghirov and Martina Schmirl. 
Their comments helped improve documentation and remove bugs.


## Version history
0.01 Initial package with everything in R 
0.1 Basic tests complete, basic optimization, and refactoring in c++. 


