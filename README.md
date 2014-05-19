# Welcome to EddyPro&reg;

EddyPro&reg; is a powerful open source software application for processing eddy
covariance data. It computes fluxes of water vapor (evapotranspiration), carbon
dioxide, methane, other trace gases, and energy with the Eddy Covariance method.

EddyPro is developed, maintained and
supported by [LI‑COR Biosciences](www.licor.com). It originates from
[ECO<sub>2</sub>S](http://gaia.agraria.unitus.it/eco2s),
the Eddy COvariance COmmunity Software project, which was developed as part
of the Infrastructure for Measurement of the European Carbon Cycle (IMECC-EU)
research project. We gratefully acknowledge the
[IMECC](http://imecc.ipsl.jussieu.fr/index.html) consortium,
the ECO<sub>2</sub>S development team, the [University of Tuscia](www.unitus.it) (Italy)
and scientists around the world who assisted with development and testing of
the original version of this software.

![EddyPro](img/app-logo-small.png)

## Overview

You can use the EddyPro installation program to install the following components:

- EddyPro Engine, the core engine, a command line software
- EddyPro GUI, the graphical user interface to configure and run the engine
- EddyPro Help, consisting in a local context-sensitive help system, a local
    html version of the online EddyPro Help system, the 'Getting Started' and
    'User Guide' documents in PDF format

Both are pre-built for a particular environment (operating system and compiler).


##Source Code Repository

EddyPro is a fully cross-platform application, which consists of a set of
command line programs and a graphical user interface (GUI).

The source code is developed using two independent Git repositories, namely:


  - git-engine
  - git-gui


## Installing EddyPro

You can download EddyPro from the LI-COR
[EddyPro website](http://www.licor.com/env/products/eddy_covariance/software.html).
The site provides download links for all supported platforms.

Start the installation program like any executable on the development platform.
Select the components that you want to install and follow the instructions of
the installation program to complete the installation.


## Building EddyPro from source

To build EddyPro follow these instructions:


### Engine

To compile the Engine use [gfortran](https://gcc.gnu.org/wiki/GFortran) (The GNU Fortran compiler) and run:

    $ make Makefile_rp
    $ make Makefile_fcc

### GUI

To compile the GUI:

1. install the [Qt framework](http://qt-project.org/)
2. open the `eddypro.pro` project file with QtCreator
3. build the project


## Utilities

To successfully run Eddypro, the program installation folder must contain the
following command line utilities under the 'bin' sub-directory:

- 7-zip
- pausep


### 7-zip

[7-Zip](http://www.7-zip.org/) is a file archiver.

The console application consists of two files:
- 7z.dll
- 7z.exe

License: [LGPL](https://www.gnu.org/licenses/lgpl.html).

### pausep

Pausep it's a Win32 process suspend/resume tool, available on
[Code Project](http://www.codeproject.com/Articles/2964/Win-process-suspend-resume-tool).

It consists of one file:
- pausep.exe

License: Code Project Open License,
[CPOL](http://www.codeproject.com/info/cpol10.aspx).


## Using EddyPro sample data

You can run EddyPro using sample data files available in the LI-COR
[EddyPro website](http://www.licor.com/env/products/eddy_covariance/software.html).


## Data Processing Options in EddyPro

+ Axis rotation for sonic anemometer tilt correction
  - Double rotation
  - Triple rotation
  - Sector-wise planar fit (Wilczak et al., 2001)
  - Sector-wise planar fit with no velocity bias (van Dijk et al., 2004)


+ Detrending of raw time series
  - Block averaging
  - Linear detrending
  - Running mean
  - Exponential running mean


+ Compensation of time lag between sonic anemometer and gas analyzer measurements
  - Automatic time lag optimization (optionally as a function of RH for H<sub>2</sub>O)
  - Maximum covariance with default (circular correlation)
  - Maximum covariance without default
  - Constant
  - None (option to not apply compensation)


+ Statistical tests for raw time series data (Vickers and Mahrt, 1997)
  - Spike count/removal
  - Amplitude resolution
  - Dropouts
  - Absolute limits
  - Skewness and kurtosis
  - Discontinuities
  - Time lags
  - Angle of attack
  - Steadiness of horizontal wind
  - Individually selectable and customizable


+ Compensation for air density fluctuations
  - Webb et al., 1980 (open path) / Ibrom et al., 2007a (closed path)
  - Use (or convert to) mixing ratio (Burba et al., 2012)
  - Optional off-season upatake correction for LI-7500 (Burba et al., 2008)
  - None (option to not apply compensation)


+ Correction for frequency response (attenuation)
  - Analytic high-pass filtering correction (Moncrieff et al., 2004)
  - Low-pass filtering, select and configure:
  - Moncrieff et al. (1997)
  - Horst (1997)
  - Ibrom et al. (2007b)
  - Horst and Lenschow (2009)
  - Fratini et. al. (2012)


+ Quality control tests for fluxes according to Foken et al. (2004)
  - Flagging according to Carbo Europe standard (Mauder and Foken, 2004)
  - Flagging according to Foken (2003)
  - Flagging after Göckede et al. (2004)


+ Random uncertainty estimation
  - Mann and Lenschow (1994)
  - Finkelstein and Sims (2001)


+ Flux footprint estimation
  - Kljun et al. (2004)
  - Kormann and Meixner (2001)
  - Hsieh et al. (2000)


+ Other options applied in both Express and/or Advanced Mode include:
  - Sonic temperature correction for humidity following van Dijk et al. (2004)
  - Spectroscopic correction for LI-7700 following McDermitt et al. (2011)
  - Angle of attack corrections for Gill anemometers following Nakai et al. (2006)
  - Angle of attack corrections for Gill anemometers following Nakai and Shimoyama (2012)
  - sInclusion of biomet data for improved flux computation/correction


+ Available outputs
  - Full (rich) output with fluxes, quality flags and much more (standard format
    or available results only)
  - Ameriflux format
  - GHG Europe format
  - Raw data statistics
  - Full length spectra and co-spectra
  - Binned spectra and co-spectra
  - Binned ogives
  - Ensemble averaged spectra
  - Ensemble averaged cospectra, fitted models and ideal (Kaimal) cospectra
  - Details of steady state and turbulence tests
  - Raw data time series after each statistical tests/correction
  - Averaged biomet data


## EddyPro Trademark and Logo Policy

In order to help users who want to cite EddyPro on posters or publications,
LI-COR provides [guidelines](docs/EddyPro_Trademark_Policy.pdf) for the proper
use of the EddyPro wordmark and logo.


## Want to Know More?

More information is available at:

  - http://www.licor.com/env/products/eddy_covariance/software.html
  - http://envsupport.licor.com/help/EddyPro5/index.htm

Be sure to check out the
'[What's new](http://envsupport.licor.com/help/EddyPro5/index.htm#Whats_New.htm)'
page, which will list any known problems or limitations of the current and
past versions.


See also the CHANGELOG available in the /docs folder.

---

We hope you will enjoy using EddyPro!
