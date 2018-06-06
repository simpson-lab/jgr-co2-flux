# Generalized additive models of climatic and metabolic controls of sub‐annual variation in pCO2 in productive hardwater lakes

* Code Repository DOI: [![Repository DOI](https://zenodo.org/badge/125887868.svg)](https://zenodo.org/badge/latestdoi/125887868)
* Paper DOI: [10.1029/2018JG004506](https://doi.org/10.1029/2018JG004506)
* Authors: E. Wiik, H. A. Haig, N. M. Hayes, K. Finlay, G. L. Simpson, R. J. Vogt, & P. R. Leavitt
* Journal of Geophysical Research: Biogeosciences

## Acces to private data

Some of the data used here are long-term monitoring data from the Qu'Appelle Long-term Ecological Research Program that are not at this time publicly available. These data are stored privately at [github.com/simpson-lab/jgr-co2-flux-private-data](https://github.com/simpson-lab/jgr-co2-flux-private-data). Access to this repository can be requested by emailing Peter Leavitt (Peter.Leavitt@uregina.ca).

## Original code repository

The original repository for all R code on the paper can be found at [https://github.com/ewiik/flux](https://github.com/ewiik/flux), where code may still be changing due to it being used for multiple papers.

## What was done for this archive version and paper

All code for data manipulation and plots for producing this paper were rerun on 10.07.2017 based on a copy of the active repository, and then again on 15.02.2018 following revisions, modified as follows:

### `docs/`

All private and open docs removed apart from a copy of `makefile-prep.txt`, renamed `metadata-code.md` for this archive and fleshed out

### `data/`

Data not relating to the paper was removed, and manually created
files were retained:

* `maunaloa2014summer.csv`,
* `weatherstations.csv`,
* `superstationz.csv`,
* `dataforallteleconnectionsuptoFeb2015.csv`

private original data was retained, R-generated data was removed

### code rerun
Code was rerun from scratch, following order outlined in `metadata-code.md`. New supplementary figures following revision are created in `flux-appendix-figs.R`. Typos of units of oxygen and DOC in one plot were corrected.

## Metadata information

### data/

`weatherstations.csv`: List of meteorological stations used for the
study sites and information on the timeline of coverage for each station

`superstationz.csv`: Which meteorological station applies to which study site

`maunaloasummer2014.csv`: Quick grab from online as other maunaloa information,
just for the summer of 2014.

`dataforallteleconnectionsuptoFeb2015.csv`: file received from colleague
	on climate indices covering years of survey data

### data/private
Private files sourced in R codes contain the following information:

* `Chl.csv`: chlorophyll data from all study sites over the study period
* `db_metadata_outliers.csv`: file which identifies true data outliers (errors)
* `fluxquery1.csv`: database query output for wind, DIC and pH for the study sites
* `fluxquery2.csv`: database query output for lake name, lake abbreviation, and lake elevation
* `fluxquery3.csv`: database query output for temperature, conductivity and salinity for
	study sites
* `IceOut.csv`: available ice-out dates for study sites
* `kerri_co2flux_results*.csv`: CO2 flux calculated by K Finlay for cross-comparison &
	reproducibility checks
* `Lakes.csv`: database query summary table of study sites (drainage area,
	mean depth, max depth, etc.)
* `Monthly_Inflow_Estimates2.csv`, `Lake_Inflows.csv`: inflow estimates for study sites
* `pCO2_vars_yearlyaverageswithclimate`: yearly summaries of previous analyses
* `Production.csv`: database query output for bottle incubation data (respiration,
	gross primary production, net primary production)
* `qDICDOCupdate.csv`: updated database query output for DIC for 2013 and 2014
* `qpco2supportdata.csv`: database query output for dissolved nutrients etc.
* `qprodsupportdata.csv`: database query output for supporting information to calculate
	bottle incubation metrics
* `qprofilesextrarecordsoxtemp.csv`, `qprofilesoxtemp.csv`: database query output for oxygen
	and temperature for study sites
* `qsecchietal.csv`: database query output for secchi and wind
* `salcalc.R`, `ysi85constants.csv`: script created to correct false salinity readings in the data using
	equations and constants obtained from YSI (sonde manufacturer)

### docs/

`metadata-code.md`: Shows the order in which to execute all R scripts from scratch,
and also outlines the dependencies and outputs of each script

## Archive of private original and private output files

All original private files used on 2017-07-10 as well as 2018-02-15 have been archived privately by Emma Wiik/Gavin Simpson. Figures generated on 2018-02-15 were used for the revised version of the manuscript.

