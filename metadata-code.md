# Instructions for working with the code and data in the repository

This file outlines the order in which scripts should be executed (from scratch), and describes dependencies and outputs of each script.

All scripts assume that the working directory is `../scripts/`

`(*)` = does not depend on parent script (but may depend on a function existing)

To produce paper objects, run this:
===================================

**Caveat: no safequards are in place for there not being a package installed...**

```r
source("gasexchange_pressuredata.R")
source("weathermanipulations.R")
source("pressuremanipulations.R")
source("precipmanipulations-routines.R")
source("co2_scenarios.R")
source("co2ExplVar.R")
source("regression_routines.R")
source("climate-weather-modeling.R")
source("regression_routines_eachlake.R")
source("eachlake-gamplots-lagged.R")
source("sensitivity.R")
```

1-4: prepare data for calculating co2 flux
===========================================
1. `gasexchange_pressuredata.R`: grabs raw hourly data from online and saves dmet.
2. `weathermanipulations.R` `(*)`: grabs dmet and calculates mean wind, pressure, temperature and humidity; spits out all data frames as separate rds's (`[*]data.rds`)
3. `pressuremanipulations.R`: grabs all rds's created by (2) and combines them with database query tables; changes wind and DIC into correct units. Saves this table as `params.rds`
	+++ saves archaic one using `gasExchange.R` as `gasFlux.rds`
4. `precipmanipulations-routines` `(*)`: archaically sourced functions/gasexchange_precipdata.R but now runs with updated bottomoftheheap scripts for DAILY precipitation data from online, saves it as precip.rds

5: calculate CO2 flux
============================================

5. this does two things
    i. `co2_scenarios.R`: reads in params. requires maunaloa, salcalc, gasExchangeFlex. Inserts also kerri's values for winds. Creates SalCalc in params and merges maunaloa pressure with params
	    +++ replaces outlier cond, pH, salinity values with NA. Saves updated params as `params-flux.rds`
    ii. `co2data_comparisons.R`: reads in archaic `gasFlux.rds` and tests similarity with Kerri Finlay's results. Requires archaic gasExchange functions. No output saved.

6: Prepare data for regressions for CO2 flux
=============================================

a. `co2ExplVar.R`:

   1. takes pdo from online and creates an annual mean column. saves as `pdo.rds`
   2. takes soi from online and saves as `soi_stand.rds` (nonstand for unstandardised data)
   3. takes nao from online and saves as `naoseasonal.rds`
   4. takes `temperature.rds` and computes monthly and annual means
   5. takes `relhum.rds` (already monthly reso)
   6. takes various supporting data from database queries. creates means for chl a and bottle production estimates
      * makes a production outlier NA
   7. merges all above excepting climate indices, and saves as co2explained.rds
      * this has all the POTENTIAL predictors too	

b. `regression_routines.R`:

   * combines co2 flux with predictors
   * subsets available predictors from co2explained
   * reads in `params-flux.rds`; `co2explained.rds`
   * deals with remaining outliers
   * saves data frame with selected predictors and `NA`s removed as `regvars.rds`
 
c. `climate-weather-modeling.R`:

   * takes most of previous data and incorporates evaporation and SPEI index into appropriate measures, saves as `weathers.rds`

7: regressions for CO2 flux
==============================================

a. `regression_routines_models.R`: ***archaic***: regresses CO2 flux against variables of interest, first stab at different ways of incorporating Lake as random effect vs factor etc. reads in regvars. no output produced
b. `regression_routines_ph.R`: ***archaic***: regresses pH against variables of interest, Year and Lake as	random effect. reads in regvars. no output.
c. `regression_routines_eachlake.R`: ***models developed and then selected for paper***: including rationale. subsets to WW for Matt's paper. reads in regvars, weathers.rds. produces models and model summaries as output. 

8: Something
=========================

a. `eachlake-gamplots.R`: ***archaic***: plots old models and saves output; should not be run unless
	old models wanted
b. `eachlake-gamplots-lagged.R`: plots paper figures from regression models, and summary plot

9: Sensitivity analyses for CO2 flux:
==============================================

a. `sensitivity.R`: requires regvars and fluxes, and function gasExchangeSensitivity; computes sensitivity analysis and saves plots and objects
b. `sensitivity-larea.R`: runs sensitivity analysis with lake area as a variable as well. Not used for paper

