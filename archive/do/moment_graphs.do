/*******************************************************************************
    Title: moment_graphs.do
    Purpose: Produce visualization of calibration on selected series.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

local outputdir $outdir/stata/calibration
tempfile model modelH data dataH

* Statistics over entire population
import delimited $dboxdir/model/moments_1yearbin.csv, clear
capture drop model
gen model = "model"
save `model', replace

* Statistics for only houseowners
import delimited $dboxdir/model/momentsH_1yearbin.csv, clear
capture drop model
gen model = "model"
save `modelH', replace

* These two datasets are from the SCF calibration.
* Renaming is needed since "wide format" makes plotting the series easier.
import delimited $dboxdatadir/../graph_data_byageind.csv, clear
keep ageind netliquid_exhousingdebt homeownershiprate
rename * (ageind avgassetsexdebt fracown)
gen model = "data"
save `data', replace

import delimited $dboxdatadir/../graph_data_byageind_homeowners.csv, clear
keep ageind netliquid_exhousingdebt housepos
rename * (ageind avgassetsexdebtowners avghouseowners)
merge 1:1 ageind using `data'
save `dataH', replace


* Merge. Convert back to age bins again.
use `model', clear
merge 1:n model agebin using `modelH'
drop _m
append using `dataH'
drop _m
replace agebin = 5*ageind + 17.5 if missing(agebin)
drop ageind
drop if agebin > 70 | agebin < 21
order model, first

* Build the graphs. Calibration is on last two, homeownership rate and housing size.
* Five series are matched, though for presentation likely  only housing size,
* homeownership rate and liquid assets (2) are shown.
twoway (scatter avgassetsexdebt agebin if model == "data") ///
       (connect avgassetsexdebt agebin if model == "model" & agebin >= 22, lp("--") m(X)), ///
        legend(lab(1 "Data") lab(2 "Model")) $graphconfig
graph export "`outputdir'/avgAssetsExDebt.pdf", as(pdf) replace
twoway (scatter avgassetsexdebtowners agebin if model == "data") ///
       (connect avgassetsexdebtowners agebin if model == "model" & agebin >= 22, lp("--") m(X)), ///
        legend(lab(1 "Data") lab(2 "Model")) $graphconfig
graph export "`outputdir'/avgAssetsExDebtOwners.pdf", as(pdf) replace
twoway (scatter fracown agebin if model == "data") ///
       (connect fracown agebin if model == "model" & agebin >= 22, lp("--") m(X)), ///
        legend(lab(1 "Data") lab(2 "Model")) $graphconfig
graph export "`outputdir'/fracOwn.pdf", as(pdf) replace
twoway (scatter avghouseowners agebin if model == "data") ///
       (connect avghouseowners agebin if model == "model" & agebin >= 22, lp("--") m(X)), ///
        legend(lab(1 "Data") lab(2 "Model")) yscale(r(0 4)) ylab(1(01)4) $graphconfig
graph export "`outputdir'/avgHouseOwners.pdf", as(pdf) replace
