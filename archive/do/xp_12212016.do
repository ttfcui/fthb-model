/*******************************************************************************
    Title: xp_12212016.do
    Purpose: Testing out auxiliary equations?
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

capture log c
log using $outdir/stata/xp_12212016.smcl, replace

* Number of linear splines for age (heterog. effects of shocks by age)
local spline_n 5 
 
* TC processed the raw output file from the model into CSV form, with other
* fixtures (see xp_12212016.py)
import delimited $dboxdatadir/fthb.csv, clear

xtset id age
by id: gen owning = (rent[_n] == 0 & rent[_n-1] == 1)
by id: gen rerenting = (rent[_n] == 1 & rent[_n-1] == 0)
order owning rerenting, after(age)

foreach var in income_val assets {
    by id: gen `var'_chg = `var' - `var'_l
}
rename income_val_chg income_chg
mkspline agebin `spline_n' = age

* Preliminary binscatters to see if a linear model fits (assuming LPM)
local depvar owning
binscatter `depvar' income_chg, control(age income_val_l) reportreg // Model 1
binscatter `depvar' income_chg, control(c.age##c.income_val_l) reportreg // Model 2
binscatter `depvar' income_chg, control(agebin? income_val_l) reportreg // Model 3
binscatter `depvar' income_chg, control(agebin? c.age##(c.income_val_l)) reportreg // Model 4
binscatter `depvar' income_chg, control(agebin? c.age##(c.income_val_l*)) reportreg // Model 5
binscatter `depvar' income_chg, control(agebin? c.agebin?##(c.income_val_l)) reportreg // Model 6

* %< Testing out different binary choice models
local model1 income_chg age income_val_l
local model2 income_chg c.age##c.income_val_l
local model3 income_chg agebin? income_val_l
local model4 income_chg agebin? income_val_l*
local model5 income_chg agebin? c.age##(c.income_val_l)
local model6 income_chg agebin? c.age##(c.income_val_l*)
local model7 income_chg agebin? c.agebin?##(c.income_val_l)
local model8 income_chg agebin? c.age##(c.assets_l)
local model9 income_chg agebin? c.age##(c.income_val_l c.assets_l)

forv k=1/9 {
    display "MODEL `k': `model`k''"
    *probit `depvar' `model`k'', vce(cluster id) difficult
    logit `depvar' `model`k'', vce(cluster id) difficult
    reg `depvar' `model`k'', vce(cluster id)
}
exit

* %>


* %< The same thing as above, but with assets instead of income
* NOTE: The ML algorithm has some problems fitting the probit models,
* so they are skipped where indicated (namely, ones with interactions w/ splines)
foreach model in probit logit reg {
    display "MODEL 1"
    `model' `depvar' c.assets##c.age 
    if "`model'" != "probit" {
        display "MODEL 2"
        `model' `depvar' c.assets##c.agebin?
    }
    display "MODEL 3"
    `model' `depvar' assets agebin? assets_l 
    if "`model'" != "probit" {
        display "MODEL 4"
        `model' `depvar' assets c.assets#c.age agebin? assets_l 
    }
    if "`model'" != "probit" {
        display "MODEL 5"
        `model' `depvar' c.assets##c.agebin? assets_l 
    }
    display "MODEL 6"
    `model' `depvar' (c.assets c.assets_l )##c.age agebin? 
    if "`model'" != "probit" {
        display "MODEL 7"
        `model' `depvar' (c.assets c.assets_l )##c.agebin? 
    }
}
log c
* %>
