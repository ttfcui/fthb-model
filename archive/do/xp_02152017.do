/*******************************************************************************
    Title: xp_02152017.do
    Purpose: Testing out auxiliary equations on PSID and model data.
    Author: TC, DB
    Date: 2017-02-23

*******************************************************************************/

capture log c
log using $outdir/stata/xp_02152017.smcl, replace
* TODO: add stuff in the regression programs that exports tables.

* %< Relevant programs
capture program drop trim_tails
program define trim_tails
    syntax varlist, level(real)
    foreach x of varlist `varlist' {
        rename `x' `x'_ut
        winsor `x'_ut, generate(`x'_w) p(`level')
        gen `x' = cond(`x'_w == `x'_ut, `x'_ut, .)
    }
end

capture program drop model_binscatter
program define model_binscatter
    syntax varname(name=depvar)

    * Preliminary binscatters to see if a linear model fits (assuming LPM)
    // Model 1
    binscatter `depvar' income_chg, control(age income_val_l) reportreg
    sleep 5000 // Model 2
    binscatter `depvar' income_chg, control(c.age##c.income_val_l) reportreg
    sleep 5000 // Model 3
    binscatter `depvar' income_chg, control(agebin? income_val_l) reportreg
    sleep 5000 // Model 5
    binscatter `depvar' income_chg, control(agebin? c.age##(c.income_val_l)) reportreg
    sleep 5000 // Model 6
    binscatter `depvar' income_chg, control(agebin? c.age##(c.income_val_l*)) reportreg
    sleep 5000 // Model 7
    binscatter `depvar' income_chg, control(agebin? c.agebin?##(c.income_val_l)) reportreg

end

* %< Testing out different binary choice models with (z, t) on RHS
capture program drop model_increg
program define model_increg
    syntax varname [pweight aweight]

    local model1 income_chg age income_val_l
    local model2 income_chg c.age##c.income_val_l
    local model3 income_chg agebin? income_val_l
    local model4 income_chg agebin? income_val_l*
    local model5 income_chg agebin? c.age##(c.income_val_l)
    local model6 income_chg agebin? c.age##(c.income_val_l*)
    local model7 income_chg agebin? c.agebin?##(c.income_val_l)
    local model8 income_chg agebin? c.age##(c.assets_val_l)
    local model9 income_chg agebin? c.age##(c.income_val_l c.assets_val_l)

    forv k=1/9 {
        display "INCOME MODEL `k': `model`k''"
        /* DISABLING PROBIT since running it on model simulation yields problems
        probit `varlist' `model`k'' [`weight' `exp'], vce(cluster id) */
        logit `varlist' `model`k'' [`weight' `exp'], vce(cluster id)
        reg `varlist' `model`k'' [`weight' `exp'], vce(cluster id)
    }
end /* %> */

* %< Testing out different binary choice models with (a, t) on RHS
capture program drop model_assetreg
program define model_assetreg
    syntax varname [pweight aweight]

    local model2 assets_chg c.age##c.assets_val_l
    local model3 assets_chg agebin? assets_val_l
    local model5 assets_chg agebin? c.age##(c.assets_val_l)
    local model7 assets_chg agebin? c.agebin?##(c.assets_val_l)
    local model8 assets_chg agebin? c.age##(c.income_val_l)
    local model9 assets_chg agebin? c.age##(c.income_val_l c.assets_val_l)

    foreach k in 2 3 5 7 8 9 {
        display "ASSETS MODEL `k': `model`k''"
        /* DISABLING PROBIT since running it on model simulation yields problems
        probit `varlist' `model`k'' [`weight' `exp'], vce(cluster id) */
        logit `varlist' `model`k'' [`weight' `exp'], vce(cluster id)
        reg `varlist' `model`k'' [`weight' `exp'], vce(cluster id)
    }
end /* %> */


* %>

* Number of linear splines for age (heterog. effects of shocks by age)
local spline_n 5 
* Cut off data to the age range specified by model
local endage = 61
tempfile data modelsim

* %< PSID data prep
use "$dboxdir\data\psid\c_TC.dta", clear
keep if age < `endage' & age > 20
tempfile orig depvar const_wgt
save `orig'
by id: egen own_c = total(own)
by id: egen all_obs = count(own)
keep if own_c > 0 & own_c != all_obs & !missing(housevalue)
by id: gen owning = (rent[_n-1] > 0 & own == 1)
by id: gen rerenting = (housevalue[_n-1] >0 & own==0)
keep id year *ing own_c all_obs
save `depvar'

use `orig', clear
keep id year weight
drop if missing(weight)
by id: gen const_wgt = weight[_N]
save `const_wgt'

use `orig', clear
merge 1:1 id year using `depvar'
drop _m
merge 1:1 id year using `const_wgt'
replace owning = 0 if missing(owning) & year != 1997
replace rerenting = 0 if missing(rerenting) & year != 1997
by id: egen const_weight = max(const_wgt)
drop const_wgt
order owning rerenting, after(year)
order const_weight, after(weight)

svyset id [pw=const_weight]

svy: mean income
matrix b = e(b)
scalar income_mean = b[1,1]

foreach var in income assets {
    by id: gen `var'_chg = (`var' - L2.`var')/income_mean
    *by id: gen `var'_chg = (`var' - L2.`var')/L2.`var'
    trim_tails `var'_chg, level(.01)
    by id: gen `var'_val_l = L2.`var'/income_mean
    by id: gen `var'_val_l2 = L4.`var'/income_mean
}
mkspline agebin `spline_n' = age
save `data'
* %>
 
* %< Processing data from model
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
rename (income_val_chg assets_l) (income_chg assets_val_l)
mkspline agebin `spline_n' = age
save `modelsim'
* %>

local depvar owning

use `data', clear
*model_binscatter `depvar'
display "`depvar'
model_increg `depvar' [pw=const_weight]
/*
CONCLUSION: 3 and 5 seems to work well.
    4 and 6 could be gotten away if you think income is AR(2), though there are
    serious significance issues + income_chg moves around..
    Interactions in 7 seem silly and coefficients tiny.
    Adding assets (8, 9) does not change things much.
    Logit + Probit > LPM, though neither of them fits that well.
*/
model_assetreg `depvar' [pw=const_weight]
/*
CONCLUSION: Effect of assets_chg is smaller and has opposite sign as
    income_chg. assets_val_l alone is more significant than income_val_l
    alone, but when split into interactions income_val_l has much larger
    magnitudes. Agebin effects don't change much.
*/

use `modelsim', clear
*model_binscatter `depvar'
model_increg `depvar'
/*
CONCLUSION: Similar sign to PSID data, though magnitudes are a lot larger?
    (could be due to income process/lack of unemployment).
    More importantly agebin coefficients are different: PSID effects
    turn negative after first spline, model simulation effects don't.
    Also, switch from ownership to rental vastly different (could jeopardize
    this second regression equation?)
*/
model_assetreg `depvar' [pw=const_weight]
/*
CONCLUSION:

*/
