/*******************************************************************************
    Title: IRS_clean.do
    Purpose: Get parameters of interest from IRS FTHB data.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

* %< Work with FTHB age distribution over years
* Most of the code is from first FTHB paper by EZ. Block comments are used
* to indicate new code for model.

* TODO: calculate this in this do file
* This is the "Time FE" that corrects for bias in the age/income-specific
* elasticities computed at the end of the do-file
* Look at Dropbox, elas_denom.csv for calculation
scalar time_fe = -0.1314


* Getting rid of any 2009-specific time FE by differencing btwn
* percentiles (so like a DID?)
capture program drop elas_diff_spec
program define elas_diff_spec
    syntax anything, var(varname) [BOOTSTRAP]
    
    sum treat_pct, meanonly
    local valmax = `r(max)'
    local maxpct = `valmax'-`1'
    keep if treat_pct <= `1' | treat_pct >= `maxpct'
    if "`bootstrap'" != "" local extravar run
    keep `var' tax_yr treat_pct elas_temp `extravar'
    reshape wide elas_temp, i(`var' tax_yr `extravar') j(treat_pct)
    egen elas_treat = rowmean(elas_temp`maxpct'-elas_temp`valmax')
    egen elas_ctrl = rowmean(elas_temp1-elas_temp`1')
    gen elas_diff = elas_treat - elas_ctrl

end
    




tempfile age_counts orig claims temp
use AGE_YR_TPCT.dta, clear
* Increase size of percentiles
replace treat_pct = floor(treat_pct/2)
collapse (sum) count, by(age_primary tax_yr treat_pct)
* To deal with weird refinance boom (?), focus on people under 65.
keep if age_primary < 60
save `age_counts'

use `age_counts', clear
collapse (sum) count, by(tax_yr treat_pct)
rename count count_tot
sort tax_yr treat_pct
save `temp'

use `age_counts', clear
merge m:1 tax_yr treat_pct using `temp', keep(3) nogen
gen share = count/count_tot
gsort tax_yr treat_pct age_primary
gen claim_data = 0
save `orig'

* Regression to measure change in distribution in each group.
tempfile coef point
gen policy = tax_yr == 2009
gen share_hat = .
gen share_se = .
sum treat_pct
forv i=1/`r(max)' {
    qui reg share age_primary##policy if treat_pct==`i'
    predict share_hat2 if e(sample), xb
    replace share_hat = share_hat2 if e(sample)
    gen share_se2 = _se[25.age_primary#1.policy] if e(sample)
    replace share_se = share_se2 if e(sample)
    drop *_*2
}
* Note that the SEs are the same for each dummy in the policy period.
keep if tax_yr == 2008
keep share_hat share_se age_primary treat_pct
save `coef'

* Bootstrap (?)
tempfile only09 bsamp bstemp
use `orig', clear
keep if tax_yr==2009
save `only09'
local j=0
qui{
forv i=1/500 {
    use `orig', clear
    bsample if tax_yr != 2009, strata(treat_pct age_primary)
    bys treat_pct age_primary: gen id = _n
    keep share treat_pct age_primary id
    rename share share_hat`i'
    if `i' > 1 merge 1:1 id treat_pct age_primary using `bstemp'
    if `i' > 1 drop _merge
    save `bstemp', replace
    local ++j
}
}
collapse (mean) share*, by(treat_pct age_primary)
reshape long share_hat, i(age_primary treat_pct) j(run)
save `bsamp'

* Convert changes to levels.
use `orig', clear
keep if tax_yr == 2009
merge m:1 age_primary treat_pct using `coef', keep(3) nogen
gen count_hat = share_hat*count_tot
gen count_diff = count - count_hat
* Special case where the se for the prediction is 
* just a linear rescaling of the se from the earlier regression.
gen count_diff_y1 = count_diff + 1.96*count_tot*share_se
gen count_diff_y2 = count_diff - 1.96*count_tot*share_se

/* Plotting extensive margin elasticities by age - these elasticities
   are identified following arguments in the 20170411 note */
gen elas_temp = log(share/share_hat)
elas_diff_spec 6 , var(age_primary)
* Take out identified FE
gen elas_corr = elas_diff - time_fe
keep elas_corr age_primary
save `point'

* Regenerate bootstrap results
use `orig', clear
keep if tax_yr== 2009
merge 1:n age_primary treat_pct using `bsamp'
gen elas_temp = log(share/share_hat)
elas_diff_spec 6 , var(age_primary) bootstrap

gen elas_corr = elas_diff - time_fe
collapse (p1) elas_corr_p1=elas_corr (p50) elas_corr_bs=elas_corr ///
    (p99) elas_corr_p99=elas_corr, by(age_primary)
merge 1:1 age_primary using `point'
drop _m
save $outdir/stata/age_extelas_v2, replace


* PLOT estimates with bootstrapped error bands
#delimit ;
twoway
    (connect elas_corr age_primary) (line elas_corr_p1 age_primary, lc(maroon) lp(dash))
    (line elas_corr_p99 age_primary, lc(maroon) lp(dash)),
    yline(0, lc(gs9) lp(dash)) legend(off) xti("Homebuyer Age")
    yti("Ext. Margin Elasticity") $graphconfig;
#delimit cr
graph export "$outdir/stata/age_extelas2.pdf", as(pdf) replace


* INCOME ELASTICITIES
tempfile age_counts orig claims temp
use INC_YR_TPCT_JNT.dta,clear
keep if joint==1
keep if inc_bin5 <= 2e5
replace inc_bin5 = floor(inc_bin5/10000)*10000
replace treat_pct = floor(treat_pct/2)
collapse (sum) count, by(tax_yr treat_pct inc_bin5)
save `age_counts'

use `age_counts', clear
collapse (sum) count, by(tax_yr treat_pct)
rename count count_tot
sort tax_yr treat_pct
save `temp'

use `age_counts', clear
merge m:1 tax_yr treat_pct using `temp', keep(3) nogen
gen share = count/count_tot
gsort tax_yr treat_pct inc_bin5
gen claim_data = 0
save `orig'

* Regression to measure change in distribution in each group.
tempfile coef point
gen policy = tax_yr == 2009
gen share_hat = .
gen share_se = .
sum treat_pct
forv i=1/`r(max)' {
    qui reg share inc_bin5##policy if treat_pct==`i' & !inlist(tax_yr, 2004, 2005, 2006)
    predict share_hat2 if e(sample), xb
    replace share_hat = share_hat2 if e(sample)
    gen share_se2 = _se[100000.inc_bin5#1.policy] if e(sample)
    replace share_se = share_se2 if e(sample)
    drop *_*2
}
* Note that the SEs are the same for each dummy in the policy period.
keep if tax_yr == 2008
keep share_hat share_se inc_bin5 treat_pct
save `coef'

* Bootstrap (?)
tempfile only09 bsamp bstemp
use `orig', clear
keep if tax_yr==2009
save `only09'
local j=0
qui{
forv i=1/500 {
    use `orig', clear
    bsample if !inlist(tax_yr, 2009, 2004, 2005, 2006), ///
        strata(treat_pct inc_bin5)
    bys treat_pct inc_bin5: gen id = _n
    keep share treat_pct inc_bin5 id
    rename share share_hat`i'
    if `i' > 1 merge 1:1 id treat_pct inc_bin5 using `bstemp'
    if `i' > 1 drop _merge
    save `bstemp', replace
    local ++j
}
}
collapse (mean) share*, by(treat_pct inc_bin5)
reshape long share_hat, i(inc_bin5 treat_pct) j(run)
save `bsamp'

* Convert changes to levels.
use `orig', clear
keep if tax_yr == 2009
merge m:1 inc_bin5 treat_pct using `coef', keep(3) nogen
gen count_hat = share_hat*count_tot
gen count_diff = count - count_hat
* Special case where the se for the prediction is 
* just a linear rescaling of the se from the earlier regression.
gen count_diff_y1 = count_diff + 1.96*count_tot*share_se
gen count_diff_y2 = count_diff - 1.96*count_tot*share_se

/* Plotting extensive margin elasticities by age - these elasticities
   are identified following arguments in the 20170411 note */
gen elas_temp = log(share/share_hat)
elas_diff_spec 6 , var(inc_bin5)
* Take out identified FE
gen elas_corr = elas_diff - time_fe
keep elas_corr inc_bin5
save `point'

* Regenerate bootstrap results
use `orig', clear
keep if tax_yr== 2009
merge 1:n inc_bin5 treat_pct using `bsamp'
gen elas_temp = log(share/share_hat)
elas_diff_spec 6 , var(inc_bin5) bootstrap

gen elas_corr = elas_diff - time_fe
collapse (p1) elas_corr_p1=elas_corr (p50) elas_corr_bs=elas_corr ///
    (p99) elas_corr_p99=elas_corr, by(inc_bin5)
merge 1:1 inc_bin5 using `point'
drop _m
save $outdir/stata/inc_extelas_v2, replace


#delimit ;
twoway
    (connect elas_corr inc_bin5) (line elas_corr_p1 inc_bin5, lc(maroon) lp(dash))
    (line elas_corr_p99 inc_bin5, lc(maroon) lp(dash)),
    yline(0, lc(gs9) lp(dash)) legend(off) xti("Homebuyer Income (Joint, 10K bins)")
    yti("Ext. Margin Elasticity") $graphconfig;
#delimit cr
graph export "$outdir/stata/inc_extelas2.pdf", as(pdf) replace



