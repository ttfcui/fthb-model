/*******************************************************************************
    Title: IRS_clean.do
    Purpose: Get parameters of interest from IRS FTHB data.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

* Most of the code is from first FTHB paper by EZ. Block comments are used
* to indicate new code for model.
args dboxdir outdir

* Get the median FTHB age moment %<
tempfile age_counts orig claims
import delimited `dboxdir'/AGE_COUNTS.csv, clear
* To deal with weird refinance boom (?), focus on people under 65.
keep if age_primary < 60
save `age_counts'

/* Generate the median FTHB age absent policy (i.e. in steady-state) */
drop if tax_yr==2009
sum age_primary [fw=count], d
scalar med_fthb = `r(p50)' - 21 // 20 is the base year
* %>

* %< Work with FTHB age distribution over years
use `age_counts', clear
collapse (sum) count, by(tax_yr)
rename count count_tot
sort tax_yr
list tax_yr count_tot
tempfile temp
save `temp'

use `age_counts', clear
merge m:1 tax_yr using `temp', keep(3) nogen
gen share = count/count_tot
gsort tax_yr age_primary
gen claim_data = 0
save `orig'


* Regression to measure change in distribution in each group.
tempfile coef
gen policy = tax_yr == 2009
reg share age_primary##policy
predict share_hat, xb
* Note that the SEs are the same for each dummy in the policy period.
gen share_se = _se[25.age_primary#1.policy]
keep if tax_yr == 2008
keep share_hat share_se age_primary
save `coef'

* Convert changes to levels.
use `orig', clear
keep if tax_yr == 2009
merge m:1 age_primary using `coef', keep(3) nogen
gen count_hat = share_hat*count_tot
gen count_diff = count - count_hat
* Special case where the se for the prediction is 
* just a linear rescaling of the se from the earlier regression.
gen count_diff_y1 = count_diff + 1.96*count_tot*share_se
gen count_diff_y2 = count_diff - 1.96*count_tot*share_se

/* Plotting only the counterfactual distribution */
#delimit;
graph bar share_hat if age_primary >= 22, over(age_primary) 
    xsize(11) ytitle("Count in 2009") 
    legend(label(2 "Actual First-Time Homeowners")
           label(1 "Counterfactual First-Time Homeowners"))
    $graphconfig;
#delimit cr

/* Get the age gradients over certain subsections of the data */
tsset age_primary
gen binmean = (share_h+ L.share_h + L2.share_h)/3
local init = 6 // so age 25
scalar agerise = log(binmean[`init'+5]) - log(binmean[`init'])
scalar agefall = log(binmean[`init'+10]) - log(binmean[`init'+5])

/* Export just the counterfactual distribution for comparison purposes */
preserve
keep age_primary share_hat

export delimited `outdir'/FTHB_dist_data.csv if age_primary >= 22, novar replace
restore

* There's more code below, but unnecessary for current analysis.
exit

/* %> */

/* %< Age-ZIP cross-section data */
* This is not actually an age-ZIP dataset: it is indexed by ZIP and year
* but gets aggregate counts for each cell, along with counts for "young buyers".
use `dboxdir'/IRS/age_zip_counts, clear

* Merge on 2005 IRS SOI data
mmerge zip using "`dboxdir'/_irszipdata_2005_noagi", ///
    type(n:n) unm(master)
gen income_bin = .
* Normalization step: Run xtile on the AGI measures for each ZIP to get
* income bins within each year
forv yr=2002/2013 {
    fastxtile incbin_temp=avgagi if tax_yr==`yr' [fw=numret], n(40)
	replace income_bin = incbin_temp if tax_yr==`yr'
	drop incbin_temp
}
* Aggregate up to an agg. FTHB income distribution by year
collapse (sum) total_buyers young_buyers (min) avgagi if ///
    income_bin > 2 & income_bin < 39, by(tax_yr income_bin)

* Proceed as before in getting counterfactual shares
tempfile orig coef
by tax_yr: egen count_tot = total(young_buyers)
gen share = young_buyers/count_tot
save `orig', replace
gen policy = tax_yr == 2009
reg share income_bin##policy if tax_yr==2009 | !inlist(tax_yr, 2006, 2007, 2008)
predict share_hat, xb
keep if tax_yr == 2008
keep share_hat income_bin
save `coef', replace

use `orig', clear
keep if tax_yr == 2009
merge m:1 income_bin using `coef', keep(3) nogen
twoway (line share income_bin) (line share_hat income_bin), $graphconfig

gen elas_temp = log(share/share_hat)
* Take out identified FE
gen elas_corr = elas_temp - time_fe
/* Remove FE
sum elas_temp, meanonly
replace elas_temp = elas_temp - `r(mean)'
scalar elas_base = elas_temp[9]
* The negative sign is because elas_base is negative as computed
replace elas_temp = -(elas_temp - elas_base)/elas_base
*/
* Plot (with absolute AGI values for reference)
local break1 = round(avgagi[6], .01)*1000
local break2 = round(avgagi[13], .01)*1000
local break3 = round(avgagi[27], .01)*1000

replace income_bin = income_bin*2.5
#delimit ;
twoway
   (connect elas_corr income_bin),
    xline(20 37.5 72.5, lp(dash) lc(gs9)) text(
    .035 26 " = \$`break1'" .035 44 " = \$`break2'" .035 79 " = \$`break3'") 
    legend(off) xti("Binned income (Percentiles)") 
    yti("Ext. Margin Elasticity, Scale Normalized") $graphconfig
    note("Dollar captions show 2009 nominal AGI at the marked percentiles.");
#delimit cr
graph export "`outdir'/inc_extelas.pdf", as(pdf) replace

/* ... %> */
