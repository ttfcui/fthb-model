
capture program drop elas_estimate
program define elas_estimate
     syntax anything(name=mnums) [fweight], binvar(varname) propvar(varname) [quantiles(numlist)]

	replace loginc = loginc/log(101/100)
	replace logast = logast/log(101/100)
	gen dummy = .
	label var dummy "Policy"


	foreach j of numlist `mnums' {
	    replace dummy = experiment`j'
	    local formula c.loginc c.loginc#c.loginc c.logast c.logast#c.logast c.loginc#c.logast
	    if "`weight'" == "" {
	    eststo: reg fthbs dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment4 == 1
	    * eststo: logit fthbs dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment4 == 1, ///
	    *     r
	    }
	    else {
	    eststo: reg fthbs dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment4 == 1
	    * eststo: glm fthbs dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment4 == 1, ///
	    *   link(logit) family(binomial) r
	    }
	    * margins dummy, predict(p) at(loginc=(`quantiles') logast=0)
	    * margins dummy, predict(p) at(loginc=(`quantiles') (mean) logast)
	    predict pred if experiment`j' == 1 | experiment4 == 1, ys(0, 1)
	    gen fthbs_exp = experiment`j'*fthbs
	    sum fthbs_exp pred if `propvar' > 0
	    bys `binvar': egen effectPred`j' = total(`propvar'*experiment`j'*pred)
	    by `binvar': egen effect`j' = total(`propvar'*fthbs_exp)
	    drop pred fthbs_exp
	}

	label var loginc "Log income"
	label var logast "Log liquid wealth"
	label var fthbs "Prob. purchase"
	order effectPred*, after(dummy)
	order effect?, after(dummy)
	capture order effect??, after(dummy)

end

/* Individual-level dataset */
tempfile init init2

foreach model in monetary_prElas monetary_little monetary monetary_nodown FTHBRB_1 monetary_prElas5 monetary_prElas5_nodown {
import delimited fthb/propensity_experiment_`model'.csv, clear
gen fthbs = (inframarg == 1 | marginal == 1)
replace inframarg = 0 if missing(inframarg)
xtile income_bin = income_val, n(15)
xtile assets_bin = assets, n(15)

gen loginc = log(income_val)
gen logast = log(assets+sqrt(assets^2+1)) // ihs transformation approx. logs

if "`model'" == "monetary_prElas" {
    pctile inc_pct = loginc, n(10)
    levelsof inc_pct, l(inc_pct)
    drop inc_pct
    save `init', replace
}
else {
    append using `init'
    save `init', replace
}

}
replace model="experiment_elas" if model=="experiment_monetary_prElas"
tab model, gen(experiment)
bys income_bin model: egen income_tot = count(id)
gen asset_prop = 1/income_tot

estimates clear
elas_estimate 4 2 6 3 7 5 1, binvar(income_bin) propvar(asset_prop) quantiles(`inc_pct')
esttab using elas_indiv.tex, b(%6.4f) tex replace ///
    mti("Stationary" "1\% Subsidy" "5\% Subsidy" "\\$8K Subsidy" "5\% Subsidy Nodown" ///
        "\\$8K Subsidy Nodown" "\\$8K Subsidy Nodown, No Trans. Fees") ///
    ti("Probability of Becoming FTHB") l nostar
    
collapse (count) id (first) effect* if model=="experiment_elas", by(income_bin)
save fthb/elas_astinc_indiv.dta, replace


/* Group-level dataset */
use fthb/pol_astinc, clear
rename (fthbs ownVal) (fthbs_elas ownVal_elas)
reshape long fthbs ownVal, i(cross_bin) j(model) string
tab model, gen(experiment)
bys income_bin model: egen income_tot = total(id)
gen asset_prop = id/income_tot

estimates clear
elas_estimate 4 2 6 3 7 5 1 [fw=id], quantiles(`inc_pct')
esttab using elas_grouped.tex, b(%6.4f) tex replace ///
    mti("Stationary" "1\% Subsidy" "5\% Subsidy" "\\$8K Subsidy" "5\% Subsidy Nodown" ///
        "\\$8K Subsidy Nodown" "\\$8K Subsidy Nodown, No Trans. Fees") ///
    ti("Probability of Becoming FTHB") l nostar

collapse (sum) id (first) effect* if model=="_elas", by(income_bin)
save fthb/elas_astinc_grouped.dta, replace

/* Group-level FTHB Size dataset */
use fthb/pol_FTHBSize_astinc, clear
rename (fthbs ownVal) (fthbs_FTHBSize_003 ownVal_FTHBSize_003)
reshape long fthbs ownVal, i(cross_bin) j(model) string
tab model, gen(experiment)
bys assets_bin model: egen asset_tot = total(id)
gen income_prop = id/asset_tot

*bys income_bin model: egen income_tot = total(id)
*gen asset_prop = id/income_tot

estimates clear
elas_estimate 12 1/11 [fw=id], binvar(assets_bin) propvar(income_prop)
collapse (sum) id (first) effect* if model=="_monetary_little", by(assets_bin)
order effect10 effect11, after(effect9)
order effectPred10 effectPred11, after(effectPred9)
save fthb/elas_FTHBSize_astinc_grouped.dta, replace

/* Group-level FTHB Size On Down dataset */
use fthb/pol_FTHBOnDownSize_astinc, clear
rename (fthbs ownVal) (fthbs_FTHBOnDownSize_003 ownVal_FTHBOnDownSize_003)
reshape long fthbs ownVal, i(cross_bin) j(model) string
tab model, gen(experiment)
bys assets_bin model: egen asset_tot = total(id)
gen income_prop = id/asset_tot
*bys income_bin model: egen income_tot = total(id)
*gen asset_prop = id/income_tot

estimates clear
elas_estimate 12 1/11 [fw=id], binvar(assets_bin) propvar(income_prop)
collapse (sum) id (first) effect* if model=="_monetary_little", by(assets_bin)
order effect10 effect11, after(effect9)
order effectPred10 effectPred11, after(effectPred9)
save fthb/elas_FTHBOnDownSize_astinc_grouped.dta, replace



/* WTP analysis combining nodown and on down datasets */
use fthb/elas_FTHBSize_astinc_grouped, clear
gen modelType = "nodown"
append using fthb/elas_FTHBOnDownSize_astinc_grouped
replace modelType = "ondown" if missing(modelType)
reshape wide effect*, i(assets_bin id) j(modelType, string)

set obs 165
egen noDownLevel = seq(), from(1) to(11) block(9)
egen assetsBin = seq(), from(1) to(9) block(1)
*egen noDownLevel = seq(), from(1) to(11) block(15)
*egen incomeBin = seq(), from(1) to(15) block(1)

forval i=2/6 {
    gen effectless`i'ondown = .
    forval j=1/11 {
        forval k=1/9 {
	    local index = (`j'-1)*9 + `k'
	    qui replace effectless`i'ondown = ///
	        effect`j'nodown[`k'] - effect`i'ondown[`k'] in `index'
        }
    }
}

* Graph test
twoway (line effectless?ondown noDownLevel if assetsBin == 1, yline(0, lc(gs9)))
twoway (line effectless?ondown noDownLevel, yline(0, lc(gs9)) xlab(1/11)) ///
    if mod(assetsBin, 2)==1, by(assetsBin, yr xr graphregion(c(white)))
