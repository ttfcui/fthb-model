/*******************************************************************************
    Title: prop_estimate.do
    Purpose: Estimate statistics of interest for comparative policy analysis.
    Author: TC, DB
    Date: 2018-06-20

*******************************************************************************/

capture program drop elas_estimate
program define elas_estimate
     syntax anything(name=mnums) [fweight], depvar(varname) binvar(varname) ///
         propvar(varname) [quantiles(numlist) noReg]

        *
	gen logincR = loginc/log(101/100)
	gen logastR = logast/log(101/100)
	gen dummy = .
	label var dummy "Policy"
	tokenize `mnums'
	local first `1'


	foreach j of numlist `mnums' {
	    replace dummy = experiment`j'
	    if "`reg'" == "" {
	    local formula c.logincR c.logincR#c.logincR c.logastR c.logastR#c.logastR c.logincR#c.logastR
	    if "`weight'" == "" {
	    eststo: reg `depvar' dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment`first' == 1
	    * eststo: logit fthbs dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment`first' == 1, ///
	    *     r
	    }
	    else {
	    eststo: reg `depvar' dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment`first' == 1
	    * eststo: glm fthbs dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment`first' == 1, ///
	    *   link(logit) family(binomial) r
	    }
	    * margins dummy, predict(p) at(loginc=(`quantiles') logast=0)
	    * margins dummy, predict(p) at(loginc=(`quantiles') (mean) logast)
	    predict pred if experiment`j' == 1 | experiment`first' == 1, ys(0, 1)
	    }
	    gen dep_exp = experiment`j'*`depvar'
	    sort `binvar'
	    if "`reg'" == "" {
	        by `binvar': egen effectPred`j' = total(`propvar'*experiment`j'*pred)
		drop pred
	    }
	    by `binvar': egen effect`j' = total(`propvar'*dep_exp)
	    drop dep_exp
	}

	label var logincR "Log income"
	label var logastR "Log liquid wealth"
	if "`reg'" == "" order effectPred*, after(dummy)
	order effect?, after(dummy)
	capture order effect??, after(dummy)

        drop dummy
	
end

/* Individual-level dataset */
capture program drop indiv_assembly
program define indiv_assembly
	syntax anything(name=mnames)

        tempfile init
	tokenize `mnames'
	local first `1'
	
        foreach model of local mnames {
	import delimited fthb/propensity_experiment_`model'.csv, clear
	destring leverage_pol?, ignore("infNa") replace
	gen fthbs = (inframarg == 1 | marginal == 1)
	replace inframarg = 0 if missing(inframarg)
	xtile income_bin = income_val, n(15)
	xtile assets_bin = assets, n(15)

	* Total consumption expenditures (excluding amortization of mortgage)
	gen expend = cons + exp(0.173)*nextdurables*(fthbs*0.20 + (1-fthbs)*(0.0548+0.00966))
	gen expend_share = exp(0.173)*nextdurables*(fthbs*0.20 + (1-fthbs)*(0.0548+0.00966))/cons
	*
	gen loginc = log(income_val)
	gen logast = log(assets+sqrt(assets^2+1)) // ihs transformation approx. logs
	gen constrained = (assets==0)
	* Debt paid during policy period by homeowners (exp(0.173) is price level of durable) 
	gen nxtast = cond(fthbs==1, nextassets + exp(0.173)*(1-0.20)*nextdurables, .)

	
	if "`model'" == "`first'" {
	    bys assets_bin: sum assets
	    bys income_bin: sum income_val
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

end


capture program drop bins_ratios_assembly
program define bins_ratios_assembly
	syntax, default(integer)

	tab model, gen(experiment)
	bys assets_bin model: egen asset_tot = count(id)
	gen income_prop = 1/asset_tot
	bys income_bin model: egen income_tot = count(id)
	gen asset_prop = 1/income_tot

	sort id age model
	* The `default' is to mark levels in the experiment_little (0 subsidies) case
	* Also note that the mean of indiv-level ratios does *not* equal
	* the weighted means of the group-level mean ratios.
	by id age: gen cons_ratio = cond(fthbs==1, consumption/consumption[`default'] - 1, .)
	by id age: gen exp_ratio = cond(fthbs==1, expend/expend[`default'] - 1, .)
	gen durVal = expend - consumption
	by id age: gen durables_ratio = cond(fthbs==1, durVal/durVal[`default'] - 1, .)
	
end

indiv_assembly monetary_prElas monetary_little monetary monetary_nodown ///
    FTHBDF_0 FTHBRB_1 monetary_prElas5 monetary_prElas5_nodown
replace model="experiment_elas" if model=="experiment_monetary_prElas"
bins_ratios_assembly, default(5)
preserve

* TBD Deleveraging variation stuff?
gen leverage0 = cond(fthbs==1, (nxtast/exp(0.173) - (1-0.20)*nextdurables)/nextdurables, .)
gen leverage_diff = leverage_pol1 - leverage0
gen leverage_diff3 = (leverage_pol3 - leverage0)/3
twoway (hist leverage_diff, color(gs12) gap(25) bin(20)) ///
    (hist leverage_diff3, fc(none) lc(navy) bin(20)) if inlist(model, ///
    "experiment_monetary_little", "experiment_monetary_nodown", "experiment_monetary", "experiment_FTHBDF_0"), by(model)
twoway (hist leverage_diff, color(gs12) gap(25) bin(40)) ///
    (hist leverage_diff3, fc(none) lc(navy) bin(40)) if inlist(model, ///
    "experiment_monetary_little", "experiment_monetary_nodown", "experiment_monetary", "experiment_FTHBDF_0") ///
    & inframarg==1, by(model)


* Aggregation of hetero. responses along one state var. axis (income/assets)
foreach axis in income assets {

local propvar asset_prop
if "`axis'" == "assets" local propvar income_prop

tempfile indivStats
estimates clear
elas_estimate 5 3 1 7 4 8 6 2, depvar(fthbs) binvar(`axis'_bin) propvar(`propvar') quantiles(`inc_pct')
label var fthbs "Prob. purchase"
esttab using elas_indiv_fthbs.tex, b(%6.4f) tex replace ///
    mti("Stationary" "1\% Subsidy" "\\$8K Transfer" "5\% Subsidy" "\\$8K Subsidy" ///
        "5\% Subsidy Nodown" "\\$8K Subsidy Nodown" "\\$8K Subsidy Nodown, No Trans. Fees") ///
    ti("Probability of Becoming FTHB") l nostar
drop log*R
collapse (count) id (first) effect* if model=="experiment_elas", by(`axis'_bin)
gen variable = "Prob. purchase"
save `indivStats'

estimates clear
restore, preserve
elas_estimate 5 3 1 7 4 8 6 2, depvar(exp_ratio) binvar(`axis'_bin) propvar(`propvar') quantiles(`inc_pct')
label var exp_ratio "Expenditures semielasticity"
esttab using elas_indiv_expend.tex, b(%6.4f) tex replace ///
    mti("Stationary" "1\% Subsidy" "\\$8K Transfer" "5\% Subsidy" "\\$8K Subsidy" ///
        "5\% Subsidy Nodown" "\\$8K Subsidy Nodown" "\\$8K Subsidy Nodown, No Trans. Fees") ///
    ti("Probability of Becoming FTHB") l nostar
drop log*R
collapse (count) id (first) effect* if model=="experiment_elas", by(`axis'_bin)
gen variable = "Expenditures semielasticity from SS"
append using `indivStats'
save `indivStats', replace

estimates clear
restore, preserve
elas_estimate 5 3 1 7 4 8 6 2, depvar(cons_ratio) binvar(`axis'_bin) propvar(`propvar') quantiles(`inc_pct')
label var cons_ratio "Consumption semielasticity"
drop log*R
collapse (count) id (first) effect* if model=="experiment_elas", by(`axis'_bin)
gen variable = "Consumption semielasticity from SS"
append using `indivStats'
save `indivStats', replace

saveold fthb/elas_`axis'_indiv.dta, replace
twoway (line effect5 effect6 effect4 effect1 `axis'_bin), name(comp_`axis', replace) ///
    legend(label(1 "Steady-State Value") label(2 "\\$8K Subsidy Inapplicable on DP") ///
           label(3 "\\$8K Subsidy Applicable on DP") label(4 "\\$8K Lump-Sum Transfer") ///
    size(small)) subti(,size(medsmall)) by(variable, yr graphregion(c(white)))
graph export "fthb/policy_decomp_`axis'.pdf", as(pdf) replace
twoway (line effect5 effect3 `axis'_bin), name(onepp_`axis', replace)  ///
    legend(label(1 "Steady-State Value") label(2 "1% Subsidy on House Value") ///
    size(small)) subti(,size(medsmall)) by(variable, yr graphregion(c(white)))
graph export "fthb/policy_elas_`axis'.pdf", as(pdf) replace
restore, preserve
}
restore

/* Group-level dataset */
use fthb/pol_astinc, clear
rename (fthbs ownVal) (fthbs_elas ownVal_elas)
reshape long fthbs ownVal, i(cross_bin) j(model) string
tab model, gen(experiment)
bys income_bin model: egen income_tot = total(id)
gen asset_prop = id/income_tot

estimates clear
elas_estimate 5 3 1 7 4 8 6 2 [fw=id], depvar(fthbs) binvar(income_bin) propvar(asset_prop) ///
    quantiles(`inc_pct')
label var fthbs "Prob. purchase"
esttab using elas_grouped.tex, b(%6.4f) tex replace ///
    mti("Stationary" "1\% Subsidy" "\\$8K Transfer" "5\% Subsidy" "\\$8K Subsidy" ///
        "5\% Subsidy Nodown" "\\$8K Subsidy Nodown" "\\$8K Subsidy Nodown, No Trans. Fees") ///
    ti("Probability of Becoming FTHB") l nostar
drop log*R

collapse (sum) id (first) effect* if model=="_elas", by(income_bin)
save fthb/elas_astinc_grouped.dta, replace

/* Individual-level FTHB Size dataset */
forv i=0.015(0.015)0.19 {
    local nodown_models `nodown_models' FTHBSize_0`i'
    local down_models `down_models' FTHBonDownSize_0`i'
}
indiv_assembly `nodown_models' monetary_little
bins_ratios_assembly, default(13)
preserve

foreach axis in income assets {

local propvar asset_prop
if "`axis'" == "assets" local propvar income_prop

tempfile indivStats
estimates clear
elas_estimate 13 1/12, depvar(fthbs) binvar(`axis'_bin) propvar(`propvar') noreg
label var fthbs "Prob. purchase"
drop log*R
collapse (count) id (first) effect* if model=="experiment_monetary_little", by(`axis'_bin)
gen variable = "Prob. purchase"
save `indivStats'

estimates clear
restore, preserve
elas_estimate 13 1/12, depvar(exp_ratio) binvar(`axis'_bin) propvar(`propvar') noreg
label var exp_ratio "Expenditures semielasticity"
drop log*R
collapse (count) id (first) effect* if model=="experiment_monetary_little", by(`axis'_bin)
gen variable = "Expenditures semielasticity from SS"
append using `indivStats'
save `indivStats', replace

estimates clear
restore, preserve
elas_estimate 13 1/12, depvar(cons_ratio) binvar(`axis'_bin) propvar(`propvar') noreg
label var cons_ratio "Consumption semielasticity"
drop log*R
collapse (count) id (first) effect* if model=="experiment_monetary_little", by(`axis'_bin)
gen variable = "Consumption semielasticity from SS"
append using `indivStats'
save `indivStats', replace

saveold fthb/elasNoDown_`axis'_indiv.dta, replace
reshape long effect, i(`axis'_bin id variable) j(noDownlevel)
replace noDownlevel = 0 if noDownlevel == 13
label var noDownlevel "Face Value of Subsidy Inapplicable to DP (000s)"
reshape wide effect, i(id variable noDownlevel) j(`axis'_bin)
if "`axis'" == "income" {
    twoway (line effect8 effect10 effect12 effect14 noDownlevel), by(variable, yr graphregion(c(white))) ///
        legend(label(1 "8th Income Bin (Median)") label(2 "10th Income Bin") ///
               label(3 "12th Income Bin") label(4 "14th Income Bin") size(small)) ///
	subti(,size(medsmall))
}
else if "`axis'" == "assets" {
    twoway (line effect1 effect8 effect11 effect14 noDownlevel), by(variable, yr graphregion(c(white))) ///
        legend(label(1 "0 Assets (1st-7th Bins)") label(2 "8th Asset Bin (Median)") ///
               label(3 "11th Asset Bin") label(4 "14th Asset Bin") size(small)) ///
	subti(,size(medsmall))
}
graph export "fthb/policy_nodownval_`axis'.pdf", as(pdf) replace
restore, preserve
}

restore
indiv_assembly `down_models' monetary_little
bins_ratios_assembly, default(13)
preserve

foreach axis in income assets {

local propvar asset_prop
if "`axis'" == "assets" local propvar income_prop

tempfile indivStats
estimates clear
elas_estimate 13 1/12, depvar(fthbs) binvar(`axis'_bin) propvar(`propvar') noreg
label var fthbs "Prob. purchase"
drop log*R
collapse (count) id (first) effect* if model=="experiment_monetary_little", by(`axis'_bin)
gen variable = "Prob. purchase"
save `indivStats'

estimates clear
restore, preserve
elas_estimate 13 1/12, depvar(exp_ratio) binvar(`axis'_bin) propvar(`propvar') noreg
label var exp_ratio "Expenditures semielasticity"
drop log*R
collapse (count) id (first) effect* if model=="experiment_monetary_little", by(`axis'_bin)
gen variable = "Expenditures semielasticity from SS"
append using `indivStats'
save `indivStats', replace

estimates clear
restore, preserve
elas_estimate 13 1/12, depvar(cons_ratio) binvar(`axis'_bin) propvar(`propvar') noreg
label var cons_ratio "Consumption semielasticity"
drop log*R
collapse (count) id (first) effect* if model=="experiment_monetary_little", by(`axis'_bin)
gen variable = "Consumption semielasticity from SS"
append using `indivStats'
save `indivStats', replace

saveold fthb/elasOnDown_`axis'_indiv.dta, replace
reshape long effect, i(`axis'_bin id variable) j(OnDownlevel)
replace OnDownlevel = 0 if OnDownlevel == 13
label var OnDownlevel "Face Value of Subsidy Applicable to DP (000s)"
reshape wide effect, i(id variable OnDownlevel) j(`axis'_bin)
if "`axis'" == "income" {
    twoway (line effect8 effect10 effect12 effect14 OnDownlevel), by(variable, yr graphregion(c(white))) ///
        legend(label(1 "8th Income Bin (Median)") label(2 "10th Income Bin") ///
               label(3 "12th Income Bin") label(4 "14th Income Bin") size(small)) ///
	subti(,size(medsmall))
}
else if "`axis'" == "assets" {
    twoway (line effect1 effect8 effect11 effect14 OnDownlevel), by(variable, yr graphregion(c(white))) ///
        legend(label(1 "0 Assets (1st-7th Bins)") label(2 "8th Asset Bin (Median)") ///
               label(3 "11th Asset Bin") label(4 "14th Asset Bin") size(small)) ///
	subti(,size(medsmall))
}
graph export "fthb/policy_ondownval_`axis'.pdf", as(pdf) replace
restore, preserve
}


/* Group-level FTHB Size dataset */

capture program drop fthb_size_clean
program define fthb_size_clean

	tempfile inc_bins
	rename (fthbs ownVal) (fthbs_FTHBOnDownSize_0015 ownVal_FTHBOnDownSize_0015)
	reshape long fthbs ownVal, i(cross_bin) j(model) string
	tab model, gen(experiment)
	bys assets_bin model: egen asset_tot = total(id)
	gen income_prop = id/asset_tot
	bys income_bin model: egen income_tot = total(id)
	gen asset_prop = id/income_tot
	preserve

	estimates clear
	elas_estimate 13 1/12 [fw=id], depvar(fthbs) binvar(income_bin) propvar(asset_prop) noreg
	collapse (sum) id (first) effect* if model=="_monetary_little", by(income_bin)
	order effect10 effect11 effect12, after(effect9)
	capture order effectPred10 effectPred11 effectPred12, after(effectPred9)
	save `inc_bins'

	restore
	estimates clear
	elas_estimate 13 1/12 [fw=id], depvar(fthbs) binvar(assets_bin) propvar(income_prop) noreg
	collapse (sum) id (first) effect* if model=="_monetary_little", by(assets_bin)
	order effect10 effect11 effect12, after(effect9)
	capture order effectPred10 effectPred11 effectPred12, after(effectPred9)
	append using `inc_bins'
	order income_bin, after(assets_bin)

end

use fthb/pol_FTHBSize_astinc, clear
fthb_size_clean
saveold fthb/elas_FTHBSize_astinc_grouped.dta, replace

/* Group-level FTHB Size On Down dataset */
use fthb/pol_FTHBOnDownSize_astinc, clear
fthb_size_clean
saveold fthb/elas_FTHBOnDownSize_astinc_grouped.dta, replace


/* WTP analysis combining nodown and on down datasets */
use fthb/elas_FTHBSize_astinc_grouped, clear
gen modelType = "nodown"
append using fthb/elas_FTHBOnDownSize_astinc_grouped
replace modelType = "ondown" if missing(modelType)
drop if missing(assets_bin)
reshape wide effect*, i(assets_bin id) j(modelType, string)


capture program drop impute_effects
program define impute_effects
    syntax anything(name=binname), nsize(integer) nbins(integer)
    local obs = `nsize'*`nbins'
    set obs `obs'
    egen noDownLevel = seq(), from(1) to(`nsize') block(`nbins')
    egen `binname' = seq(), from(1) to(`nbins') block(1)

    forval i=0/6 {
    if `i' == 0 {
        gen effectnodown = .
	gen effectondown = .
    }
    else {
        gen effectless`i'ondown = .
	label var effectless`i'ondown "Net of Effects of \$`i'K on DP"
    }
    forval j=1/`nsize' {
        forval k=1/`nbins' {
	    local index = (`j'-1)*`nbins' + `k'
	    if `i' == 0 {
	        qui replace effectnodown = effect`j'nodown[`k'] in `index'
		qui replace effectondown = effect`j'ondown[`k'] in `index'
	    }
	    else {
	        qui replace effectless`i'ondown = ///
	            effect`j'nodown[`k'] - effect`i'ondown[`k'] in `index'
            }
        }
    }
    }
end

impute_effects assetsBin, nsize(12) nbins(9)
replace assetsBin = assetsBin + 6 if assetsBin > 1
* Graph test
label var noDownLevel "Face Value of Subsidy (000s)"
twoway (line effectless?ondown noDownLevel if assetsBin == 1, yline(0, lc(gs9)))
twoway (line effectless?ondown noDownLevel, yline(0, lc(gs9)) xlab(1/12)) ///
    if mod(assetsBin, 2)==1, by(assetsBin, yr xr graphregion(c(white))) ///
    legend(size(small)) yti("Net Chg in Purchase Propensity, pp") 
graph export "fthb/ondown_wtp.pdf", as(pdf) replace
twoway (line effect??down noDownLevel, xlab(1/12)) ///
    if mod(assetsBin, 2)==1, by(assetsBin, xr graphregion(c(white))) ///
    legend(label(1 "Inapplicable on DP") label(2 "Applicable on DP") size(small))
graph export "fthb/novsondow.pdf", as(pdf) replace
