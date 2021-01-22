
capture program drop elas_assembly
program define elas_assembly
	syntax anything(name=mnames), file1(string) file2(string) [genheat]

	tempfile init init2
	tokenize `mnames'
	local first `1'
        display "`2'"
        if "`2'" == "" {  // if only one model is called
            local mnames `mnames' `mnames'
            local infraFlag TRUE
        }

	foreach model of local mnames {
	import delimited propensity_experiment_`model'.csv, clear
	destring leverage_pol?, ignore("infNa") replace
        if "`infraFlag'" != "" {
            replace marginal = . // copying prop file to get SS data
            local infraFlag
        }
	gen fthbs = (inframarg == 1 | marginal == 1)
	rename (owner_pol1 owner_pol9 leverage_pol1 leverage_pol3 consumption) ///
	    (ownP1 ownP9 levP1 levP3 cons)
	replace inframarg = 0 if missing(inframarg)
	xtile income_bin = income_val, n(15)
	xtile assets_bin = assets, n(15)
        destring cons, force replace // risky but this gets rid of pathological agents
	* House value decisions over policy takers (hence fthbs)
	* The coefficient is DP on house divided by rental rate
	gen FTHVal = cond(fthbs==1, (.20)/(0.0548+0.00966)*nextdurables/durables, .)
	* Total consumption expenditures (excluding amortization of mortgage)
	gen expend = cons + exp(0.173)*nextdurables*(fthbs*0.20 + (1-fthbs)*(0.0548+0.00966))
	*
	gen loginc = log(income_val)
	gen logast = log(assets+sqrt(assets^2+1)) // ihs transformation approx. logs
	gen constrained = (assets==0)
	* Debt paid during policy period by homeowners (exp(0.173) is price level of durable) 
	gen nxtast = cond(fthbs==1, nextassets + exp(0.173)*(1-0.20)*nextdurables, .)
	if "`model'" == "`first'" {
	    bys assets_bin: sum assets
	    bys income_bin: sum income_val
	}
	gen age_bin = ceil(age/4)*4 - 2
	gen cross_bin = 1e3*income_bin + assets_bin
	gen cross_bin2 = 1e3*age_bin + assets_bin
	* House size decisions over unvarying cohort (inframarginal buyers)
	gen ownVal = cond(inframarg == 1, nextdurables, .)

	unique cross_bin
	preserve

	collapse (mean) fthbs ownP? inframarg ownVal expend cons nxtast constrained ///
	    (min) income_bin assets_bin (mean) loginc logast (count) id, by(cross_bin)
	if "`genheat'" != "" {
	heatmap fthbs income_bin assets_bin, zti("Purchase Propensity") ///
	    xti("Renter Asset Dist. Quantile") yti("Renter Income Dist. Quantile") ///
	    save(heatmap_`model'_astinc.pdf)
	heatmap inframarg income_bin assets_bin, zti("Purchase Propensity")  ///
	    xti("Renter Asset Dist. Quantile") yti("Renter Income Dist. Quantile") /// 
	    save(heatmap_`model'_astinc_ss.pdf)
	}
	if "`model'" == "`first'" & "`firstFlag'" == "" {
	    save `init', replace
            local firstFlag FALSE
	}
	else {
	    local modelfix = subinstr("`model'", ".", "", .)
            local modelfix = substr("`modelfix'", 1, 25)  // varname restrictions
	    rename (fthbs ownP1 ownP9) (fthbs_`modelfix' ownP1_`modelfix' ownP9_`modelfix')
	    rename (cons expend ownVal) ///
	        (cons_`modelfix' expend_`modelfix' ownVal_`modelfix')
	    merge 1:1 cross_bin using `init'
	    drop _m
	    save `init', replace
	}
	restore

	collapse (mean) fthbs ownP? inframarg ownVal expend cons nxtast constrained ///
	    (first) age_bin assets_bin (mean) age logast (count) id, by(cross_bin2)
	tostring age_bin, replace
	if "`genheat'" != "" {
	heatmap fthbs age_bin assets_bin, zti("Purchase Propensity") ///
	    xti("Renter Asset Dist. Quantile") yti("Renter Age Range") ///
	    save(heatmap_`model'_astage.pdf)
	heatmap inframarg age_bin assets_bin, zti("Purchase Propensity") ///
	    xti("Renter Asset Dist. Quantile") yti("Renter Age Range") ///
	    save(heatmap_`model'_astage_ss.pdf)
	}
	if "`model'" == "`first'" {
	    save `init2', replace
	}
	else {
	    local modelfix = subinstr("`model'", ".", "",.)
            local modelfix = substr("`modelfix'", 1, 25)
	    rename (fthbs ownP1 ownP9) (fthbs_`modelfix' ownP1_`modelfix' ownP9_`modelfix')
	    rename (cons expend ownVal) ///
	        (cons_`modelfix' expend_`modelfix' ownVal_`modelfix')
	    merge 1:1 cross_bin using `init2'
	    drop _m
	    save `init2', replace
	}
	}
	
	use `init', clear
	saveold `file1', replace
	use `init2', clear
	saveold `file2', replace
	
end

capture program drop elas_prep
program define elas_prep
	order fthbs*, after(cross_b*)
	foreach vartype in ownVal ownP9 ownP1 expend cons nxtast {
	    order `vartype'*, after(fthbs)
	}
	order id*, after(ownVal)
	gen elas_v1 = fthbs/inframarg - 1
	gen elas_v2 = (fthbs_monetary_prElas5/inframarg - 1)/5
	gen elas_v3 = (fthbs_monetary_prElas5_nodown/inframarg - 1)/5

	capture label var income_bin "Renter Income Dist. Quantile"
	capture label var age_bin "Renter Age Range"
	label var assets_bin "Renter Asset Dist. Quantile"
	label var inframarg "Proportion of FTHBs absent policy (inframarginal)"
	label var id "Total renters within bin"
	foreach vartype in fthbs ownP1 ownP9 ownVal expend cons nxtast {
	    if "`vartype'" == "fthbs" local vardesc "Proportion of FTHBs"
	    else if "`vartype'" == "ownP1" local vardesc "Proportion owners 1 year post policy"
	    else if "`vartype'" == "ownP9" local vardesc "Proportion owners 9 years post policy"
	    else if "`vartype'" == "ownVal" local vardesc "Mean normalized house value"
	    else if "`vartype'" == "expend" local vardesc "All consumption expenditures"
	    else if "`vartype'" == "cons" local vardesc "Mean non-housing consumption"
	    else if "`vartype'" == "nxtast" local vardesc "Leverage ratio upon homeownership"

	    label var `vartype' "`vardesc', 1% price subsidy experiment"
	    label var `vartype'_monetary_prElas5 "`vardesc', 5% price subsidy experiment"
	    label var `vartype'_monetary_prElas5_nodown "`vardesc', 5% price subsidy not on down experiment"
	    label var `vartype'_monetary "`vardesc', $8000 FTHB subsidy experiment"
	    label var `vartype'_monetary_nodown "`vardesc', $8000 FTHB subsidy not on down experiment"
        }
	label var elas_v1 "Purchase Propensity Elasticity wrt Price, 1% subsidy"
	label var elas_v2 "Purchase Propensity Elasticity wrt Price, 5% subsidy"
	label var elas_v3 "Purchase Propensity Elasticity wrt Price, 5% subsidy nodown"
	
end

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
	    local formula c.logincR##c.logastR c.logincR#c.logincR c.logastR#c.logastR
	    if "`weight'" == "" {
	    eststo: reg `depvar' dummy##(`formula' constrained c.logincR#constrained) [`weight' `exp'] if experiment`j' == 1 | experiment`first' == 1
	    * eststo: logit fthbs dummy##(`formula') [`weight' `exp'] if experiment`j' == 1 | experiment`first' == 1, ///
	    *     r
	    }
	    else {
	    eststo: reg `depvar' dummy##(`formula' c.constrained c.logincR#c.constrained) [`weight' `exp'] if experiment`j' == 1 | experiment`first' == 1
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


*TODO: add models to be regressed on as output
elas_assembly `1', file1(pol_astinc) file2(pol_astage)
    
*use pol_astinc, clear
*elas_prep
*saveold pol_astinc, replace

/* Group-level dataset */
use pol_astinc, clear
rename (fthbs ownVal) (fthbs_elas ownVal_elas)
reshape long fthbs ownVal, i(cross_bin) j(model) string
tab model, gen(experiment)
bys income_bin model: egen income_tot = total(id)
gen asset_prop = id/income_tot
* TODO:
local modelorder 1 2
estimates clear
* TODO:
elas_estimate `modelorder' [fw=id], depvar(fthbs) binvar(income_bin) propvar(asset_prop) ///
    quantiles(`inc_pct')
label var fthbs "Prob. purchase"
drop log*R

collapse (sum) id (first) effect* if model=="_elas", by(income_bin)
save output/elas_astinc_`2'_grouped.dta, replace

exit, STATA
