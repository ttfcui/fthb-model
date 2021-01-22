
capture program drop elas_assembly
program define elas_assembly
	syntax anything(name=mnames), file1(string) file2(string) [genheat]

	tempfile init init2
	tokenize `mnames'
	local first `1'

	foreach model of local mnames {
	import delimited fthb/propensity_experiment_`model'.csv, clear
	gen fthbs = (inframarg == 1 | marginal == 1)
	replace inframarg = 0 if missing(inframarg)
	xtile income_bin = income_val, n(15)
	xtile assets_bin = assets, n(15)
	gen loginc = log(income_val)
	gen logast = log(assets+sqrt(assets^2+1)) // ihs transformation approx. logs
	if "`model'" == "`first'" {
	    bys assets_bin: sum assets
	    bys income_bin: sum income_val
	}
	gen age_bin = ceil(age/4)*4 - 2
	gen cross_bin = 1e3*income_bin + assets_bin
	gen cross_bin2 = 1e3*age_bin + assets_bin
	* House size decisions over unvarying cohort (inframarginal buyers)
	gen ownVal = inframarg*nextdurables
	replace ownVal = . if ownVal == 0  // These are renters
	unique cross_bin
	preserve

	collapse (mean) fthbs inframarg ownVal (min) income_bin assets_bin ///
	    (mean) loginc logast (count) id, by(cross_bin)
	if "`genheat'" != "" {
	heatmap fthbs income_bin assets_bin, zti("Purchase Propensity") ///
	    xti("Renter Asset Dist. Quantile") yti("Renter Income Dist. Quantile") ///
	    save(fthb/heatmap_`model'_astinc.pdf)
	heatmap inframarg income_bin assets_bin, zti("Purchase Propensity")  ///
	    xti("Renter Asset Dist. Quantile") yti("Renter Income Dist. Quantile") /// 
	    save(fthb/heatmap_`model'_astinc_ss.pdf)
	}
	if "`model'" == "`first'" {
	    save `init', replace
	}
	else {
	    local modelfix = subinstr("`model'", ".", "", .)
	    rename fthbs fthbs_`modelfix'
	    rename ownVal ownVal_`modelfix'
	    merge 1:1 cross_bin using `init'
	    drop _m
	    save `init', replace
	}
	restore

	collapse (mean) fthbs inframarg ownVal (first) age_bin assets_bin ///
	    (mean) age logast (count) id, by(cross_bin2)
	tostring age_bin, replace
	if "`genheat'" != "" {
	heatmap fthbs age_bin assets_bin, zti("Purchase Propensity") ///
	    xti("Renter Asset Dist. Quantile") yti("Renter Age Range") ///
	    save(fthb/heatmap_`model'_astage.pdf)
	heatmap inframarg age_bin assets_bin, zti("Purchase Propensity") ///
	    xti("Renter Asset Dist. Quantile") yti("Renter Age Range") ///
	    save(fthb/heatmap_`model'_astage_ss.pdf)
	}
	if "`model'" == "`first'" {
	    save `init2', replace
	}
	else {
	    local modelfix = subinstr("`model'", ".", "",.)
	    rename fthbs fthbs_`modelfix'
	    rename ownVal ownVal_`modelfix'
	    merge 1:1 cross_bin using `init2'
	    drop _m
	    save `init2', replace
	    sum
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
	order ownVal*, after(fthbs)
	order id*, after(ownVal)
	gen elas_v1 = fthbs/inframarg - 1
	gen elas_v2 = (fthbs_monetary_prElas5/inframarg - 1)/5
	gen elas_v3 = (fthbs_monetary_prElas5_nodown/inframarg - 1)/5

	capture label var income_bin "Renter Income Dist. Quantile"
	capture label var age_bin "Renter Age Range"
	label var assets_bin "Renter Asset Dist. Quantile"
	label var inframarg "Proportion of FTHBs absent policy (inframarginal)"
	label var id "Total renters within bin"
	foreach vartype in fthbs ownVal {
	    if "`vartype'" == "fthbs" local vardesc "Proportion of FTHBs"
	    if "`vartype'" == "ownVal" local vardesc "Mean Normalized House Value"

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
/*
elas_assembly monetary_prElas monetary_little monetary monetary_nodown ///
    FTHBRB_1 monetary_prElas5 monetary_prElas5_nodown, ///
    file1(fthb/pol_astinc) file2(fthb/pol_astage)
  
use fthb/pol_astinc, clear
elas_prep
saveold fthb/pol_astinc, replace

use fthb/pol_astage, clear
elas_prep
saveold fthb/pol_astage, replace
*/
forv i=0.03(0.015)0.19 {
    local nodown_models `nodown_models' FTHBSize_0`i'
    local down_models `down_models' FTHBonDownSize_0`i'
}
display "`nodown_models'"
display "`down_models'"

elas_assembly `nodown_models' monetary_little, ///
    file1(fthb/pol_FTHBSize_astinc) file2(fthb/pol_FTHBSize_astage)
elas_assembly `down_models' monetary_little, ///
    file1(fthb/pol_FTHBOnDownSize_astinc) file2(fthb/pol_FTHBOnDownSize_astage)

