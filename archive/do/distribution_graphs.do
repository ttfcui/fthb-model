/*******************************************************************************
    Title: distribution_graphs.do
    Purpose: Produce visualization of housing, asset wealth distributions.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

local outputdir $outdir/stata/calibration

/* %< LORENZ CURVES FOR WEALTH, DATA VS. MODEL */
tempfile scf_filter
use "$dboxdatadir/../lorenz_test.dta", clear
keep if mod(_n, 20) == 0
twoway (line houses_hCum houses_hPct if !missing(houses_h)) ///
    (line labIncWork_hCum labIncWork_hPct if !missing(labIncWork_h)), legend ( ///
	 label(1 "Housing, homeowners") label(2 "Income, homeowners")) $graphconfig	
graph export "`outputdir'/wealthDist_data.pdf", as(pdf) replace
keep *Cum *Pct
rename (labIncWorkCum labIncWorkPct netWorthncCum netWorthncPct ///
        netWorthnc_rCum netWorthnc_rPct housesCum) ///
    (income incomePct networth networthPct networth_r networth_rPct houses)
gen type = "data"
save `scf_filter'

import delimited $dboxdir/model/lorenzTest.csv, clear
foreach var of varlist income* houses* netw* {
        sum pct if !missing(`var')
        gen `var'Pct = pct/`r(max)'*100
}

keep if mod(_n, 50) == 0 // Clears up graph, fewer points

twoway (line houses_h houses_hPct if !missing(houses_h)) ///
    (line income_h income_hPct if !missing(income_h)), legend( ///
	 label(1 "Housing, homeowners") label(2 "Income, homeowners")) $graphconfig
graph export "`outputdir'/wealthDist_model.pdf", as(pdf) replace

append using `scf_filter'
twoway (line houses housesPct if type != "data", lc(maroon)) ///
       (line houses housesPct if type == "data", lp("--.")), ///
        xti("") yti("Housing") legend(label(1 "Model") label(2 "Data")) $graphconfig
graph export "`outputdir'/houses_dist_model.pdf", as(pdf) replace

twoway (line income incomePct if type != "data", lc(maroon)) ///
       (line income incomePct if type == "data" & incomePct <=100, lp("--.")), ///
        xti("") yti("Income") legend(label(1 "Model") label(2 "Data")) $graphconfig
graph export "`outputdir'/income_dist_model.pdf", as(pdf) replace

twoway (line networth networthPct if type != "data", lc(maroon)) ///
       (line networth networthPct if type == "data", lp("--.")), ///
        xti("") yti("Net Worth") legend(label(1 "Model") label(2 "Data")) $graphconfig
graph export "`outputdir'/networth_dist_model.pdf", as(pdf) replace
/* %> */

/* %< INDIVIDUAL ASSETS DATA PROCESSING */
import delimited $dboxdir/model/indivAssets.csv, clear
gen leverage = -nextassets/nextdurables*(1-rent)
replace leverage = . if leverage <= 0
sort id age

* Normalize to mean income simulated, not base value
foreach var of varlist nextassets networth {
    sum income_val, meanonly
    replace `var' = `var'/`r(mean)'
}
* Produces tables similar to BGLV Table 2.
sum leverage if age <= 38, d
bys rent: sum nextassets networth if age <= 38, d

* Lorenz curve for renters' financial assets (data TBD)
preserve
keep if rent==1
sum nextassets
sort nextassets
gen nextassetscum = sum(nextassets)/`r(sum)'*100
egen Pct = rank(nextassets), u
gen nextassetsPct = Pct/_N*100
sort nextassetsPct
twoway line nextassetscum nextassetsPct, yti("Assets for Renters") xti("") $graphconfig
/* %> */

/* %< DISTRIBUTION OF ASSETS FOR RENTERS: IS THERE BUNCHING AT 0?
import delimited $outdir/rentAssets.csv, clear
label var nextassets "Assets among renters"
hist nextassets if nextassets < 3, width(0.025) $graphconfig
graph export "`outputdir'/rentAssets_model.pdf", as(pdf) replace
hist nextassets if nextassets < 3 & age < 40, width(0.025) $graphconfig
graph export "`outputdir'/rentAssets_workers_model.pdf", as(pdf) replace
%> */

/* %< DISTRIBUTION OF H/W RATIOS, MODEL VS. DATA */
/*
import delimited $outdir/ownerWealth.csv, clear
preserve
tempfile all
xtile hw_ratio_pct=housetonetw, n(100) altdef
collapse (max) hw_ratio_model=housetonetw, by(hw_ratio_pct)
save `all', replace

restore
xtile hw_ratio_pct=housetonetw if age < 40, n(100) altdef
collapse (max) hw_ratio_model_work=housetonetw, by(hw_ratio_pct)
merge 1:1 hw_ratio_pct using `all'
drop _m
merge 1:1 hw_ratio_pct using $outdir/hw_ratios
drop _m
qqplot hw_ratio_model hw_ratio, $graphconfig m(Th) msize(small) ysc(log) ///
    xti("Quantiles of H/W ratio, data") yti("Quantiles of H/W ratio, model") ///
    ysc(log) xsc(log)
graph export "`outputdir'/ownerWealth_model.pdf", as(pdf) replace
qqplot hw_ratio_model_work hw_ratio_work, $graphconfig m(Th) msize(small) ///
    xti("Quantiles of H/W ratio, data") yti("Quantiles of H/W ratio, model") ///
    ysc(log) xsc(log)
graph export "`outputdir'/ownerWealth_workers_model.pdf", as(pdf) replace

/*
gen netw_inv = 1/networth
label var netw_inv "Reciprocal of net worth, owners"
binscatter housetonetw netw_inv, n(30) m(Th) xscale(rev) reportreg
graph export "$outdir/stata/hw_relation.pdf", as(pdf) replace
*/ */
/* %> */

