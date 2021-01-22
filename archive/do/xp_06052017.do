/*******************************************************************************
    Title: xp_06072017.do
    Purpose: Diagnostics for counterfactual policies.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/
capture program drop plot_agggraph
program define plot_agggraph
    syntax anything(name=fname), baseval(real) XTItle(string) [WITHSS]

    use fthb_`fname'comp_final.dta, clear
    sort agebin compval
    label define ages 0 "All ages" 20 "Ages 22-29" 30 "Ages 30-39" 40 "Ages 40-49"
    label val agebin ages
    
    unique compval
    local comp_no = `r(sum)'
    sum compval, meanonly
    local min = `r(min)'
    local comp_int = (`r(max)' - `min')/(`comp_no' - 1)
    local baseid = round(`baseval'/`comp_int', 1)
    
    by agebin: gen buychar_pol = Buyer_characteristics__m/Buyer_characteristics__m[`baseid']*100
    by agebin: gen buychar_ss = Buyer_characteristics/Buyer_characteristics[`baseid']*100
    label var buychar_pol "Marginal Buyers"
    label var buychar_ss "Inframarginal Buyers"

    if "`withss'" != "" local yvars buychar_*
    else local yvars buychar_pol    
    twoway line `yvars' compval, by(agebin) yti("Policy Magnitude (100=base)") ///
        xti(`xtitle')
        
end

plot_agggraph Size, baseval(.12) xti("Value of Subsidy (in increments of $1000)")
plot_agggraph F, baseval(.06) xti("Proportional Fixed Cost") withss
plot_agggraph Down, baseval(0.2) xti("Minimum Down Payment Proportion") withss
plot_agggraph House, baseval(.25) xti("Minimum House Size (Model Units)") withss
plot_agggraph I, baseval(2.5) xti("1/EIS") withss

capture program drop plot_elasHeatmap
program define plot_elasHeatmap
    syntax namelist [if], ytitle(string) xtitle(string) ztitle(string) save(string)

    preserve
    local normalizer = 67.250 // From initial calibration in SCF
    tokenize `namelist'
    local zvar `1'
    macro shift
    local varlist `*'
    
    foreach var of local varlist {
        capture gen `var' = string(`var'bin*67.250, "%5.1f")
    }
    local varlist `varlist' compval
    keep if !missing(Proportion_of_buyers)
    replace Proportion_of_buyers = 1.0 if Proportion_of_buyers > 1.0
    replace compval=round(compval/.015, 1)
    
    list in 1/10
    tokenize `varlist'
    heatmap `zvar' `1' `2' `if', polbr(8) ytitle("`ytitle'") ///
        out customf(0.15 0.30 0.50 0.70 0.85) save("`save'") ///
        xtitle("`xtitle'") zti("`ztitle'")

end

use fthb_SizeElas_final, clear
plot_elasHeatmap income_val
plot_elasHeatmap assets
forv i=1/12 {
*plot_elasHeatmap assets income if compval==`i', ytitle() xtitle() save()
}

/*
gen incomes = string(income_valbin*67.250, "%5.1f")
replace compval=compval/.015
forv i=1/5 {
replace prop_buyers_policy`i' = 1.0 if prop_buyers_policy`i' > 1.0 & !missing(prop_buyers_policy`i')
heatmap prop_buyers_policy`i' incomes compval if !missing(prop_buyers_policy`i'), ///
    polbr(8) ytitle("Income (000s)") save(heatmap`i'.pdf) ///
    out customf(0.15 0.30 0.50 0.70 0.85) ///
    xtitle("Value of Subsidy (in increments of $1000)") zti("Prob. becoming FTHB")
}

gen assets = string(assetsbin*67.250, "%5.1f")
replace compval=round(compval/.015, 1)
keep if !missing(Proportion_of_buyers)
replace Proportion_of_buyers = 1.0 if Proportion_of_buyers > 1.0
heatmap Proportion_of_buyers assets compval, polbr(8) ytitle("Liquid Assets (000s)") ///
    out customf(0.15 0.30 0.50 0.70 0.85) zlab(2) save(heatmap_assets.pdf) ///
    xtitle("Value of Subsidy (in increments of $1000)") zti("Prob. becoming FTHB")

forv i=1/12 {
heatmap Proportion_of_buyers assets income if compval==`i', ///
    ytitle("Liquid Assets (000s)") save(output/stata/heatmap_compval`i'.pdf) ///
    out zlab(2) customf(0.15 0.30 0.50 0.70 0.85) ///
    xtitle("Permanent Income State") zti("Prob. becoming FTHB")
}
*/
