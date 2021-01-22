/*******************************************************************************
    Title: xp_03152017.do
    Purpose: Outputting aggregate statistics from "fthb_stats" file. 
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/


local IRSdiv = 0.00654 // K-L divergence on 2009 dist. vs. 10-year counterfactual (share_hat)

/* %< FTHB_stats file processing */
use $outdir/fthb_stats_final.dta, clear
label def agebin 0 "All agents" 20 "Ages 22-29" 30 "Ages 30-39" 40 "Ages 40-49" ///
    50 "Ages 50-59" 60 "Ages 60-69" 70 "Ages 70-79"
label val agebin agebin

* Creating select statistics to compare different policies
tempfile trunc
scalar fthb_trans = 11725 // Captured from running experiment.py - UPDATE IF NECESSARY
scalar cars_trans = 115111 // Captured from running experiment.py too
gen comp_stats = cond( ///
    regexm(desc, "1 period fwd, total") == 1 & subtype == "count", 1, cond( ///
    regexm(desc, "2 periods fwd, total") == 1 & subtype == "count", 2, cond( ///
    regexm(desc, "2\+ periods fwd, total") == 1 & subtype == "count", 3, cond( ///
    regexm(desc, "More .+ over lifecycle$") == 1, 4, cond( ///
    regexm(desc, "Fewer .+ over lifecycle$") == 1, 5, cond( ///
    desc == "Buyer characteristics" & subtype == "id", 99, .))))))
preserve
collapse (first) desc subtype (sum) valuemonetary valueCARS valueFTHBD-valueFTHBS ///
    if agebin > 0 & agebin < 60 & !inlist(comp_stats, 4, 5), by(comp_stats)
save `trunc'
restore, preserve
collapse (first) desc subtype (sum) valuemonetary valueCARS valueFTHBD-valueFTHBS ///
    if agebin == 0 & inlist(comp_stats, 4, 5), by(comp_stats)
append using `trunc'
drop if missing(comp_stats)

foreach var of varlist valuem-valueFTHBS {
    gen Fac`var' = `var'/`var'[_N] if comp_stats < 4
    if regexm("`var'", "CARS") == 0 replace Fac`var' = ///
        (`var'[1]-`var'[2])/fthb_trans if comp_stats >= 4
    else replace Fac`var' = (`var'[1]-`var'[2])/cars_trans if comp_stats >= 4
}
drop if comp_stats > 4
label var Facvaluem "Main policy experiment"
label var FacvalueC "CARS policy experiment"
label var FacvalueFTHBD "Lump-sum tax credit"
label var FacvalueFTHBT "Two-period tax credit"
label var FacvalueFTHBS "Small value tax credit"

eststo tabA: estpost tabstat Facvaluem-FacvalueCARS, by(comp_stats) ///
    statistics(mean) columns(statistics) not
eststo tabB: estpost tabstat Facvaluem FacvalueFTHBD-FacvalueFTHBS, by(comp_stats) ///
    statistics(mean) columns(statistics) not

local labels eql("Total timing margin, 1 period" "Total timing margin, 2 periods" ///
    "Total timing margin, 2+ periods" "Net extensive margin")
esttab tabA using $outdir/stata/comp_table1.tex, unstack main(mean) l `labels' ///
    nostar not noobs tex replace
esttab tabB using $outdir/stata/comp_table2.tex, unstack main(mean) l `labels' ///
    nostar not noobs tex replace


* Creating a table of positive extensive margin for slides
gen posext = subtype if regexm(desc, "Extensive from") == 1
destring posext, replace
replace posext = 0 if missing(posext)
gen ext_dummy = cond(regexm(desc, "from 0") == 1 & posext == 1, 1, cond( ///
                     regexm(desc, "from 0") == 1 & posext > 1, 2, cond( ///
                     regexm(desc, "from 1") == 1 & posext == 2, 3, cond( ///
                     regexm(desc, "from 1") == 1 & posext > 2, 4, cond( ///
                     regexm(desc, "from 2") == 1 & posext == 3, 5, cond( ///
                     regexm(desc, "from 2") == 1 & posext < 2, 6, .))))))
replace ext_dummy = 99 if desc == "Buyer characteristics" & subtype == "id"

preserve
collapse (first) desc-subtype (sum) valuemonetary, by(agebin ext_dummy)
drop if missing(ext_dummy)
by agebin: gen value_prop = value/value[_N]
drop if ext_dummy == 99
estpost tab agebin ext_dummy [iw=value_prop] if agebin <= 50, notot
esttab . using $outdir/stata/extensive_table.tex, noobs not nostar unstack ///
    eql("0-1 Extensive" "0-1+ Extensive" "1-2 Extensive" "1-2+ Extensive" ///
        "2-3 Extensive" "2-1 Extensive") tex replace

* Creating graph of CARS durable adjusters by capital vintage
restore, preserve
keep if regexm(desc, "Proportion of buyers by capital vintage") == 1
collapse (max) valueCARS*, by(agebin desc subtype)
destring subtype, gen(vintage)
sort agebin desc subtype
encode desc, gen(vintType)

twoway (scatter valueCARS vintage if vintType==2 & vintage > 0, m(Oh)) ///
       (scatter valueCARS vintage if vintType==1 & vintage > 0) if agebin <= 0, ///
       xline(5, lc(gs9)) by(agebin, $graphconfig) yti("Proportion of buyers") xti("Capital vintage") ///
       legend(label(1 "Stationary equilibrium") label(2 "Policy simulation"))
*graph export $outdir/stata/CARSvintage.pdf, replace
/* %> */

/* %< FTHB_EIS file processing (select statistics over many models) */
use $outdir/fthb_EIScomp_final.dta, clear
sort agebin desc subtype
* Reorganizing all the different measures.
foreach var of varlist value* {
    * Difference of means btwn policy/SS distribution
    by agebin: gen `var'_meandiff = `var'[5] - `var'[2]
    * Nonparametric skew of just the policy distribution
    by agebin: gen `var'_skew = (`var'[5] - `var'[4])/`var'[6]
    * Difference of nonparametric skews
    by agebin: gen `var'_skewdiff = ///
        (`var'[5] - `var'[4])/`var'[6] - (`var'[2] - `var'[1])/`var'[3]
    * Hellinger distance (not in same row as KL divergence so have to recast it)
    by agebin: gen `var'_hell = `var'[7]
    order `var'_h `var'_m `var'_s*, after(`var')
}
by agebin: keep if _n == _N
* By transposing, each column becomes a set of statistics for a given age
* group. Then data is organized on EIS-measure level.
xpose, clear

drop in 1/3
rename (v*) (ageAll age22_29 age30_39 age40_49 age50_59)
gen EIS = 1.5+((ceil(_n/5)-1)/4)
replace EIS = ceil(_n/5) - 10 if EIS > 5
bys EIS: gen metric = _n

label def metric 1 "KL divergence" 2 "Hellinger distance" 3 ///
    "Nonpara. skew of policy dist." 4 "Difference of means" ///
    5 "Difference of skews"
label val metric metric

* Graphs are outputted below
twoway (connect ageAll EIS, lp("-.")), by(metric, yr)
graph export $outdir/stata/EISplots.pdf, as(pdf) replace
twoway (connect age22_29 EIS, lp("-.")), by(metric, yr)
graph export $outdir/stata/EISplots_20s.pdf, as(pdf) replace
twoway (connect age30_39 EIS, lp("-.")), by(metric, yr)
graph export $outdir/stata/EISplots_30s.pdf, as(pdf) replace

twoway (connect ageAll EIS if metric==1, lp("-.") m(S) msize(small)), ///
    graphregion(c(white)) yti("KL Divergence Btwn Distributions") ///
    yline(`IRSdiv', lp(dash))
graph export $outdir/stata/EISsurface.pdf, as(pdf) replace
/* %> */

