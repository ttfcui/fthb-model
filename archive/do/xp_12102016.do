/*******************************************************************************
    Title: xp_12102016.do
    Purpose: Looking closer at weird trends in FTHB model.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

local outputdir $outdir/stata/xp_12102016
/* %< Look at asset dynamics of renters before they buy during working age.
    The idea here is that we look at most renters who buy before retirement
    (see xp_12102016.py for the procedure) and then take the means of
    savings (a) for each age cohort, and the initial income state from which
    their simulation began. Of course the agent could have gotten higher
    incomes by the time they buy. */
import delimited $outdir/asset_panel.csv, clear
keep if mod(incomes_init+1, 3) == 0
reshape wide assets_*, i(age) j(incomes_init)

foreach pct in 75 50 {
    forv s = 2(3)11 {
        label var assets_`pct'`s' "`pct'th pctile, init. state `s'"
    }
}

* Two lines on the graph: 0.12 is about the size of the credit, normalized
* to median household income. 0.6 is about the average downpayment on a house
* (20% of 3 units of housing expenditure)
twoway (connect assets_50* age, lp(l "-" "_-" "_..") ///
        color(blue*1.4 blue*0.3 navy purple)), yline(0.12 0.6, lp("-") lc(gs9)) ///
        $graphconfig xti("Age") yti("Assets a, units of household income")
graph export `outputdir'/asset_race_median.pdf, as(pdf) replace
twoway (connect assets_75* age, lp(l "-" "_-" "_..") ///
        color(blue*1.4 blue*0.3 navy purple)), yline(0.12 0.6, lp("-") lc(gs9)) ///
        $graphconfig xti("Age") yti("Assets a, units of household income")
graph export `outputdir'/asset_race_75p.pdf, as(pdf) replace

/* %> */

/* %< Looking at actual simulated income in the model - do they match
  the BPP data which the simulation should replicate? */
tempfile model
import delimited $outdir/model_incomes.csv, clear
replace age = age + 20
save `model', replace

make_BPP // This program is from BPP_calibration.do, run it first!
merge 1:1 age using `model'
gen age1 = age+0.2
gen age2 = age-0.2
twoway (scatter uy age1) (scatter logy_sim age2) ///
    (rcap uy_upper uy_lower age1, lc(navy) lp("-........")) ///
    (rcap logy_sim75 logy_sim25 age2, lc(maroon) lp("-........")), $graphconfig
*twoway scatter uy logy_model_wo logy_sim age, title("After taxes")
graph export `outputdir'/model_incomes.pdf, as(pdf) replace
/* Conclusion: after decreasing s.d. below unconditional value and adding
   a slightly negative mean, mean income matches better with data.
   The income distribution in model, assuming Gaussian shocks, doesn't
   seem to capture the top end well. */

/* %> */

/* %< Housing wealth-consumption ratios and asset levels for FTHBs 
      in steady-state, during policy and post-policy.
      Timing: period 1 is steady-state. Policy only activates in period 2.
      The transition continues up to period 38 just to see how fast
      the return to steady-state is. */
import delimited $outdir/transition_collapse_12112016.csv, clear
replace period = period + 2 // Notational change
keep if agebought < 40 // Working agents only

* This plots heatmaps of experiment 1) in the model slides, where the x axis
* is time period, y axis is age cohort and variable plotted is
* housing wealth-consumption ratios and financial wealth, respectively.
* (because we only track FTHBs, financial wealth always nonnegative)
heatmap h_c_ratio agebought period if model != "experiment_monetary_little", ///
    polbreak(2, 4) save(`outputdir'/heatmap_h.pdf) out zti("Housing -Consumption Ratio")
heatmap assets agebought period if model != "experiment_monetary_little", ///
    count polbreak(2, 4) save(`outputdir'/heatmap_a.pdf) out zti(Financial Assets)
heatmap consumption_chg agebought period if model != "experiment_monetary_little", ///
    polbreak(2, 4) customf(0.10 0.20 0.25 0.50 0.75) save(`outputdir'/heatmap_cchg.pdf) out zti("YoY Consumption change")
heatmap housing_chg agebought period if model != "experiment_monetary_little", ///
    polbreak(2, 4) customf(0.05 0.10 0.20 0.50 0.75) save(`outputdir'/heatmap_hchg.pdf) out zti("YoY Housing change")
* This plots heatmaps of experiment 3) in the model slides
heatmap h_c_ratio agebought period if model != "experiment_monetary_nodown", ///
    polbreak(2, 4) save(`outputdir'/heatmap_h_little.pdf) zti("Housing -Consumption Ratio")
heatmap assets agebought period if model != "experiment_monetary_nodown", ///
    count polbreak(2, 4) save(`outputdir'/heatmap_a_little.pdf) zti(Financial Assets)
