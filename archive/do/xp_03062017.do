/*******************************************************************************
    Title: xp_03062017.do
    Purpose: Playing with "master" individual-level model simulations.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

* Read in master data and select canonical policy
do $repodir/do/masterdata.do
local outputdir $outdir/stata/xp_03072017
*keep if model == "experiment_monetary_nodown"

* %< Labels
gen extensive = cond(!missing(posext), 1, cond(!missing(negext), 2, .))
gen asset_filter = (assets <= 0.01)
* Decomposing policy response
gen decomp = cond(missing(pullforward), 1, cond( ///
    pullforward == 1 & adjDiff == 0, 2, cond( ///
    pullforward == 2 & adjDiff == 0, 3, cond( ///
    pullforward >= 3 & pullforward <= 5 & adjDiff == 0, 4, cond( ///
    pullforward > 5 & adjDiff == 0, 5, cond( ///
    adjustment_pol == 1 & adjustment_ss == 0, 6, cond( ///
    adjDiff > 0 & !(adjustment_pol == 1 & adjustment_ss == 0), 7, cond( ///
    adjDiff < 0, 8, .))))))))
* Indicators for coming in/out of income shocks
gen posIncShock = (income > income_l)
gen negIncShock = (income < income_l)
gen posIncShock2 = (income < income_f1)
gen negIncShock2 = (income > income_f1)
egen groups = group(marginal extensive)

* Graph labelling
label def groups 1 "Infra/Positive" 2 "Infra/Negative" ///
    3 "Marginal/Positive" 4 "Marginal/Negative"
label def decomp 1 "Inframarginal" 2 "Timing margin, 1 period" ///
    3 "Timing margin, 2 periods" 4 "Timing margin, 3-5 periods" ///
    5 "Timing margin, 5+ periods" 6 "0-1 Extensive Margin" ///
    7 "Other Pos. Extensive Margin" 8 "Neg. Extensive Margin"
label def assetSplit 0 "Assets > 1e-2" 1 "Assets <= 1e-2"
label val decomp decomp
label val groups groups
label val asset_filter assetSplit
label var posIncShock "Proportion w/ Positive income shock, policy period"
label var negIncShock "Proportion w/ Negative income shock, policy period"
label var posIncShock2 "Proportion w/ Positive income shock, period after policy"
label var negIncShock2 "Proportion w/ Negative income shock, period after policy"
*preserve
* %>

* %< Decompose policy-induced transactions and housing investment
collapse (count) id (sum) nextDurables ss=nextDurables_ss, by(decomp)
scalar inv_ss = ss[1]
gen idprop = id/id[1]
gen durprop = nextD/inv_ss
eststo counts: estpost tab decomp [iw=idprop]
eststo inv: estpost tab decomp [iw=durprop]
esttab counts inv using `outputdir'/decomp_table.tex, noobs not nostar ///
    varlabels(`e(labels)') mlab("FTHBs" "Investment") tex replace
* %>

* %< 2x2 histograms of different groups on timing/extensive margin
restore, preserve
hist income, by(groups, $graphconfig) width(1) lcolor(purple) fcolor(none)
graph export `outputdir'/incomes_comp.pdf, as(pdf) replace

replace age = age + 21
label var age "Age"
hist age, by(groups, $graphconfig) width(1) lcolor(purple) fcolor(none)
graph export `outputdir'/ages_comp.pdf, as(pdf) replace

hist assets if assets < 3, by(groups, $graphconfig) lcolor(purple) fcolor(none)
graph export `outputdir'/assets_comp.pdf, as(pdf) replace

tempfile rent_stationary
import delimited using $outdir/rentAssets.csv, clear
rename nextassets nextAssets
save `rent_stationary'
restore, preserve
append using `rent_stationary'
replace marginal = -1 if missing(marginal)
twoway (hist nextAssets if marginal==-1 & nextAssets < 3, ///
        width(0.05) fcolor(blue*0.2) lcolor(white)) ///
       (hist assets if marginal==0, width(0.05) fcolor(none) lcolor(purple)), ///
       legend(label(1 "Stationary equilibrium") label(2 "Inframarginal")) $graphconfig
graph export `outputdir'/assets_planned.pdf, as(pdf) replace
twoway (hist nextAssets if marginal==-1 & nextAssets < 3, ///
        width(0.05) fcolor(blue*0.1) lcolor(white)) ///
       (hist assets if marginal==1, width(0.05) fcolor(none) lcolor(purple)), ///
       legend(label(1 "Stationary equilibrium") label(2 "Marginal")) $graphconfig
graph export `outputdir'/assets_marginal.pdf, as(pdf) replace
* %>

* %< Summary statistics for variables by above groups. The interesting
* thing here is really the income shock indicator proportions
local varlist nextAssets* *IncSho*
bys groups: sum `varlist'
by groups: sum `varlist' if income==7
by groups: sum `varlist' if income==9
eststo groups: estpost tabstat nextAssets posIncShock posIncShock2 negIncShock2, ///
    by(groups) statistics(mean) columns(statistics) nototal
esttab groups using `outputdir'/extstat_table.tex, main(mean) unstack ///
    noobs not nostar l replace
* %>

* %< Heatmaps of income history for policy takers, in policy period and period
* before. Attempt at showing both distribution of income and if people
* purchase coming in or out of income shocks
clear
tempfile income_temp // first, build a square grid of income values
set obs 9
egen income_l = fill(3/11)
expand 9
bys income_l: gen income = _n + 2
expand 2
bys income_l income: gen marginal = _n - 1
save `income_temp'
restore, preserve

collapse (count) id, by(marginal income income_l)
by marginal: egen id_total = total(id)
gen id_prop = id/id_total*100
merge 1:1 marginal income income_l using `income_temp'
heatmap id_prop income income_l if marginal==1 & income > 3 & income_l > 3, ///
    save(`outputdir'/income_m.pdf) count yti("Income state, policy period") ///
    xti("Income state, policy -1") zti("%age of group") zlab(3)
heatmap id_prop income income_l if marginal==0 & income > 3 & income_l > 3, ///
    save(`outputdir'/income.pdf) count yti("Income state, policy period") ///
    xti("Income state, policy -1") zti("%age of group")zlab(3)
* %>

* %< From how many periods forward are those with zero assets pulling
* purchase decision? Everyone else?
restore, preserve
twoway (hist pullforward if pullforward < 40, width(1) fc(none) lc(red)), ///
    by(asset_filter, r(2) $graphconfig)
graph export `outputdir'/assets_pull.pdf, as(pdf) replace
* %>

* %< Correlates between extra transactions pulled forward
binscatter adjDiff incF_mean if adjDiff > 0, control(income) n(30) ///
    xti("Mean, future income state realizations") ///
    yti("Extra transactions due to policy") reportreg
graph export `outputdir'/adjdiff_income.pdf, as(pdf) replace
binscatter adjDiff assets if adjDiff > 0 & assets < 3, control(income) n(30) ///
    xti("Financial assets before policy") yti("Extra transactions due to policy") ///
    reportreg
graph export `outputdir'/adjdiff_assets.pdf, as(pdf) replace

* Correlates between magnitude of negative ext. margin (fewer purchases
* over lifecycle)
binscatter negext incF_std, n(30) control(income) line(qfit) ///
    xti("S.D., future income state realizations") ///
    yti("First period on neg. extensive margin") reportreg
graph export `outputdir'/negext_income.pdf, as(pdf) replace

binscatter negext assets if assets < 3, n(30) control(income) line(qfit) ///
    xti("Financial assets before policy") ///
    yti("First period on neg. extensive margin") reportreg
graph export `outputdir'/negext_assets.pdf, as(pdf) replace
* %>

* %<
* Generating partial plot of income state by periods housing decision
* pulled forward, as well as people going from 0 to multiple
* purchases (the "40")
* These plots are also divided further according to extensive margin:
* aggregate, more purchases, same #, fewer.
restore, preserve
tempfile pullf_inc temp
collapse (mean) income if pullforward <= 10 | pullforward==40 , ///
    by(pullforward)
save `pullf_inc'
restore, preserve
collapse (mean) income if pullforward <= 10 | pullforward==40, ///
    by(pullforward extensive)
replace extensive = 0 if missing(extensive)
save `temp'
use `pullf_inc', clear
append using `temp'
save `pullf_inc', replace

* Generating income means for those pulling 2+ periods, 5+ periods, 10+ periods
foreach val in 2 5 10 {
    restore, preserve
    collapse (mean) income_plus=income if ///
        pullforward > `val' & pullforward < 40
    gen pullf_bar = `val' + 1
    save `temp', replace
    local m`val' = income_plus[1]
    use `pullf_inc', clear
    append using `temp'
    save `pullf_inc', replace

    restore, preserve
    collapse (mean) income_plus=income if ///
        pullforward > `val' & pullforward < 40, by(extensive)
    replace extensive = 0 if missing(extensive)
    gen pullf_bar = `val' + 1
    save `temp', replace
    local m`val'pos = income_plus[1]
    local m`val'neg = income_plus[2]
    use `pullf_inc', clear
    append using `temp'
    save `pullf_inc', replace

}
restore
use `pullf_inc', clear

replace pullf_bar = pullforward if pullforward <= 2 | pullforward == 40
replace pullf_bar = 4 if pullf_bar == 6
replace pullf_bar = 5 if pullf_bar == 11
replace pullf_bar = 6 if pullf_bar == 40
label def pullf 1 "1" 2 "2" 3 "2+" 4 "5+" 5 "10+" 6 "Pos. Extensive"
label val pullf_ pullf

#delimit ;
twoway (bar income pullf_bar, barw(0.7) color(blue*0.1))
    (bar income_plus pullf_bar, barw(0.7) color(blue*0.1))
    if missing(extensive), ysc(r(7 7.8)) ylab(7(0.2)7.8) legend(off)
    xti("Imputed periods durable purchased moved forward") yti("Income state, mean")
    xlab(,valuelabel) $graphconfig 
    note("Note: Income is an AR(1) process discretized into 13 states.");
graph export `outputdir'/pulledInc_1.pdf, as(pdf) replace;

twoway (bar income pullf_bar if missing(extensive), barw(0.7) color(blue*0.1))
    (bar income_plus pullf_bar if missing(extensive), barw(0.7) color(blue*0.1))
    (bar income pullf_bar if extensive == 1, barw(0.7) fcolor(none) lcolor(purple) lw(thick))
    (bar income_plus pullf_bar if extensive == 1, barw(0.7) fcolor(none) lcolor(purple) lw(thick))
    , ysc(r(7 7.8)) ylab(7(0.2)7.8) xti("Imputed periods durable purchased moved forward")
    yti("Income state, mean") xlab(,valuelabel) $graphconfig
    legend(order(1 3) label(1 "Aggregate") label(3 "Positive Ext Margin"))
    note("Note: Income is an AR(1) process discretized into 13 states.");
graph export `outputdir'/pulledInc_2.pdf, as(pdf) replace;
#delimit cr
exit

/* OLD SCATTERPLOTS: TBD
* Replacement for axis readability issues
replace pullforward = 13 if pullforward == 40
* Plot aggregate and more purchases income levels together
twoway (scatter income pullforward if missing(extensive)), ///
    yline(`m2' `m5' `m10', lc(gs9) lp(dash)) ysc(r(`m10' 7)) ///
    name(all, replace) $graphconfig
graph export `outputdir'/pulledInc_1.pdf, as(pdf) replace
twoway (scatter income pullforward if missing(extensive)) ///
    (scatter income pullforward if extensive==1), ///
    yline(`m2pos' `m5pos' `m10pos', lc(gs9) lp(dash)) ///
    ysc(r(`m10pos' 7)) name(posext, replace) legend(off) $graphconfig
graph export `outputdir'/pulledInc_2.pdf, as(pdf) replace
* %> */
