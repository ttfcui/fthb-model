/*******************************************************************************
    Title: xp_12202016.do
    Purpose: Event study plots for policy-induced buyers.
    Author: TC, DB
    Date: 2016-03-31

*******************************************************************************/

local outputdir $outdir/stata/xp_12202016

/* %< Steady-state */
import delimited $dboxdir/model/stationary_eventstudy.csv, clear
local lowinc = 5
local highinc = 9
tab v1
sort v1 age income time
order time, after(income)
egen group = group(income age), l

#delimit ;
twoway (line assets housing time if v1=="negshock", lc(navy maroon) lp(solid))
    (line assets housing time if v1=="nof1shock", lc(navy maroon) lp(dash dash))
    if income >= `lowinc' & income <= `highinc', by(group, colf $graphconfig) xline(0, lc(gs9))
    legend(label(1 "Liquid Assets, Pos. Shock") label(2 "Housing, Pos. Shock")
           label(3 "Liquid Assets, No Shock") label(4 "Housing, No Shock")) name(pos, replace);
twoway (line assets housing time if v1=="negshock", lc(navy maroon) lp(solid))
    (line assets housing time if v1=="nof1shock", lc(navy maroon) lp(dash dash))
    if income >= `lowinc' & income <=  `highinc', by(group, colf $graphconfig) xline(0, lc(gs9))
    legend(label(1 "Liquid Assets, Neg. Shock") label(2 "Housing, Neg. Shock")
           label(3 "Liquid Assets, No Shock") label(4 "Housing, No Shock")) name(neg, replace);
#delimit cr
/* %> */

capture program drop get_labels
program define get_labels

    label define var_pol -2 "var_l2" -1 "var_l1" 0 "var_pol0" 1 "var_pol1" ///
        2 "var_pol2" 3 "var_pol3" 4 "var_pol9", replace
    label define var_cf -2 "var_l2" -1 "var_l1" 0 "var_ss" 1 "var_pol1_cf" ///
        2 "var_pol2_cf" 3 "var_pol3_cf" 4 "var_pol9_cf", replace
    encode variable if regexm(variable, "cf|ss") == 0, gen(var1) label(var_pol)
    encode variable if regexm(variable, "_l|cf|ss") == 1, gen(var2) label(var_cf)
    label drop var_pol var_cf
    label define period_split 4 "9"
    label val var? period_split
end

* %< Load the data (generated from xp_12202016.py)
import delimited $dboxdir/model/transition_eventstudy.csv, clear
* 10 periods are included in the data, 6 for world with policy and 6
* for the no-policy counterfactual (so the 2 periods before policy enactment
* is identical for both)
tostring income_state, replace
replace income_state = income_state + "/13"
replace age = age + 20
get_labels

* The actual plots are made here. Adjust the list of age values in i to
* create plots for other age groups, or change it around so there are separate
* plots by age, not by income state during policy takeup
* VARIABLES: C is the consumption good, H_own is owned housing (=0 for renters)
* H_rent is rental housing (=0 for owners). Plotted points are medians.
foreach i in 24 27 30 35 40 /* 50 */ {
foreach variable in C H_own H_rent {
twoway (connect value var1 if var1 < 4) (scatter value var1 if var1 == 4, mc(navy)) ///
    (connect value var2 if var2 < 4, lp("-..") color(gs9) m(Dh)) ///
    (scatter value var2 if var2 == 4, color(gs9) m(Dh)) ///
    if age == `i' & var == "`variable'", legend(order(1 3) label(1 "Policy response") ///
    label(3 "Steady-state counterfactual") rows(2)) xline(0, lp("-...")) ///
    xline(3.2, lc(gs7)) xlab(-2(1)4, valuelabel) ///
    ylab(, angle(0)) by(income_state, title("`variable', age=`i'") $graphconfig) ///
    name(`variable', replace) yti(Variable value) xti(Period around policy takeup)
}
graph export "`outputdir'/eventstudy_age_`i'_var_C.pdf", as(pdf) name(C) replace
* The housing plots are combined together to more easily compare upsizing/downsizing
graph combine H_own H_rent, graphregion(c(white)) ycom xcom
graph export "`outputdir'/eventstudy_age_`i'_var_housing.pdf", as(pdf) replace
}
* %>

* %< Load the data indexed by wealth bins instead of age (also from same py)
import delimited $dboxdir/model/transition_eventstudy_wealth.csv, clear
tostring income_state, replace
gen wealth = string(finwealthbin)
replace wealth = ".0" if wealth=="0"
replace income_state = income_state + "/13"
get_labels

foreach z in 5 7 9 /* 50 */ {
foreach variable in C H_own H_rent {
twoway (connect value var1 if var1 < 4) (scatter value var1 if var1 == 4, mc(navy)) ///
    (connect value var2 if var2 < 4, lp("-..") color(gs9) m(Dh)) ///
    (scatter value var2 if var2 == 4, color(gs9) m(Dh)) ///
    if income_state == "`z'/13" & var == "`variable'", legend(order(1 3) label(1 "Policy response") ///
    label(3 "Steady-state counterfactual") rows(2)) xline(0, lp("-...")) ///
    xline(3.2, lc(gs7)) xlab(-2(1)4, valuelabel) ///
    ylab(, angle(0)) by(wealth, title("`variable', income=`z'/13") $graphconfig) ///
    name(`variable', replace) yti(Variable value) xti(Period around policy takeup)
}
graph export "`outputdir'/wealthstudy_inc_`z'_var_C.pdf", as(pdf) name(C) replace
graph export "`outputdir'/wealthstudy_inc_`z'_var_rental.pdf", as(pdf) name(H_rent) replace
graph export "`outputdir'/wealthstudy_inc_`z'_var_own.pdf", as(pdf) name(H_own) replace
}
* %>

* %< Load data but filtered to agents receiving positive shock
import delimited $dboxdir/model/transition_eventstudy_negshock.csv, clear
* 10 periods are included in the data, 6 for world with policy and 6
* for the no-policy counterfactual (so the 2 periods before policy enactment
* is identical for both)
tostring income_state, replace
replace income_state = income_state + "/13"
replace age = age + 20
get_labels

* The actual plots are made here. Adjust the list of age values in i to
* create plots for other age groups, or change it around so there are separate
* plots by age, not by income state during policy takeup
* VARIABLES: C is the consumption good, H_own is owned housing (=0 for renters)
* H_rent is rental housing (=0 for owners). Plotted points are medians.
foreach i in 24 27 30 35 40 /* 50 */ {
foreach variable in C H_own H_rent {
twoway (connect value var1 if var1 < 4) (scatter value var1 if var1 == 4, mc(navy)) ///
    (connect value var2 if var2 < 4, lp("-..") color(gs9) m(Dh)) ///
    (scatter value var2 if var2 == 4, color(gs9) m(Dh)) ///
    if age == `i' & var == "`variable'", legend(order(1 3) label(1 "Policy response") ///
    label(3 "Steady-state counterfactual") rows(2)) xline(0, lp("-...")) ///
    xline(3.2, lc(gs7)) xlab(-2(1)4, valuelabel) ///
    ylab(, angle(0)) by(income_state, title("`variable', age=`i'") $graphconfig) ///
    name(`variable', replace) yti(Variable value) xti(Period around policy takeup)
}
graph export "`outputdir'/eventstudy_negshock_age_`i'_var_C.pdf", as(pdf) name(C) replace
* The housing plots are combined together to more easily compare upsizing/downsizing
graph combine H_own H_rent, graphregion(c(white)) ycom xcom
graph export "`outputdir'/eventstudy_negshock_age_`i'_var_housing.pdf", as(pdf) replace
}
* %>

* %<
import delimited $dboxdir/model/transition_eventstudy_wealth_negshock.csv, clear
tostring income_state, replace
gen wealth = string(finwealthbin)
replace wealth = ".0" if wealth=="0"
replace income_state = income_state + "/13"
get_labels

foreach z in 5 7 9 /* 50 */ {
foreach variable in C H_own H_rent {
twoway (connect value var1 if var1 < 4) (scatter value var1 if var1 == 4, mc(navy)) ///
    (connect value var2 if var2 < 4, lp("-..") color(gs9) m(Dh)) ///
    (scatter value var2 if var2 == 4, color(gs9) m(Dh)) ///
    if income_state == "`z'/13" & var == "`variable'", legend(order(1 3) label(1 "Policy response") ///
    label(3 "Steady-state counterfactual") rows(2)) xline(0, lp("-...")) ///
    xline(3.2, lc(gs7)) xlab(-2(1)4, valuelabel) ///
    ylab(, angle(0)) by(wealth, title("`variable', income=`z'/13") $graphconfig) ///
    name(`variable', replace) yti(Variable value) xti(Period around policy takeup)
}
graph export "`outputdir'/wealthstudy_negshock_inc_`z'_var_C.pdf", as(pdf) name(C) replace
graph export "`outputdir'/wealthstudy_negshock_inc_`z'_var_rental.pdf", as(pdf) name(H_rent) replace
graph export "`outputdir'/wealthstudy_negshock_inc_`z'_var_own.pdf", as(pdf) name(H_own) replace
}
* %>


capture program drop trans_labels
program define trans_labels

    tostring income_state, replace
    gen wealth = string(finwealthbin)
    replace wealth = ".0" if wealth=="0"
    replace income_state = income_state + "/13"

    label define var_pol -2 "var_l2" -1 "var_l1" 0 "var_pol0" 1 "var_pol1" ///
        2 "var_pol2" 3 "var_pol3" 4 "var_pol9", replace
    label define var_cf -2 "var_l2" -1 "var_l1" 0 "var_ss" 1 "var_pol1_cf" ///
        2 "var_pol2_cf" 3 "var_pol3_cf" 4 "var_pol9_cf", replace
    encode variable if regexm(variable, "cf|ss") == 0, gen(var1) label(var_pol)
    encode variable if regexm(variable, "_l|cf|ss") == 1, gen(var2) label(var_cf)
    label drop var_pol var_cf
    label define period_split 4 "9"
    label val var? period_split
end

capture program drop event_plots
program define event_plots
    syntax, datadir(string)

    tempfile noshock posshock negshock
    import delimited `datadir'/transition_eventstudy_wealth.csv, clear
    trans_labels
    gen fullLab = income_state + "  " + var
    keep if regexm(fullLab, "(5|7|9).+(C|H_)") == 1
    save `noshock'

    import delimited `datadir'/transition_eventstudy_wealth_posshock.csv, clear
    trans_labels
    gen fullLab = income_state + " Pos  " + var
    keep if regexm(fullLab, "(7).+(C|H_)") == 1
    save `posshock'

    import delimited `datadir'/transition_eventstudy_wealth_negshock.csv, clear
    trans_labels
    gen fullLab = income_state + " Neg  " + var
    keep if regexm(fullLab, "(7|9).+(C|H_)") == 1
    save `negshock'

    use `noshock', clear
    append using `posshock'
    append using `negshock'

    twoway (connect value var1 if var1 < 4) (scatter value var1 if var1 == 4, mc(navy)) ///
        (connect value var2 if var2 < 4, lp("-..") color(gs9) m(Dh)) ///
        (scatter value var2 if var2 == 4, color(gs9) m(Dh)) ///
        if wealth ==".0" , legend(order(1 3) label(1 "Policy response") ///
        label(3 "Steady-state counterfactual") rows(2)) xline(0, lp("-...")) ///
        xline(3.2, lc(gs7)) xlab(-2(1)4, valuelabel) ///
        ylab(, angle(0)) by(fullLab, row(3) yrescale colfirst $graphconfig) ///
        name(event, replace) yti(Variable value) xti(Period around policy takeup)

end
