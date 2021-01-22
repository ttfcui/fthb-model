*use $dboxdir/model/masterdata_test.dta, clear

label var id "ID within age cohort"
label var age "Age Minus 20"

order nextDurables nextAssets consumption nextDurables_ss ///
    nextAssets_ss consumption_ss, after(age)
label var nextDurables "Durable choice given policy"
label var nextAssets "Assets/Loan choice given policy"
label var consumption "Consumption choice given policy"
label var nextDurables_ss "Durable choice absent policy"
label var nextAssets_ss "Assets/Loan choice absent policy"
label var consumption_ss "Consumption choice absent policy"

order durables assets rent marginal durTime durablesOrig durablesGap, ///
    after(consumption_ss)
label var durables "Existing owned/rented durable stock"
label var assets "Existing Assets/Loans"
label var rent "1 == Durable is rented"
label var marginal "1 == Policy induced buyer"
label var durTime "Age of last durable adjustment"
label var durablesOrig "Imputed value of last durable adjustment"
label var durablesGap "Imputed difference btwn current and last durable"

order income_l2 income_l, after(durablesGap)
label var income_l2 "Income state, T-2"
label var income_l "Income state, T-1"
label var income "Income state, T (policy period)"
label var income_f1 "Income state, T+1"
label var income_f2 "Income state, T+2"
label var income_f3 "Income state, T+3"
label var income_f4 "Income state, T+4"
label var incF_mean "Mean over future income history"
label var incF_std "Standard deviation over future income history"

label var adjustment_pol "Lifetime Tracked Adjustments, Policy"
label var adjustment_ss "Lifetime Tracked Adjustments, Steady State"
label var adjDiff "adj*_pol - adj*_ss"

label var pullforward "Imputed periods durable purchased moved forward"
label var posext "First purchase on pos. extensive margin"
label var negext "First purchase on neg. extensive margin"

order C*, after(negext)
label var C_pol0 "Identical to consumption"
label var C_pol1 "Consumption in policy simulation, T+1"
label var C_pol2 "Consumption in policy simulation, T+2"
label var C_pol3 "Consumption in policy simulation, T+3"
label var C_pol1_cf "Consumption in stationary counterfactual, T+1"
label var C_pol2_cf "Consumption in stationary counterfactual, T+2"
label var C_pol3_cf "Consumption in stationary counterfactual, T+3"

*saveold masterdata_test.dta, replace
