/*******************************************************************************
    Title: mpcTest.do
    Purpose: Test the distribution of MPCs simulated from model.
    Author: TC, DB
    Date: 2020/10

*******************************************************************************/

import delimited indivMPCs.csv, clear
foreach var of varlist consumption* {
	sum `var', d
	replace `var' = . if `var' <= `r(p5)' | `var' >= `r(p95)' 
}

#delimit ;
twoway (kdensity consumption_mpc if rentalmpc==0)
    (kdensity consumption_mpc if rentalmpc==1, lp(dash)),
    legend(label(1 "Owners after Transfer") label(2 "Renters after Transfer"))
    xti("MPCs, Consumption Only") yti("") name(mpc, replace);
graph export "mpcs.pdf", as(pdf) replace;
twoway (kdensity consumption_mpcwdur if rentalmpc==0)
    (kdensity consumption_mpcwdur if rentalmpc==1, lp(dash)),
    legend(label(1 "Owners after Transfer") label(2 "Renters after Transfer"))
    xti("MPCs, Consumption & Owned Housing Value") yti("") name(mpcDur, replace);
graph export "mpcsDur.pdf", as(pdf) replace;
twoway (kdensity consumption_mpcwusecosts if rentalmpc==0)
    (kdensity consumption_mpcwusecosts if rentalmpc==1, lp(dash)),
    legend(label(1 "Owners after Transfer") label(2 "Renters after Transfer"))
    xti("MPCs, Consumption & Housing User Costs") yti("") name(mpcUseCosts, replace);
graph export "mpcsUserCost.pdf", as(pdf) replace;
twoway (kdensity consumption_mpc if rentalmpc==0)
    (kdensity consumption_mpc if rentalmpc==1, lp(dash)) if pastrent==1,
    legend(label(1 "Owners after Transfer") label(2 "Renters after Transfer"))
    xti("MPCs, Consumption Only") yti("") name(mpcRenters, replace);
graph export "mpcs_renters.pdf", as(pdf) replace;
twoway (kdensity consumption_mpcwdur if rentalmpc==0)
    (kdensity consumption_mpcwdur if rentalmpc==1, lp(dash)) if pastrent==1,
    legend(label(1 "Owners after Transfer") label(2 "Renters after Transfer"))
    xti("MPCs, Consumption & Owned Housing Value") yti("") name(mpcDurRenters, replace);
graph export "mpcsDur_renters.pdf", as(pdf) replace;
twoway (kdensity consumption_mpcwusecosts if rentalmpc==0)
    (kdensity consumption_mpcwusecosts if rentalmpc==1, lp(dash)) if pastrent==1,
    legend(label(1 "Owners after Transfer") label(2 "Renters after Transfer"))
    xti("MPCs, Consumption & Housing User Costs") yti("") name(mpcUseCostsRenters, replace);
graph export "mpcsUserCosts_renters.pdf", as(pdf) replace;
#delimit cr

binscatter consumption_mpc age, discrete by(pastrent) line(none) legend(label(1 "Homeowners before Transfer") label(2 "Renters before Transfer") cols(2)) yti("Mean MPCs, Nondurable Consumption")
binscatter consumption_mpc age, discrete by(pastrent) line(none) legend(label(1 "Homeowners before Transfer") label(2 "Renters before Transfer") cols(1)) yti("Mean MPCs, Nondurable Consumption")
