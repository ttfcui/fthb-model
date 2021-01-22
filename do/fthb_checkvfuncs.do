/*******************************************************************************
    Title: fthb_checkvfuncs.do
    Purpose: Visualizes vfuncs and policies for a certain model simulation.
    Author: TC, DB
    Date: 2019-06-19

*******************************************************************************/

args incomestate hetero

capture program drop vfunc_plotting
program define vfunc_plotting
    syntax, value(integer) incomestate(integer) hetero(string)
    
	import delimited fthb/vfuncs/vfunc`value'`hetero'.txt, delimiter(" ", collapse) clear
	drop v1
	
	egen incstate = rank(v3), track
	replace v6 = . if v6 < -1e4
	replace v7 = . if v7 < -1e4
	replace v8 = . if v8 < -1e4
	egen incval = group(incstate)
        replace incval = (incval-1)*2 + 1
	tab incval
	sum v13 if incval==`incomestate'

	egen housestate = rank(v4), track
	egen hval = group(housestate)
	keep if hval == 1
	tab v4

	rename (v6-v9 v5 v10-v12) (EV Vadjust Vrent Vnoadjust assets NextAssets NextOwnDurables NextDurables)
	gen Vdiff = Vadjust - Vrent

	twoway (line Vadjust Vrent assets if incval==`incomestate' & assets >= 0 & assets <= 1), ///
	    name(vdiff_age`value', replace)
	graph export "fthb/vfuncs/Vfuncs`value'_`incomestate'`hetero'.pdf", as(pdf) replace

        twoway (line NextDur NextOwn NextAssets assets if incval==`incomestate' & assets >= 0 & assets <= 1) ///
	    (function y=x, range(0 1) lc(gs9) lp(dash)), name(vdiff_policies`value', replace)
        graph export "fthb/vfuncs/policies`value'_`incomestate'`hetero'.pdf", as(pdf) replace

end

capture program drop vfunc_lifecycle
program define vfunc_lifecycle
    syntax, incomestate(integer) hetero(string)
    
	import delimited fthb/vfuncs/vfuncLC`hetero'.txt, delimiter(" ", collapse) clear
	drop v1
	keep if v4 == 0.0  // now this is assets, not durables
	egen incstate = rank(v3), track
        egen incval = group(incstate)
        replace incval = (incval-1)*2 + 1
	replace v6 = . if v6 < -1e7
	replace v7 = . if v7 < -1e7
	replace v8 = . if v8 < -1e7
	tab incval
	sum v13 if incval==`incomestate'

	rename (v6-v9 v5 v10-v12) (EV Vadjust Vrent Vnoadjust year NextAssets NextOwnDurables NextDurables)
	gen Vdiff = Vadjust - Vrent

	twoway (line Vdiff year if incval==`incomestate'), ///
	    name(vdiff_LC, replace)

        twoway (line NextDur NextOwn NextAssets year if incval==`incomestate'), ///
	    name(vpol_LC, replace)

end

* Individual calls to function calls here
vfunc_plotting, value(1) incomestate(`incomestate') hetero(`hetero') // age 5/27
vfunc_plotting, value(2) incomestate(`incomestate') hetero(`hetero') // age 37/59
vfunc_plotting, value(3) incomestate(`incomestate') hetero(`hetero') // age 37/59
vfunc_plotting, value(4) incomestate(`incomestate') hetero(`hetero') // age 38/60
vfunc_plotting, value(5) incomestate(`incomestate') hetero(`hetero') // age 62/84
vfunc_plotting, value(6) incomestate(`incomestate') hetero(`hetero') // age 63/85*
