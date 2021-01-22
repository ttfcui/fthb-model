local IRSdiv = 0.00654 // K-L divergence on 2009 dist. vs. 10-year counterfactual (share_hat)

sort agebin desc subtype
foreach var of varlist value* {
    by agebin: gen `var'_meandiff = `var'[5] - `var'[2]
	by agebin: gen `var'_skew = (`var'[5] - `var'[4])/`var'[6]
	by agebin: gen `var'_skewdiff = ///
	    (`var'[5] - `var'[4])/`var'[6] - (`var'[2] - `var'[1])/`var'[3]
	by agebin: gen `var'_hell = `var'[7]
	order `var'_h `var'_m `var'_s*, after(`var')
}
by agebin: keep if _n == _N
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

twoway (connect ageAll EIS, lp("-.")), by(metric, yr)
graph export "C:\Users\ttfcui\Documents\student-debt-fthb\code\model\EISplots.pdf", as(pdf) replace
twoway (connect age22_29 EIS, lp("-.")), by(metric, yr)
graph export "C:\Users\ttfcui\Documents\student-debt-fthb\code\model\EISplots_20s.pdf", as(pdf) replace
twoway (connect age30_39 EIS, lp("-.")), by(metric, yr)
graph export "C:\Users\ttfcui\Documents\student-debt-fthb\code\model\EISplots_30s.pdf", as(pdf) replace

twoway (connect ageAll EIS if metric==1, lp("-.") m(S) msize(small)), ///
    graphregion(c(white)) yti("KL Divergence Btwn Distributions") ///
	yline(`IRSdiv', lp(dash))
graph export "C:\Users\ttfcui\Documents\student-debt-fthb\code\model\output\stata\EISsurface.pdf", as(pdf) replace

