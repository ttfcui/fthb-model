/*******************************************************************************
    Title: SCF_calibration.do
    Purpose: Clean SCF data for model moment calibration.
    Author: TC, DB
    Date: 2020/11

*******************************************************************************/

args datadir initdir

************************************************************
* %< 0. Data prep programs.
************************************************************

* %< SORT VARIABLES FROM FULL SURVEY DATA 

tempfile more1998 more2001 more2004 more2007
local pensions X4206 X4219 X4806 X4819
scalar CPIbase = 3438 // 2013 index
* Inflation matrix. First row are denominators for conversion to 2013 dollars,
* Second and third rows form ratio for income stats from last year
matrix inflAdj = (2405, 2618, 2788, 3062 \ 2397, 2600, 2774, 3045 \ 2364, 2529, 2701, 2961)
local j = 1

forval years = 1998(3)2007 {

    # delimit ;
    * New pension questions from the 04 survey onward;
    if `years' >= 2004 local pensions X11044 X11049 X11244 X11249;
    local vars Y1 X410 X7973 X7976 X413 X421 X424 X427 X430
        X414 X7132 X816 X432 X407 X409 X5704 X5724 X5725 X3015 X3016
        X3017 X7510 X7509 X4112 X4712 X4511 X5111 X4113 X4713 X4131 X4731 X4132 X4732 `pensions';
    use `vars' using `datadir'/raw_data/scf`years'.dta, clear;
    gen year = `years';
    order year `vars';
    rename (`vars') (y1 hasCC hasMcVisa hasAmex revbalance1 revbalance2 revbalance3
        revbalance4 revbalance5 maxCredit intrateCC intrateMM payFreq turnedDown
        notApplied hhSelfY othInc sourceOthInc noSaveBor noSaveZero saveWhatTa
        spendMoreY buyHome labEarn1 labEarn2 fullEarn1 fullEarn2 freqLe1 freqLe2
	selfEarn1 selfEarn2 freqSe1 freqSe2 hContrib hEmpContrib wContrib wEmpContrib);
    # delimit cr

    * Credit card dummies, sum debt for those paying it off
    replace hasCC = 0 if hasCC == 5 | (hasMcVisa == 5 & hasAmex == 5)
    foreach var of varlist revbalance1-intrateMM hContrib-wEmpContrib{
        replace `var' = 0 if `var' == -1 | `var' == -2
    }
    gen ccDebt = 0
    replace ccDebt = (revbalance1 + revbalance2 + revbalance3 + revbalance4 ///
        + revbalance5) if payFreq > 1
    gen credDenied = (inlist(turnedDown, 1, 3))

    * Not sure what these are?
    gen htm1 = (noSaveBor == 1 | noSaveZero == 1)
    gen htm2 = (htm1 == 1 | saveWhatTa == 1)
    gen htm4 = (spendMoreY == 1 | spendMoreY == 2) & buyHome == 5
    gen htm3 = (htm4 == 1 & buyHome == 5)

    * Replace codes in variables with their meaning, nothing or missing
    replace othInc = 0 if othInc == -1
    foreach var of varlist *Earn1 *Earn2 {
        replace `var' = 0 if `var' < 0
    }
    foreach var of varlist freq?e? {
        replace `var' = . if `var' < 0
    }
    gen labIncAlt = labEarn1 + labEarn2 + selfEarn1 + selfEarn2

    * Adjust for inflation to sync with summary stats
    foreach var of varlist ccDebt labIncAlt hhSelfY othInc labIncAlt *Contrib {
        if inlist("`var'", "hhSelfY", "othInc") {
            replace `var' = `var' * inflAdj[2, `j'] / inflAdj[3, `j']
        }
        replace `var' =  `var' * CPIbase / inflAdj[1, `j']
    }

    save `more`years'', replace
    local ++j
}

* %>

* %< LOAD AHS SUMMARY STATISTICS TO RECOVER MOMENT USED IN CALIBRATION
/* I (TC) transcribed median house prices and median housing costs for rentors
   from the AHS, tables 3-14 & 4-13 respectively. Values are recorded for
   every other year between 1999 and 2005.
   Right now we take the time average first, then take the ratio of the two.
   Does it make that big a difference? */
import delimited `datadir'/raw_data/ahs_stats.csv, clear
gen pi = 3438/cpi
collapse (mean) housevalue-rentcost [iw=pi]
scalar hp_ratio = housevalue[1]/(rentcost[1]*12)
* %>

* %< LOAD SUMMARY EXTRACT DATA AND JOIN WITH EARLIER TABLES
clear
cd `datadir'/SDA_output
do "scf-9807dict.txt"
rename *, lower

**DROP TOP PERCENTILES FOR NETWORTH
** THIS IS AN IMPORT SAMPLE RESTRICTION
foreach t of numlist 1989(3)2007{
	qui _pctile networth if year ==`t' [pweight = wgt], p(1 95)
	drop if (networth>=r(r2) | networth <= r(r1)) & year == `t' // What about lower tail?
}
*AGE RESTRICTIONS
keep if age>=18 & age<=79

preserve
tempfile demogs
clear
do "scf-9807Demo.txt"
rename *, lower
save `demogs'
restore
merge 1:1 caseid wgt using `demogs', keep(3)

*MERGE, KEEP ONLY YEARS DESIRED
keep if year >= 1998 & year <= 2004
forval years = 1998(3)2004 {
    mmerge year y1 using `more`years'', type(1:1) unm(master) update
    drop _m
}

*%>

* %< FURTHER VARIABLE DEFINITIONS

* GENERATE LABOR INCOME FOR WORKERS
gen labInc = wageinc + transfothinc - othInc
replace othInc = 0 if inlist(sourceOthInc, 11, 14, 30, 36)
replace labInc = labInc + othInc

gen labIncPlus = labInc + hhSelfY
gen incomeNoSelfY = income - hhSelfY // all income minus self-employment y
drop if labInc < 0
gen labIncWork = labInc if age < 60 // Model agents stop working at 60

* AGE BINS
gen ageind = (ceil(age / 5) - 4)*5 + 20 // Base category is 20-24
gen ageyoung = (ceil(age / 2) - 10)*2 + 20 if age <= 70 // Base category is 20-22
replace ageind = 70 if ageind > 70 // Effective cap at 70

*-----------------
*DEFINITIONS
*-----------------
* $70 from 2008 Survey of Consumer Payment Choice Table 10, adjusted for inflation to
* 2013 dollars, multiplied by two people in household, divided by
* median liq in dataset
_pctile liq [pw=wgt], n(2)
local cashFrac = (70*(3438/3213)*1.07*2/r(r1))

gen liqPos      = liq*(1+`cashFrac')
gen liqPosNoret = call + checking + `cashFrac'*liq
gen direct      = nmmf + stocks + bond
gen directEq    = deq
gen directBond  = direct - deq
gen housePos    = houses // + oresre + nnresre
gen houseNeg    = mrthel + resdbt
gen netHouse    = housePos - houseNeg
gen netCars     = vehic
gen sb          = savbnd
gen certdep     = cds
gen retAcc      = retqliq
gen retAccEq    = reteq
gen lifeIns     = cashli
gen netbus      = 0 // bus

gen othAssets   = othfin + othma + othnfin //EXCLUDE because we also exclude other debt

gen netLiq      = liqPos - ccDebt
gen illiqPos    = housePos + netCars + direct + sb + certdep + retAcc + lifeIns + netbus 
gen netIlliq    = illiqPos - houseNeg
gen netWorthScf = networth - netbus

gen brLiqPos      = liqPos + direct
gen netBrLiq      = brLiqPos - ccDebt
gen brIlliqPos    = illiqPos - direct
gen netBrIlliq    = brIlliqPos - houseNeg
gen netWorthNew   = netBrIlliq + netBrLiq
gen netWorthnc    = netWorthNew - netCars

gen netliquid = netBrLiq
gen netliquid_exhousingdebt = netBrLiq - houseNeg
gen networth_exhousingvalue = netWorthnc - housePos
gen homeownershiprate = (houses > 0)

gen LTV = houseNeg/housePos
gen LTVcond = houseNeg/housePos if houseNeg > 0
gen H_W = housePos/(brLiqPos +brIlliqPos ) // retirement accounts excluded (matters here)
gen W = brLiqPos +brIlliqPos // // retirement accounts excluded (matters here)

* %>
* %>

************************************************************
* %< 1. Model initialization programs.
************************************************************

capture program drop initial_grid_calibration // %<
program define initial_grid_calibration
    syntax anything [pweight], initdir(string) theta(real) normed(real) hoyears(real)
    /* Output young adult observations to initialize Fortran simulation.

       Args:
           anything: Two integers that mark the age range of the subsample.
           initdir: output folder for the calibration file.
           theta: Calibrated down payment percentage.
           normed: Mean HH income in sample over which we normalize incomes.
           hoyears: Separate homeowners' age to further pare down sample.
    */
    preserve
    tempfile quantiles
    local outvars age homeownershiprate labIncWork housePos houseNeg netWorthnc ///
        netliquid netliquid_ex networth_ex wgt 
    keep `outvars'
    order `outvars'
    tokenize `anything'
    keep if age >= `1' & age <= `2'

    * "a" in Fortran model actually "voluntary equity," 
    * assets go in "model_l" instead
    gen model_a = netliquid_ex + (1-`theta')*housePos
    gen model_l = netliquid_ex
    gen model_h = housePos - min(model_a, 0)/(1-`theta')
    replace model_a = max(model_a, 1e-5)

    sum labIncWork, d
    replace labIncWork = labIncWork/`normed'
    pctile labPct = labIncWork [`weight' `exp'], n(4)
    foreach var of varlist model_a-model_h {
         replace `var' = `var'/`normed'
    }
    * Correlation between young income and assets for renters is
    * low, so initialize using the unconditional ECDF of the variables
    _pctile model_a if homeownershiprate == 0 [`weight' `exp'], per(1 90)
    replace model_a = . if homeownershiprate == 0 & (model_a <= `r(r1)' | model_a >= `r(r2)')
    pctile model_a_ecdf = model_a if homeownershiprate == 0 [`weight' `exp'], n(20)
    * generating separate income quantiles
    xtile labQtile = labIncWork [`weight' `exp'], n(4)
    save `quantiles'
    
    collapse (count) N=labIncWork (mean) model_a (sd) modelsd_h = model_h modelsd_a=model_a ///
        (p25) model_lo=model_h (p75) model_hi=model_h,  by(labQtile homeownershiprate)
    list
    
    use `quantiles', clear
    * Trim the home value distribution somewhat
    _pctile model_h if homeownershiprate == 1 [`weight' `exp'], per(2.5 97.5)
    replace model_h = . if homeownershiprate == 1 & (model_h <= `r(r1)' | model_h >= `r(r2)')
    * So three categories: "Below median" and top 2 quartiles
    replace labQtile = 2 if labQtile == 1
    matrix drop _all
    * Instead of using binned means, the initial state distributions are fit
    * using gamma distributions over different income bins (so like a mixture)
    forv i=2/4 {
    	svy: mean homeownershiprate if labQtile == `i' & age <= `hoyears'
	matrix l_coefs`i' = (labPct[`i'] \ e(b))
	matrix l_coefs`i'own = (labPct[`i'] \ e(b))
	sum model_a if homeownershiprate==1 & labQtile==`i' [aw`exp'], d
	sum model_h if homeownershiprate==1 & labQtile==`i' [aw`exp'], d
	gammafit model_a if labQtile==`i' & homeownershiprate==0 [aw`exp'], r
	matrix l_coefs`i' = l_coefs`i' \ e(alpha) \ e(beta)
	gammafit model_a if labQtile==`i' & homeownershiprate==1 [aw`exp'], r
	matrix l_coefs`i'own = l_coefs`i'own \e(alpha) \ e(beta)
	svmat double l_coefs`i'
	_pctile model_h if labQtile==`i' & homeownershiprate == 1 [`weight' `exp'], per(25 75)
	matrix l_coefs`i'own = l_coefs`i'own \ r(r1) \ r(r2)
	svmat double l_coefs`i'own

    }    
    
    * Export, convert so each row represents each variable
    keep l_coefs*
    keep if !missing(l_coefs2own)
    list in 1
    foreach var of varlist l_coefs21-l_coefs4own1 {
        replace `var' = 0.0 if missing(`var')
    }
    xpose, clear
    format * %-16.10f
    list

    export delimited using `initdir'/initial_grid_calibration.txt, ///
        delimiter(tab) novar datafmt replace

end
* %>

capture program drop make_BPP /* %< */
program define make_BPP
    syntax anything(name=initdir)
    /* Using SCF data, generate a smooth function of age effects in the
    income process (similar to KV 2014) */

    preserve
    tab lf fullEarn1
    keep if fullEarn1 == 1 // full-time employment
    gen educlvl = cond(educ < 12, 1, cond(educ==12, 2, cond(educ > 12 & educ < 16, 3, 4)))
    gen logLabInc = log(labIncWork)
    qui reg logLab i.year i.educlvl i.race i.occat1 i.kids i.married i.hhsex ///
        i.educlvl#i.year i.lf#i.year i.occat1#i.year i.race#i.year
    predict uy, r
    keep if age >= 22 & age < 60 // Retirement happens ON YEAR 60

    binscatter uy age, n(40) title("BEFORE taxes") name(raw, replace)
    sleep 3000
    mkspline ages = age, cubic nk(5)
    
    * Fit age effects on a spline of age.
    * Note the if condition: see later note
    reg uy ages? if age > 30
    predict logy_model_wo, xb

    svyset [pw=wgt]
    svy: mean uy
    matrix b = e(b)
    scalar means = b[1,1]
    display exp(means)

    collapse (first) logy_model (mean) uy (p75) uy_upper=uy (p25) uy_lower=uy, ///
        by(age)
    * Diagnostic: how well does the polynomial fit?
    twoway scatter uy logy_model_wo age, title("BEFORE taxes") $graphconfig ///
        legend(label(1 "Data Residual") label(2 "Fitted means")) name(fitted, replace)

    format logy_model_wo %13.10f
    keep logy_model_wo
    export delimited `initdir'/ageearnings.txt, novar dataf replace
    
end
* %>


* Initializing the income distribution.
capture program drop init_incomes
program define init_incomes
        syntax anything(name=initdir), maxgrid(real) normed(real) agefe(real)
        /* Args:
            maxgrid: the max value the income shocks grid supports, found by
                checking the model log file.
            normed: Mean HH income in sample over which we normalize incomes.
            agefe: Term netting out the cohort fixed effect in log income over the
            subsample of young adults.
        */
	
	preserve
	keep if age >= 22 & age <= 25 & fullEarn1 == 1
	sum labIncWork, d
	gen labIncFrac = labIncWork/(exp(`maxgrid')*exp(`agefe')*`normed')
	sum labIncFrac, d
	replace labIncFrac = . if labIncFrac > 1 | labIncFrac <= 0
	cumul labIncWork, gen(IncECDF)
	twoway (scatter IncECDF labIncFrac) 
	hist labIncFrac, name(init_incomeshist, replace)
	betafit labIncFrac [aw=wgt]
	
	local alphaval = e(alpha)
	local betaval = e(beta)
	!echo "`alphaval',`betaval'" > `initdir'/init_incomes.csv

end

capture program drop make_moments /* %< */
program define make_moments
    syntax anything(name=initdir), normed(real)
    /* Generate a one-row output of moments that we match the model to.
        Args:
            initdir: directory to which the file is outputted.
            normed: Mean HH income in sample over which we normalize incomes.
    */

    preserve
    matrix mom = J(9,1,.)

    * SCF data
    capture gen worth_earn_ratio = netWorthnc/labIncWork
    _pctile worth_earn_ratio [pw=wgt], n(2)
    matrix mom[1,1] = r(r1) // Median NW/Earnings ratio
    svy: mean homeownershiprate if age < 60
    matrix out = e(b)
    matrix mom[2,1] = out[1,1] // Average homeownership rate
    svy: mean homeownershiprate if age < 30
    matrix out = e(b)
    matrix mom[3,1] = out[1,1] // Homeownership rate, < 30 y.o.
    svy: mean homeownershiprate if age > 65 & age <= 75
    matrix out = e(b)
    matrix mom[4,1] = out[1,1] // Homeownership rate, > 65 y.o.
    
    * These are just recalled from IRS_clean.do
    matrix mom[5,1] = med_fthb     // Median FTHB age from IRS data
    matrix mom[6,1] = agerise      // (28-30)/(22-25) FTHB gradient from IRS data
    matrix mom[7,1] = agefall      // (35-33)/(28-30) FTHB gradient from IRS data
    matrix mom[8,1] = medinc_fthb/`normed'  // 50th pctile FTHB income from IRS data
    matrix mom[9,1] = inc75_fthb /`normed'  // 75th pctile FTHB income from IRS data

    * The *10 scales some moments that do not vary as much in absolute
    * value as others, which skews calibration toward certain moments.
    matrix weights = (1 \ 10 \ 10 \ 10 \ 1 \ 10 \ 10 \ 5 \ 5)
    matrix mom = mom, weights

    clear
    svmat mom
    xpose, clear
    export delimited `initdir'/moments.csv, novar replace

end /* %> */


************************************************************
* %< 2. Additional calibration moment programs.
************************************************************

capture program drop bin_moments /* %< */
program define bin_moments
    syntax varname [pweight], normed(real) datadir(string)

    local outvars labIncWork netliquid* housePos houseNeg networth_ex ///
        netWorthnc retAcc homeownershiprate
    preserve
    collapse `outvars' [`weight' `exp'], by(`varlist')
    order `outvars', after(`varlist')

    foreach var of varlist labIncWork-retAcc {
        replace `var' = `var'/`normed'
    }
    outsheet using "`datadir'/raw_data/../graph_data_by`varlist'.csv", ///
        comma names replace

    restore, preserve
    keep if homeownershiprate == 1 // owners 
    collapse `outvars' [`weight' `exp'], by(`varlist')
    order `outvars', after(`varlist')

    foreach var of varlist labIncWork-retAcc {
        replace `var' = `var'/`normed'
    }
    outsheet using "`datadir'/raw_data/../graph_data_by`varlist'_homeowners.csv", ///
        comma names replace

    restore, preserve
    keep if homeownershiprate == 0 // renters 
    collapse `outvars' [`weight' `exp'], by(`varlist')
    order `outvars', after(`varlist')

    foreach var of varlist labIncWork-retAcc {
        replace `var' = `var'/`normed'
    }
    outsheet using "`datadir'/raw_data/../graph_data_by`varlist'_renters.csv", ///
        comma names replace

end
* %>

capture program drop hw_ratios /* %< */
program define hw_ratios
    syntax anything(name=datadir)

    gen hw_ratio = housePos/netWorthnc if housePos > 0
    qui _pctile hw_ratio [pweight = wgt], p(5 95)
    replace hw_ratio = . if (hw_ratio>=r(r2) | hw_ratio <= r(r1))
    preserve

    tempfile all
    xtile hw_ratio_pct=hw_ratio [pw=wgt], n(100)
    collapse (max) hw_ratio, by(hw_ratio_pct)
    save `all', replace

    restore, preserve
    xtile hw_ratio_pct=hw_ratio if !missing(labIncWork) [pw=wgt], n(100)
    collapse (max) hw_ratio_work=hw_ratio, by(hw_ratio_pct)
    merge 1:1 hw_ratio_pct using `all'
    drop _m
    keep if !missing(hw_ratio)
    saveold `datadir'/raw_data/../hw_ratios, replace

end /* %> */

capture program drop var_distrib /* %< */
program define var_distrib
    syntax varlist [if], datadir(string) normed(real) append(string) [propvar(varlist)]

    preserve
    if "`if'" != "" keep `if'
    local collapse_query
    foreach pctile in 10 25 50 75 90 {
        local collapse_query `collapse_query' (p`pctile')
        foreach var in `varlist' `propvar' {
            local collapse_query `collapse_query' `var'p`pctile' = `var'
        }
    }
    collapse `collapse_query' [pw=wgt]
    gen t = "moments"
    reshape long `varlist' `propvar', i(t) j(type, string)
    foreach var in `varlist' {
        capture replace `var' = `var'/`normed'
    }
    list
    outsheet using "`datadir'/raw_data/../var_distrib_by`append'.csv", ///
        comma names replace

end /* %> */

capture program drop wealth_lorenz /* %< */
program define wealth_lorenz
    syntax anything(name=datadir)

    * Processes the SCF to build a dataset containing the cumulative distributions
    * of total housing wealth and income, which is then shown to produce
    * something like figure 4 in Diaz, Luongo-Prado (2010).

    preserve
    keep if age < 60
    tempfile orig
    gen labIncWork_h = labIncWork*homeownershiprate
    gen houses_h = houses*homeownershiprate
    gen netliquid_r = max(netliquid_ex*(1-homeownershiprate), 0.0)
    save `orig', replace
    
    local mergevars
    foreach var of varlist labIncWork* houses* netWorthnc netliquid_r {
        tempfile `var'_file
        keep `var' homeownershiprate
        drop if missing(`var')
	* Trimming the top (because our income process isn't exactly like reality)
        sort `var'
	_pctile `var', p(97.5)
	replace `var' = . if !missing(`var') & `var' >= `r(r1)'
	* Creating share of wealth variable
        sum `var'
        gen `var'Cum = sum(`var')/`r(sum)'*100
	* Creating percentile variable
        if regexm("`var'", "_h$") == 1 {
            egen Pct = rank(`var') if homeownershiprate == 1, u
        }
        else if regexm("`var'", "_r$") == 1 {
            egen Pct = rank(`var') if homeownershiprate == 0, u
        }
        else {
            egen Pct = rank(`var'), u
        }
        sum `var'Cum Pct
        drop if missing(Pct)
        save ``var'_file', replace
        use `orig', clear
	local mergevars `mergevars' `var'
    }
    sum homeownershiprate

    use `labIncWork_file', clear
    foreach var of local mergevars {
        if "`var'" != "labIncWork" {
        merge 1:1 Pct using ``var'_file', keep(1 2 3)
        drop _m
	}
    }
    sum
    foreach var of varlist labIncWork houses netWorthnc *_h *_r {
        sum Pct if !missing(`var')
        gen `var'Pct = Pct/`r(max)'*100
    }
    twoway (line housesCum housesPct) (line labIncWorkCum labIncWorkPct if !missing(labIncWork))
    sleep 3000
    twoway (line houses_hCum houses_hPct if !missing(houses_h)) ///
        (line labIncWork_hCum labIncWork_hPct if !missing(labIncWork_h))
    saveold `datadir'/raw_data/../lorenz_test, replace

end
* %>

* %>

/* Creating initialization and target moment files */
capture program drop main_SCF
program define main_SCF
    syntax, initdir(string) datadir(string)

* normalize by the labor income or working age people including renters and homeowners 	
* NOTE MEANS WEIGHED BY IPWS USED FROM HERE
svyset [pweight=wgt]
svy: mean labIncWork if age>=22 & age< 60
matrix out = e(b)
local incNormed = out[1,1]
_pctile netWorthnc [pw=wgt], p(5 95)
local netw1 = r(r1)
local netw2 = r(r2)
display "Normalizing household income factor: `incNormed'"

* SUMMARY STATISTICS (NOTE SVY COMMANDS USED)
* overall including renters and homeowners
svy: mean labIncWork netliquid* housePos houseNeg networth_ex netWorthnc

* conditional on owning a house
svy: mean labIncWork netliquid* housePos houseNeg networth_ex netWorthnc if houses > 0

/* Initialize deterministic parts of income process */
make_BPP `initdir'
initial_grid_calibration 22 25 [pw=wgt], initdir(`initdir') ///
    theta(0.20) normed(`incNormed') hoyears(24)
make_moments `initdir', normed(`incNormed')
init_incomes `initdir', normed(`incNormed') agefe(-0.30)

/* Creating other non-fitted calibration moments here */
qui bin_moments ageind [pw=wgt], datadir(`datadir') normed(`incNormed')
qui bin_moments ageyoung [pw=wgt], datadir(`datadir') normed(`incNormed')

* Distribution of assets across 
replace netliquid_exhousingdebt = netliquid_exhousingdebt + retAcc if age >= 60
/* hw_ratios */
var_distrib labIncWork netliquid netliquid_exhousingdebt housePos houseNeg ///
    netWorthnc W if netWorthnc >= `netw1' & netWorthnc <= `netw2' & ///
    homeownershiprate == 0 & age < 60, datadir(`datadir') normed(`incNormed') append(renters)
var_distrib labIncWork netliquid netliquid_exhousingdebt housePos houseNeg ///
    netWorthnc W if netWorthnc >= `netw1' & netWorthnc <= `netw2' & ///
    homeownershiprate == 1 & age < 60, datadir(`datadir') propvar(LTV LTVcond) ///
    normed(`incNormed') append(owners)
wealth_lorenz `datadir'

end

* Main function call
main_SCF, initdir(`initdir') datadir(`datadir')
cd `initdir'
