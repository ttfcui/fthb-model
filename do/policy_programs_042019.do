/***********************************************************
   Name: policy_programs_042019.do
   Purpose: Collection of programs used on model simulations
       for analysis purposes.
   Date: 2021/01

***********************************************************/

************************************************************
* %< 0. Data prep programs.
************************************************************

/* %< 0.1. DATA LOADING PROGRAMS */


capture program drop load_analysis_data /* %< */
program define load_analysis_data
    syntax anything(name=fname), datadir(string) [modelname(string)]
    /* Wrapper to load the specific dta files used in this code
       (note the dta files are generated from the simulations elsewhere
       Args of note:
           datadir: directory where the dta files are stored.
    */

    if "`fname'" == "stats" {
        use `datadir'/fthb_stats_final.dta, clear
        polstats_prep
    }
    else if "`fname'" == "microdata" {
        use `datadir'/masterdata_finalanalysis.dta, clear
        do `datadir'/masterdata.do
        if "`modelname'" != "" keep if model == "`modelname'"
        saveold `datadir'/masterdata_main.dta, replace
        microdata_prep
    }
    else if "`fname'" == "facevalue" {
        use `datadir'/masterdata_facevalues.dta, clear
        do `datadir'/masterdata.do
        microdata_prep
    }
    else if "`fname'" == "lives" {
        use `datadir'/fthb_trans_lives.dta, clear
        if "`modelname'" != "" keep if model == "`modelname'"
    }

end /* %> */


capture program drop indiv_assembly /* %< */
program define indiv_assembly
	syntax anything(name=mnames)
       /* Using the source propensity files for different model iterations,
          generate a panel of simulations at the individual-age cohort-model
          level.
          Args of note:
              mnames: list of model names as defined in the CSV names,
                  delineated by spaces.
       */

        tempfile init
	tokenize `mnames'
	local first `1'
	
        foreach model of local mnames {
	import delimited fthb/propensity_experiment_`model'.csv, clear
	destring leverage_pol?, ignore("infNa") replace
	gen fthbs = (inframarg == 1 | marginal == 1)
	replace inframarg = 0 if missing(inframarg)
	xtile income_bin = income_val, n(15)
	xtile assets_bin = assets, n(15)
	*
	gen loginc = log(income_val)
	gen logast = log(assets+sqrt(assets^2+1)) // ihs transformation approx. logs
	gen constrained = (assets==0)
    
        gen H_own_pol = nextdurables*fthbs
	gen H_rent_pol = nextdurables*(1-fthbs)
	local theta = 1-(1-0.20)*(1-0.0548)
        gen_stimulus consump H_own_pol H_rent_pol nextassets, usercost(0.0548) ///
	    theta(`theta') rentprem(0.00966) pricecoef(0.173)
        
	
	if "`model'" == "`first'" {
	    bys assets_bin: sum assets
	    bys income_bin: sum income_val
	    pctile inc_pct = loginc, n(10)
	    levelsof inc_pct, l(inc_pct)
	    drop inc_pct
	    save `init', replace
	}
	else {
	    append using `init'
	    save `init', replace
	}
        }

end /* %> */


capture program drop elas_assembly /* %< */
program define elas_assembly
	syntax anything(name=mnames), file1(string) file2(string)
       /* Using the source propensity files for different model iterations,
          Generated aggregated purchase propensities for heterogeneous agent
          bins before and during the policy.
          Args of note:
              mnames: list of model names as defined in the CSV names,
                  delineated by spaces.
           file1/file2: output file names for the disaggregated file by income
               tiers and in aggregate, respectively.
       */


	tempfile init init2
	tokenize `mnames'
	local first `1'

	foreach model of local mnames {
	import delimited fthb/propensity_experiment_`model'.csv, clear
	destring leverage_pol?, ignore("infNa") replace
	gen fthbs = (inframarg == 1 | marginal == 1)
	rename (owner_pol1 owner_pol9 leverage_pol1 leverage_pol3 consumption) ///
	    (ownP1 ownP9 levP1 levP3 cons)
	replace inframarg = 0 if missing(inframarg)
	
	gen_groupings
	gen age_bin = ceil(age/4)*4 - 2
        egen infra_groups = group(infra_inc infra_ass)
	* House size decisions over unvarying cohort (inframarginal buyers)
	gen ownVal = cond(inframarg == 1, nextdurables, .)
	preserve

	collapse (mean) fthbs ownP? ///
	    (min) infra_inc infra_ass (count) id, by(infra_groups)

	if "`model'" == "`first'" {
	    save `init', replace
	}
	else {
	    local modelfix = subinstr("`model'", ".", "", .)
	    rename (fthbs ownP1 ownP9) (fthbs_`modelfix' ownP1_`modelfix' ownP9_`modelfix')
	    merge 1:1 infra_groups using `init'
	    drop _m
	    save `init', replace
	}
	restore

	collapse (mean) fthbs ownP? ///
	    (min) infra_inc infra_ass (count) id
	gen agg = "agg"

	if "`model'" == "`first'" {
	    save `init2', replace
	}
	else {
	    local modelfix = subinstr("`model'", ".", "",.)
	    rename (fthbs ownP1 ownP9) (fthbs_`modelfix' ownP1_`modelfix' ownP9_`modelfix')
	    merge 1:1 agg using `init2'
	    drop _m
	    save `init2', replace
	}
	}
	
	use `init', clear
	saveold `file1', replace
	use `init2', clear
	saveold `file2', replace
	
end /* %> */


capture program drop gen_facevalues /* %< */
program define gen_facevalues

forv i=0.015(0.015)0.19 {
    local fnum = string(`i', "%05.4f")
    local nodown_models `nodown_models' FTHBSize_`fnum'
    local down_models `down_models' FTHBonDownSize_`fnum'
}
display "`nodown_models'"
display "`down_models'"

elas_assembly `nodown_models' monetary_little, ///
    file1(fthb/pol_FTHBSize_groups) file2(fthb/pol_FTHBSize_agg)
elas_assembly `down_models' monetary_little, ///
    file1(fthb/pol_FTHBOnDownSize_groups) file2(fthb/pol_FTHBOnDownSize_agg)

end /* %> */


/* %> */


/* %< 0.2. DATA CLEANING PROGRAMS */


capture program drop microdata_prep /* %< */
program define microdata_prep

    * Labels
    gen extensive = cond(!missing(posext), 1, cond(!missing(negext), 2, .))
    gen asset_filter = (assets <= 0.01)
    * Decomposing policy response
    gen decomp = cond(marginal==0, 1, cond( ///
        pullforward == 1 & adjDiff == 0, 2, cond( ///
        pullforward == 2 & adjDiff == 0, 3, cond( ///
        pullforward >= 3 & pullforward <= 5 & adjDiff == 0, 4, cond( ///
        pullforward > 5 & adjDiff == 0, 5, cond( ///
        adjustment_pol == 1 & adjustment_ss == 0, 6, cond( ///
        adjDiff > 0 & !(adjustment_pol == 1 & adjustment_ss == 0), 7, cond( ///
        adjDiff < 0, 8, .))))))))
    * Decomposing agents with timing + extensive margin
    gen decomp_ext = cond(missing(pullforward), 1, cond( ///
        pullforward == 1, 2, cond( ///
        pullforward == 2, 3, cond( ///
        pullforward >= 3 & pullforward <= 5, 4, cond( ///
        pullforward > 5, 5,.)))))
    * Indicators for coming in/out of income shocks
    gen posIncShock = (income > income_l)
    gen negIncShock = (income < income_l)
    gen posIncShock2 = (income < income_f1)
    gen negIncShock2 = (income > income_f1)
    egen groups = group(marginal extensive)

    * Graph labelling
    label def groups 1 "Infra/Positive" 2 "Infra/Negative" ///
        3 "Marginal/Positive" 4 "Marginal/Negative"
    label def decomp 1 "Inframarginal" 2 "Timing Margin, 1 period" ///
        3 "Timing Margin, 2 periods" 4 "Timing Margin, 3-5 periods" ///
        5 "Timing Margin, 5+ periods" 6 "Extensive Margin, 0-1" ///
        7 "Extensive Margin, Other +ve" 8 "Extensive Margin, -ve"
    label def assetSplit 0 "Assets > 1e-2" 1 "Assets <= 1e-2"
    label val decomp decomp
    label val decomp_ext decomp
    label val groups groups
    label val asset_filter assetSplit
    label var posIncShock "Proportion w/ Positive income shock, policy period"
    label var negIncShock "Proportion w/ Negative income shock, policy period"
    label var posIncShock2 "Proportion w/ Positive income shock, period after policy"
    label var negIncShock2 "Proportion w/ Negative income shock, period after policy"

end /* %> */


capture program drop gen_groupings /* %< */
program define gen_groupings

        count if inlist(model, "FTHB Tax Credit, $8K", "experiment_monetary_nodown")
	if `r(N)' == 0 {
	    preserve
	    qui indiv_assembly monetary_nodown
	    replace marginal = 0 if inframarg == 1
	    drop inframarg
	}
	_pctile income_val if marginal==0 & inlist(model, "FTHB Tax Credit, $8K", "experiment_monetary_nodown"), pe(1)
	local lev1 = `r(r1)'

	* All working-age income dist
	capture preserve
	import delimited "fthb/fthb_experiment_monetary_nodown.txt", delim(" ", collapse) clear
	rename (v3 v11 v16) (age assets income_val)
	sum income_val if age >= 1 & age <= 39, d
	local lev2 =  `r(p50)'
	display "1: `lev1';   2: `lev2'"

	gen infra_asset_tiers = (assets == 0) // because we don't want to count people with mortgages
	gen infra_incval_tiers = cond(income_val > `lev2', 3, cond(income_val > `lev1', 2, 1))
	tab infra_asset_tiers if age >= 1 & age <= 39
	tab infra_incval_tiers if age >= 1 & age <= 39

        restore
	gen infra_asset_tiers = (assets <= 0)
	gen infra_incval_tiers = cond(income_val > `lev2', 3, cond(income_val > `lev1', 2, 1))
	
end /* %> */


capture program drop relabel_models /* %< */
program define relabel_models

	replace model = "FTHB Tax Credit, $8K" if model=="experiment_monetary_nodown"
	replace model = "CARS Program" if model=="experiment_CARS"
	replace model = "CARS Program, No Borrowing" if model=="experiment_CARS_nocoll"
	replace model = "FTHB Bridge Loan, $5K" if model=="experiment_monetary"
	replace model = "FTHB Bridge Loan, $5K w/out Forced Depreciation" if model=="experiment_monetary_dep"

end /* %> */


/* %> */


/* %<  0.3. CALCULATION PROGRAMS */


capture program drop gen_stimulus /* %< */
program define gen_stimulus
    syntax varlist(max=4), ///
        usercost(real) theta(real) rentprem(real) pricecoef(real) [hetero(string)]
    
    tokenize `varlist'
    local consvar `1'
    local ownvar `2'
    local rentvar `3'
    local assetvar `4'
    
    display "Consumption Variable: `consvar'"
    display "Owned House Variable: `ownvar'"
    display "Rental House Variable: `rentvar'"
    display "Assets Variable: `assetvar'"
    
    *exp(*) term is price. h_own coefficient is down payment percentage.
    *H_rent coefficient is rental price (user cost + premium)    
    gen expend`hetero' = `consvar' + exp(`pricecoef')*(`ownvar'*`theta' + `rentvar'*(`usercost'+`rentprem'))
    gen expend_share`hetero' = exp(`pricecoef')*(`ownvar'*`theta' + `rentvar'*(`usercost'+`rentprem'))/expend`hetero'
    gen durexpend`hetero' = max(`ownvar'*`theta', `rentvar'*(`usercost'+`rentprem'))
    gen durlevels`hetero' = max(`ownvar', `rentvar')
    capture gen nxtast`hetero' = cond(fthbs==1, `assetvar' + exp(`pricecoef')*(1-`theta')*`ownvar', .)

end /* %> */


capture program drop stimulus_prep /* %< */
program define stimulus_prep

local usercost = 0.0548
local theta = 1-(1-0.20)*(1-`usercost')
local rentprem = 0.0234
local cardep = 0.125
local housedep = 0.022
local pricecoef = -0.0696

local prepvars pol0
forv t=1/5 {
    local prepvars `prepvars' pol`t' pol`t'_cf
}
foreach period of local prepvars {

    if regexm("`period'", "_cf") == 0 {
        * Correction: counterfactual counts level before dep, policy sim doesn't
	* So H_own divided by (1-dep. rate)
        replace H_own_`period' = H_own_`period'/(1-`cardep') if regexm(model, "CARS") == 1
    }
    gen_stimulus C_`period' H_own_`period' H_rent_`period', usercost(`usercost') ///
        theta(`theta') rentprem(`rentprem') pricecoef(`pricecoef') hetero(_`period')
}
rename (nextDurables_ss) (durlevels_pol0_cf)
gen H_own_pol0_cf = durlevels_pol0_cf*(1-marginal)
gen H_rent_pol0_cf = durlevels_pol0_cf*marginal
gen durexpend_pol0_cf = max(H_own_pol0_cf*`theta', H_rent_pol0_cf*(`usercost'+`rentprem'))

foreach i of numlist 0/5 {
    *investment coefficient is depreciation rate (TODO separate rate for Cars)
    gen durables_diff`i' = durexpend_pol`i' - durexpend_pol`i'_cf
    gen durlevels_diff`i' = durlevels_pol`i' - durlevels_pol`i'_cf
    if `i' > 0 {
        local j = `i' - 1
        gen investment_diff`i' = `housedep'*H_own_pol`j' + (durexpend_pol`i' - durexpend_pol`j') ///
	- (durexpend_pol`i'_cf + `housedep'*H_own_pol`j'_cf - durexpend_pol`j'_cf)
        * CARS depreciates so cannot add maintenace costs in benefit
	replace investment_diff`i' = investment_diff`i' - (`housedep'*H_own_pol`j' - `housedep'*H_own_pol`j'_cf) ///
	    if regexm(model, "CARS") == 1
    }
}

rename (welfare) (V_pol0_cf)
foreach i of numlist 0/5 {
    * Welfare measurements
    gen welfare_diff`i' = -(V_pol`i' - V_pol`i'_cf)/V_pol0_cf*100 // pctage terms
}

end /* %> */


capture program drop cost_benefit_durables /* %< */
program define cost_benefit_durables
    syntax anything [if/], discount(real) [noNORMalize]
    
    tokenize `anything'
    local basevar `1'
    local diffvars `2'
    
	preserve
	local ifquery "if fthbs==1"
	if `"`if'"' != "" local ifquery `"`ifquery' & `if'"'
	collapse (mean) `basevar' `diffvars'* (count) id `ifquery', by(model)
	forv i=1/5 {
	    gen `diffvars'`i'disc = `discount'*`diffvars'`i'
	    order `diffvars'`i'disc, before(`diffvars'1)
	    egen npv_`i' = rowtotal(`basevar'-`diffvars'`i'disc)
	    if "`normalize'" == "" {
	    replace npv_`i' = npv_`i'/0.12 if inlist(model, "FTHB Tax Credit, $8K")
	    replace npv_`i' = npv_`i'/0.075 if inlist(model, "FTHB Bridge Loan, $5K")
	    replace npv_`i' = npv_`i'/0.12*(id/2.4e5) if inlist(model, "Unconditional $8K Transfer")
	    replace npv_`i' = npv_`i'/0.06 if inlist(model, "CARS Program")
	    replace npv_`i' = npv_`i'/0.03 if regexm(model, "2K") == 1
	    }
	}
	
	if "`normalize'" == "" {
	replace `basevar' = `basevar'/0.12 if inlist(model, "FTHB Tax Credit, $8K")
	replace `basevar' = `basevar'/0.075 if inlist(model, "FTHB Bridge Loan, $5K")
	replace `basevar' = `basevar'/0.12*(id/2.4e5) if inlist(model, "Unconditional $8K Transfer")
	replace `basevar' = `basevar'/0.06 if inlist(model, "CARS Program")
	replace `basevar' = `basevar'/0.03 if regexm(model, "2K") == 1
	}

	list model `basevar' npv_* id

end /* %> */


capture program drop impute_effects /* %< */
program define impute_effects
    syntax anything(name=binname), nsize(integer) nbins(integer)
    local obs = `nsize'*`nbins'
    set obs `obs'
    egen noDownLevel = seq(), from(1) to(`nsize') block(`nbins')
    egen `binname' = seq(), from(1) to(`nbins') block(1)

    forval i=0/8 {
    if `i' == 0 {
        gen effectnodown = .
	gen effectondown = .
    }
    else {
        gen effectless`i'ondown = .
	label var effectless`i'ondown "Net of Effects of \$`i'K on DP"
    }
    forval j=1/`nsize' {
        forval k=1/`nbins' {
	    local index = (`j'-1)*`nbins' + `k'
	    if `i' == 0 {
	        qui replace effectnodown = fthbs`j'nodown[`k'] in `index'
		qui replace effectondown = fthbs`j'ondown[`k'] in `index'
	    }
	    else {
	        qui replace effectless`i'ondown = ///
	            100*(fthbs`j'nodown[`k'] - fthbs`i'ondown[`k']) in `index'
            }
        }
    }
    }
end /* %> */


/* %> */


********** %>


************************************************************
* %< 1. Data Analysis Programs.
************************************************************


/* %<  1.1. SIMPLE TABLE+FIGURE PROGRAMS */


capture program drop decomp_table /* %< */
program define decomp_table
    syntax anything(name=fname) [if/], tabvar(varname) [PCT]
    
    if `"`if'"' != "" local ifcond `"& `if'"'
    if "`pct'" != "" local cval cells(b(fmt(4)) rowpct(fmt(2)))
    
    estpost tab decomp `tabvar' [iw=decomped] if decomp != 1 `ifcond'
    esttab using `fname', `cval' unstack label noobs nonum tex replace
    
end /* %> */


capture program drop pullf_dist /* %< */
program define pullf_dist
    syntax varname [if], compgroup(numlist max=2) label1(string) label2(string) name(string)
    
    tokenize `compgroup'
    local redgrp `1'
    local graygrp `2'
    
twoway (hist pullforward if `varlist'==`redgrp' & pullforward < 40, width(1) fc(none) lc(red)) ///
    (hist pullforward if `varlist'==`graygrp' & pullforward < 40 , width(1) fc(none) lc(gs10) lp(dash_dot)) ///
    `if', xsize(7) ///
    legend(label(1 "`label1'") label(2 "`label2'")) by(model, $graphopts note("")) name("`name'", replace)
    graph export fthb/pullf_`name'.pdf, as(pdf) replace

end /* %> */


/* %> */


/* %<  1.2. DATA MERGING PROGRAMS */


capture program drop dataprep_takers /* %< */
program define dataprep_takers
    /* Relabel and clean dataset with stacked simulated data across simulations. */

    load_analysis_data microdata, datadir(fthb)
    * Where are these people coming from? Seem like HH that were not
    * policy takers but also changed their extensive margin
    drop if model=="" | missing(nextDurables)

    relabel_models
    
    preserve
    tempfile denoms
    count if model=="FTHB Tax Credit, $8K" & marginal == 0
    local denom = `r(N)'
    keep if decomp == 1
    collapse (count) id, by(model decomp)
    rename id denoms
    keep model denoms
    replace denoms = `denom' if regexm(model, "FTHB") == 1 & denoms != `denom'
    save `denoms'

    restore
    merge n:1 model using `denoms'
    drop _m
    
end /* %> */


capture program drop dataprep_renters /* %< */
program define dataprep_renters
        /* Merging the propensity file with policy period decisions on
           the fthb file with counterfactual statistics. */

	tempfile master_main
	save `master_main'

	indiv_assembly monetary monetary_nodown monetary_dep CARS CARS_nocoll
	relabel_models
	replace marginal = 0 if inframarg == 1
	drop inframarg
	merge 1:1 id age model using `master_main', update

	gen_groupings
	tab model infra_asset_tiers, row
	tab model infra_incval_tiers, row
	tab model infra_incval_tiers if marginal==0, row
	egen infra_groups = group(infra_inc infra_ass)

	label def infra_groups 1 "Poorest, Positive Assets" 3 "Poor, Positive Assets" ///
	    4 "Poor, No Assets" 5 "Non-poor, Positive Assets" ///
	    6 "Non-poor HH, No Assets"
	label val infra_groups infra_groups
	tab infra_groups if model=="FTHB Tax Credit, $8K"
	tab infra_groups if model=="FTHB Tax Credit, $8K" & marginal==0

end /* %> */


capture program drop gen_faceveffects /* %< */
program define gen_faceveffects
        /* Clean dataset collecting aggregate effects of policy over
           various face value sizes. */

	tempfile orig
	use fthb/pol_FTHBSize_groups, clear
	rename (fthbs ownP?) (fthbs_FTHBSize_0015 ownP?_FTHBSize_0015)
	reshape long fthbs ownP1 ownP9, i(infra_groups) j(model) string
	encode model, gen(Model)
	save `orig'

	use fthb/pol_FTHBOnDownSize_groups, clear
	rename (fthbs ownP?) (fthbs_FTHBOnDownSize_0015 ownP?_FTHBOnDownSize_0015)
	reshape long fthbs ownP1 ownP9, i(infra_groups) j(model) string
	encode model, gen(Model)
	gen modelType = "ondown"
	append using `orig'
	replace modelType = "nodown" if missing(modelType)
	drop if missing(infra_groups)
	drop model

	drop ownP?
	reshape wide fthbs, i(infra_groups modelType) j(Model)
	reshape wide fthbs*, i(infra_groups) j(modelType, string)


	impute_effects infraGroups, nsize(12) nbins(6)
	saveold fthb/effectlessondown_table, replace

end /* %> */


/* %> */


********** %>


************************************************************
* %< 2. Table programs.
************************************************************


capture program drop timing_tables /* %< */
program define timing_tables
    syntax [if/], tabvar(varname)
    /* Replicates Table XX which decomposes the timing margins
       for all marginal buyers with FTHB v CARS */
    
    preserve
    if `"`if'"' != "" local ifcond `"& `if'"'

    keep if decomp <= 5 `ifcond'
    replace decomp = 2 if decomp == 3 // join "2 period" with "1 period"
    replace decomp = 5 if decomp == 4 // join "3-5 period" with "5+ periods"
    collapse (count) id (first) denoms, by(`tabvar' decomp)
    bys `tabvar': gen decomped = id/denoms

    label def decomp_tab3 2 "1-2 Periods" 5 "3+ Periods"
    label val decomp decomp_tab3
    
    if "`tabvar'" == "model" { 
    
    decomp_table fthb/fthbvcars_timing.tex ///
        if inlist(model, "FTHB Tax Credit, $8K", "CARS Program"), tabvar(`tabvar')

    decomp_table fthb/fthbalt_timing.tex ///
        if regexm(model, "CARS") == 0, tabvar(`tabvar')

    decomp_table fthb/carsalt_timing.tex ///
        if regexm(model, "CARS") == 1, tabvar(`tabvar')
	
    }
    else {
        decomp_table fthb/fthbhetero_timing.tex, tabvar(`tabvar') pct
    }

end /* %> */


capture program drop extensive_tables /* %< */
program define extensive_tables
    syntax [if/], tabvar(varname) fthbtrans(real) carstrans(real)
    /* Replicates Table XX which decomposes the extensive margins for HH
       with lifetime purchases different from counterfactual */

    preserve
    if `"`if'"' != "" local ifcond `"& `if'"'
    
    keep if (decomp > 5 | decomp == 1) `ifcond'
    * Reallocating inframarg people into ext. margin categories
    replace decomp = 6 if decomp == 1 & adjustment_pol == 1 & adjustment_ss == 0
    replace decomp = 7 if decomp == 1 & adjDiff > 0 & !(adjustment_pol == 1 & adjustment_ss == 0)
    replace decomp = 8 if decomp == 1 & adjDiff < 0  
    gen extDiff = adjDiff/`carstrans' if regexm(model, "CARS") == 1
    replace extDiff = adjDiff/`fthbtrans' if regexm(model, "CARS") == 0    
    
    collapse (count) id (sum) adjDiff extDiff (first) denoms if decomp != 1, by(`tabvar' decomp)
    bys `tabvar': egen extTot = total(extDiff)
    by `tabvar': gen decomped = id/denoms

    label def decomp_tab3 6 "0-1 Lifetime Transactions" ///
        7 "Other Positive Margin" 8 "Negative Margin", replace
    label val decomp decomp_tab3
    
    if "`tabvar'" == "model" { 
    
    decomp_table fthb/fthbvcars_extensive.tex ///
        if inlist(model, "FTHB Tax Credit, $8K", "CARS Program"), tabvar(`tabvar')

    decomp_table fthb/fthbalt_extensive.tex ///
        if regexm(model, "CARS") == 0, tabvar(`tabvar')

    decomp_table fthb/carsalt_extensive.tex ///
        if regexm(model, "CARS") == 1, tabvar(`tabvar')  
    }
     else {
        decomp_table fthb/fthbhetero_extensive.tex, tabvar(`tabvar') pct
    }   
    
    collapse (first) extTot, by(`tabvar')
    list

end /* %> */


capture program drop cost_benefit_table /* %< */
program cost_benefit_table
        /* Bang-for-buck estimates for policy based on stimulated spending/cost */

	cost_benefit_durables durables_diff0 investment_diff, discount(1)
	local discount = 1/1.03
	cost_benefit_durables durables_diff0 investment_diff, discount(`discount')
	local discount = 1/1.1
	cost_benefit_durables durables_diff0 investment_diff, discount(`discount')

	cost_benefit_durables durables_diff0 investment_diff if infra_groups <= 4, discount(1)
	local discount = 1/1.03
	cost_benefit_durables durables_diff0 investment_diff if infra_groups <= 4, discount(`discount')
	local discount = 1/1.1
	cost_benefit_durables durables_diff0 investment_diff if infra_groups <= 4, discount(`discount')
end /* %> */


capture program drop welfare_table /* %< */
program welfare_table
        /* Pctage change in welfare estimates based on counterfactual utility estimates */

        rename welfare_diff0 weldiff_base
	cost_benefit_durables weldiff_base welfare_diff, discount(1) nonorm
	local discount = 1/1.03
	cost_benefit_durables weldiff_base welfare_diff, discount(`discount') nonorm
	local discount = 1/1.1
	cost_benefit_durables weldiff_base welfare_diff, discount(`discount') nonorm

	cost_benefit_durables weldiff_base welfare_diff if infra_groups <= 4, discount(1) nonorm
	local discount = 1/1.03
	cost_benefit_durables weldiff_base welfare_diff if infra_groups <= 4, discount(`discount') nonorm
	local discount = 1/1.1
	cost_benefit_durables weldiff_base welfare_diff if infra_groups <= 4, discount(`discount') nonorm
end /* %> */


capture program drop hetero_tables /* %< */
program define hetero_tables
* Timing and Extensive Tables By Group

timing_tables if model=="FTHB Bridge Loan, $5K", tabvar(infra_groups)
extensive_tables if model=="FTHB Bridge Loan, $5K", tabvar(infra_groups) fthbtrans(14671) carstrans(29699)
!mv fthb/fthbhetero_timing.tex fthb/fthbhetero_timing_ondown.tex
!mv fthb/fthbhetero_extensive.tex fthb/fthbhetero_extensive_ondown.tex
timing_tables if model=="FTHB Tax Credit, $8K", tabvar(infra_groups)
extensive_tables if model=="FTHB Tax Credit, $8K", tabvar(infra_groups) fthbtrans(14671) carstrans(29699)

end /* %> */


capture program drop extensive_sumstats /* %< */
program define extensive_sumstats
* Extensive margin means
        preserve
	gen income_valdiff = income_val - income_val_l
	forv i=1/4 {
	    gen income_difff`i' = income_f`i' - income
	}

	collapse (mean) age income_val income*diff income_difff? durlevels_diff1 durlevels_diff3 (count) id, ///
	    by(adjDiff model)
	sort model adjDiff
	list model adjDiff age income* id
end /* %> */


capture program drop wtp_calcs /* %< */
program define wtp_calcs
        preserve
	* Down Payment WTP calculations
	foreach val in 1 2 4 {
	restore, preserve
	sort infraGroups noDownLevel
	keep if (effectless`val'ondown[_n] >= 0 & effectless`val'ondown[_n-1] < 0) | ///
	    (effectless`val'ondown[_n] < 0 & effectless`val'ondown[_n+1] >= 0)
	bys infraGroups: gen intercept = noDownLevel[_n-1] + ///
	    abs(effectless`val'ondown[_n-1])/(abs(effectless`val'ondown[_n-1]) + effectless`val'ondown[_n])
	gen wtp = (intercept - `val')*1e3
	list infraGroups noDownLevel effectless`val'ondown intercept wtp
	}
end /* %> */


********** %>

************************************************************
* %< 3. Figure programs.
************************************************************


capture program drop fig_barplots /* %< */
program define fig_barplots
        syntax [if]
* Bar plot of policy group decompositions
	preserve
	gen infra_count = 1/denoms
	capture keep `if' 
	
	collapse (count) id (first) infra*tiers if marginal == 1 [iw=infra_count], ///
	    by(model infra_groups)
	reshape wide id infra_ass infra_inc, i(model) j(infra_groups)
	list model id?
	gen poor_saver = id1 + id3
	capture gen id2 = .  // No poorest/no assets stimulated
	replace id2 = . // Can there be poorest/no assets? Maybe but small.
	replace id1=0 if missing(id1)
	gen id1min = 0
	local rbar_call "rbar id1min id1 model_lab, horiz barw(0.67)"
	forval i =2/6 {
	    replace id`i'=0 if missing(id`i')
	    local j = `i' - 1
	    gen id`i'min = id`j'
	    replace id`i' = id`i' + id`j'
	    if `i' != 2 local rbar_call "`rbar_call' || rbar id`i'min id`i' model_lab, horiz barw(0.67)"
	}
	tab model

	encode model, gen(model_lab)
	label def model_lab 1 "CARS Program" 2 "$5K Bridge Loan" 3 "$8K FTHBTC", replace
	display "`rbar_call'"
	list
	#delimit ;
	twoway `rbar_call' ||, ysc(reverse) ylab(,angle(0) val labs(small)) yti("") 
	    ylab(1(1)3) xti("Ratio of HH to Inframarginal Policy Takers") $graphopts
	    legend(col(3) order(1 2 3 4 5 -) 
	    label(1 "Poorest, Positive {it:a}{sub:-1}") label(2 "Poor, Positive {it:a}{sub:-1}")
	    label(3 "Poor, No {it:a}{sub:-1}") label(4 "Non-poor, Positive {it:a}{sub:-1}")
	    label(5 "Non-poor, No {it:a}{sub:-1}") label(6 "Poor+Poorest, Positive {it:a}{sub:-1}")
	    rowgap(*1.75) symx(*0.75) size(small) region(c(none))) xsize(6.5);
	#delimit cr

	graph export fthb/bargrp_plot.pdf, as(pdf) replace

end /* %> */


capture program drop fig_propplots /* %< */
program define fig_propplots
    syntax [if]
  /* Propensity Plots */
  
    tempfile modelprob main
    bys infra_groups model: egen renter_count = count(id)
    replace renter_count = 1/renter_count
    save `main'
    
    capture keep `if' 
    collapse (count) id (first) infra*tiers if fthbs == 1 [iw=renter_count], ///
        by(model infra_groups)
    save `modelprob'
    use `main', clear
    keep if model=="FTHB Bridge Loan, $5K"
    replace model=" No Policy"
    collapse (count) id (first) infra*tiers if marginal == 0 [iw=renter_count], ///
        by(model infra_groups)
    append using `modelprob'

    encode model, gen(modelv)
    drop model
    reshape wide id, i(infra_*) j(modelv)
    #delimit ;
    graph bar id*, over(infra_groups,  ///
       relabel(1 `""Poorest" "Positive {it:a}{sub:-1}""' 2 `""Poor" "Positive {it:a}{sub:-1}""'
       3 `""Poor" "No {it:a}{sub:-1}""' 4 `""Non-poor" "Positive {it:a}{sub:-1}""' 
       5 `""Non-poor" "No {it:a}{sub:-1}""'))
       legend(label(1 "No Policy") label(2 "$5K Bridge Loan")
       label(3 "$8K FTHBTC") label(4 "$8K Uncond. Transfer") region(c(none)))
       yti("Propensity of Purchase Within Group", ali(top) m(small)) xsize(6.5) $graphopts;
    #delimit cr
    graph export fthb/propgrp_plot.pdf, as(pdf) replace
    use `main', clear
    drop renter_count
    
end /* %> */


capture program drop pullf_distall /* %< */
program define pullf_distall
/* Pullforward distribution degeneracy */

    preserve
    keep if pullforward < 40 & adjDiff == 0 
    replace pullforward = 15 if pullforward > 15 & !missing(pullforward)
    
    pullf_dist infra_incval_tiers if inlist(model, "FTHB Bridge Loan, $5K", "FTHB Tax Credit, $8K"), ///
        compgroup(2 3) label1("Poor HH") label2("Non-Poor HH") name(hist_1)

    pullf_dist infra_incval_tiers if inlist(model, "CARS Program", "FTHB Tax Credit, $8K"), ///
        compgroup(2 3) label1("Poor HH") label2("Non-Poor HH") name(hist_1cars)

    pullf_dist infra_incval_tiers if inlist(model, "FTHB Bridge Loan, $5K", "FTHB Tax Credit, $8K"), ///
        compgroup(2 1) label1("Poor HH") label2("Poorest HH") name(hist_2)

    pullf_dist infra_asset_tiers if inlist(model, "FTHB Bridge Loan, $5K", "FTHB Tax Credit, $8K"), ///
        compgroup(0 1) label1("Positive Assets HH") label2("No Assets HH") name(hist_3)
	
end /* %> */


capture program drop extensive_binscatter /* %< */
program define extensive_binscatter
/* Extensive margin binscatters */

preserve
local legendlab `"xtitle("Income Level (000s)") ytitle("Change in Durable Levels," "3 Years Out of Policy (000s)") legend(label(1 "No Extensive") label(2 "Positive Extensive") label(3 "Negative Extensive"))"'
gen adjdummy = cond(adjDiff==0, 1, cond(adjDiff > 0, 2, 3)) if !missing(adjDiff)
* REQUIRES EGENMORE
bys model adjdummy: egen adjbins = xtile(income_val), n(11)
foreach var in income_val durlevels_diff3 {
    replace `var' = `var'*67.250 // unnormalizing
}
local thres1 = .4188*67.250 // inframarginal tail
local thres2 = .7704*67.250 // median
* The income_val restriction is like a winsorization routine
binscatter durlevels_diff3 income_val if model=="FTHB Tax Credit, $8K" & income_val < 110, ///
    by(adjdummy) xline(`thres1' `thres2') xsc(r(`thres1' `thres2')) rd(`thres2') xq(adjbins) `legendlab' name(fthb_bs, replace)
binscatter durlevels_diff3 income_val if model=="CARS Program", ///
    by(adjdummy) xline(`thres1' `thres2') xsc(r(`thres1' `thres2')) rd(`thres2') xq(adjbins) `legendlab' name(cars_bs, replace)
binscatter durlevels_diff3 income_val if model=="FTHB Bridge Loan, $5K" & income_val < 110, ///
    by(adjdummy) xline(`thres1' `thres2') xsc(r(`thres1' `thres2')) rd(`thres2') xq(adjbins) `legendlab' name(fthbbl_bs, replace)

end /* %> */


capture program drop ondown_plots /* %< */
program define ondown_plots
/* Down Payment Trend Plots */

	preserve
	drop effectless1ondown effectless3ondown effectless5on-effectless7on
	* Graph test
	label def infra_groups 1 "Poorest, Positive Assets" 3 "Poor, Positive Assets" ///
	    4 "Poor, No Assets" 5 "Non-poor, Positive Assets" ///
	    6 "Non-poor HH, No Assets"
	label val infraGroups infra_groups
	label var noDownLevel "Face Value of Subsidy (000s)"
	twoway (line effectless?ondown noDownLevel if infraGroups == 1, yline(0, lc(gs9)))
	twoway (line effectless?ondown noDownLevel, yline(0, lc(gs9)) xtick(1/12) xlab(1(3)11)) ///
	    if infraGroups != 2, by(infraGroups, xr $graphopts ) ///
	    legend(size(small)) yti("Net Chg in Purchase Propensity, pp") xsize(7)
	graph export "fthb/ondown_wtp.pdf", as(pdf) replace
	twoway (line effect??down noDownLevel, xtick(1/12) xlab(1(3)11)) ///
	    if infraGroups != 2, by(infraGroups, xr $graphopts note("")) ///
	    legend(label(1 "Inapplicable on DP") label(2 "Applicable on DP") size(small)) xsize(7)
	graph export "fthb/novsondow.pdf", as(pdf) replace

end /* %> */


********** %>


************************************************************
* %< MAIN PROGRAM (And Misc.)
************************************************************

capture program drop main
program define main

capture log c
log using fthb/policy_redo.log, replace
global graphopts graphregion(c(white))

dataprep_takers
* Timing/Extensive Crosstabs
timing_tables, tabvar(model)
extensive_tables, tabvar(model) fthbtrans(15082) carstrans(29869)

dataprep_renters
fig_barplots if !inlist(model, "FTHB Bridge Loan, $5K w/out Forced Depreciation", "CARS Program, No Borrowing")
fig_propplots if inlist(model, "FTHB Tax Credit, $8K", "FTHB Bridge Loan, $5K")
stimulus_prep

pullf_distall
hetero_tables
cost_benefit_table
exit
welfare_table
extensive_sumstats
extensive_binscatter

gen_facevalues
gen_faceveffects
ondown_plots
wtp_calcs

end

main


capture program drop cost_benefit_durables
program define cost_benefit_durables
    syntax [if], discount(real)
/* Bang-for-buck calculations over face value */
    
	preserve
	collapse (mean) durables_diff0 investment_diff* (count) id `if', by(model)
	gen denom = real(subinstr(substr(model, -5, .), "_", "", 1))
	forv i=1/5 {
	    gen investment_diff`i'disc = `discount'*investment_diff`i'
	    order investment_diff`i'disc, before(investment_diff1)
	    egen npv_`i' = rowtotal(durables_diff0-investment_diff`i'disc)
	    replace npv_`i' = npv_`i'/denom
	}

	replace durables_diff0 = durables_diff0/denom
	list model durables_diff0 npv_* id
	
	replace denom = denom/.015
	twoway (line npv_5 denom if regexm(model, "FTHBSize") == 1) ///
	    (line npv_5 denom if regexm(model, "FTHBonDown") == 1)

end

capture program drop cost_benefit_facevalues
program define cost_benefit_facevalues

load_analysis_data facevalue, datadir(fthb)
gen_groupings
egen infra_groups = group(infra_inc infra_ass)
stimulus_prep

cost_benefit_durables, discount(1)
local discount = 1/1.03
cost_benefit_durables, discount(`discount')
local discount = 1/1.1
cost_benefit_durables, discount(`discount')
end

cost_benefit_facevalues

/* %> */
