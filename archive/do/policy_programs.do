
************************************************************
* %< 0. Data prep programs.
************************************************************
capture program drop microdata_prep /* %< */
program define microdata_prep

    * Labels
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
        5 "Timing Margin, 5+ periods" 6 "0-1 Extensive Margin" ///
        7 "Other Pos. Extensive Margin" 8 "Neg. Extensive Margin"
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

capture program drop polstats_prep /* %< */
program define polstats_prep

    label def agebin 0 "All agents" 20 "Ages 22-29" 30 "Ages 30-39" 40 "Ages 40-49" ///
        50 "Ages 50-59" 60 "Ages 60-69" 70 "Ages 70-79"
    label val agebin agebin

    gen comp_stats = cond( ///
        regexm(desc, "1 period fwd, total") == 1 & subtype == "count", 1, cond( ///
        regexm(desc, "2 periods fwd, total") == 1 & subtype == "count", 2, cond( ///
        regexm(desc, "2\+ periods fwd, total") == 1 & subtype == "count", 3, cond( ///
        regexm(desc, "More .+ over lifecycle$") == 1, 4, cond( ///
        regexm(desc, "Fewer .+ over lifecycle$") == 1, 5, cond( ///
        desc == "Buyer characteristics" & subtype == "id", 99, .))))))

    label var value "FTHB Credit, delivered period following purchase"
    label var valuemonetary "FTHB Credit, applicable on down payment"
    label var valuemonetary_dep "FTHB Credit w/ mandatory depreciation"
    label var valuemonetary_dep_nodown "FTHB Credit w/ mandatory depreciation, applicable on down payment"
    label var valueCARS "CARS Repeat Buyer policy on down payment"
    label var valueCARS_nodown "CARS Repeat Buyer policy on whole budget"
    label var valueCARS_nocoll "CARS Repeat Buyer policy w/ no leverage"
    label var valueCARS_altScrap "CARS Repeat Buyer policy w/ zero scrap value" 
    label var valueFTHBTP_2 "FTHB Credit, two periods"
    label var valueHIMD_1_0 "Mortgage Interest Deduction Steady State"
    label var valueFTHBDep_0_03 "FTHB Credit w/ higher depreciation"
    label var valueFTHBSize_0_03 "FTHB Credit w/ \$2,000 face value"
    label var valueFTHBSize_0_165 "FTHB Credit w/ \$11,000 face value"
    label var valueFTHBPerm_99_0 "Permament FTHB Credit, inapplicable on down"
    label var valueFTHBPermonDown_99_0 "Permanent FTHB Credit, applicable on down"

end /* %> */

capture program drop trans_labels /* %< */
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
end /* %> */

capture program drop load_analysis_data /* %< */
program define load_analysis_data
    syntax anything(name=fname), datadir(string) [modelname(string)]

    if "`fname'" == "stats" {
        use `datadir'/fthb_stats_final.dta, clear
        polstats_prep
    }
    else if "`fname'" == "microdata" {
        use `datadir'/masterdata_final.dta, clear
        do $repodir/do/masterdata.do
        if "`modelname'" != "" keep if model == "`modelname'"
        saveold `datadir'/masterdata_main.dta, replace
        microdata_prep
    }
    else if "`fname'" == "lives" {
        use `datadir'/fthb_trans_lives.dta, clear
        if "`modelname'" != "" keep if model == "`modelname'"
    }
    else if "`fname'" == "comp" {
        use `datadir'/fthb_`modelname'comp_final.dta, clear
    }
    else if "`fname'" == "Elas" {
        use `datadir'/fthb_`modelname'Elas_final.dta, clear
    }

end /* %> */

* %>

************************************************************
* %< 1. Table programs.
************************************************************
capture program drop table_aggdecomp /* %< */
program define table_aggdecomp
    syntax anything(name=outputdir), [hetero(string)]
    * TODO: Automatically generate extra rows, syntax in final table
    * (like summing up pure timing margin)

    preserve
    collapse (count) id (sum) nextDurables ss=nextDurables_ss, by(decomp)
    scalar inv_ss = ss[1]
    gen idprop = id/id[1]
    gen durprop = nextD/inv_ss
    eststo counts: estpost tab decomp [iw=idprop]
    eststo inv: estpost tab decomp [iw=durprop]
    esttab counts inv using `outputdir'/decomp_table`hetero'.tex, noobs not nostar ///
        varlabels(`e(labels)') mlab("FTHBs" "Investment") nonumber tex replace
end /* %> */

capture program drop table_extmeans /* %< */
program define table_extmeans
    syntax anything(name=outputdir), [vars(varlist)]

    if "`vars'" == "" local vars nextAssets posIncShock posIncShock2 negIncShock2
    eststo groups: estpost tabstat `vars', by(groups) statistics(mean) ///
        columns(statistics) nototal
    esttab groups using `outputdir'/extstat_table.tex, main(mean) unstack ///
        noobs not nostar nonumber l tex replace
end /* %> */

capture program drop table_exttiming /* %< */
program define table_exttiming
    syntax anything(name=outputdir), [trunc(integer 3) hetero(string)]

    preserve
    gen adjDiff_trunc = sign(adjDiff)*min(abs(adjDiff), `trunc')
    * TODO
    label def trunc -3 "2+ Negative" -2 "2 Negative" -1 "1 Negative" ///
        0 "No Extensive" 1 "1 Positive" 2 "2 Positive" 3 "2+ Positive", replace
    label val adjDiff_trunc trunc

    eststo tabular: estpost tab decomp_ext adjDiff_trunc
    esttab tabular using `outputdir'/exttiming_table`hetero'.tex, c(rowpct(fmt(2))) ///
        unstack noobs collabels(none) nonumber not nostar l tex replace
end /* %> */

capture program drop table_polmargins /* %< */
program define table_polmargins
    syntax anything(name=outputdir), fthb_trans(integer) cars_trans(integer)

    tempfile trunc
    preserve
    collapse (first) desc subtype (sum) value valuemonetary valueCARS valueFTHBT-valueFTHBSize_0_165 ///
        if agebin == 0 & !inlist(comp_stats, 4, 5), by(comp_stats)
    save `trunc'
    restore, preserve
    collapse (first) desc subtype (sum) value valuemonetary valueCARS valueFTHBT-valueFTHBSize_0_165 ///
        if agebin == 0 & inlist(comp_stats, 4, 5), by(comp_stats)
    append using `trunc'
    drop if missing(comp_stats)

    foreach var of varlist value* {
        gen Fac`var' = `var'/`var'[_N] if comp_stats < 4
        if regexm("`var'", "CARS") == 0 replace Fac`var' = ///
            (`var'[1]-`var'[2])/`fthb_trans' if comp_stats >= 4
        else replace Fac`var' = (`var'[1]-`var'[2])/`cars_trans' if comp_stats >= 4
    }
    drop if comp_stats > 4
    label var Facvalue "Main policy experiment"
    label var FacvalueC "CARS policy experiment"
    label var Facvaluem "Credit applicable at purchase"
    label var FacvalueFTHBT "Two-period tax credit"
    label var FacvalueFTHBSize_0_03 "Small value tax credit (\\\$2,000)"
    label var FacvalueFTHBSize_0_165 "Large value tax credit (\\\$11,000)"

    eststo tabA: estpost tabstat Facvalue FacvalueCARS, by(comp_stats) ///
        statistics(mean) columns(statistics) not
    eststo tabB: estpost tabstat Facvalue Facvaluem FacvalueFTHBT-FacvalueFTHBSize_0_165, by(comp_stats) ///
        statistics(mean) columns(statistics) not

    local labels eql("Total timing margin, 1 period" "Total timing margin, 2 periods" ///
        "Total timing margin, 2+ periods" "Net extensive margin")
    esttab tabA using `outputdir'/comp_table1.tex, unstack main(mean) l `labels' ///
        nostar not noobs tex nonumber replace
    esttab tabB using `outputdir'/comp_table2.tex, unstack main(mean) l `labels' ///
        nostar not noobs tex nonumber replace
end /* %> */

capture program drop table_extdecomp /* %< */
program define table_extdecomp
    syntax anything(name=outputdir), mvar(varname)

    preserve
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

    collapse (first) desc-subtype (sum) `mvar', by(agebin ext_dummy)
    list
    drop if missing(ext_dummy)
    by agebin: gen value_prop = `mvar'/`mvar'[_N]
    drop if ext_dummy == 99
    estpost tab agebin ext_dummy [iw=value_prop] if agebin <= 50, notot
    esttab . using `outputdir'/extensive_table_`mvar'.tex, noobs not nostar nonumber unstack ///
        eql("0-1 Extensive" "0-1+ Extensive" "1-2 Extensive" "1-2+ Extensive" ///
            "2-3 Extensive" "2-1 Extensive") tex replace

end /* %> */

* %>

************************************************************
* %< 2. Figure programs.
************************************************************
capture program drop figure_heatmaptrans /* %< */
program define figure_heatmaptrans
    syntax anything(name=outputdir)

    tempfile income_temp // first, build a square grid of income values
    preserve
    clear
    set obs 9
    egen income_l = fill(4/12)
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
end /* %> */

capture program drop figure_rentassets /* %< */
program define figure_rentassets
    syntax anything(name=outputdir), datadir(string) [hetero(string)]

    tempfile rent_stationary
    preserve
    import delimited using `datadir'/indivAssets.csv, clear
    keep if rent == 1 & age < 40
    keep id-nextassets
    save `rent_stationary'

    restore
    append using `rent_stationary'
    replace marginal = -1 if missing(marginal)
    twoway (hist nextassets if marginal==-1 & nextassets < 3, ///
            width(0.05) fcolor(blue*0.2) lcolor(white)) ///
           (hist assets if marginal==0, width(0.05) fcolor(none) lcolor(purple)), ///
           legend(label(1 "Stationary equilibrium") label(2 "Inframarginal")) $graphconfig
    graph export `outputdir'/assets_planned.pdf, as(pdf) replace
    twoway (hist nextassets if marginal==-1 & nextassets < 3, ///
            width(0.05) fcolor(blue*0.1) lcolor(white)) ///
           (hist assets if marginal==1, width(0.05) fcolor(none) lcolor(purple)), ///
           legend(label(1 "Stationary equilibrium") label(2 "Marginal")) $graphconfig
    graph export `outputdir'/assets_marginal`hetero'.pdf, as(pdf) replace
end /* %> */

capture program drop figure_infra /* %< */
program define figure_infra
    syntax anything(name=outputdir), [hetero(string)]

    preserve
    keep if marginal==0
    keep if mod(_n, 50) == 0
    gen leverage = -nextAssets/nextDurables
    gen leverage_ss = -nextAssets_ss/nextDurables_ss
    label var leverage_ss "Leverage absent policy"
    label var leverage "Leverage given policy"

    qqplot nextDurables nextDurables_ss if marginal==0 & leverage > 0, ///
        msize(small) $graphconfig
    graph export `outputdir'/infra_durables`hetero'.pdf, as(pdf) replace
    qqplot consumption consumption_ss if marginal==0 & leverage > 0, ///
        msize(small) $graphconfig
    graph export `outputdir'/infra_consumption`hetero'.pdf, as(pdf) replace
    qqplot leverage leverage_ss if marginal==0 & leverage > 0, msize(small) ///
        $graphconfig
    graph export `outputdir'/infra_leverage`hetero'.pdf, as(pdf) replace

end /* %> */

capture program drop figure_incpull /* %< */
program define figure_incpull
    syntax anything(name=outputdir)

    tempfile pullf_inc temp
    local normalizer = 67.250 // From initial calibration in SCF
    preserve

    collapse (mean) income_val if pullforward <= 10 | pullforward==40 , ///
        by(pullforward)
    save `pullf_inc'
    restore, preserve
    collapse (mean) income_val if pullforward <= 10 | pullforward==40, ///
        by(pullforward extensive)
    replace extensive = 0 if missing(extensive)
    save `temp'
    use `pullf_inc', clear
    append using `temp'
    save `pullf_inc', replace

    * Generating income means for those pulling 2+ periods, 5+ periods, 10+ periods
    foreach val in 2 5 10 {
        restore, preserve
        collapse (mean) income_plus=income_val if ///
            pullforward > `val' & pullforward < 40
        gen pullf_bar = `val' + 1
        save `temp', replace
        local m`val' = income_plus[1]
        use `pullf_inc', clear
        append using `temp'
        save `pullf_inc', replace

        restore, preserve
        collapse (mean) income_plus=income_val if ///
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
    restore, preserve
    use `pullf_inc', clear
    list

    replace pullf_bar = pullforward if pullforward <= 2 | pullforward == 40
    replace pullf_bar = 4 if pullf_bar == 6
    replace pullf_bar = 5 if pullf_bar == 11
    replace pullf_bar = 6 if pullf_bar == 40
    label def pullf 1 "1" 2 "2" 3 "2+" 4 "5+" 5 "10+" 6 "Pos. Extensive"
    label val pullf_ pullf
    replace income_val = income_val*`normalizer'
    replace income_plus = income_plus*`normalizer'
    list

    #delimit ;
    twoway (bar income_val pullf_bar, barw(0.7) color(blue*0.1))
        (bar income_plus pullf_bar, barw(0.7) color(blue*0.1))
        if missing(extensive), legend(off)
        xti("Imputed periods durable purchased moved forward") yti("Mean Income Level")
        xlab(,valuelabel) $graphconfig 
        note("Note: Income is normalized in the model by average SCF Income "
              "(approx. $67,250 2013 USD).");
    graph export `outputdir'/pulledInc_1.pdf, as(pdf) replace;

    twoway (bar income_val pullf_bar if missing(extensive), barw(0.7) color(blue*0.1))
        (bar income_plus pullf_bar if missing(extensive), barw(0.7) color(blue*0.1))
        (bar income_val pullf_bar if extensive == 1, barw(0.7) fcolor(none) lcolor(purple) lw(thick))
        (bar income_plus pullf_bar if extensive == 1, barw(0.7) fcolor(none) lcolor(purple) lw(thick))
        , xti("Imputed periods durable purchased moved forward")
        yti("Mean Income Level") xlab(,valuelabel) $graphconfig
        legend(order(1 3) label(1 "Aggregate") label(3 "Positive Ext Margin"))
        note("Note: Income is normalized in the model by average SCF Income "
              "(approx. $67,250 2013 USD).");
    graph export `outputdir'/pulledInc_2.pdf, as(pdf) replace;
    #delimit cr

end /* %> */

capture program drop figure_margevent /* %< */
program define figure_margevent
    syntax anything(name=outputdir), datadir(string)

    local datadir $dboxdir/model
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

    #delimit ;
    local plotquery1 (connect value var1 if var1 < 4) 
        (scatter value var1 if var1 == 4, mc(navy))
        (connect value var2 if var2 < 4, lp("-..") color(gs9) m(Dh))
        (scatter value var2 if var2 == 4, color(gs9) m(Dh));
    local plotquery2 legend(order(1 3) label(1 "Policy response")
        label(3 "Steady-state counterfactual") rows(2)) xline(0, lp("-..."))
        xline(3.2, lc(gs7)) xlab(-2(1)4, valuelabel)
        ylab(, angle(0)) yti(Variable value) xti(Period around policy takeup);

    twoway `plotquery1' if wealth ==".0" , `plotquery2'
        by(fullLab, row(3) yrescale colfirst $graphconfig) name(event, replace);
    graph export `outputdir'/margEvent_0.pdf, as(pdf) replace;

    preserve;
    keep if regexm(fullLab, "Pos|Neg") == 0;
    encode fullLab, gen(facet);
    label def event_cap 1 "Poor Agent, Consumption" 2 "Poor Agent, Owned Housing"
        3 "Poor Agent, Rental Housing" 4 "Middle Agent, Consumption"
        5 "Middle Agent, Owned Housing" 6 "Middle Agent, Rental Housing"
        7 "Rich Agent, Consumption" 8 "Rich Agent, Owned Housing"
        9 "Rich Agent, Rental Housing";
    label val facet event_cap;

    twoway `plotquery1' if wealth ==".0" & inlist(var, "H_rent", "H_own") ,
        `plotquery2' by(facet, row(2) yrescale colfirst $graphconfig)
        subtitle(,size(small)) name(event, replace);
    graph export `outputdir'/margEvent_1.pdf, as(pdf) replace;

    twoway `plotquery1' if wealth ==".0" & inlist(income_state, "5/13", "9/13") ,
        `plotquery2' by(facet, row(3) yrescale colfirst $graphconfig)
        subtitle(,size(small)) name(event, replace);
    graph export `outputdir'/margEvent_2.pdf, as(pdf) replace;

    restore;
    keep if regexm(fullLab, "Pos|Neg") == 1;
    encode fullLab, gen(facet);
    label def event_cap 4 "Pos Shock to Middle, Consumption" 5 "Pos Shock to Middle, Owned Housing"
        6 "Pos Shock to Middle, Rental Housing" 7 "Neg Shock to Rich, Consumption"
        8 "Neg Shock to Rich, Owned Housing" 9 "Neg Shock to Rich, Rental Housing";
    label val facet event_cap;

    twoway `plotquery1' if wealth ==".0", `plotquery2'
       by(facet, row(3) yrescale colfirst $graphconfig) subtitle(,size(small)) name(event, replace);
    graph export `outputdir'/margEvent_shk.pdf, as(pdf) replace;
    #delimit cr

end /* %> */

capture program drop figure_translives /* %< */
program define figure_translives
    syntax anything(name=identifier), datadir(string) outdir(string)

    preserve
    bys id cohort: egen adj_ss = total(adjust_ss)
    bys id cohort: egen adj_pol = total(adjust_pol)
    replace adjust_ss = policy*adjust_ss
    replace adjust_pol = policy*adjust_pol

    gen marklow = -1.5e-1
    gen markup = -5e-2
    local ret_time = 39-`2'
    tokenize `identifier'
    keep if id==`1' & cohort==`2' & cohort + policy < 50

    twoway (line income policy, lc(navy*0.2) lw(thick)) ///
       (line h_ss rent_ss policy, c(J) lp(solid dash) lw(medthick medthick)) ///
       (rcap markup marklow adjust_ss, lc(gs2)), legend(label(1 "Income") ///
       label(2 "Owned Durable") label(3 "Rented Durable") order(1 - 2 3)) ///
       title("Without Policy") saving(ss, replace) $graphconfig
    twoway (line income policy, lc(navy*0.2) lw(thick)) ///
       (line h_pol rent_pol policy, c(J) lp(solid dash) lw(medthick medthick)) ///
       (rcap markup marklow adjust_pol, lc(gs2)), legend(off) ///
       title("With Policy") yti("Model Normalized Values") ///
       xti("Years Since 1-Period Policy") saving(pol, replace) $graphconfig
    graph combine ss.gph pol.gph, r(2) graphregion(c(white))
    erase ss.gph
    erase pol.gph
    graph export `outdir'/transLives_`1'_`2'.pdf, as(pdf) replace
end /* %> */

capture program drop figure_agggraph /* %< */
program define figure_agggraph
    syntax anything(name=outputdir), baseval(real) XTItle(string) [hetero(string) WITHSS]

    sort agebin compval
    label define ages 0 "All ages" 20 "Ages 22-29" 30 "Ages 30-39" 40 "Ages 40-49"
    label val agebin ages
    
    preserve
    duplicates drop compval, force
    local comp_no = _N
    restore
    sum compval, meanonly
    local min = `r(min)'
    local comp_int = (`r(max)' - `min')/(`comp_no' - 1)
    local baseid = round(`baseval'/`comp_int', 1)
    replace compval=round(compval/.015, 1)
    local baseval = `baseval'/.015
    
    by agebin: gen buychar_pol = Buyer_characteristics__m/Buyer_characteristics__m[`baseid']*100
    by agebin: gen buychar_ss = Buyer_characteristics/Buyer_characteristics[`baseid']*100
    label var buychar_pol "Marginal Buyers"
    label var buychar_ss "Inframarginal Buyers"

    if "`withss'" != "" local yvars buychar_*
    else local yvars buychar_pol    
    twoway line `yvars' compval, by(agebin, $graphconfig) ///
        yti("Policy Magnitude (100=base)") xti(`xtitle') xline(`baseval')
    graph export `outputdir'/aggLPlots`hetero'.pdf, as(pdf) replace
        
end /* %> */

capture program drop figure_elasHeatmap /* %< */
program define figure_elasHeatmap
    syntax namelist [if], ytitle(string) xtitle(string) ztitle(string) save(string)

    preserve
    local normalizer = 67.250 // From initial calibration in SCF
    tokenize `namelist'
    local zvar `1'
    macro shift
    local varlist `*'
    
    foreach var of local varlist {
        capture gen `var' = string(`var'bin*`normalizer', "%5.1f")
        gsort `var' -compval
        * If takeup is zero, do not just code as "N/A"
        replace Proportion_of_buyers = 0.0 if !missing(Proportion_of_buyers[_n-1])
    }
    local varlist `varlist' compval
    keep if !missing(Proportion_of_buyers)
    replace Proportion_of_buyers = 1.0 if Proportion_of_buyers > 1.0
    replace compval=round(compval/.015, 1)
    
    *list in 1/10
    tokenize `varlist'
    heatmap `zvar' `1' `2' `if', polbr(8) ytitle("`ytitle'") ///
        out customf(0.15 0.30 0.50 0.70 0.85) save(`save') ///
        xtitle("`xtitle'") zti("`ztitle'")

end /* %> */

capture program drop figure_hotrends /* %< */
program define figure_hotrends
    syntax anything(name=outputdir)

    preserve

    keep if regexm(desc, "Per-period homeownership rate")
    drop valueCARS-valueFTHBTP_2
   
    destring subtype, replace force
    sort agebin subtype
   
    replace valueHIMD = . if missing(subtype)
    bys agebin: replace valueHIMD = valueHIMD[_n+1]
    twoway (line value valuemonetary valueFTHBPerm_ valueFTHBPermon valueHIMD subtype) ///
        if agebin==0, xtitle("Periods Since Policy Start") ///
        legend(rows(5) size(small)) ytitle("Homeowner Proportion") $graphconfig
    graph export `outputdir'/compPol_horates.pdf, as(pdf) replace
   
    twoway (line value valuemonetary valueFTHBPerm_ valueFTHBPermon valueHIMD subtype), ///
        by(agebin, $graphconfig) xtitle("Periods Since Policy Start") ///
        legend(rows(5) size(small)) ytitle("Homeowner Proportion")
    graph export `outputdir'/compPol_horates_age.pdf, as(pdf) replace

end
 

capture program drop figure_carsvintage /* %< */
program define figure_carsvintage
    syntax anything(name=outputdir), [elig(integer 3)]

    preserve
    keep if regexm(desc, "Proportion of buyers by capital vintage") == 1
    collapse (max) valueCARS*, by(agebin desc subtype)
    destring subtype, gen(vintage)
    sort agebin desc subtype
    encode desc, gen(vintType)

    twoway (scatter valueCARS vintage if vintType==2 & vintage > 0 & vintage <= 15, m(Oh)) ///
           (scatter valueCARS vintage if vintType==1 & vintage > 0 & vintage <= 15) if agebin <= 0, ///
           xline(`elig', lc(gs9)) by(agebin, $graphconfig) yti("Proportion of buyers") xti("Capital vintage") ///
           legend(label(1 "Stationary equilibrium") label(2 "Policy simulation"))
    graph export `outputdir'/CARSvintage.pdf, replace

end /* %> */

* %>


* %< MICRODATA PROGRAM
capture program drop gen_microdata
program define gen_microdata
    syntax anything(name=datadir)

    load_analysis_data microdata, datadir(`datadir') ///
        modelname(experiment_monetary_nodown)

    table_aggdecomp $outdir/stata, hetero(_nodown)
    table_exttiming $outdir/stata, hetero(_nodown)
    figure_heatmaptrans $outdir/stata
    figure_rentassets $outdir/stata, datadir(`datadir') hetero(_nodown)
    figure_infra $outdir/stata, hetero(_nodown)
    table_extmeans $outdir/stata
    figure_incpull $outdir/stata

    figure_margevent $outdir/stata, datadir(`datadir')

    load_analysis_data microdata, datadir(`datadir') ///
        modelname(experiment_monetary)
    table_aggdecomp $outdir/stata, hetero(_ondown)
    table_exttiming $outdir/stata, hetero(_ondown)
    figure_infra $outdir/stata, hetero(_ondown)
    figure_rentassets $outdir/stata, datadir(`datadir') hetero(_ondown)

    load_analysis_data microdata, datadir(`datadir') ///
        modelname(experiment_CARS)
    table_aggdecomp $outdir/stata, hetero(_CARS)
    figure_infra $outdir/stata, hetero(_CARS)

    load_analysis_data microdata, datadir(`datadir') ///
        modelname(experiment_CARS_nodown)
    table_aggdecomp $outdir/stata, hetero(_CARS_nodown)
    figure_infra $outdir/stata, hetero(_CARS_nodown)

end /* %> */

* %< TRANS. LIVES PROGRAM
capture program drop gen_lives
program define gen_lives
    syntax anything(name=datadir)

    load_analysis_data lives, datadir(`datadir') ///
        modelname(experiment_monetary_nodown)
    figure_translives 2399 13 , datadir(`datadir') outdir($outdir/stata)
    copy "$outdir/stata/transLives_2399_13.pdf" "$outdir/stata/transLives_posExt.pdf", replace
    figure_translives 2449 23 , datadir(`datadir') outdir($outdir/stata)
    copy "$outdir/stata/transLives_2449_23.pdf" "$outdir/stata/transLives_negExt.pdf", replace
    load_analysis_data lives, datadir(`datadir') ///
        modelname(experiment_CARS)
    figure_translives 1309 8 , datadir(`datadir') outdir($outdir/stata) // CARS example
    copy "$outdir/stata/transLives_1309_8.pdf" "$outdir/stata/transLives_CARS.pdf", replace
end /* %> */

* %< MODEL STATS PROGRAM
capture program drop gen_stats
program define gen_stats
    syntax anything(name=datadir)

    load_analysis_data stats, datadir(`datadir')
    figure_carsvintage $outdir/stata
    figure_hotrends $outdir/stata
    // Arguments are captured from running experiment.py - UPDATE IF NECESSARY
    table_polmargins $outdir/stata, fthb_trans(14671) cars_trans(30303)
    table_extdecomp $outdir/stata, mvar(value)
end /* %> */

* %< ALT. SUBSIDIES COMP PROGRAM
capture program drop gen_comps
program define gen_comps
    syntax anything(name=datadir)

    load_analysis_data comp, datadir(`datadir') modelname(FTHBSize)
    figure_agggraph $outdir/stata, baseval(.12) ///
        xti("Value of Subsidy (in increments of $1000)") hetero(_nodown)
    load_analysis_data comp, datadir(`datadir') modelname(FTHBonDownSize)
    figure_agggraph $outdir/stata, baseval(.12) ///
        xti("Value of Subsidy (in increments of $1000)") hetero(_ondown)
    load_analysis_data Elas, datadir(`datadir') modelname(FTHBSize)
    figure_elasHeatmap Proportion_of_buyers income_val, ytitle("Income (000s)") ///
        xtitle("Value of Subsidy (in increments of $1000)") ///
        ztitle("Prob. becoming FTHB") save($outdir/stata/heatmap_incomes.pdf)
    * AD HOC call to heatmap where policy applicable on down
    load_analysis_data Elas, datadir(`datadir') modelname(FTHBSizeondown)
    figure_elasHeatmap Proportion_of_buyers income_val, ytitle("Income (000s)") ///
        xtitle("Value of Subsidy (in increments of $1000)") ///
        ztitle("Prob. becoming FTHB") save($outdir/stata/heatmap_incomesOD.pdf)

end /* %> */


capture program drop main
program define main

    gen_ssfit $outdir/stata/calibration
    gen_microdata $dboxdir/model
    gen_lives $dboxdir/model
    gen_stats $dboxdir/model
    gen_comps $dboxdir/model
end 
