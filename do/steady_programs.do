
/* %< LORENZ CURVES FOR WEALTH, DATA VS. MODEL */
capture program drop steady_lorenz
program define steady_lorenz
    syntax anything(name=varquery), datadir(string) simuldir(string) outdir(string)

    tempfile scf_filter
    use "`datadir'/lorenz_test.dta", clear
    keep if mod(_n, 20) == 0
    twoway (line houses_hCum houses_hPct if !missing(houses_h)) ///
        (line labIncWork_hCum labIncWork_hPct if !missing(labIncWork_h)), legend ( ///
             label(1 "Housing, homeowners") label(2 "Income, homeowners")) $graphconfig	
    graph export "`outdir'/wealthDist_data.pdf", as(pdf) replace
    keep *Cum *Pct
    rename (labIncWorkCum labIncWorkPct netWorthncCum netWorthncPct ///
            netliquid_rCum netliquid_rPct housesCum) ///
        (income incomePct networth networthPct netliquid_r netliquid_rPct houses)
    gen type = "data"
    save `scf_filter'

    import delimited `simuldir'/lorenzTest.csv, clear
    foreach var of varlist income* houses* netw* {
            sum pct if !missing(`var')
            gen `var'Pct = pct/`r(max)'*100
    }

    keep if mod(_n, 50) == 0 // Clears up graph, fewer points

    if regexm("`varquery'", "DLP") == 1 {
    twoway (line houses_h houses_hPct if !missing(houses_h)) ///
        (line income_h income_hPct if !missing(income_h)), legend( ///
             label(1 "Housing, homeowners") label(2 "Income, homeowners")) $graphconfig
    graph export "`outdir'/wealthDist_model.pdf", as(pdf) replace
    }

    append using `scf_filter'
    if regexm("`varquery'", "houses") == 1 {
    twoway (line houses housesPct if type != "data", lc(maroon)) ///
           (line houses housesPct if type == "data", lp("--.")), ///
            xti("") yti("Housing") legend(label(1 "Model") label(2 "Data")) $graphconfig
    graph export "`outdir'/houses_dist_model.pdf", as(pdf) replace
    }

    if regexm("`varquery'", "income") == 1 {
    twoway (line income incomePct if type != "data", lc(maroon)) ///
           (line income incomePct if type == "data" & incomePct <=100, lp("--.")), ///
            xti("") yti("Income") legend(label(1 "Model") label(2 "Data")) $graphconfig
    graph export "`outdir'/income_dist_model.pdf", as(pdf) replace
    }

    if regexm("`varquery'", "networth") == 1 {
    twoway (line networth networthPct if type != "data", lc(maroon)) ///
           (line networth networthPct if type == "data", lp("--.")), ///
            xti("") yti("Net Worth") legend(label(1 "Model") label(2 "Data")) $graphconfig
    graph export "`outdir'/networth_dist_model.pdf", as(pdf) replace
    }

end
/* %> */

/* %< INDIVIDUAL ASSETS DATA PROCESSING */
capture program drop steady_sumstats
program define steady_sumstats
    syntax, datadir(string) simuldir(string) outdir(string)

    import delimited `simuldir'/indivAssets.csv, clear
    gen leverage = -nextassets/nextdurables*(1-rent)
    replace leverage = . if leverage <= 0
    sort id age

    * Normalize to mean income simulated, not base value
    foreach var of varlist nextassets networth {
        sum income_val, meanonly
        replace `var' = `var'/`r(mean)'
    }
    * Produces tables similar to BGLV Table 2.
    eststo lev: estpost sum leverage if age <= 38, d
    eststo ownstats: estpost sum nextassets networth if rent == 0 & age <= 38, d
    eststo rentstats: estpost sum nextassets networth if rent == 1 & age <= 38, d
    esttab ownstats rentstats using `outdir'/ss_asset_stats.tex, ///
        c(p10 p25 p50 p75 p90) noobs not nostar nonumber tex replace

    * Lorenz curve for renters' financial assets (data TBD)
    keep if rent==1 & age <= 38
    sum nextassets
    sort nextassets
    gen nextassetscum = sum(nextassets)/`r(sum)'*100
    egen Pct = rank(nextassets), u
    gen nextassetsPct = Pct/_N*100
    sort nextassetsPct
    keep if mod(_n, 50) == 0 // Clears up graph, fewer points
    append using `datadir'/lorenz_test.dta

    twoway (line nextassetscum nextassetsPct) ///
        (line netliquid_rCum netliquid_rPct if !missing(netliquid_rCum)), ///
        yti("Assets for Renters") xti("") legend( ///
        label(1 "Model") label(2 "Data")) $graphconfig
    graph export "`outdir'/networth_rent_model.pdf", as(pdf) replace

end
/* %> */

capture program drop steady_moments /* %< */
program define steady_moments
    syntax anything(name=varquery), datadir(string) simuldir(string) outdir(string)

    tempfile model modelH data dataH
    * Statistics over entire population
    import delimited `simuldir'/moments_1yearbin.csv, clear
    capture drop model
    gen model = "model"
    save `model', replace

    * Statistics for only houseowners
    import delimited `simuldir'/momentsH_1yearbin.csv, clear
    capture drop model
    gen model = "model"
    save `modelH', replace

    * These two datasets are from the SCF calibration.
    * Renaming is needed since "wide format" makes plotting the series easier.
    import delimited `datadir'/graph_data_byageind.csv, clear
    keep ageind netliquid_exhousingdebt homeownershiprate
    rename * (ageind avgassetsexdebt fracown)
    gen model = "data"
    save `data', replace

    import delimited `datadir'/graph_data_byageind_homeowners.csv, clear
    keep ageind netliquid_exhousingdebt housepos
    rename * (ageind avgassetsexdebtowners avghouseowners)
    merge 1:1 ageind using `data'
    save `dataH', replace


    * Merge. Convert back to age bins again.
    use `model', clear
    merge 1:n model agebin using `modelH'
    drop _m
    append using `dataH'
    drop _m
    *replace agebin = 5*ageind + 17.5 if missing(agebin)
    replace agebin = ageind if missing(agebin)
    drop ageind
    drop if agebin > 70 | agebin < 21
    order model, first
    
    * Build the graphs. Calibration is on last two, homeownership rate and housing size.
    * Five series are matched, though for presentation likely  only housing size,
    * homeownership rate and liquid assets (2) are shown.
    local graphcap legend(lab(1 "Data") lab(2 "Model")) xti("Agent Age")
    if regexm("`varquery'", "assetsexdebt") == 1 {
    twoway (scatter avgassetsexdebt agebin if model == "data") ///
           (connect avgassetsexdebt agebin if model == "model" & agebin >= 22, lp("--") m(X)), ///
            `graphcap' yti("Net Financial Assets") $graphconfig
    graph export "`outdir'/avgAssetsExDebt.pdf", as(pdf) replace
    }
    if regexm("`varquery'", "assetsexdebtowners") == 1 {
    twoway (scatter avgassetsexdebtowners agebin if model == "data") ///
           (connect avgassetsexdebtowners agebin if model == "model" & agebin >= 22, lp("--") m(X)), ///
            `graphcap' yti("Net Financial Assets, Homeowners") $graphconfig
    graph export "`outdir'/avgAssetsExDebtOwners.pdf", as(pdf) replace
    }
    if regexm("`varquery'", "own") == 1 {
    twoway (scatter fracown agebin if model == "data") ///
           (connect fracown agebin if model == "model" & agebin >= 22, lp("--") m(X)), ///
            `graphcap' yti("Homeownership Rate") $graphconfig
    graph export "`outdir'/fracOwn.pdf", as(pdf) replace
    }
    if regexm("`varquery'", "houses") == 1 {
    twoway (scatter avghouseowners agebin if model == "data") ///
           (connect avghouseowners agebin if model == "model" & agebin >= 22, lp("--") m(X)), ///
            `graphcap' yti("House Size") yscale(r(0 4)) ylab(1(01)4) $graphconfig
    graph export "`outdir'/avgHouseOwners.pdf", as(pdf) replace
    }

end /* %> */

capture program drop gen_ssfit
program define gen_ssfit
    syntax anything(name=outdir), dboxdir(string)

    steady_moments own houses assetsexdebtowners, datadir(`dboxdir'/SCF) ///
        simuldir(`dboxdir'/model) outdir(`outdir')
    steady_lorenz DLP houses income networth, datadir(`dboxdir'/SCF) ///
        simuldir(`dboxdir'/model) outdir(`outdir')
    steady_sumstats, datadir(`dboxdir'/SCF) simuldir(`dboxdir'/model) outdir(`outdir')
    * OTHER STATS TO COME
end
