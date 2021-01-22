/*******************************************************************************
    Title: elasticity_graphs.do
    Purpose: Produce visualization of individual elasticities (explain?)
    Author: TC, DB
    Date: 2016-03-31
*******************************************************************************/


/* %< TEMP: extensive margin elasticities */
import delimited $outdir/stats_elasticities_experiment_monetary_dep.csv, clear
drop if age >= 38 & !missing(age)
encode desc, gen(model)
preserve
foreach var of varlist age assetsbin income_valbin {
    keep if !missing(`var') & missing(income)
	keep `var' value model
    reshape wide value, i(`var') j(model)
    gen elas_temp = log(value1/value2)
	list
    twoway line elas_temp `var', yti("Model Elasticity Estimate, Logs")
	sleep 5000
	restore, preserve
}
/* %> */


/*
Output elas.csv with this Pandas query on server
(merging steady-state and policy outputs)

for model1, model2, title in [('experiment_CARSnoF', 'experiment_CARS', 'elas')
                             ,('experiment_CARS_nocollnoF',
                               'experiment_CARS_nocoll', 'elas_nocoll'),
                              ('experiment_little','experiment_theta','elas_fthb')]:
    policy = (modelCal.appendModels('transition_fthb', model=[model1, model2])
              .appended.drop_duplicates(['id', 'age', 'model']))
    if title == 'elas_fthb': # Looking at ext and int margin buyers
        ss = modelCal.appendModels('fthb', model=[model1, model2]).appended
        elas = pd.merge(ss, policy[policy['age'] == policy['ageBought']],
                        on=['id', 'age', 'model'], suffixes=('', 'Pol'))
    else: # looking at intensive margin buyers
        ss = modelCal.appendModels('dist_fthb', model=[model1, model2]).appended
        elas = pd.merge(ss, policy[policy['age'] == policy['ageBought']],
                        left_on=['id', 'ageBought', 'model'],
                        right_on=['id', 'age', 'model'], suffixes=('', 'Pol'))
    elas.describe()
    elas.to_csv('%s.csv' % title, index=False)


*/

import delimited $outdir/elas.csv, clear /* %< */
drop if poltaken == "First-time"
foreach var of varlist leverage* {
    replace `var' = -`var'
    replace `var' = 0 if `var' < 0
}

gen Dchg = nextdurablespol - nextdurables
sum agebought* durables* leverage* wealth* // Ensure identical between merged data

* Caption locals, change depending on what kind of model is simulated
local titlesL ytitle(Log difference in housing stock) title("No adjustment costs") replace
local titlesF ytitle(Log difference in housing stock) title("Baseline") replace
local explanation caption("Simulation of a 5% house price subsidy lasting for" ///
   "1 year, applicable to repeat buyers.")
   
  
binscatter Dchg agebought if model=="experiment_CARSnoF", n(80) ///
    reportreg xtitle("Age of House Purchase") `titlesL' savegraph(1.gph)
binscatter Dchg agebought if model=="experiment_CARS", n(80) ///
    reportreg xtitle("Age of House Purchase") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasAge.pdf, replace

binscatter Dchg durables if model=="experiment_CARSnoF", n(40) ///
    reportreg xtitle("Current-period durable stock") `titlesL' savegraph(1.gph)
binscatter Dchg durables if model=="experiment_CARS", n(40) ///
    reportreg xtitle("Current-period durable stock") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasDurables.pdf, replace

binscatter Dchg leverage if model=="experiment_CARSnoF", n(40) ///
    reportreg xtitle("Leverage (= 0 if assets +ve)") `titlesL' savegraph(1.gph)
binscatter Dchg leverage if model=="experiment_CARS", n(40) ///
    reportreg xtitle("Leverage (= 0 if assets +ve)") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasLeverage.pdf, replace

binscatter Dchg wealth if model=="experiment_CARSnoF", n(40) ///
    reportreg xtitle("Wealth (assets + income shock)") `titlesL' savegraph(1.gph)
binscatter Dchg wealth if model=="experiment_CARS", n(40) ///
    reportreg xtitle("Wealth (assets + income shock)") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasWealth.pdf, replace

* %>

import delimited $outdir/elas_nocoll.csv, clear /* %< */
drop if poltaken == "First-time"

gen Dchg = nextdurablespol - nextdurables
sum agebought* durables* leverage* wealth* // Ensure identical between merged data

* Caption locals, change depending on what kind of model is simulated
local titlesL ytitle(Log difference in housing stock) title("No collateral, no adj. costs") replace
local titlesF ytitle(Log difference in housing stock) title("No collateral") replace
local explanation caption("Simulation of a 5% house price subsidy lasting for" ///
   "1 year, applicable to repeat buyers.")
   
  
binscatter Dchg agebought if model=="experiment_CARS_nocollnoF", n(80) ///
    reportreg xtitle("Age of House Purchase") `titlesL' savegraph(1.gph)
binscatter Dchg agebought if model=="experiment_CARS_nocoll", n(80) ///
    reportreg xtitle("Age of House Purchase") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasAge_nocoll.pdf, replace

binscatter Dchg durables if model=="experiment_CARS_nocollnoF", n(40) ///
    reportreg xtitle("Current-period durable stock") `titlesL' savegraph(1.gph)
binscatter Dchg durables if model=="experiment_CARS_nocoll", n(40) ///
    reportreg xtitle("Current-period durable stock") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasDurables_nocoll.pdf, replace

binscatter Dchg wealth if model=="experiment_CARS_nocollnoF", n(40) ///
    reportreg xtitle("Wealth (assets + income shock)") `titlesL' savegraph(1.gph)
binscatter Dchg wealth if model=="experiment_CARS_nocoll", n(40) ///
    reportreg xtitle("Wealth (assets + income shock)") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasWealth_nocoll.pdf, replace

* %>

import delimited $outdir/elas_fthb.csv, clear /* %< */


gen Dchg = nextdurablespol - nextdurables


* Intensive margin: same as before
preserve
drop if renting == 1

* Caption locals, change depending on what kind of model is simulated
local titlesL yline(0, lc(gs3)) ytitle(Log difference in housing stock) ///
    title("Equivalent discount on down") replace
local titlesF yline(0, lc(gs3)) ytitle(Log difference in housing stock) ///
    title("Subsidy unapplicable to down") replace
local explanation caption("Simulation of a 5% house price subsidy lasting for" ///
   "1 year, applicable to first-time homebuyers.")

binscatter Dchg agebought if model=="experiment_theta", n(80) ///
    reportreg xtitle("Age of House Purchase") `titlesL' savegraph(1.gph)
binscatter Dchg agebought if model=="experiment_little", n(80) ///
    reportreg xtitle("Age of House Purchase") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasAge_fthb.pdf, replace

binscatter Dchg wealth if model=="experiment_theta", n(40) ///
    reportreg xtitle("Wealth (assets + income shock)") `titlesL' savegraph(1.gph)
binscatter Dchg wealth if model=="experiment_little", n(40) ///
    reportreg xtitle("Wealth (assets + income shock)") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasWealth_fthb.pdf, replace

binscatter Dchg durablespol if model=="experiment_theta", n(40) ///
    reportreg xtitle("Current-period rental housing") `titlesL' savegraph(1.gph)
binscatter Dchg durablespol if model=="experiment_little", n(40) ///
    reportreg xtitle("Current-period rental housing") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing elasticities relative to steady-state counterfactual")
graph export $outdir/stata/HouseElasDurables_fthb.pdf, replace

* Extensive margin: looking at those not found in dist_fthb files.
restore, preserve
drop if renting == 0

local titlesL yline(0, lc(gs3)) ytitle(Log difference in housing vs. rental stock) ///
    title("Equivalent discount on down") replace
local titlesF yline(0, lc(gs3)) ytitle(Log difference in housing vs. rental stock) ///
    title("Subsidy unapplicable to down") replace
local explanation caption("Comparison of housing levels relative to planned" ///
     "rental purchase, for extensive margin buyers.")

binscatter Dchg wealth if model=="experiment_theta", n(40) ///
    reportreg xtitle("Wealth (assets + income shock)") `titlesL' savegraph(1.gph)
binscatter Dchg wealth if model=="experiment_little", n(40) ///
    reportreg xtitle("Wealth (assets + income shock)") `titlesF' savegraph(2.gph)
graph combine 1.gph 2.gph, graphregion(c(white)) `explanation' ///
    title("Housing up/downsizing from rental value")
graph export $outdir/stata/HouseChg_fthb.pdf, replace

* %>
erase 1.gph
erase 2.gph

