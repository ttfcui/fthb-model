#!/usr/bin/python

""" Title: experiment.py
    Purpose: Master file for running FTHB model subsidy experiments.
    Author: TC, DB
    Date: 2019-05-12

"""

# CHECK WHETHER PE_SUBSIDIES OR GE_SUBSIDIES SET BEFORE RUNNING.
# FOR NOW PE_SUBSIDIES IS ASSUMED.

import experiment_programs as run
import numpy as np

from model.model_iterate import lifecycle_iterate
import pandas as pd

# %< Principal variations of temporary policy

# FTHB, required maintenance of durable
run.simFTHB([('experiment_monetary_little', 0.00, 0.00),
             ('experiment_monetary', 0.075, 1.00),
             ('experiment_monetary_nodown', 0.075, 0.00)], 1.0,
             gen_micro=True, gen_props=True)

# FTHB, no maintenance allowed
run.simFTHB([('experiment_monetary_dep_nodown', 0.12, 0.00),
             ('experiment_monetary_dep', 0.075, 1.00)], 0.0,
             gen_micro=True, gen_props=True)
# CARS, required maintenance
# TODO: There are some serious problems with this now - debug.
# run.simCARS([('experiment_repeat', 0.20, 1)], 1.0)

# CARS, no maintenance allowed
run.simCARS([('experiment_CARS', 0.20, 0, 0.035),
             ('experiment_CARS_nocoll', 1, 0, 0.035),
             ('experiment_CARS_altScrap', 0.20, 0, 0.01)], 0.0,
             gen_micro=True, gen_props=True)
# %>

# %< FTHB, general equilibrium (TBD?)
run.simFTHB_GE([('experiment_monetary_GE', 0.12, 1.00),
                ('experiment_monetary_nodown_GE', 0.12, 0.00)], 1.0)
# %>

# %< Simple comparative statistics checks

# Tax credit counted as 1% price deduction on housing
run.simFTHB([('experiment_monetary_prElas', 0.01, 1.00)], 1.0,
             percentages=True, gen_props=True, gen_micro=True)
# Tax credit counted as 5% price deduction on housing
run.simFTHB([('experiment_monetary_prElas5', 0.05, 1.00)], 1.0,
             percentages=True, gen_props=True)
# Tax credit counted as 5% price deduction on housing, next period
run.simFTHB([('experiment_monetary_prElas5_nodown', 0.05, 0.00)], 1.0,
             percentages=True, gen_props=True)

# Kaplan-Violante (2014) lump-sum transfer
run.simFTHB([('experiment_monetary_KV', 0.00, 0.00)], 1.0,
             gen_micro=True, gen_props=True, aggPol=True, stop_model=True)

# CARS experiment, fixed car size
# Note this has to be run SEPARATELY from lump-sum transfer,
# the polParams.txt file must be changed in between.
run.simCARS([('experiment_CARS_fixed', 0.20, 1, 0.035)], 0.0,
             gen_micro=True, gen_props=True, aggPol=True)

# Nodown policy but with DP shock (FHA like)
run.simFTHB([('experiment_monetary_FHAonly', 0.00, 0.00), 
             ('experiment_monetary_FHA_nodown', 0.12, 0.00), 
             ('experiment_monetary_FHA', 0.12, 1.00)], 1.0,
             gen_micro=True, gen_props=True, aggPol=True)

# Tax credit available as lump-sum
run.compPols(('experiment_FTHBDF', 'downflag', [0]), 'First-time', 1.0,
             share=1.00, gen_props=True, gen_micro=True)
# Tax credit with transaction cost rebate
run.compPols(('experiment_FTHBRB', 'rebatedflag', [1]), 'First-time', 1.0,
             gen_props=True)

# Tax credit, two periods
run.compPols(('experiment_FTHBTP', 'PolEnd', [2]), 'First-time', 1.0)

# Different down payment
run.compPols(('experiment_FTHBDown', 'thetamatlab',
              np.linspace(0.05, 0.35, 9)), 'First-time', 1.0, delModel=True)
# Different depreciation
run.compPols(('experiment_FTHBDep', 'delta',
              [0.025, 0.03, 0.05]), 'First-time', 1.0, delModel=True)
# Different fixed costs
run.compPols(('experiment_FTHBF', 'F', np.linspace(0.02, 0.12, 11)),
             'First-time', 1.0, delModel=True)

# Different quadratic fixed costs
run.compPols(('experiment_FTHBF2', 'F2', [0.005, 0.025, 0.05, 0.1]),
             'First-time', 1.0, delModel=True)

# Tax credit size variation
run.compPols(('experiment_FTHBSize', 'adjTransfer',
              np.linspace(0.015, 0.18, 12)), 'First-time', 1.0, delModel=True,
              gen_props=True, gen_micro=True)

run.compPols(('experiment_FTHBonDownSize', 'adjTransfer',
              np.linspace(0.015, 0.18, 12)), 'First-time', 1.0, delModel=True,
              gen_props=True, gen_micro=True, share=1.00)

run.compPols(('experiment_CARSSize', 'adjTransfer',
              np.linspace(0.015, 0.18, 12)), 'Repeat', 0.0, delModel=True)
# Different minimum house size
run.compPols(('experiment_FTHBHouse', 'Dmin', np.linspace(0.0, 1.0, 9)),
             'First-time', 1.0, delModel=True)


# EIS, CARS
EISgrid = np.concatenate([np.linspace(1.5, 3, 7), np.linspace(4, 8, 5)])
print EISgrid

run.compPols(('experiment_FTHBI', 'elasticity', EISgrid), 'First-time',
              1.0, delModel=True)
# delta, CARS
run.compPols(('experiment_CARSd', 'delta', [0.2, 0.075, 0.025]), 'Repeat', 0.0)
# minimum size, CARS
run.compPols(('experiment_CARSLevel', 'Dmin', [0.27]), 'Repeat', 0.0)
# eligibility, CARS
run.compPols(('experiment_CARSE', 'EligYrsR', [1, 10]), 'Repeat', 0.0)

# HIMD counterfactual simulation
run.compPols(('experiment_HIMD', 'himdflag', [1.0]), 'First-time', 1.0,
              tvalue=0.0)
# Permanent FTHB counterfactual simulation
run.compPols(('experiment_FTHBPerm', 'PolEnd', [99.0]), 'First-time', 1.0)
run.compPols(('experiment_FTHBPermonDown', 'PolEnd', [99.0]), 'First-time', 1.0,
             share=1.00)
"""
# period considered as quarter and not year
# NOTICE THAT, IN ORDER FOR THIS TO WORK, DETERMINISTIC INCOME
# PROCESSES, DEATH PROBS AND UTIL SCALES MUST ALSO BE TURNED
# OFF (SEE SHARE.F90)
call(['matlab', '-nodesktop', '-nodisplay', '-r',
      '"cd(\'matlab\');Tretire=152;simulateearningsprocess;quit;"'])
modelDict = genParams(carsParamsDict, {'thetamatlab': 0.20,
                      'Eta_transfer': 1, 'maint': 0.0})
modelDict['beta2'] = modelDict['beta2']**(1.0/4)
modelDict['delta'] = modelDict['delta']/4
modelDict['EligYrsF'] = 999
modelDict['EligYrsR'] = 28
modelDict.update({'Tretire': 152, 'Tdie': 153})

modelPol = lifecycle_iterate(modelDict)
modelPol.matlabPlot(pol='Return', life=153, model=['qtrtest'])
call(['matlab', '-nodesktop', '-nodisplay', '-r',
      '"cd(\'matlab\');simulateearningsprocess;quit;"'])

# %>
"""

# %< Output relevant statistics for each policy in a joint Stata file
facevaluelist = []
for value in np.linspace(0.015, 0.18, 12):
    facevaluelist += ['FTHBSize_%5.3f' % value, 'FTHBonDownSize_%5.3f' % value,
                      'FTHBSize_%4.2f' % value, 'FTHBonDownSize_%4.2f' % value]
    
run.polMerge(['monetary_nodown', 'monetary', 'monetary_dep', 'monetary_KV',
              'monetary_dep_nodown', 'CARS', 'CARS_fixed', 'CARS_nocoll',
              'CARS_altScrap', 'FTHBTP_2', 'HIMD_1.0', 'monetary_FHA_nodown',
              'monetary_FHA', 'FTHBPerm_99.0', 'FTHBPermonDown_99.0']
              + facevaluelist)
# %>

# %< Clean file containing life cycle profiles of policy takers
livesModels = ['experiment_monetary_nodown', 'experiment_monetary',
               'experiment_CARS']
lives = lifecycle_iterate.appendModels('transition_lives',
                                       model=livesModels).appended
lives.loc[:, 'rent_ss'] = lives['rent']*lives['h_ss']
lives.loc[:, 'h_ss'] = (1 - lives['rent'])*lives['h_ss']
lives.loc[lives['adjust_ss'] == 0, 'adjust_ss'] = np.nan
lives.loc[lives['adjust_pol'] == 0, 'adjust_pol'] = np.nan
lives.to_stata('fthb_trans_lives.dta')
# %>

# %< Select partial agg. statistics from comp. statics models
run.compMerge([('FTHBSize', np.linspace(0.015, 0.18, 12)),
               ('FTHBonDownSize', np.linspace(0.015, 0.18, 12)),
               ('FTHBHouse', np.linspace(0.0, 1.0, 9)),
               ('FTHBDown', np.linspace(0.05, 0.35, 9)),
               ('FTHBF', np.linspace(0.02, 0.12, 11)),
               ('FTHBI', EISgrid)], 'income_valbin')
run.compMerge([('CARSSize', np.linspace(0.03, 0.18, 11))], 'assetsbin')
run.compMerge([('monetary', ['', 'nodown', 'little'])], 'income_valbin')
run.compMerge([('CARS', ['', 'nocoll'])], 'assetsbin')
# %>

# %< Collate partial statistics of model stats by varying EIS
total = pd.read_csv('stats_final_experiment_FTHBI_1.5.csv')
for pol in EISgrid[1:]:
    out = pd.read_csv('stats_final_experiment_FTHBI_%s.csv' % pol)
    total = pd.merge(total, out, on=['agebin', 'desc', 'subtype'],
                     how='outer', suffixes=('', pol))
    total = total[(total['desc'].str.endswith("distributions")) |
                  (total['desc'].str.startswith("Buyer chara") &
                   total['subtype'].str.startswith("age"))]
total = total.set_index(['agebin', 'desc']).replace([np.inf, -np.inf], np.nan)
print total
total.sort_index().to_stata('fthb_EIScomp_final.dta')
# %>
