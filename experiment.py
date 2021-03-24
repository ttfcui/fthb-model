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
run.simCARS([('experiment_CARS', 0.10, 0, 0.020),
             ('experiment_CARS_lowcoll', 0.90, 0, 0.020),
             ('experiment_CARS_altScrap', 0.10, 0, 0.01)], 0.0,
             gen_micro=True, gen_props=True)
# %>

# %< FTHB, general equilibrium (TBD?)
run.simFTHB([('experiment_monetary_GE', 0.075, 1.00)], 1.0,
             gen_micro=True, gen_props=True, transLen=3,
             loadPrices=True)
run.simFTHB([('experiment_monetary_nodown_GE', 0.075, 0.00)], 1.0,
             gen_micro=True, gen_props=True, transLen=3,
             loadPrices=True)
# %>

# %< Comparative statistics checks

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

# HIMD counterfactual simulation
run.compPols(('experiment_HIMD', 'himdflag', [1.0]), 'First-time', 1.0,
              tvalue=0.0)
# Permanent FTHB counterfactual simulation
run.compPols(('experiment_FTHBPerm', 'PolEnd', [99.0]), 'First-time', 1.0)
run.compPols(('experiment_FTHBPermonDown', 'PolEnd', [99.0]), 'First-time', 1.0,
             share=1.00)

# %>

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
