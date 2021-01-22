#!/usr/bin/python

""" Title: experiment_programs.py
    Purpose: Master file for running FTHB model subsidy experiments.
    Author: TC, DB
    Date: 2019-02-22

"""

import pandas as pd
from model.model_iterate import lifecycle_iterate
from subprocess import call
from os import remove
from shutil import rmtree

# Parameters of interest in model (preset or calibrated)
algDict = {'agridsize': 120, 'Dgridsize': 55, 'zgridsize': 13,
           'hpnodes(1)': 0.173}
commonParamsDict = {'beta2': 0.9175, 'elasticity': 2.25, 'F': 0.06,
                    'F2': 0.00, 'psi': 5, 'ret_wealth': 0.05, 'tax1': 0.225,
                    'r': 0.024, 'rho_z': 0.91, 'sigma_z': 0.21,
                    'poismean': 5.8, 'balancer': -0.010, 'numhouseholds': 6000}
fthbParamsDict = {'EligYrsR': 99, 'EligYrsF': 3, 'delta': 0.022,
                  'thetamatlab': 0.20, 'rentUtil': 0.95, 'elasticity2': 0.859,
                  'rentPrem': 0.966e-2, 'Dmin': 1.10, 'polEnd': 1}
carsParamsDict = {'EligYrsR': 3, 'EligYrsF': 99, 'rentPrem': 9.9, 'Dmin': 0.45,
                  'elasticity2': 0.95, 'delta': 0.125, 'rentPremRetire': 9.9, 
                  'scrapped': 1.0, 'adjTransfer': 0.06, 'Eta_transfer': 1.00,
                  'polEnd': 1, 'custom_start': '.TRUE.'}


# Convenience parameter dict editor
def genParams(polParamsDict, polDict, **kwargs):
    modelDict = algDict.copy()
    if 'percentages' in kwargs:
        modelDict.update({'pctageflag': '.TRUE.'})
    if 'aggPol' in kwargs:
        modelDict.update({'agg_policies': '.TRUE.'})
    modelDict.update(commonParamsDict)
    modelDict.update(polParamsDict)
    modelDict.update(polDict)
    return modelDict


# FTHB simulation template
def simFTHB(polSpecList, maint, **kwargs):
    for model, value, prop in polSpecList:
        modelDict = genParams(fthbParamsDict, {'adjTransfer': value,
                              'Eta_transfer': prop, 'maint': maint}, **kwargs)
        statsList = ['%f' % modelDict['delta'], '%f' % modelDict['Dmin'],
		     '%f' % modelDict['polEnd']]
        propsList = ['%f' % modelDict['thetamatlab'], '%f' % modelDict['hpnodes(1)']]
        modelPol = lifecycle_iterate(modelDict)
        if 'stop_model' not in kwargs:
        modelPol.execSh(model=model)
        modelPol.matlabPlot(model=[model])
        call(['./gen_stats.py', model, 'First-time'] + statsList)
        if 'gen_elas' in kwargs:
            call(['./elasticities.py', model, 'First-time', statsList[0]])
        if 'gen_micro' in kwargs:
            call(['./gen_microdata.py', model, 'First-time', statsList[0]])
        if 'gen_props' in kwargs:
            call(['./gen_propensity.py', model, 'First-time'] + propsList)


# CARS simulation template
def simCARS(polSpecList, maint, **kwargs):
    for model, valueT, valueD, scrapV in polSpecList:
        modelDict = genParams(carsParamsDict, {'thetamatlab': valueT,
                              'downflag': valueD, 'scrapvalue': scrapV,
                              'maint': maint}, **kwargs)
        statsList = ['%f' % modelDict['delta'], '%f' % modelDict['Dmin'],
		     '%f' % modelDict['polEnd']]
        propsList = ['%f' % modelDict['thetamatlab'], '%f' % modelDict['hpnodes(1)']]
        modelPol = lifecycle_iterate(modelDict)
        if 'stop_model' not in kwargs:
        modelPol.execSh(model=model)
        modelPol.matlabPlot(pol='Return', model=[model])
        call(['./gen_stats.py', model, 'Repeat'] + statsList)
        if 'gen_elas' in kwargs:
            call(['./elasticities.py', model, 'Repeat', statsList[0]])
        if 'gen_micro' in kwargs:
            call(['./gen_microdata.py', model, 'Repeat', statsList[0]])
        if 'gen_props' in kwargs:
            call(['./gen_propensity.py', model, 'Repeat'] + propsList)
        # Disk space management (this transition file is massive)
        remove('output/%s/transition_fthb_%s.txt' % (model, model))


# FTHB, general equilibrium (TBD?)
def simFTHB_GE(polSpecList, maint, transLen=12, **kwargs):
    for model, value, prop in polSpecList:
        modelDict = genParams(fthbParamsDict, {'adjTransfer': value,
                              'Eta_transfer': prop, 'maint': maint,
                                               'hptransLength': transLen})
        modelPol = lifecycle_iterate(modelDict)
        if 'stop_model' not in kwargs:
        modelPol.execSh(model=model)  # WARNING: THIS TAKES A LONG TIME
        modelPol.matlabPlot(model=[model], life=transLen)
        aggs = pd.read_table('model/transition_aggs.txt', sep='\s+',
                             header=None, names=['price', 'hstock', 'rents',
                                                 'resales', 'newhomes'])
        aggs.to_csv('output/%s/transition_aggs.csv' % model, index=False)
        call(['./gen_stats.py', model, 'First-time',
              '0.022', '%f' % fthbParamsDict['Dmin']])
        if 'gen_micro' in kwargs:
            call(['./gen_microdata.py', model, 'First-time', '0.022'])


# Simple comparative statistics checks
def compPols(polSpec, polType, maint, share=0.00, tvalue=0.12, **kwargs):
    if polType == 'First-time':
        polDict = fthbParamsDict
        polParams = {'adjTransfer': tvalue, 'Eta_transfer': share, 'maint': maint}
    elif polType == 'Repeat':
        polDict = carsParamsDict
        polParams = {'thetamatlab': 0.20, 'downflag': 1, 'maint': maint}
    statsList = ['%f' % polDict['delta'], '%f' % polDict['Dmin'],
		 '%f' % polDict['polEnd']]
    modelB, stat, vlist = polSpec
    for val in vlist:
        model = '{}_{:6.4f}'.format(modelB, val)
        modelDict = genParams(polDict, polParams)
        propsList = ['%f' % modelDict['thetamatlab'], '%f' % modelDict['hpnodes(1)']]
        modelDict[stat] = val
        modelPol = lifecycle_iterate(modelDict)
        if 'stop_model' not in kwargs:
        modelPol.execSh(model=model)
        if polType == 'First-time':
            modelPol.matlabPlot(model=[model])
        elif polType == 'Repeat':
            modelPol.matlabPlot(pol='Return', model=[model])
        call(['./elasticities.py', model, polType] + [statsList[0]])
        call(['./gen_stats.py', model, polType] + statsList)
        if 'gen_micro' in kwargs:
            call(['./gen_microdata.py', model, polType] + [statsList[0]])
        if 'gen_props' in kwargs:
            call(['./gen_propensity.py', model, polType] + propsList)
        if 'delModel' in kwargs:
            rmtree(modelPol.mainD + 'output/%s/' % model)


# Output relevant statistics for each policy in a joint Stata file
def polMerge(models):
    initial = models[0]
    total = pd.read_csv('stats_final_experiment_%s.csv' % initial)
    aggVal = pd.read_csv('output/experiment_%s/PolSeries.txt'% initial,
                         names=['CumTrans', 'CumFTHB', 'PeriodInv'],
                         header=None)
    aggVal = aggVal.reset_index()
    aggVal['model'] = 'experiment_%s' % initial
    micro = pd.read_csv('masterdata_experiment_%s.csv' % initial)
    for pol in models[1:]:
        try:
        out = pd.read_csv('stats_final_experiment_%s.csv' % pol)
        total = pd.merge(total, out, on=['agebin', 'desc', 'subtype'],
                         how='outer', suffixes=('', pol))
        out2 = pd.read_csv('output/PolSeries_experiment_%s.txt'% pol,
                           names=['CumTrans', 'CumFTHB', 'PeriodInv'],
                           header=None)
        out2 = out2.reset_index()
        out2['model'] = 'experiment_%s' % pol
        aggVal = aggVal.append(out2)
        except:
            print('Stat Files for Model  %s not found!' % pol)
        try:
            out3 = pd.read_csv('masterdata_experiment_%s.csv' % pol)
            micro = micro.append(out3)
        except:
            print('Microdata for Model  %s not found!' % pol)

    total = total.set_index(['agebin', 'desc', 'subtype'])
    aggVal = aggVal.set_index(['model', 'index'])
    micro = micro.set_index(['model', 'id', 'age'])
    total.sort_index().to_stata('fthb_stats_final.dta')
    aggVal.sort_index().to_csv('fthb_aggVal_final.csv')
    micro.sort_index().to_stata('masterdata_final.dta')


# Collate partial statistics of model stats
def compMerge(compSpecList, binVar):
    from numpy import inf, nan
    for param, par_list in compSpecList:
        compStats = []
        compElas = []
        for par_val in par_list:
            par_suf = '_%s' % par_val if par_val != '' else par_val
            test = pd.read_csv('stats_final_experiment_%s%s.csv'
                               % (param, par_suf))
            stats = (test[(test['agebin'].isin([0, 20, 30, 40])) &
                          (test['subtype'] == 'id')].groupby('agebin').head(2))
            stats['compval'] = par_val
            stats = stats.pivot_table('value', ['agebin', 'compval'], 'desc')
            compStats.append(stats)

            try:
                test = pd.read_csv('stats_elasticities_experiment_%s%s.csv'
                                   % (param, par_suf))
                if isinstance(binVar, str):
                    stats = test.loc[(test[binVar].notnull()) &
                                     (test.notnull().sum(axis=1) == 4),
                                     ['desc', binVar, 'value']]
                    idVar = [binVar, 'compval']
                # Stopgap for assets/income observations
                elif isinstance(binVar, list):
                    stats = test.loc[(test.notnull().sum(axis=1) == 3 +
                                     len(binVar)) & (test[binVar[0]].notnull()),
                                     ['desc'] + binVar + ['value']]
                    idVar = binVar + ['compval']
                stats['compval'] = par_val
                stats = stats.pivot_table('value', idVar, 'desc')
                compElas.append(stats)
            except:
                pass
        final = (pd.concat(compStats).sort_index(0, ['agebin', 'compval'])
                 .replace([inf, -inf], nan))
        final.to_stata('fthb_%scomp_final.dta' % param)
	final = (pd.concat(compElas).sort_index(0, idVar)
                 .replace([inf, -inf], nan))
	final.to_stata('fthb_%sElas_final.dta' % param)
        try:
            (pd.concat(compElas).sort_index(0, idVar)
             .to_stata('fthb_%sElas_final.dta' % param))
        except:
            pass
