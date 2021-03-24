#!/usr/bin/python

""" Title: experiment_programs.py
    Purpose: Master file for running FTHB model subsidy experiments.
    Author: TC, DB
    Date: 2018-06-20

"""

import pandas as pd
from model.model_iterate import lifecycle_iterate
from subprocess import call
from os import remove
from shutil import rmtree

# Parameters of interest in model (preset or calibrated)

algDict = {'agridsize': 135, 'Dgridsize': 35, 'zgridsize': 17,
            'hpnodes(1)': 0.0599}
commonParamsDict = {'beta2': 0.925, 'elasticity': 3.0, 'unempprob': 0.000,
                    'psi': 0.5, 'ret_wealth': 0.00, 'tax1': 0.175, 'F': 0.06,
                     'r': 0.03,  'sigma_temp': 0.041, 'rentUtil': 1.00,
                     'balancer': -0.0000, 'numhouseholds': 30000,
                     'ss_only': '.FALSE.', 'print_micro': '.TRUE.'}
fthbParamsDict = {'EligYrsR': 99, 'EligYrsF': 3, 'delta': 0.022, 'movprob': 0.015,
                   'F2': 0.43, 'movprobR': 0.000, 'thetamatlab': 0.2, 'elasticity2': 0.759,
                   'rentPrem': 0.0065, 'Dmin': 1.1, 'polEnd': 1, 'dtau': 0.01, 
                   'numhouseholds': 30000}
carsParamsDict = {'EligYrsR': 3, 'EligYrsF': 99, 'rentPrem': 0.00, 'Dmin': 0.05,
                   'F2': 0.0, 'elasticity2': 0.90, 'delta': 0.125, 'rentPremRetire': 0.00,
                   'scrapped': 1.0, 'adjTransfer': 0.067, 'Eta_transfer': 1.00,
                   'polEnd': 1, 'hptransLength': 1,
                   'numhouseholds': 6000, 'custom_start': '.TRUE.'}

def genParams(polParamsDict, polDict, **kwargs):
    """ Parameter dict wrapper - appends calibrated parameters to the
        inputted policy parameters, plus additional settings.

    Args:
        polParamsDict: A dictionary of calibrated parameters set for
            each model environment (i.e. CARS vs. FTHB)
        polDict: the inputted policy parameters, governing the size
            of the subsidy, application to the down, etc.
        **percentages: If true, the values in polDict is in percentage
            change from house value, not levels.
        **aggPol: If true, turns on an aggregate shock that also kicks in
            when the temporary policy does.
    """

    modelDict = algDict.copy()
    if 'percentages' in kwargs:
        modelDict.update({'pctageflag': '.TRUE.'})
    if 'aggPol' in kwargs:
        modelDict.update({'agg_policies': '.TRUE.'})
    if 'loadPrices' in kwargs:
        modelDict.update({'loadpricepath': '.TRUE.'})
    modelDict.update(commonParamsDict)
    modelDict.update(polParamsDict)
    modelDict.update(polDict)
    return modelDict


def simFTHB(polSpecList, maint, **kwargs):
    """ A wrapper for FTHB-like simulated policies, where temporary subsidies
        stimulate current renters to accelerate their homeownership
        decisions.

    Args:
        polSpecList: A list of triple tuples, each corresponding to one
            policy simulation. The tuple elements respectively mark the name
            associated with the sim, the subsidy level and the share of subsidy
            applicable to the down (in [0,1]).
        maint: Dummy variable, taking 1 if the model imposes automatic
            maintenance of the durable and 0 if the durable good depreciates.
        **stop_model: If True, do not run the model again and instead
            just run the data generation scripts on the available output.
        **gen_elas: Runs the elasticity generation script given
            in elasticities.py.
        **gen_micro: Runs the detailed microdata script given
            in gen_microdata.py.
        **gen_props: Runs the takeup propensity script given in
            gen_propensity.py.
    """
    for model, value, prop in polSpecList:
        modelDict = genParams(fthbParamsDict, {'adjTransfer': value,
                              'Eta_transfer': prop, 'maint': maint,}, **kwargs)

        if 'transLen' in kwargs:
            modelDict.update({'hptransLength': kwargs['transLen'], 'pe_start': '.FALSE.', 'ge_start': '.TRUE.'})

        statsList = ['%f' % modelDict['delta'], '%f' % modelDict['Dmin'],
		     '%f' % modelDict['polEnd']]
        pyArgs = [model, 'First-time', statsList[0]]
        propsList = ['%f' % modelDict['thetamatlab'], '%f' % modelDict['hpnodes(1)']]
        if 'transLen' in kwargs:
            statsList.append('%f' % modelDict['hptransLength'])
            pyArgs.append('%f' % modelDict['hptransLength'])
        modelPol = lifecycle_iterate(modelDict)
        if 'stop_model' not in kwargs:
            modelPol.execSh(model=model)
        else:
            modelPol.execSh(model=model, stop_model=True)

        try:
            modelPol.matlabPlot(model=[model], end=kwargs['transLen'])
        except:
            modelPol.matlabPlot(model=[model])

        if value == 0.0 and 'override' not in kwargs:
            print('No Policy Simulated')
            continue

        call(['python', 'gen_stats.py', model, 'First-time'] + statsList)
        if 'gen_elas' in kwargs:
            call(['python', 'elasticities.py'] + pyArgs)
        if 'gen_micro' in kwargs:
            call(['python', 'gen_microdata.py'] + pyArgs)
        if 'gen_props' in kwargs:
            call(['python', 'gen_propensity.py', model, 'First-time'] + propsList)


# CARS simulation template
def simCARS(polSpecList, maint, keepTransition=False, **kwargs):
    """ A wrapper for CARS-like simulated policies, where temporary subsidies
        stimulate existing car owners to accelerate their car replacement
        decisions.

    Args:
        polSpecList: A list of quadruple tuples, each corresponding to one
            policy simulation. The tuple elements respectively mark the name
            associated with the sim; the subsidy level; the down payment
            requirements for a new car; and the common scrap value of selling
            the car, as in Adda and Cooper (2000).
        maint: Dummy variable, taking 1 if the model imposes automatic
            maintenance of the durable and 0 if the durable good depreciates.
        **kwargs: Identical to the kwargs in simFTHB.
    """

    for model, valueT, valueD, scrapV in polSpecList:
        modelDict = genParams(carsParamsDict, {'thetamatlab': valueT,
                              'downflag': valueD, 'scrapvalue': scrapV,
                              'maint': maint}, **kwargs)

        if 'transLen' in kwargs:
            modelDict.update({'hptransLength': kwargs['transLen'], 'pe_start': '.FALSE.', 'ge_start': '.TRUE.'})
        statsList = ['%f' % modelDict['delta'], '%f' % modelDict['Dmin'],
		     '%f' % modelDict['polEnd']]
        pyArgs = [model, 'Repeat', statsList[0]]
        propsList = ['%f' % modelDict['thetamatlab'], '%f' % modelDict['hpnodes(1)']]
        if 'transLen' in kwargs:
            statsList.append('%f' % modelDict['hptransLength'])
            pyArgs.append('%f' % modelDict['hptransLength'])
        modelPol = lifecycle_iterate(modelDict)

        if 'stop_model' not in kwargs:
            modelPol.execSh(model=model)
        else:
            modelPol.execSh(model=model, stop_model=True)

        try:
            modelPol.matlabPlot(pol='Return', model=[model], end=kwargs['transLen'])
        except:
            modelPol.matlabPlot(pol='Return', model=[model])

        if modelDict['adjTransfer'] == 0.0:
            print('No Policy Simulated')
            continue

        call(['python', 'gen_stats.py', model, 'Repeat'] + statsList)
        if 'gen_elas' in kwargs:
            call(['python', 'elasticities.py'] + pyArgs)
        if 'gen_micro' in kwargs:
            call(['python', 'gen_microdata.py'] + pyArgs)
        if 'gen_props' in kwargs:
            call(['python', 'gen_propensity.py', model, 'Repeat'] + propsList)
        # Disk space management (this transition file is massive)
        if not keepTransition:
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
        else:
            modelPol.execSh(model=model, stop_model=True)
        modelPol.matlabPlot(model=[model], life=transLen)
        aggs = pd.read_table('model/transition_aggs.txt', sep='\s+',
                             header=None, names=['price', 'hstock', 'rents',
                                                 'resales', 'newhomes'])
        aggs.to_csv('output/%s/transition_aggs.csv' % model, index=False)
        call(['python', 'gen_stats.py', model, 'First-time',
              '0.022', '%f' % fthbParamsDict['Dmin']])
        if 'gen_micro' in kwargs:
            call(['python', 'gen_microdata.py', model, 'First-time', '0.022'])


# Simple comparative statistics checks
def compPols(polSpec, polType, maint, share=0.00, tvalue=0.075, **kwargs):
    """ Modified function allowing for iteration of FTHB or CARS-like policies
        over various values of a single parameter.

    Args:
        polSpec: A triple tuple. The elements respectively denote the prefix
            common to all model iterations; the parameter being iterated; and
            the list of values in the loop.
        polType: Accepts 'First-time' for FTHB models and 'Repeat' for
            CARS policies.
        maint: Dummy variable, taking 1 if the model imposes automatic
            maintenance of the durable and 0 if the durable good depreciates.
        share: The share of subsidy applicable to the down (in [0,1]).
        tvalue: The subsidy level.
        **kwargs: Identical to the kwargs in simFTHB.
    """

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
        call(['python', 'elasticities.py', model, polType] + [statsList[0]])
        call(['python', 'gen_stats.py', model, polType] + statsList)
        if 'gen_micro' in kwargs:
            call(['python', 'gen_microdata.py', model, polType] + [statsList[0]])
        if 'gen_props' in kwargs:
            call(['python', 'gen_propensity.py', model, polType] + propsList)
        if 'delModel' in kwargs:
            rmtree(modelPol.mainD + 'output/%s/' % model)


def polMerge(models):
    """ Output relevant statistics for each policy in a joint Stata file.

    Args:
        models: List of strings marking model names to be appended
            together. Use the names on the directories in the output folder.
    """
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
