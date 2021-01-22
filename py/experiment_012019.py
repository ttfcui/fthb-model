#!/usr/bin/python

""" Title: experiment.py
    Purpose: Master file for running FTHB model subsidy experiments.
    Author: TC, DB
    Date: 2018-06-20

"""
import experiment_programs as run
import numpy as np
import subprocess

from model.model_iterate import lifecycle_iterate
from main_calibration import momOut_9mom, momOut_extension
import pandas as pd

# Baseline policy


def fthbMomentOut(mTitle, mDesc, outFile='altParams_6mom_032019.txt', **kwargs):
    moments, target = momOut_9mom('model/')
    momentsH, momentsExtra = momOut_extension('experiment_{}'.format(mTitle))
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    if 'extra' in kwargs:
        momentOut += ['{:6.4f},'.format(valH) for valH in momentsExtra.tolist()]
    subprocess.call(['echo', mDesc] + momentOut, stdout=open(outFile, 'a'))
    try:
        elasTitle = kwargs['elasTitle']
    except:
        elasTitle = mTitle
    subprocess.call(['stata', 'do', 'propensity_bins.do',
                     mTitle, elasTitle])


def fthbCalibCycle(param1, param1Desc, **kwargs):

    param1Lab = param1.keys()[0]
    try:
       param1Var = run.fthbParamsDict[param1Lab] 
       paramDict = run.fthbParamsDict
    except:
       param1Var = run.commonParamsDict[param1Lab] 
       paramDict = run.commonParamsDict
    param1Copy = param1Var
    for p1Val in param1.values()[0]:
            if 'param2' in kwargs:
	        param2Lab = kwargs['param2'].keys()[0]
                paramDict[param1Lab] = p1Val
                for p2Val in kwargs['param2'].values()[0]:
                    mFname = 'FTHB{}{:6.4f}{}'.format(param1Lab, p1Val, param2Lab)
                    try:
                        mFname += kwargs['hetero']
                    except:
                        pass
                    mFDesc = ('\'Ondown, {}={:6.4f}, {}={:6.4f}\','.format(
			      param1Desc, p1Val, kwargs['param2Desc'], p2Val))
	            run.compPols(('experiment_{}'.format(mFname),
				  param2Lab, [p2Val]), 'First-time', 1.0,
                                  share=1.00, tvalue=kwargs['tvalue'],
                                  gen_props=True, delModel=True)
                    fthbMomentOut('{}_{:6.4f}'.format(mFname, p2Val), mFDesc, **kwargs)
            else:
                mFname = 'FTHB{}'.format(param1Lab)
                try:
                    mFname += kwargs['hetero']
                except:
                    pass
                mFDesc = '\'Ondown, {}={:6.4f}\','.format(param1Desc, p1Val)
	        run.compPols(('experiment_{}'.format(mFname),
                              param1Lab, [p1Val]), 'First-time', 1.0,
                              share=1.00, tvalue=kwargs['tvalue'],
                              gen_props=True, delModel=True)
                fthbMomentOut('{}_{:6.4f}'.format(mFname, p1Val), mFDesc, **kwargs)

    paramDict[param1Lab] = param1Copy



run.simFTHB([('experiment_monetary_nodown', 0.12, 0.00)], 1.0,
            gen_micro=True, gen_props=True)
fthbMomentOut('monetary_nodown', '\'baseline (nodown)\',',
	      elasTitle='nodown')

run.simFTHB([('experiment_FTHBSize_0.075', 0.075, 0.00)], 1.0,
            gen_micro=True, gen_props=True)
fthbMomentOut('FTHBSize_0.075', '\'baseline (nodown) 5K\',',
	      elasTitle='FTHB5K')

run.simFTHB([('experiment_FTHBOnDownBase', 0.075, 1.00)], 1.0,
            gen_micro=True, gen_props=True)
fthbMomentOut('FTHBOnDownBase', '\'5K Bridge Loan (new baseline)\',',
	      elasTitle='ondown', extra=True)
exit()
"""
moments, target = momOut_6mom('model/', False)
momentsH, momentsExtra = momOut_extension('experiment_FTHBOnDownBase')
print(moments)
momentOut = ['{:6.4f},'.format(val) for val in moments.tolist()[0]]
momentOut += ['{:6.4f},'.format(val) for val in momentsH.tolist()]
momentOut += ['{:6.4f},'.format(valH) for valH in momentsExtra.tolist()]
subprocess.call(['echo', '\'5K Bridge Loan (new baseline)\','] + momentOut,
                stdout=open('altParams_6mom_022019.txt', 'a'))

subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBOnDownBase',
                 'ondown'])
exit()


# Different quadratic fixed costs
#F2params = [0.025, 0.05, 0.1, 0.25]
#fthbCalibCycle({'F2': F2params}, 'Quad FC', tvalue=0.075)

for val in F2params:

    run.compPols(('experiment_FTHBF2', 'F2', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True, share=1.00, tvalue=0.03)
    moments, target = momOut_6mom('model/', False)
    momentsH = momOut_extension('experiment_FTHBF2_{:6.4f}'.format(val))
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    momentOut += ['{:6.4f},'.format(valH) for valH in momentsH.tolist()]
    subprocess.call(['echo', '\'Bridge Loan, Quad FC={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom_022019.txt', 'a'))

    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBF2_{:6.4f}'.format(val),
                     'FTHBF2_{:6.4f}'.format(val)])
"""
# Different minimum house sizes
Dparams = [1.35, 1.5]
#fthbCalibCycle({'Dmin': Dparams}, 'Minimum D', tvalue=0.075, hetero='NoEndow')
run.fthbParamsDict['Dmin'] = 1.3500
run.fthbParamsDict['movprobR'] = 0.025
fthbCalibCycle({'elasticity': [2.0, 2.25, 2.5, 3, 5]}, '2.5% Renter moves, 1/EIS ', param2={'poismean': [5.0, 4.6, 4.2, 5.5]},
               param2Desc='Poissson Mean', tvalue=0.075, hetero='MovSm')
exit()
fthbCalibCycle({'Dmin': Dparams}, 'Minimum D', param2={'F2': [0.10, 0.20, 0.35]},
               param2Desc='Quadratic F', tvalue=0.075)
"""
for movval in [0.25, 0.30]:  # [0.02, 0.05, 0.10]:
        run.fthbParamsDict['movprobR'] = movval
	for val in Dparams:
	    run.compPols(('experiment_FTHBmov{:5.3f}D'.format(movval), 'Dmin', [val]),
                          'First-time', 1.0, share=1.00, tvalue=0.075,
                          gen_props=True, delModel=True)
	    moments, target = momOut_6mom('model/', False)
	    momentsH = momOut_extension('experiment_FTHBmov{:5.3f}D_{:6.4f}'.format(movval, val))
	    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
	    momentOut += ['{:6.4f},'.format(valH) for valH in momentsH.tolist()]
	    subprocess.call(['echo', '\'Ondown 5K, Renter Move={:5.3f}, Minimum D={}\','.format(movval, val)] + momentOut,
			    stdout=open('altParams_6mom_022019.txt', 'a'))
	    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBmov{:5.3f}D_{:6.4f}'.format(movval, val),
			     'FTHBmov{:5.3f}D_{:6.4f}'.format(movval, val)])
"""

run.fthbParamsDict['Dmin'] = 1.500
run.fthbParamsDict['movprobR'] = 0.1000
fthbCalibCycle({'rentUtil': [0.90, 0.93, 0.85]}, 'Min D=1.5, 20% Renter Move, rentutil',
               param2={'poismean': [0.0]}, param2Desc='Poisson Mean',
               tvalue=0.075, hetero='')
exit()
"""
run.fthbParamsDict['movprobR'] = 0.300
for F2val in [0.025, 0.04, 0.05]:
        run.fthbParamsDict['F2'] = F2val
	for val in Dparams:
	    run.compPols(('experiment_FTHBMovF2{:5.3f}D'.format(F2val), 'Dmin', [val]),
                          'First-time', 1.0, share=1.00, tvalue=0.075,
                          gen_props=True, delModel=True)
	    moments, target = momOut_6mom('model/', False)
	    momentsH = momOut_extension('experiment_FTHBMovF2{:5.3f}D_{:6.4f}'.format(F2val, val))
	    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
	    momentOut += ['{:6.4f},'.format(valH) for valH in momentsH.tolist()]
	    subprocess.call(['echo', '\'Ondown, Renter move =30\%, Quad Fixed={:5.3f}, Minimum D={}\','.format(F2val, val)] + momentOut,
			    stdout=open('altParams_6mom_022019.txt', 'a'))
	    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBMovF2{:5.3f}D_{:6.4f}'.format(F2val, val),
			     'FTHBMovF2{:5.3f}D_{:6.4f}'.format(F2val, val)])
exit()

run.fthbParamsDict['rentUtil'] = 0.85
run.fthbParamsDict['F2'] = 0.1
for movval in [0.02, 0.05, 0.10]:  # [0.06, 0.08, 0.10]:
        run.fthbParamsDict['movprobR'] = movval
	for val in Dparams:
	    run.compPols(('experiment_FTHBUDmov{:5.3f}D'.format(movval), 'Dmin', [val]),
                          'First-time', 1.0, share=1.00, tvalue=0.075,
                          gen_props=True, delModel=True)
	    moments, target = momOut_6mom('model/', False)
	    momentsH = momOut_extension('experiment_FTHBUDmov{:5.3f}D_{:6.4f}'.format(movval, val))
	    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
	    momentOut += ['{:6.4f},'.format(valH) for valH in momentsH.tolist()]
	    subprocess.call(['echo', '\'Ondown 5K, Quad cost=0.1, Rent Disutil=0.85, Renter Move={:5.3f}, Minimum D={}\','.format(movval, val)] + momentOut,
			    stdout=open('altParams_6mom_022019.txt', 'a'))
	    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBUDmov{:5.3f}D_{:6.4f}'.format(movval, val),
			     'FTHBUDmov{:5.3f}D_{:6.4f}'.format(movval, val)])

exit()

# Different minimum house sizes, with different fcost share
run.fthbParamsDict['pr'] = 0.20
Dparams = [1.1, 1.2, 1.35, 1.5]
for val in Dparams:
    run.compPols(('experiment_FTHBAltPrD', 'Dmin', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True, tvalue=0.03)
    moments, target = momOut_6mom('model/', False)
    momentsH = momOut_extension('experiment_FTHBAltPrD_{:6.4f}'.format(val))
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    momentOut += ['{:6.4f},'.format(valH) for valH in momentsH.tolist()]
    subprocess.call(['echo', '\'Buyer=80\% Fcost share, Minimum D={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom_022019.txt', 'a'))
    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBAltPrD_{:6.4f}'.format(val),
                     'FTHBAltPrD_{:6.4f}'.format(val)])
run.fthbParamsDict['Dmin'] = 1.1
exit()
"""



## OLD COUNTERFACTUAL EXPERIMENTS

"""
# Different quadratic FC + linear cost
run.fthbParamsDict['F'] = 0.01
for val in F2params:
    run.compPols(('experiment_FTHBFLowF2', 'F2', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True)
    moments, target = momOut_6mom('model/', False)
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    subprocess.call(['echo', '\'F=0, Quad FC={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom.txt', 'a'))
    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBFLowF2_{:6.4f}'.format(val),
                     'FTHBFLowF2_{:6.4f}'.format(val)])

run.fthbParamsDict['F'] = 0.03
for val in F2params:
    run.compPols(('experiment_FTHBFHighF2', 'F2', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True)
    moments, target = momOut_6mom('model/', False)
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    subprocess.call(['echo', '\'F=0.03, Quad FC={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom.txt', 'a'))
    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBFHighF2_{:6.4f}'.format(val),
                     'FTHBFHighF2_{:6.4f}'.format(val)])

# Different quadratic FC + linear cost
run.fthbParamsDict['F'] = 0.01
run.fthbParamsDict['Dmin'] = 1.5
for val in F2params:
    run.compPols(('experiment_FTHBFLowDHiF2', 'F2', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True)
    moments, target = momOut_6mom('model/', False)
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    subprocess.call(['echo', '\'F=0, D=1.2 Quad FC={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom.txt', 'a'))
    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBFLowDHiF2_{:6.4f}'.format(val),
                     'FTHBFLowDHiF2_{:6.4f}'.format(val)])

run.fthbParamsDict['F'] = 0.06
for val in F2params:
    run.compPols(('experiment_FTHBDHiF2', 'F2', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True)
    moments, target = momOut_6mom('model/', False)
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    subprocess.call(['echo', '\'D=1.2, Quad FC={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom.txt', 'a'))
    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBDHiF2_{:6.4f}'.format(val),
                     'FTHBDHiF2_{:6.4f}'.format(val)])
run.fthbParamsDict['Dmin'] = 1.1

# Different quadratic FC + EIS
run.fthbParamsDict['elasticity'] = 3

for val in F2params:
    run.compPols(('experiment_FTHBGammaHiF2', 'F2', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True)
    moments, target = momOut_6mom('model/', False)
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    subprocess.call(['echo', '\'Gamma=3, Quad FC={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom.txt', 'a'))
    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBGammaHiF2_{:6.4f}'.format(val),
                     'FTHBGammaHiF2_{:6.4f}'.format(val)])

run.fthbParamsDict['elasticity'] = 5
for val in F2params:
    run.compPols(('experiment_FTHBGammaVHiF2', 'F2', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True)
    moments, target = momOut_6mom('model/', False)
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    subprocess.call(['echo', '\'Gamma=5, Quad FC={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom.txt', 'a'))
    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBGammaVHiF2_{:6.4f}'.format(val),
                     'FTHBGammaVHiF2_{:6.4f}'.format(val)])

run.fthbParamsDict['elasticity'] = 2.25

# Different property taxes
run.fthbParamsDict['dtau'] = 0.025

for val in F2params:
    run.compPols(('experiment_FTHBTaxHiF2', 'F2', [val]), 'First-time', 1.0,
                 gen_props=True, delModel=True)
    moments, target = momOut_6mom('model/', False)
    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
    subprocess.call(['echo', '\'dtau=2.5%, Quad FC={}\','.format(val)] + momentOut,
                    stdout=open('altParams_6mom.txt', 'a'))
    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBTaxHiF2_{:6.4f}'.format(val),
                     'FTHBTaxHiF2_{:6.4f}'.format(val)])

run.fthbParamsDict['dtau'] = 0.01


# Much higher share of adjustments due to exogenous moving shocks

run.fthbParamsDict['F'] = 0.045
run.fthbParamsDict['pr'] = 0.0
for mval in [0.25]:  # [0.06, 0.08, 0.10]:
        run.fthbParamsDict['movprob'] = mval
	for val in F2params:
	    run.compPols(('experiment_FTHBmov{:5.3f}v2F2'.format(mval), 'F2', [val]),
                          'First-time', 1.0, gen_props=True, delModel=True)
	    moments, target = momOut_6mom('model/', False)
	    momentOut = ['{:6.4f},'.format(vals) for vals in moments.tolist()[0]]
	    subprocess.call(['echo', '\'F=4.5%, movprob={:5.3f}, Quad FC={}\','.format(mval, val)] + momentOut,
			    stdout=open('altParams_6mom.txt', 'a'))
	    subprocess.call(['stata', 'do', 'propensity_bins.do', 'FTHBmov{:5.3f}v2F2_{:6.4f}'.format(mval, val),
			     'FTHBmov{:5.3f}v2F2_{:6.4f}'.format(mval, val)])


"""
