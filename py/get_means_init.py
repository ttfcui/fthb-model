
from model_iterate import lifecycle_iterate
from numpy import log, linspace, power
from subprocess import call
import pandas as pd
import os
baseFile = 'lifecycle_TC05022016.f90'

grid1 = linspace(0.0530, 0.0710, 13)
grid2 = linspace(0.80, 0.87, 8)
for x, y in [(x, y) for x in grid1 for y in grid2]:
    model = lifecycle_iterate({'r_rentalO': x, 'elasticity2O': y,
                               'numhouseholds': 5000, 'EligYrsR': 66.0},
                              baseFile)
    name = 'r%.4fa%.2f' % (x, y) 
    model.execSh(model=name)
    # Removing the largest output files
    for fname in ['fthb', 'transition_fthb', 'householdresults']:  
        os.remove('output/%s/%s_%s.txt' % (name, fname, name))


model = lifecycle_iterate({}, baseFile)
model.appendModels('lifecycleprofiles', model=['r%.4fa%.2f' % (x, y)
	for x in grid1 for y in grid2])


momH = model.genMoments(5, ['avgHouseOwners', 'avgAssetsOwners',
        'avgAssetsExDebtOwners'], 'NOwn')
mom  = model.genMoments(5, ['avgAssets', 'avgAssetsExDebt', 'fracOwn'])

scf1 = pd.read_csv('output/graph_data_byageind.csv')[['netliquid_exhousingdebt', 'homeownershiprate']]
scf2 = pd.read_csv('output/graph_data_byageind_homeowners.csv')[['housePos']]

scf1.index.names = ['ageBin']; scf1.index = 15 + 5*scf1.index
scf2.index.names = ['ageBin']; scf2.index = 15 + 5*scf2.index

quad = lambda x: power(x, 2)
calibrate = (mom[['avgAssetsExDebt', 'fracOwn']].merge
	(momH[['avgHouseOwners']], left_index=True, right_index=True).merge
	(scf1, left_index=True, right_index=True).merge
	(scf2, left_index=True, right_index=True))
calibrate['cal1'] = ((
	calibrate['avgAssetsExDebt'] - calibrate['netliquid_exhousingdebt']).
	map(quad))
calibrate['cal2'] = 10*((
	calibrate['fracOwn'] - calibrate['homeownershiprate']).map(quad))
calibrate['cal3'] = ((
	calibrate['avgHouseOwners'] - calibrate['housePos']).map(quad))

idx = pd.IndexSlice
opt = calibrate.loc[idx[:,20:50],:].groupby(level=0)[['cal2', 'cal3']
	].sum().sum(axis=1).idxmin()
print(opt) # Result from grid search

momH.loc[opt].to_csv('momentsH_5yearbin.csv')
mom.loc[opt].to_csv('moments_5yearbin.csv')