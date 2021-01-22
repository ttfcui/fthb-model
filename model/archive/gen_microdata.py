#!/usr/bin/python

""" Title: xp_03012017.py
    Purpose: Creating a panel simulation dataset with many variables.
    Author: TC, DB
    Date: 2016-12-10

"""

""" Summary: This file creates statistics grouped into "three margins,"
    as well as homeownership rates during the transition and
    substitution effects during policy takeup. It copies code
    from other py files.
"""
from model.model_iterate import lifecycle_iterate
from sys import argv
import numpy as np
from pandas import merge, DataFrame, read_table, concat, to_numeric
# Useful functions for data processing
def ageTrim(dat, end=39):
    try:
        return dat.appended.loc[dat.appended['age'] <= end, :]
    except:
        return dat.appended.loc[dat.appended['ageBought'] <= end, :]


readfunc = lifecycle_iterate.appendModels
# Call relevant model output files
view1 = ageTrim(readfunc('transition_adjust', model=[argv[1]]))
view2 = ageTrim(readfunc('lifecycle_adjust', model=[argv[1]]))
view3 = (ageTrim(readfunc('fthb', model=[argv[1]]))
         [['id', 'age', 'nextAssets', 'assets', 'consumption', 'durTime', 'rent',
         'income_val', 'income', 'income_l', 'income_val_l',
         'durables', 'nextDurables', 'welfare']])
view3 = view3.sort_values(['id', 'age'])
view4 = ageTrim(readfunc('transition_fthb', model=[argv[1]]))
view5 = ageTrim(readfunc('dist_fthb', model=[argv[1]]))

def merge1():
    print('Dataset #1: \n \n')
    merged = merge(view1, view2, on=['id', 'age', 'model'],
                   suffixes=('_pol', '_ss'))
    view1.info()
    view2.info()
    merged.info()
    merged['adjDiff'] = merged['adjustment_pol'] - merged['adjustment_ss']
    print(merged.loc[merged['adjDiff'] != 0].describe())
    print(merged.loc[merged['adjDiff'] == 0].describe())
    return merged

def merge2(mdir='model', tretire=39, end=39):
    print('Dataset #2: \n \n')
    data = read_table('%s/transition_difffull.txt' % mdir,
                      sep='\s+', engine='python', header=None)
    cols = data.columns.tolist()
    cols[0:2] = ['id', 'age']
    data.columns = cols
    data['age'] += 1  # to reconcile age with other datasets
    data = data.loc[(data['age'] <= end)]

    pullforward = data.loc[data[2] == 1]
    pullforward.info()
    agent, dates = np.where(pullforward == -1)
    pullforward.loc[:, 'pullforward'] = np.nan
    pullforward.loc[:, 'pullforward'].iloc[np.unique(agent)] = (
        pullforward.iloc[np.unique(agent), :-1].idxmin(axis=1) - 2)

    pullforward_ext = (pullforward.loc[pullforward['pullforward'].isnull()]
                       .drop([2], axis=1))
    pullforward.loc[pullforward['pullforward'].isnull(),
                    'pullforward'] = tretire + 1

    for i in xrange(tretire, 2, -1):
        data.iloc[:, i] = data.iloc[:, 2:i+1].sum(axis=1)

    posext = data.loc[data[tretire] > 0]
    posext = merge(posext, pullforward[['pullforward']], left_index=True,
                   right_index=True, how='left')
    posext.loc[:, 'posext'] = np.nan
    for row, rowval in posext.iterrows():
       indices = np.where(rowval[2:] == 1)[0]
       posext.loc[row, 'posext'] = indices[np.argmax(indices > posext.loc
                                                     [row, 'pullforward'])]
    negext = data.loc[data[tretire] < 0]
    negext = merge(negext, pullforward[['pullforward']], left_index=True,
                   right_index=True, how='left')
    negext.loc[:, 'negext'] = ((negext == -1).idxmax(axis=1)) - 2

    merged = merge(pullforward[['pullforward']], posext[['posext']],
                   how='outer', left_index=True, right_index=True)
    merged = merge(merged, negext[['negext']], how='outer',
                   left_index=True, right_index=True)
    merged = merge(data[['id', 'age']], merged,
                   left_index=True, right_index=True)
    print(merged.describe())
    return merged


def merge3(polType, dep):
    print('Dataset #3: \n \n')
    # Future income shocks identical between SS and policy, for each id
    for t in xrange(1, 5):
        view3['income_valf%d' % t] = view3.groupby('id')['income_val'].shift(-1*t)
        view3['income_f%d' % t] = view3.groupby('id')['income'].shift(-1*t)
    incHolder = DataFrame()
    for j in xrange(2, 40):
        grouped = (view3.loc[(view3['age'] < 40) & (view3['age'] >= j)]
                   .groupby('id'))
        inc_data = concat([grouped['age'].first(),
                           grouped['income_val'].mean(),
                           grouped['income_val'].std()], axis=1).reset_index()
        inc_data.columns = ['id', 'age', 'incF_mean', 'incF_std']
        incHolder = incHolder.append(inc_data)

    view4['policy'] = view4['ageBought'] - view4['age']
    # Impute original value of owned durable using vintage stats
    if polType == 'Repeat':
        view4['durablesOrig'] = ((1.0/(1.0-dep))**(
                                 view4['ageBought'] - view4['durTime'])*(
                                 view4['durables']))
    # Some imputation issues means unreasonably high values filtered out
    # (but this only applies to the repeat buyer case)
        cap = (view4.loc[view4['PolTaken'] == polType, ['nextDurables']]
               .describe().loc['max', 'nextDurables'])
        durBool = (view4['durablesOrig'] <= cap) & (view4['policy'] == 0)
        durCols = ['durablesOrig', 'durTime']
    else:
        durBool = (view4['policy'] == 0)
        durCols = ['durTime']

    durShift = view4.loc[durBool & (view4['PolTaken'] == polType),
                         ['id', 'ageBought', 'nextDurables', 'nextAssets',
                          'consumption', 'welfare'] + durCols]
    durSS = view5.loc[view5['PolTaken'] == polType, ['id', 'ageBought',
                      'consumption']]

    # Merge and collapse
    merged = merge(durShift, durSS, on=['id', 'ageBought'], how='left',
                   suffixes=('', '_planned'))
    merged = merge(merged, view3, suffixes=('', '_ss'),
                   left_on=['id', 'ageBought'], right_on=['id', 'age'])
    merged = merge(merged, incHolder, on=['id', 'age'])
    merged['marginal'] = 0
    merged.loc[merged['consumption_planned'].isnull(), 'marginal'] = 1
    merged = merged.drop(['ageBought', 'consumption_planned', 'durTime_ss'],
                         axis=1)
    if polType == 'Repeat':
        merged['durablesGap'] = (merged['durablesOrig']
                                 - merged['nextDurables'])
    merged['model'] = argv[1]
    return merged


def merge4(view, futVar):
    print('Dataset #4, variable %s: \n \n' % futVar)
    view = view.loc[view['variable'] == futVar]
    view['ageNext'] = view['age'] + 5
    merged = merge(view, view, left_on=['id', 'ageNext', 'variable', 'model'],
                   right_on=['id', 'age', 'variable', 'model'], how='left',
                   suffixes=('', '_cf'))
    merged = merged.iloc[:, 0:23]  # We want lags, not leads, of counterfactual
    colNames = merged.columns.tolist()
    colNames[-5:] = ['var_pol1_cf', 'var_pol2_cf', 'var_pol3_cf', 'var_pol4_cf', 'var_pol5_cf']
    colNames = [col.replace('var_', '%s_' % futVar) for col in colNames]
    merged.columns = colNames
    futCols = (['%s_pol%s' % (futVar, i) for i in xrange(0, 6)] +
               ['%s_pol%s_cf' % (futVar, i) for i in xrange(1, 6)])
    merged = merged[['id', 'age'] + futCols]
    return merged

if __name__ == '__main__':
    print('Creating master data file: \n \n')
    merged1 = merge1()
    merged2 = merge2()
    merged3 = merge3(argv[2], float(argv[3]))

    merged_adj = merge(merged1, merged2, how='left', on=['id', 'age'])
    merged = merge(merged3, merged_adj, how='left', on=['id', 'age', 'model'])
    for varname in ['C', 'H_own', 'H_rent', 'V']:
         lagsleads = (lifecycle_iterate.readModel('transition_lagsleads')
                      .appended.sort_values(['id', 'age', 'variable']))
         floatCol = lagsleads.columns[~lagsleads.columns.isin(['variable', 'model'])].tolist()
         lagsleads.loc[:,floatCol] = lagsleads.loc[:,floatCol].apply(to_numeric, errors='coerce')
         merged4 = merge4(lagsleads, varname)
         merged = merge(merged, merged4, how='left', on=['id', 'age'])
    merged = merged.set_index(['id', 'age'])
    merged.info()
    merged.to_csv('masterdata_%s.csv' % argv[1])
