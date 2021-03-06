#!/usr/bin/python

""" Title: xp_12202016.py
    Purpose: Data cleaning for analysis in xp_12202016.do.
    Author: TC, DB
    Date: 2016-12-10

"""

""" Summary: This file is split into two main tasks:
    1) Process a file that tracks variable realizations around
    policy take-up;
    2) Create a dataset for policy-induced agents that has their
    house purchased with the subsidy, their counterfactual
    housing decision without subsidy, and other covariates.
"""
from model.model_iterate import lifecycle_iterate
import numpy as np
from pandas import merge, melt, DataFrame, concat
from sys import argv


# Useful functions for data processing
def genbins(data, col='age', bval=40):
    try:
        prange = np.linspace(0.0, 100.0, bval)
        bins = np.percentile(data.loc[:, col], prange)
    except:
        bins = bval
    data.loc[:, '%sbin' % col] = np.digitize(data[col], bins)
    j = 1
    for val in bins:
        data.loc[data['%sbin' % col] == j, '%sbin' % col] = val
        j += 1

# Steady-state response to income shocks
view = lifecycle_iterate.readModel('fthb').appended.sort_values(['id', 'age'])
# "nextDurables" in original file also records rental housing size
view['ownDur'] = (1.0 - view['rent'])*view['nextDurables']
# Rent in original file records decision to rent in current period, not
# existing status
view['renter'] = view.groupby('id')['rent'].shift(1)
view['time'] = 0

# Lags and leads are iteratively created for the variables.
lagleadv = {'assets': 'nextAssets', 'cons': 'consumption', 'hh': 'ownDur'}
lagleads = {'_f1': -1, '_f2': -2, '_f3': -3, '_f4': -4, '_l1': 1}

for suf, item in lagleads.iteritems():
    for key in lagleadv:
        view.loc[:, '%s%s' % (key, suf)] = (
            view.groupby('id')[lagleadv[key]].shift(item))
    view.loc[:, 'time%s' % suf] = -1*item

view = view.loc[view['age'].isin([5, 15, 25, 30, 35])]

# Unlike policy response output, can afford to stack groups getting both
# positive and negative shock here. Counterfactual are agents who did not
# get shocks in the event period (though they could get more in the future)
out = {}
consistency_bool = ((view['income_l'] == view['income_l2']) &
                    (view['renter'] == 1))
out['posshock'] = (view[(view['income'] > view['income_l']) &
                        consistency_bool].groupby(['age', 'income_l']).median())
out['negshock'] = (view[(view['income'] < view['income_l']) &
                        consistency_bool].groupby(['age', 'income_l']).median())
out['nof1shock'] = (view[(view['income'] == view['income_l']) &
                         consistency_bool].groupby(['age', 'income_l']).median())

# Unlike the output with the policy, some reshaping here to convert the
# variables created in wide form into long form
final = concat(out).filter(regex=r'[hs]_[fl]|cons|nextA|ownD|time', axis=1)
# The tuples are a fixed family of column values, so if number of variables
# examined changes this must be changed as well. Logic of the construction
# should be the same - just need to change values.
tuples = [tuple(range(25, 21, -1)), (5, 0, 1, 4), tuple(range(9, 5, -1)),
          tuple(range(13, 9, -1)), tuple(range(17, 13, -1)),
          tuple(range(21, 17, -1))]
final = concat(DataFrame({'time': final.iloc[:, h], 'assets': final.iloc[:, i],
                          'cons': final.iloc[:, j], 'housing': final.iloc[:, k]
                          }) for h, i, j, k in tuples)
# Output
final.to_csv('stationary_eventstudy.csv')
# %>

# Task 1) %<
""" Step 1: Read in an output file tracking lags and leads of variables.
    Separate columns indicate observations of the same variable at
    separate times.
"""  # %<
view = (lifecycle_iterate.readModel('transition_lagsleads')
        .appended.sort(['id', 'age', 'variable']))
# %>

""" Step 2: Clean the output file so
    1) Only look at agents who were incentivized to buy with the policy,
    without any additional income shocks in the short term;
    2) statistics can be binned on the agent's income state when
    receiving the policy;
    3) Merge counterfactual statistics for agents in a world without
    the policy (i.e. their profile iterated forward when
    creating the steady-state)
"""  # %<
# Boolean for 1); "incentivized to buy" means "still not owning house
# in the steady-state counterfactual"
view = view.set_index(['id', 'age'])
incs = view.loc[view['variable'] == 'Y', ['var_ss', 'var_l1', 'var_pol1']]
incs.columns = ['income_state', 'income_stateL', 'income_stateF']
wealth = DataFrame(view.loc[view['variable'] == 'Q', 'var_l1'] - 0.8*(
                   view.loc[view['variable'] == 'H_own', 'var_l1']))
wealth.columns = ['finwealth']
view = merge(view, incs, left_index=True, right_index=True)
view = merge(view, wealth, left_index=True, right_index=True).reset_index()

# The third condition can be adjusted: use == (constant income state)
# or > (positive income shock when policy takeup occurs)
try:
    if argv[1].lower().startswith('pos'):
        income_bool = (view['income_state'] > view['income_stateL'])
        affix = '_posshock'
    elif argv[1].lower().startswith('neg'):
        income_bool = (view['income_state'] < view['income_stateL'])
        affix = '_negshock'
    else:
        income_bool = (view['income_state'] == view['income_stateL'])
        affix = ''
except:
    income_bool = (view['income_state'] == view['income_stateL'])
    affix = ''

try:
    if argv[2].lower().startswith('infra'):
        marg_bool = (view['var_ss'] > 0.0)
        affix += '_inframarg'
    else:
        marg_bool = (view['var_ss'] == 0.0)
except:
    marg_bool = (view['var_ss'] == 0.0)

switch_bool = ((view['var_ss'] != view['var_pol0']) & marg_bool
               & income_bool & (view['income_state'] == view['income_stateF'])
               & (view['variable'] == 'H_own'))
switch = view.loc[switch_bool, ['id', 'age']]
view2 = merge(view, switch, on=['id', 'age'])
view2 = view2.drop(['income_stateL', 'income_stateF'], axis=1)
# %>

""" Step 3: Merge counterfactual statistics, relabel to prevent confusion.
    The logic here is that the counterfactual for an agent 3 years out is
    the agent with the same id, but 3 years older. The lags of that
    agent who is 3 years older fill in the gap.
    The same cannot be said of the counterfactual 9 years out, so that
    is dealt with separately.
"""  # %<

view2['ageNext'] = view2['age'] + 3
view3 = merge(view2, view, left_on=['id', 'ageNext', 'variable', 'model'],
              right_on=['id', 'age', 'variable', 'model'], how='left',
              suffixes=('', '_cf'))
view3 = view3.iloc[:, 0:19]  # We want lags, not leads, of counterfactual
view2['ageNext'] = view2['age'] + 9
view3Alt = merge(view2, view, left_on=['id', 'ageNext', 'variable', 'model'],
                 right_on=['id', 'age', 'variable', 'model'], how='left',
                 suffixes=('', '_cf')).loc[:, ['id', 'var_ss_cf', 'variable',
                                               'model', 'age_cf']]
view3Alt['age'] = view3Alt['age_cf'] - 9
view3 = merge(view3, view3Alt, on=['id', 'age', 'variable', 'model'],
              how='left', suffixes=('', '_cfAlt'))
view3 = view3.iloc[:, 0:-1]
# Sort fin. wealth into bins of equal width
genbins(view3, 'finwealth', 20)
colNames = view3.columns.tolist()
colNames[-5:-1] = ['var_pol1_cf', 'var_pol2_cf', 'var_pol3_cf', 'var_pol9_cf']
view3.columns = colNames
# %>

""" Step 4: Aggregate dataset into age-income bins.
    Melt dataset so it's a pure long dataset, each observation being
    defined on the age-income-variable-period level.
"""  # %<
view4 = view3.groupby(['age', 'income_state', 'variable']).mean()
values = view3.columns[view3.columns.str.contains('var_')].tolist()
view4.index.names = ['age', 'income_state', 'var']
melt(view4.reset_index(), ['age', 'income_state', 'var'], values).to_csv(
    'transition_eventstudy%s.csv' % affix, index=False)

view4 = view3.groupby(['finwealthbin', 'income_state', 'variable']).mean()
values = view3.columns[view3.columns.str.contains('var_')].tolist()
view4.index.names = ['finwealthbin', 'income_state', 'var']
(melt(view4.reset_index(), ['finwealthbin', 'income_state', 'var'], values)
 .to_csv('transition_eventstudy_wealth%s.csv' % affix, index=False))
# %>

# %>

""" Task 2) %<
Output elas.csv with this Pandas query on server
(merging steady-state and policy outputs)
TODO: Revamp this for latest subsidy experiments


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
"""  # %>
