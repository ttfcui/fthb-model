#!/usr/bin/python

""" Title: xp_12102016.py
    Purpose: Data cleaning for analysis in xp_12102016.do.
    Author: TC, DB
    Date: 2016-12-10

"""

""" Summary: This file contains code for numerous tasks, working somewhat
    off of each other but mostly independent. Purpose is to check
    calibration accurate and algorithms work as they should.
"""
from model.model_iterate import lifecycle_iterate
import numpy as np
from pandas import merge, concat, read_table
test = lifecycle_iterate({})
price = np.exp(0.38)  # Approximate steady-state price

# %< Asset race prep
view = lifecycle_iterate.readModel('fthb').appended.sort(['id', 'age'])
view.loc[view['age'] == 2, 'income_init'] = view['income']
view['income_init'] = view['income_init'].fillna(method='ffill')
view['income_last'] = view.groupby('id')['income'].shift(1)
view['rental_yrs'] = view.groupby('id')['rent'].transform(sum)

asset_group = (view.loc[(view['rent'] == 1) & (view['age'] < 40) &
                        (view['rental_yrs'] < 30)
                        ].groupby(['incomeS_init', 'age'])['nextAssets'])
asset_panel = concat([asset_group.agg(lambda x: np.percentile(x, 75)),
                      asset_group.agg(lambda x: np.percentile(x, 50))],
                     axis=1)
asset_panel.columns = ['assets_75', 'assets_50']
asset_panel.to_csv('asset_panel.csv')
# %>

# %< Counting up anomalous behavior at model start
vslice = view.loc[(view['age'] <= 6) & (view['rent'] == 0),
                  ['id', 'age', 'nextDurables', 'income',
                  'income_init', 'nextAssets', 'consumption']]
vslice['nextDurables'] = vslice['nextDurables']*price
data_h = vslice.groupby(['income', 'age'])['nextDurables'].agg(
    {'housing': np.mean, 'N': len})
data_a = vslice.groupby(['income', 'age'])['nextAssets'].agg(
    {'assets': np.mean})
data_c = vslice.groupby(['income', 'age'])['consumption'].agg(
    {'consumption': np.mean})
data_d = vslice.groupby(['income', 'age'])['income_init'].agg(
    {'last_shock': np.mean})

print(data_h.join(data_a).join(data_c).join(data_d))
# %>

# %< Conditional on starting point, when do people actually buy?
# ANSWER: many at the same income state - so they're waiting
# to get rich versus saving at the beginning of life
incomeCompare = view.loc[(view['rent'] == 0) & (view['age'] < 40)
                         ].groupby(['income_init', 'age'])['income']
print(incomeCompare.median())
# %>

# %< Income means over the lifecycle
incomes = lifecycle_iterate.readModel('income_val').appended
inc_view = merge(view, incomes, how='left', on=['age', 'income'])
# Using Pandas to generate mean, 75th, 25th percentiles
# (then taking the log of the values with numpy)
income_dist = np.log(inc_view.groupby('age')['income_val'].agg(
                     {'logy_sim': np.mean,
                      'logy_sim75': lambda x: np.percentile(x, 75),
                      'logy_sim25': lambda x: np.percentile(x, 25)}))
income_dist[0:38].to_csv('model_incomes.csv')
# %>

# %< Collapse steady-state and transition FTHB data, then append them
# together with relevant aggregated statistics
# Steady-state statistics
view_ss = lifecycle_iterate.readModel('dist_fthb').appended
view_ss['nextDurables'] = view_ss['nextDurables']*price
view_ss['h_c_ratio'] = view_ss['nextDurables']/view_ss['consumption']
view_ss['consumption_chg'] = view_ss['consumption']/view_ss['consumption_l']
view_ss['housing_chg'] = view_ss['nextDurables']/price/view_ss['durables']

# Policy transition statistics
view = lifecycle_iterate.appendModels('transition_fthb', model=[
    'experiment_monetary_nodown', 'experiment_monetary_little']).appended
view['period'] = view['ageBought'] - view['age']
view['nextDurables'] = view['nextDurables']*price
view['h_c_ratio'] = view['nextDurables']/view['consumption']
view['consumption_chg'] = view['consumption']/view['consumption_l']
view['housing_chg'] = view['nextDurables']/price/view['durables']

ss_stats = view_ss.groupby('ageBought')[[
    'h_c_ratio', 'assets', 'consumption_chg', 'housing_chg']].median()
ss_stats = ss_stats.reset_index()
ss_stats['period'] = -1
ss_stats['model'] = 'Steady-state'
ss_stats = ss_stats.set_index(['model', 'period', 'ageBought'])

transition_stats = (view.groupby(['model', 'period', 'ageBought'])
                    [['h_c_ratio', 'assets', 'consumption_chg', 'housing_chg'
                      ]].median())
ss_stats.append(transition_stats).replace([np.inf], np.nan).to_csv(
    'transition_collapse_12112016.csv')
# %>

# %< Are people induced by the policy spending within their budget constraint?
merged = merge(view[(view['model'] == 'experiment_monetary_nodown') &
                    (view['period'] == 0)], view_ss, on=['id', 'ageBought'],
               how='left', suffixes=('', '_ss'))
# Looking at policy-induced people not receiving an income shock
merged2 = merged[(merged['nextDurables_ss'].isnull()) &
                 (merged['income_l'] == merged['income']) &
                 (merged['ageBought'] > 2)].iloc[:, 1:16]
# Accounting for changes in the deterministic income process
ageEarn = read_table('model/input_data/ageearnings.txt',
                     header=None).reset_index()
ageEarn['index'] += 2
ageEarn[1] = ageEarn[0].shift(1)
ageEarn['change'] = np.exp(ageEarn[0])/np.exp(ageEarn[1]) - 1
merged2 = merge(merged2, ageEarn[['index', 'change']], left_on='ageBought',
                right_on='index')

# Cashflow compares shifting composition and changes in expenditures
# between period before taking policy and period taking policy.
# By Walras's law cashflow should equal to 0.0
merged2['rentcost'] = merged2['durables']*((price*0.045)+0.0027)
merged2['dpcost'] = merged2['nextDurables']*0.236*1.01
merged2['cashflow'] = (merged2['rentcost'] - merged2['dpcost'] +
                       (merged2['consumption_l'] - merged2['consumption']) +
                       (merged2['assets'] - merged2['nextAssets']) +
                       merged2['change'] + .024*(
                       merged2['assets_l'] - merged2['assets']))
merged2['cashflow'].describe()
# %>
