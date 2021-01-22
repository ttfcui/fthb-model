#!/usr/bin/python

""" Title: xp_01262017.py
    Purpose: Breaking down marginal buyer behaviour.
    Author: TC, DB
    Date: 2016-12-10

"""

""" Summary: This file contains code for two tasks: calculating an
    "extensive margin" over the long run of people changing their durable
    adjustment plans and collapsing some marginal buyer statistics.
"""
from model.model_iterate import lifecycle_iterate
import numpy as np
from pandas import merge

# %< Calculating extensive margin
view1 = lifecycle_iterate.readModel('transition_adjust').appended
view2 = lifecycle_iterate.readModel('lifecycle_adjust').appended
merged = pd.merge(view1, view2, on=['id', 'age', 'model'],
                  suffixes=('_pol','_ss'))
# Sanity check: policy incentive shouldn't drive people to make 0 transactions
merged_bool = (merged['adjustments_pol'] == 0) & (merged['adjustments_ss'] > 0)
merged[merged_bool].groupby('age').count()

# How many people went from x adjustments in SS to y adjustments in POL?
grouped = merged.groupby(['adjustments_ss', 'adjustments_pol')
grouped['id'].count().to_csv('margins_test.csv')
# What's the long-run aggregate count on transactions induced by policy?
merged['diff'] = merged['adjustments_pol'] - merged['adjustments_ss']
merged.groupby('adjustments_ss').agg(
    {'adjustments_pol': np.sum, 'diff': np.sum})
# %>

# %< Marginal buyer statistics (mostly intended for repeat buyer policy)
view = lifecycle_iterate.readModel('transition_fthb').appended
view['policy'] = view['ageBought'] - view['age']
# Some imputation issues means unreasonably high values filtered out
cap = (view.loc[view['PolTaken'] == 'Repeat', ['nextDurables']]
       .describe().loc['max','nextDurables'])
durShift = view.loc[(view['durablesOrig'] < b) & (view['PolTaken'] == 'Repeat')
                    & (view['policy'] == 0), ['id', 'ageBought',
                    'nextDurables', 'durables','durablesOrig']]
view2 = lifecycle_iterate.readModel('dist_fthb').appended
durSS = view2.loc[(view2['PolTaken']=='Repeat'),
                  ['id', 'ageBought', 'nextDurables', 'durables']]
# Steady-state data at time of takeup
view3 = (lifecycle_iterate.readModel('fthb')
         .appended[['id', 'age', 'assets', 'income', 'income_l', 'income_l2']])

# Merge and collapse
merged = pd.merge(durShift, durSS, on=['id', 'ageBought'], how='left')
merged = pd.merge(merged, view3,
                  left_on=['id', 'ageBought'], right_on=['id', 'age'])
merged['marginal'] = 0
merged.loc[merged['durables_y'].isnull(), 'marginal'] = 1
merged['durablesGap'] = merged['durablesOrig'] - merged['nextDurables_x']
collapsed = merged.groupby(['ageBought', 'marginal']).agg({'id': len,
    'nextDurables_x': np.mean, 'nextDurables_y': np.mean,
    'durablesOrig': np.mean, 'durablesGap': np.mean,
    'income': np.mean, 'assets': np.mean})
collapsed.to_csv('fthb_intensive_collapse.csv')
