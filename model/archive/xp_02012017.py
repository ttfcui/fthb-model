#!/usr/bin/python

""" Title: xp_02012017.py
    Purpose: Creating a common set of statistics for comparing policies.
    Author: TC, DB
    Date: 2016-12-10

"""

""" Summary: This file creates statistics grouped into "three margins,"
    as well as homeownership rates during the transition and
    substitution effects during policy takeup. It copies code
    from other py files.
"""
from model.model_iterate import lifecycle_iterate
import numpy as np
from pandas import merge, IndexSlice, read_table


def genbins(data, col='age'):
    data['agebin'] = np.floor_divide(data[col], 10)*10 + 20
idx = IndexSlice

view1 = lifecycle_iterate.readModel('transition_adjust').appended
view2 = lifecycle_iterate.readModel('lifecycle_adjust').appended
view3 = (lifecycle_iterate.readModel('fthb').appended
         [['id', 'age', 'nextAssets', 'assets', 'consumption',
           'income', 'income_l', 'income_l2']])
view4 = lifecycle_iterate.readModel('transition_fthb').appended
view5 = lifecycle_iterate.readModel('dist_fthb').appended
view6 = lifecycle_iterate.readModel('housing_transit').appended

# %< PART 1: PURE EXTENSIVE MARGIN STATS
merged = merge(view1, view2, on=['id', 'age', 'model'],
               suffixes=('_pol', '_ss'))
genbins(merged)

grouped = merged.groupby(['adjustments_ss', 'adjustments_pol', 'agebin'])
grouped_a = merged.groupby(['adjustments_ss', 'adjustments_pol', 'agebin'])
print(grouped['id'].count().loc[idx[0, 1:]])
print(grouped_a['id'].count().loc[idx[0, 1:, :]])
grouped_slice = grouped['id'].count().reset_index()
grouped_a_slice = grouped_a['id'].count().reset_index()
print(grouped_slice.loc[grouped_slice['adjustments_ss'] <
                        grouped_slice['adjustments_pol'], 'id'].sum())
print(grouped_slice.loc[grouped_slice['adjustments_ss'] >
                        grouped_slice['adjustments_pol'], 'id'].sum())
print(grouped_a_slice.loc[grouped_a_slice['adjustments_ss'] <
                          grouped_a_slice['adjustments_pol'], 'id']
      .groupby('agebin').sum())
print(grouped_a_slice.loc[grouped_a_slice['adjustments_ss'] >
                          grouped_a_slice['adjustments_pol'], 'id']
      .groupby('agebin').sum())
grouped = merged.groupby(['adjustments_ss', 'adjustments_pol'])
grouped['id'].count().to_csv('margins_test.csv')
# %>

# %< PART 2: PSEUDO-TIMING MARGINS
data = read_table('model/transition_difffull.txt', sep='\s+', header=None)
cols = data.columns.tolist()
cols[0:2] = ['id', 'age']
data.columns = cols
for i in xrange(39, 2, -1):
    data.iloc[:, i] = data.iloc[:, 2:i+1].sum(axis=1)

data = merge(data, view3, on=['id', 'age'])
genbins(data)
init = (data[2] == 1)
print(data.loc[init, 'income'].count())
print(data.loc[init, ['income', 'agebin']].groupby('agebin').count())
print(data.loc[init & (data[3] == 0), 'income'].describe())
print(data.loc[init & (data[3] == 0), ['income', 'agebin']]
      .groupby('agebin').describe())
print(data.loc[init & (data[3] == 1) & (data[4] == 0), 'income'].describe())
print(data.loc[init & (data[3] == 1) & (data[4] == 0), ['income', 'agebin']]
      .groupby('agebin').describe())
print(data.loc[init & (data[3] == 1) & (data[4] == 1), 'income'].describe())
print(data.loc[init & (data[3] == 1) & (data[4] == 1), ['income', 'agebin']]
      .groupby('agebin').describe())
# %>

# %< PART 3: INTENSIVE MARGINS (also see outputted file)
view4['policy'] = view4['ageBought'] - view4['age']
# Some imputation issues means unreasonably high values filtered out
cap = (view4.loc[view4['PolTaken'] == 'Repeat', ['nextDurables']]
       .describe().loc['max', ' nextDurables'])
durShift = view4.loc[(view4['durablesOrig'] <= cap) & (view4['policy'] == 0)
                     & (view4['PolTaken'] == 'Repeat'), ['id', 'ageBought',
                     'nextDurables', 'durables', 'durablesOrig']]
durSS = view5.loc[(view5['PolTaken'] == 'Repeat'),
                  ['id', 'ageBought', 'nextDurables', 'durables']]
# Steady-state data at time of takeup

# Merge and collapse
merged = merge(durShift, durSS, on=['id', 'ageBought'], how='left')
merged = merge(merged, view3,
               left_on=['id', 'ageBought'], right_on=['id', 'age'])
merged['marginal'] = 0
merged.loc[merged['durables_y'].isnull(), 'marginal'] = 1
merged['durablesGap'] = merged['durablesOrig'] - merged['nextDurables_x']
genbins(merged, col='ageBought')

minsize_bool = (abs(merged['nextDurables_y'] - 1.3) < 0.1)
print(merged.loc[minsize_bool].groupby('agebin')['id'].count())
minsize_bool = (abs(merged['nextDurables_x'] - 1.3) < 0.1)
print(merged.loc[minsize_bool].groupby('agebin')['id'].count())

aggdict = {'id': len, 'nextDurables_x': np.mean, 'nextDurables_y': np.mean,
           'durablesOrig': np.median, 'durablesGap': np.median,
           'income': np.mean, 'assets': np.mean}
collapsed = merged.groupby(['agebin', 'marginal']).agg(aggdict)
collapsed2 = merged.groupby(['marginal']).agg(aggdict)
print(collapsed2)
collapsed.to_csv('fthb_intensive_collapse.csv')
# %>

# %< PART 4: HOMEOWNERSHIP RATES OVER POLICY (FTHB only?)

# %>

# %< PART 5: TYPES OF MARGINAL BUYERS
durStats = view4.loc[(view4['durablesOrig'] <= cap) & (view4['policy'] == 0)
                     & (view4['PolTaken'] == 'Repeat'), ['id', 'ageBought',
                     'nextAssets', 'consumption']]
merged = merge(durStats, durSS, on=['id', 'ageBought'], how='left')
merged = merge(merged, view3,
               left_on=['id', 'ageBought'], right_on=['id', 'age'])
merged['marginal'] = 0
merged.loc[merged['durables'].isnull(), 'marginal'] = 1
genbins(merged, col='ageBought')

aggdict = {'id': len, 'consumption_x': np.mean, 'consumption_y': np.mean,
           'nextAssets_x': np.mean, 'nextAssets_y': np.mean,
           'income': np.mean, 'assets': np.mean}
collapsed = merged.groupby(['agebin', 'marginal']).agg(aggdict)
collapsed2 = merged.groupby(['marginal']).agg(aggdict)
# %>
