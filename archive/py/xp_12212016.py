#!/usr/bin/python

""" Title: xp_12212016.py
    Purpose: Data cleaning for analysis in xp_12212016.do
    (i.e. auxiliary equations)
    Author: TC, DB
    Date: 2016-12-21

"""

from model.model_iterate import lifecycle_iterate
from pandas import merge, IndexSlice
from numpy import nan
import time

start = time.time()
model = lifecycle_iterate({})

# Merge FTHB hitting times onto full panel
view = model.readModel('fthb').appended
model.readModel('dist_fthb')
# Cast a dummy variable
model.appended['fthb_flag'] = 1
view_fthb = model.appended[['id', 'ageBought', 'fthb_flag']]
merged = merge(view, view_fthb, how='left',
               left_on=['id', 'age'], right_on=['id', 'ageBought'])
# Left unmatched -> missing values, so replace with 0
merged['fthb_flag'] = merged['fthb_flag'].fillna(0)
# Delete superfluous values where lags are undefined
merged.loc[merged['age'] <= 2, ['consumption_l', 'assets_l', 'income_l']
           ] = nan
merged.loc[merged['age'] <= 3, 'income_l2'] = nan

# Merge actual income flows and lagged flows onto full panel
view_incs = model.readModel('income_val').appended.iloc[:,:-1]
# To ensure the age cohort component of income matches correctly,
# have to transform ages for lagged income values
merged['age_l'] = merged['age'] - 1
merged['age_l2'] = merged['age'] - 2
merged = merge(merged, view_incs, how='left', on=['age', 'income'])
merged = merge(merged, view_incs, how='left', left_on=['age_l', 'income_l'],
               right_on=['age', 'income'], suffixes=('','_l'))
merged = merge(merged, view_incs, how='left', left_on=['age_l2', 'income_l2'],
               right_on=['age', 'income'], suffixes=('','_l2'))
merged = merged.drop(['model', 'ageBought', 'age_l', 'age_l2', 'income',
                      'income_l', 'income_l2'], axis=1)
merged = merged.set_index(['id', 'age']).sort_index()

# TODO: Is there a need to throw out homeowner data
# (since they're not FTHBs but unlike renters either?)

# Outputting for exploratory analysis in Stata. When actually
# calibrating, should be simple enough to keep all code in
# Python and save time wasted on transferring data.
merged.loc[IndexSlice[:,2:39],:].to_csv('fthb.csv')
end = time.time()
print(end - start)