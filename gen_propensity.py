#!/opt/apps/rhel7/Python-2.7.11/bin/python

""" Title: gen_propensity.py
    Purpose: Create a high N, low-variable dataset used to calculate
    conditional purchase probability.
    Author: TC, DB
    Date: 2018-11-13
"""

from model.model_iterate import lifecycle_iterate
from sys import argv
import numpy as np
from pandas import merge, IndexSlice, isnull, concat, pivot, to_numeric
from scipy.stats import mode, entropy

# Useful functions for data processing
def ageTrim(dat, end=39):
    try:
        return dat.appended.loc[dat.appended['age'] <= end, :]
    except:
        return dat.appended.loc[dat.appended['ageBought'] <= end, :]


readfunc = lifecycle_iterate.appendModels
# Call relevant model output files
# TODO: the policy functions for the non-policy takers in this file
# are kind of a ruse, since if they choose differently in a policy
# sim while not buying a house they're not recorded in any file.
orig_col = ['id', 'age', 'nextDurables', 'nextAssets', 'consumption', 'durTime',
            'durables', 'assets', 'income', 'income_l', 'income_val', 'rent', 'moved', 'model']
view1 = (ageTrim(readfunc('fthb', model=[argv[1]]))[orig_col])
view1 = view1.sort_values(['id', 'age'])
view1['renter'] = view1.groupby('id')['rent'].shift(1)
if (argv[2] == 'First-time'):
    view1 = (view1.loc[(view1['renter'] == 1) & (view1['durTime'] >= 3)]
	     .drop(['rent', 'renter', 'durTime'], axis=1))
else:
    view1['vintage'] = view1['age'] - view1['durTime'] - 1
    view1['vintageShift'] = view1.groupby('id')['vintage'].shift(-1)
    view1 = (view1.loc[(view1['rent'] == 0) & (view1['vintage'] >= 3)]
             .drop(['rent', 'renter', 'durTime'], axis=1))

view2 = (ageTrim(readfunc('transition_fthb', model=[argv[1]]))
         [['id', 'age', 'ageBought', 'consumption', 'nextDurables',
           'nextAssets', 'model']])
view3 = (ageTrim(readfunc('dist_fthb', model=[argv[1]]))
         [['id', 'ageBought', 'model']])
view3.columns = ['id', 'age', 'model']

def prop1():
    """ For all policy-eligible HH (renters or repeat buyers), merge together
        a dataset of their outcomes during the policy period as well
        as in the counterfactual. Then classify the policy takers in this sample.

    """

    print("Part 1 \n \n")
    view2['policy'] = view2['ageBought'] - view2['age']
    view2r = view2[(view2['policy'] == 0)].drop(['age', 'policy'], axis=1)
    view2r = view2r.rename(columns={'ageBought': 'age'}) 
    view3['infraMarg'] = 1

    merged = merge(view1, view2r, how='right', on=['id', 'age', 'model'],
                   suffixes=('', '_pol'))
    merged['income_diff'] = merged['income'] - merged['income_l']
    merged = (merged.drop(['nextDurables', 'nextAssets', 'consumption'], axis=1)
              .rename(columns={'nextDurables_pol': 'nextDurables',
                               'nextAssets_pol': 'nextAssets',
                               'consumption_pol': 'consumption'})) 
    view2_merge = merge(view3, merged, how='right', on=['id', 'age', 'model'])
    view2_merge['purchPol'] = 1
    try:
        view2_merge.loc[isnull(view2_merge['infraMarg']), 'marginal'] = 1
    except:
        raise ValueError('Variable assignment failed. This is likely because '
            'the model output file has zero observations.')

    merged_nonbuyer = merge(view1, view2_merge[['id', 'age', 'infraMarg',
                            'purchPol', 'marginal']], how='left', on=['id', 'age'])
    merged_nonbuyer = merged_nonbuyer[(isnull(merged_nonbuyer['purchPol']))]
    # merged_nonbuyer.info()
    return concat([view2_merge, merged_nonbuyer])

def prop2(theta, hprice):
    """ For all policy-eligible HH, compute some outcomes in the post-policy
        period. This includes future homeownership status, deleveraging
        variables, and other possible variables up to 5/9 periods in future.

    Args:
        theta: Proportion of house value needed for down payment.
        hprice: Price level of housing (either PE or determined in GE).

    """

    print("Part 2 \n \n")
    view = (lifecycle_iterate.readModel('transition_lagsleads')
            .appended.sort_values(['id', 'age', 'variable']))
    floatCol = view.columns[~view.columns.isin(['variable', 'model'])].tolist()
    view.loc[:,floatCol] = view.loc[:,floatCol].apply(to_numeric, errors='coerce')
    viewVar = view.loc[view['variable'].isin(['H_own', 'Q'])]
    try:
        viewVar = viewVar.pivot_table(index=['id', 'age'], columns='variable')
    except:
        raise ValueError('Variable assignment failed. This is likely because '
            'the model output file has zero observations.')
    idx = IndexSlice

    # Short-term homeownership reversion
    viewVar.loc[:, idx['owner_pol1', 'H_own']] = (viewVar.loc
        [:, idx['var_pol1', 'H_own']] > 0).astype(int) 
    viewVar.loc[:, idx['owner_pol3', 'H_own']] = (viewVar.loc
        [:, idx['var_pol3', 'H_own']] > 0).astype(int)
    # Long-term homeownership reversion (9 years after)
    viewVar.loc[:, idx['owner_pol9', 'H_own']] = (viewVar.loc
        [:, idx['var_pol9', 'H_own']] > 0).astype(int) 

    # Short, medium term leveraging?
    # TODO: missing the price level in the house value term.
    viewVar.loc[:, idx['leverage_pol1', 'H_own']] = (viewVar.loc[:,
        idx['var_pol1', 'Q']] - (1.0-theta)*np.exp(hprice)*viewVar.loc[:,
        idx['var_pol1', 'H_own']])/viewVar.loc[:,idx['var_pol1', 'H_own']]
    viewVar.loc[:, idx['leverage_pol3', 'H_own']] = (viewVar.loc[:,
        idx['var_pol3', 'Q']] - (1.0-theta)*np.exp(hprice)*viewVar.loc[:,
        idx['var_pol3', 'H_own']])/viewVar.loc[:,idx['var_pol3', 'H_own']]
    viewOut = viewVar.loc[:, idx[['id', 'age', 'owner_pol1', 'owner_pol3', 'owner_pol9',
                          'leverage_pol1', 'leverage_pol3'], 'H_own']]
    viewOut.columns = viewOut.columns.droplevel(1)
    return viewOut

    
if __name__ == '__main__':
    print('Creating purchase propensity data: \n \n')
    merged = prop1().set_index(['id', 'age'])
    try:
        futVars = prop2(float(argv[3]),float(argv[4]))
        merged = merged.join(futVars, how='left')
    except:
        pass
    merged = merged.drop(['purchPol', 'income_l'], axis=1)
    merged.to_csv('propensity_%s.csv' % argv[1])

