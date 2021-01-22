#!/usr/bin/python

from model.model_iterate import lifecycle_iterate
from sys import argv
import numpy as np
from pandas import merge, DataFrame, IndexSlice, concat


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


def ageTrim(dat, end=39):
    try:
        return dat.appended.loc[dat.appended['age'] <= end, :]
    except:
        return dat.appended.loc[dat.appended['ageBought'] <= end, :]


def cleaning(data, desc, *args):
    data = DataFrame(data)
    try:
        data.reset_index()['agebin'].describe()
    except:
        data['agebin'] = 0  # All ages (before retirement)
    data.loc[:, 'desc'] = desc
    data = data.reset_index().set_index(['desc', 'agebin'])
    try:
        data = data.rename(columns={args[0]: 'value'})
        data = data.rename(columns={args[1]: 'subtype'})
    except:
        pass
    return data
idx = IndexSlice


def elas_prep(polType, dep):
    print("Data cleaning... \n \n")
    viewT['policy'] = viewT['ageBought'] - viewT['age']
    # Impute original value of owned durable using vintage stats
    if polType == 'Repeat':
        viewT['durablesOrig'] = ((1.0/(1.0-dep))**(
                                 viewT['ageBought'] - viewT['durTime'])*(
                                 viewT['durables']))
        # Some imputation issues means unreasonably high values filtered out
        # (but this only applies to the repeat buyer case)
        cap = (viewT.loc[viewT['PolTaken'] == polType, ['nextDurables']]
               .describe().loc['max', 'nextDurables'])
        durBool = (viewT['durablesOrig'] <= cap) & (viewT['policy'] == 0)
        durCols = ['durables', 'durablesOrig', 'consumption', 'durTime',
                   'income', 'assets']
    else:
        durBool = (viewT['policy'] == 0)
        durCols = ['durables', 'consumption', 'durTime', 'income', 'assets']

    durShift = viewT.loc[durBool & (viewT['PolTaken'] == polType),
                         ['id', 'ageBought', 'nextDurables'] + durCols]
    durShift = merge(durShift, incomes, left_on=['ageBought', 'income'],
                     right_on=['age', 'income'])
    durSS = viewS.loc[viewS['PolTaken'] == polType, ['id', 'ageBought',
                      'nextDurables'] + durCols]
    durSS = merge(durSS, incomes, left_on=['ageBought', 'income'],
                  right_on=['age', 'income'])
    return durShift, durSS

# %>

readfunc = lifecycle_iterate.appendModels
# Call relevant model output files
view = ((ageTrim(readfunc('fthb', model=[argv[1]]))
        [['id', 'age', 'nextDurables', 'durables', 'nextAssets', 'assets',
          'consumption', 'durTime', 'rent',
          'income', 'income_l', 'income_l2']]))
view = view.sort_values(['id', 'age'])
view['renter'] = view.groupby('id')['rent'].shift(1)
viewT = ageTrim(readfunc('transition_fthb', model=[argv[1]]))
viewS = ageTrim(readfunc('dist_fthb', model=[argv[1]]))
incomes = ageTrim(readfunc('income_val', model=[argv[1]]))
view = merge(view, incomes, on=['age', 'income'])


def elasticities(polType, durShift, durSS):
    print("Generate elasticities... \n \n")

    def prop_final(denomdata, data, columns, grpfunc):
        try:
            denom = denomdata.groupby(columns).aggregate(grpfunc)
            left = data.groupby(columns).aggregate(grpfunc)
        except:
            denom = denomdata.groupby(lambda idx: 0).aggregate(grpfunc)
            left = data.groupby(lambda idx: 0).aggregate(grpfunc)

        out = merge(left, denom, how='left', left_index=True,
                    right_index=True, suffixes=('', '_ss'))
        return out
    contAggFunc = lambda x: np.mean(np.log(x))

    if polType == 'Repeat':

        # PROPORTION OF BUYERS BY OWNED DURABLE VINTAGE (CARS only?)
        # The -1 is because ageBought is model age + 1
        view['vintage'] = view['age'] - view['durTime'] - 1
        view['vintageShift'] = view.groupby('id')['vintage'].shift(-1)
        repeatBool = ((view['rent'] == 0) & (view['vintage'] >= 3))
        view2 = view[repeatBool]

        durShift['vintage'] = durShift['ageBought'] - durShift['durTime'] - 1
        durSS['vintage'] = durSS['ageBought'] - durSS['durTime'] - 1
        for colname in ['age', 'vintage', 'assets']:
            if colname == 'age':
                agedivs = [9, 19, 29, 39]
                genbins(view2, colname, agedivs)
                genbins(durShift, colname, agedivs)
                genbins(durSS, colname, agedivs)
            else:
                genbins(view2, colname)
                divs = (view2['%sbin' % colname].sort_values()
                        .unique())
                genbins(durShift, colname, divs)
                genbins(durSS, colname, divs)

        for dims in ['total', 'vintage', 'age', 'assetsbin',
                     ['agebin', 'vintage']]:
            type_prop = prop_final(view2, durSS, dims,
                                   {'id': lambda x: x.shape[0],
                                    'consumption': contAggFunc,
                                    'nextDurables': contAggFunc,
                                    'durables': contAggFunc})
            all_stats.append(cleaning(type_prop['id']/type_prop['id_ss'],
                                      'Proportion of buyers, steady state', 0))
            all_stats.append(cleaning(type_prop['consumption'],
                                      'Log consumption units, steady state',
                                      'consumption'))
            all_stats.append(cleaning(type_prop['nextDurables'] -
                                      type_prop['durables'], 'Log chg '
                                      'housing services, steady state', 0))
            type_prop = prop_final(view2, durShift, dims,
                                   {'id': lambda x: x.shape[0],
                                    'consumption': contAggFunc,
                                    'nextDurables': contAggFunc,
                                    'durables': contAggFunc})
            all_stats.append(cleaning(type_prop['id']/type_prop['id_ss'],
                                      'Proportion of buyers', 0))
            all_stats.append(cleaning(type_prop['consumption'],
                                      'Log consumption units', 'consumption'))
            all_stats.append(cleaning(type_prop['nextDurables'] -
                                      type_prop['durables'],
                                      'Log chg housing services', 0))
            all_stats.append(cleaning(type_prop['id_ss'], '# Buyers in Bin',
                                      'id_ss'))
    elif polType == 'First-time':
        fthbBool = ((view['renter'] == 1) & (view['durTime'] >= 3))
        view2 = view.loc[fthbBool]

        for colname in ['age', 'income_val', 'assets']:
            if colname == 'age':
                agedivs = [9, 19, 29, 39]
                genbins(view2, colname, agedivs)
                genbins(durShift, colname, agedivs)
                genbins(durSS, colname, agedivs)
            else:
                genbins(view2, colname)
                divs = (view2['%sbin' % colname].sort_values()
                        .unique())
                genbins(durShift, colname, divs)
                genbins(durSS, colname, divs)
        durShift.info()
        durSS.info()

        view['vintage'] = view['age'] - view['durTime'] - 1

        for dims in ['total', 'income_valbin', 'age', 'assetsbin',
                     ['income', 'agebin'], ['income', 'assetsbin']]:
            type_prop = prop_final(view2, durSS, dims,
                                   {'id': lambda x: x.shape[0],
                                    'consumption': contAggFunc,
                                    'nextDurables': contAggFunc,
                                    'durables': contAggFunc})
            all_stats.append(cleaning(type_prop['id']/type_prop['id_ss'],
                                      'Proportion of buyers, steady state', 0))
            all_stats.append(cleaning(type_prop['consumption'],
                                      'Log consumption units, steady state',
                                      'consumption'))
            all_stats.append(cleaning(type_prop['nextDurables'] -
                                      type_prop['durables'], 'Log chg '
                                      'housing services, steady state', 0))
            type_prop = prop_final(view2, durShift, dims,
                                   {'id': lambda x: x.shape[0],
                                    'consumption': contAggFunc,
                                    'nextDurables': contAggFunc,
                                    'durables': contAggFunc})
            all_stats.append(cleaning(type_prop['id']/type_prop['id_ss'],
                                      'Proportion of buyers', 0))
            all_stats.append(cleaning(type_prop['consumption'],
                                      'Log consumption units', 'consumption'))
            all_stats.append(cleaning(type_prop['nextDurables'] -
                                      type_prop['durables'],
                                      'Log chg housing services', 0))
            all_stats.append(cleaning(type_prop['id_ss'], '# Buyers in Bin',
                                      'id_ss'))


# %>


if __name__ == '__main__':
    print('Generating elasticities: \n \n')
    all_stats = []
    durShift, durSS = elas_prep(argv[2], float(argv[3]))
    elasticities(argv[2], durShift, durSS)
    final = concat(all_stats)
    try:
        final = final.drop(['adjustment_ss'], axis=1)
    except:
        pass
    final = final.swaplevel(-2, -1).sort_index()
    final.to_csv('stats_elasticities_%s.csv' % argv[1])
