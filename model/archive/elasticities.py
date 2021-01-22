#!/usr/bin/python

from model.model_iterate import lifecycle_iterate
from sys import argv
import numpy as np
from pandas import merge, DataFrame, IndexSlice, concat

readfunc = lifecycle_iterate.appendModels
# Call relevant model output files
view = (readfunc('fthb', model=[argv[1]]).appended
        [['id', 'age', 'nextAssets', 'assets', 'consumption', 'durTime',
          'rent', 'income', 'income_l', 'income_l2']])
viewT = readfunc('transition_fthb', model=[argv[1]]).appended
viewS = readfunc('dist_fthb', model=[argv[1]]).appended
incomes = readfunc('income_val', model=[argv[1]]).appended
view = merge(view, incomes, on=['age', 'income'])


# Useful functions for data processing
def genbins(data, col='age', bval=40):
    try:
        vrange = {'max': data.loc[:, col].max(), 'min': data.loc[:, col].min()}
        bins = np.linspace(vrange['min'], vrange['max'], bval)
    except:
        bins = bval
    data.loc[:, '%sbin' % col] = np.digitize(data[col], bins)
    j = 1
    for val in bins:
        data.loc[data['%sbin' % col] == j, '%sbin' % col] = val
        j += 1


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
        durBool = ((viewT['durablesOrig'] <= cap) & (viewT['policy'] == 0) &
                   (viewT['ageBought'] <= 39))  # Working age only
        durCols = ['durables', 'durablesOrig', 'durTime', 'income', 'assets']
    else:
        durBool = ((viewT['policy'] == 0) & (viewT['ageBought'] <= 39))
        durCols = ['durables', 'durTime', 'income', 'assets']

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


def elasticities(polType, durShift, durSS):
    print("Generate elasticities... \n \n")

    def prop_final(data, columns):
        try:
            denom = (view2.groupby(columns)
                     [['id']].count())
            left = data.groupby(columns)[['id']].count()
        except:
            denom = view.loc[:, ['id']].count()
            denom = DataFrame(denom).transpose()
            left = data[['id']].count()
            left = DataFrame(left).transpose()

        out = merge(left, denom, how='left', left_index=True,
                    right_index=True, suffixes=('', '_ss'))
        return out

    if polType == 'Repeat':
        repeatBool = ((view['rent'] == 0) & (view['vintage'] >= 5))
        view2 = view[repeatBool]

        # PROPORTION OF BUYERS BY OWNED DURABLE VINTAGE (CARS only?)
        # The -1 is because ageBought is model age + 1
        view2['vintage'] = view2['age'] - view2['durTime'] - 1
        view2['vintageShift'] = view2.groupby('id')['vintage'].shift(-1)

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
                genbins(durShift, colname)
                genbins(durSS, colname)

        for dims in ['total', 'vintage', 'age', 'assetsbin',
                     ['agebin', 'vintage']]:
            type_prop = prop_final(durSS, dims, repeatBool)
            all_stats.append(cleaning(type_prop['id']/type_prop['id_ss'],
                                      'Proportion of buyers, steady state', 0))
            type_prop = prop_final(durShift, dims, repeatBool)
            all_stats.append(cleaning(type_prop['id']/type_prop['id_ss'],
                                      'Proportion of buyers', 0))
    elif polType == 'First-time':
        fthbBool = ((view['rent'] == 1) & (view['durTime'] >= 3))
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
        print durShift.describe()
        print durSS.describe()

        view['vintage'] = view['age'] - view['durTime'] - 1

        for dims in ['total', 'income_valbin', 'age', 'assetsbin',
                     ['income', 'agebin'], ['income', 'assetsbin']]:
            type_prop = prop_final(durSS, dims)
            all_stats.append(cleaning(type_prop['id']/type_prop['id_ss'],
                                      'Proportion of buyers, steady state', 0))
            type_prop = prop_final(durShift, dims)
            all_stats.append(cleaning(type_prop['id']/type_prop['id_ss'],
                                      'Proportion of buyers', 0))


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
