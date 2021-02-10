#!/usr/bin/python

""" Title: gen_stats.py
    Purpose: Creating a common set of statistics for comparing policies.
    Author: TC, DB
    Date: 2018-11-13

"""

""" Summary: This file creates statistics grouped into "three margins,"
    as well as homeownership rates during the transition and
    substitution effects during policy takeup. It copies code
    from other py files.
"""
from model.model_iterate import lifecycle_iterate
from sys import argv
import numpy as np
from pandas import merge, IndexSlice, read_table, melt, concat
from scipy.stats import mode, entropy


# Definition of Hellinger distance (numpy script found online)
def hellinger2(p, q):
    from scipy.spatial.distance import euclidean
    return euclidean(np.sqrt(p), np.sqrt(q)) / np.sqrt(2)


# Useful functions for data processing
def genbins(data, col='age'):
    data['agebin'] = np.floor_divide(data[col], 10)*10 + 20


def ageTrim(dat, end=39):
    try:
        return dat.appended.loc[dat.appended['age'] <= end, :]
    except:
        return dat.appended.loc[dat.appended['ageBought'] <= end, :]


def cleaning(data, desc, *args):
    from pandas import DataFrame
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

readfunc = lifecycle_iterate.appendModels
# Call relevant model output files
view1 = ageTrim(readfunc('transition_adjust', model=[argv[1]]))
view2 = ageTrim(readfunc('lifecycle_adjust', model=[argv[1]]))
view3 = (ageTrim(readfunc('fthb', model=[argv[1]]))
         [['id', 'age', 'nextAssets', 'assets', 'consumption', 'durTime',
           'rent', 'adjust', 'income_val', 'income', 'income_l', 'income_l2']])
view3 = view3.sort_values(['id', 'age'])
view4 = ageTrim(readfunc('transition_fthb', model=[argv[1]]))
view5 = ageTrim(readfunc('dist_fthb', model=[argv[1]]))
view6 = ageTrim(readfunc('housing_transit', model=[argv[1]]))
view7 = readfunc('lifecycleprofiles', model=[argv[1]])
all_stats = []


# %< PART 1: PURE EXTENSIVE MARGIN STATS
def stats1():
    print("Part 1 \n \n")
    merged = merge(view1, view2, on=['id', 'age', 'model'],
                   suffixes=('_pol', '_ss'))
    genbins(merged)

    grouped = merged.groupby(['adjustment_ss', 'adjustment_pol'])
    grouped_a = merged.groupby(['adjustment_ss', 'adjustment_pol', 'agebin'])
    grpSlice = grouped['id'].count().reset_index()
    grpSlice['diff'] = ((grpSlice['adjustment_pol'] -
                         grpSlice['adjustment_ss'])*(grpSlice['id']))
    grpASlice = grouped_a['id'].count().reset_index()
    grpASlice['diff'] = ((grpASlice['adjustment_pol']
                          - grpASlice['adjustment_ss'])*(grpASlice['id']))

    for num in xrange(0, 3):
        try:
            all_stats.append(cleaning(grouped['id'].count().loc[idx[num, 1:]],
                                      'Extensive from {}'.format(num),
                                      'id', 'adjustment_pol'))
            all_stats.append(cleaning(grouped_a['id'].count().loc[idx[num, 1:,
                                      :]], 'Extensive from {}'.format(num),
                                      'id', 'adjustment_pol'))
        except:
            all_stats.append(cleaning([0], 'Extensive from {}'.format(num), 0))

    # TODO: May have to put try statements around agebin collapses,
    # because will fail out due to groupby clause
    for cond, sub in [('grpSlice[\'diff\'] > 0', 'More transactions'),
                      ('grpSlice[\'diff\'] < 0', 'Fewer transactions')]:
        try:
            all_stats.append(cleaning([grpSlice.loc[eval(cond), 'diff'].sum()],
                '{} over lifecycle'.format(sub), 0))
            all_stats.append(cleaning(grpASlice.loc[eval(cond), ['diff', 'agebin']]
                .groupby('agebin').sum(), '{} over lifecycle'.format(sub), 'diff'))
        except:
            all_stats.append(cleaning([0], '{} over lifecycle'.format(sub), 0))
            

# %>


# %< PART 2: PSEUDO-TIMING MARGINS
def stats2(mdir='model', tretire=39, end=39):
    print("Part 2 \n \n")
    data = read_table('%s/transition_difffull.txt' % mdir,
                      sep='\s+', engine='python', header=None)
    cols = data.columns.tolist()
    cols[0:2] = ['id', 'age']
    data.columns = cols
    data['age'] += 1  # to reconcile age with other datasets
    try:
        for i in xrange(end, 2, -1):
            data.iloc[:, i] = data.iloc[:, 2:i+1].sum(axis=1)
    except:
        raise ValueError('Variable assignment failed. This is likely because '
            'the model output file has zero observations.')

    data = merge(data, view3, on=['id', 'age'])
    data = data[data['age'] <= tretire]
    genbins(data)
    init = (data[2] == 1)
    statFunc = lambda x: x.describe().loc[['count', 'mean']]
    statGrpFunc = lambda x: x.describe().loc[idx[:, ['count', 'mean']]]
    ivar = 'income_val'

    for cond, sub in [('data[end] > 0', 'more transactions'),
            ('data[end] < 0', 'fewer transactions')]:
        try:
            all_stats.append(cleaning(statFunc(data.loc[eval(cond), ivar]),
                                      'All agents w/ {} + mean income'.format(sub),
                                      ivar, 'index'))
            all_stats.append(cleaning(statGrpFunc(data.loc[eval(cond)]
                                                  .groupby('agebin')[ivar]),
                                      'All agents w/ {} + mean income'.format(sub),
                                      ivar, 'level_1'))
        except:
            pass


    for cond, sub in [('init', 'total'),
                      ('init & (data[end] > 0)', 'more transactions'),
                      ('init & (data[end] < 0)', 'fewer transactions')]:
        try:
            all_stats.append(cleaning(statFunc(data.loc[eval(cond), ivar]),
                                      ('All marginal transactions + mean '
                                       'income, %s' % sub), ivar, 'index'))
            all_stats.append(cleaning(statGrpFunc(data.loc[eval(cond)]
                                                  .groupby('agebin')[ivar]),
                                      ('All marginal transactions + mean inco'
                                       'me, %s' % sub), ivar, 'level_1'))

            all_stats.append(cleaning(statFunc(data.loc[eval(cond) &
                                                        (data[3] == 0),
                                                        ivar]),
                                      ('Marginal buyers shifting transaction '
                                       '1 period fwd, %s' % sub),
                                      ivar, 'index'))
            all_stats.append(cleaning(statGrpFunc(data.loc[eval(cond)
                                                           & (data[3] == 0)].
                                                  groupby('agebin')[ivar]),
                                      ('Marginal buyers shifting transaction '
                                       '1 period fwd, %s' % sub),
                                      ivar, 'level_1'))

            all_stats.append(cleaning(statFunc(data.loc[eval(cond) &
                                                        (data[3] != 0) &
                                                        (data[4] - data[3]
                                                         == -1), ivar]),
                                      ('Marginal buyers shifting transaction '
                                       '2 periods fwd, %s' % sub),
                                      ivar, 'index'))
            all_stats.append(cleaning(statGrpFunc(data.loc[eval(cond) &
                                                           (data[3] != 0) &
                                                           (data[4] - data[3]
                                                            == -1)].groupby
                                                  ('agebin')[ivar]),
                                      ('Marginal buyers shifting transaction '
                                       '2 periods fwd, %s' % sub),
                                      ivar, 'level_1'))

            all_stats.append(cleaning(statFunc(data.loc[eval(cond) &
                                                        (data[3] != 0) &
                                                        (data[4] - data[3]
                                                         != -1), ivar]),
                                      ('Marginal buyers shifting transaction '
                                       '2+ periods fwd, %s' % sub),
                                      ivar, 'index'))
            all_stats.append(cleaning(statGrpFunc(data.loc[eval(cond) &
                                                           (data[3] != 0) &
                                                           (data[4] - data[3]
                                                            != -1)].groupby
                                                  ('agebin')[ivar]),
                                      ('Marginal buyers shifting transaction '
                                       '2+ periods fwd, %s' % sub),
                                      ivar, 'level_1'))
        except:
            pass

# %>


# %< PART 3: INTENSIVE MARGINS
def stats3(polType, dep, dmin, polEnd):
    print("Part 3 \n \n")
    view4['policy'] = view4['ageBought'] - view4['age']

    all_stats.append(cleaning([view3['adjust'].value_counts()[2]],
                              'Total transactions', 0))
    agentCount = float(view3.shape[0])
    all_stats.append(cleaning([view3['adjust'].value_counts()[2]/agentCount],
                              'Turnover rate', 0))
    # Impute original value of owned durable using vintage stats
    if polType == 'Repeat':
        view4['durablesOrig'] = ((1.0/(1.0-dep))**(
                                 view4['ageBought'] - view4['durTime'])*(
                                 view4['durables']))
        # Some imputation issues means unreasonably high values filtered out
        # (but this only applies to the repeat buyer case)
        cap = (view4.loc[view4['PolTaken'] == polType, ['nextDurables']]
               .describe().loc['max', 'nextDurables'])
        durBool = ((view4['durablesOrig'] <= cap) & (view4['policy'] < polEnd))
        durCols = ['durables', 'income_val', 'durablesOrig', 'durTime']
    else:
        durBool = (view4['policy'] == 0)
        durCols = ['durables', 'durTime']

    durShift = view4.loc[durBool & (view4['PolTaken'] == polType),
                         ['id', 'ageBought', 'nextDurables'] + durCols]
    durSS = view5.loc[view5['PolTaken'] == polType, ['id', 'ageBought',
                      'durTime', 'nextDurables', 'durables', 'income_val']]

    # Merge and collapse
    merged = merge(durShift, durSS, on=['id', 'ageBought'], how='left',
                   suffixes=('', '_ss'))
    merged = merge(merged, view3, suffixes=('', '_ss'),
                   left_on=['id', 'ageBought'], right_on=['id', 'age'])
    merged['marginal'] = 0
    merged.loc[merged['durables_ss'].isnull(), 'marginal'] = 1
    if polType == 'Repeat':
        merged['durablesGap'] = merged['durablesOrig'] - merged['nextDurables']
    genbins(merged, col='ageBought')

    minsize_bool = (abs(merged['nextDurables_ss'] - dmin) < 0.1)
    try:
        all_stats.append(cleaning([merged.loc[minsize_bool, 'id'].count()],
                                  ('Agents buying minimum sized durable, '
                                   'steady-state'), 0))
        all_stats.append(cleaning(merged.loc[minsize_bool]
                                  .groupby('agebin')['id'].count(),
                                  ('Agents buying minimum sized durable, '
                                   'steady-state'), 'id'))
    except:
        pass

    minsize_bool = (abs(merged['nextDurables'] - dmin) < 0.1)
    try:
        all_stats.append(cleaning([merged.loc[minsize_bool, 'id'].count()],
                                  ('Agents buying minimum sized durable, '
                                   'policy'), 0))
        all_stats.append(cleaning(merged.loc[minsize_bool]
                                  .groupby('agebin')['id'].count(),
                                  ('Agents buying minimum sized durable, '
                                   'policy'), 'id'))
    except:
        pass

    aggdict = {'id': len, 'nextDurables': np.mean, 'nextDurables_ss': np.mean,
               'durables_ss': np.median}
    if polType == 'Repeat':
        aggdict.update({'durablesOrig': np.median, 'durablesGap': np.median})
    collapsed = merged.groupby(['marginal']).agg(aggdict)
    collapsed2 = merged.groupby(['agebin', 'marginal']).agg(aggdict)
    all_stats.append(cleaning(collapsed.loc[idx[0], :],
                              'Intensive margin magnitudes', 0, 'index'))
    all_stats.append(cleaning(melt(collapsed2.loc[idx[:, 0], :].reset_index(),
                                   'agebin'), 'Intensive margin magnitudes',
                              'value', 'variable'))

    all_stats.append(cleaning(collapsed.loc[idx[1], :],
                              'Intensive margin magnitudes, marginal buyers',
                              1, 'index'))
    all_stats.append(cleaning(melt(collapsed2.loc[idx[:, 1], :].reset_index(),
                                   'agebin'),
                              'Intensive margin magnitudes, marginal buyers',
                              'value', 'variable'))
    return durShift, durSS

# %>


# %< PART 4: POLICY-SPECIFIC TIME SERIES
def stats4(polType, durShift, durSS):
    print("Part 4 \n \n")
    if polType == 'Repeat':
        # PROPORTION OF BUYERS BY OWNED DURABLE VINTAGE (CARS only?)
        # The -1 is because ageBought is model age + 1
        view3['vintage'] = view3['age'] - view3['durTime'] - 1
        view3['vintageShift'] = view3.groupby('id')['vintage'].shift(-1)
        genbins(view3)
        durShift['vintage'] = durShift['ageBought'] - durShift['durTime'] - 1
        genbins(durShift, col='ageBought')
        durSS['vintage'] = durSS['ageBought'] - durSS['durTime'] - 1
        genbins(durSS, col='ageBought')

        # "Control" sample by limiting only to agents close to
        # median income (10 pctile band)
        ss_incBool = view3['income_val'].quantile([0.4, 0.5, 0.6]).tolist()
        vintageBool = ((view3['rent'] == 0) & (view3['vintage'] > 0) &
                       (view3['income_val'] >= ss_incBool[0]) &
                       (view3['income_val'] <= ss_incBool[-1]))

        def vintage_final(data, columns):
            vintage_count = (view3.loc[vintageBool].groupby(columns)
                             [['id']].count())
            vintage_inelig = (view3.loc[vintageBool & (view3['vintage'] < 3)
                                        & (view3['vintageShift'] == 1)]
                              .groupby(columns)[['id']].count())
            vintage_left = data.groupby(columns)[['id']].count()
            out = merge(vintage_left, vintage_count, how='left',
                        left_index=True, right_index=True,
                        suffixes=('', '_ss'))
            out = out.append(merge(vintage_inelig, vintage_count, how='left',
                                   left_index=True, right_index=True,
                                   suffixes=('', '_ss')))
            return out

        durSS_input = durSS.loc[(durSS['income_val'] >= ss_incBool[0]) &
                                (durSS['income_val'] <= ss_incBool[-1])]
        durShift_input = durShift.loc[
            (durShift['income_val'] >= ss_incBool[0]) &
            (durShift['income_val'] <= ss_incBool[-1])]

        vintage_prop = vintage_final(durSS_input, 'vintage')
        all_stats.append(cleaning(vintage_prop['id']/vintage_prop['id_ss'],
                                  ('Proportion of buyers by capital vintage'
                                   ', steady state'), 0, 'vintage'))
        vintage_prop = vintage_final(durShift_input, 'vintage')
        all_stats.append(cleaning(vintage_prop['id']/vintage_prop['id_ss'],
                                  'Proportion of buyers by capital vintage',
                                  0, 'vintage'))

        vintage_prop = vintage_final(durSS_input, ['agebin', 'vintage'])
        all_stats.append(cleaning(vintage_prop['id']/vintage_prop['id_ss'],
                                  ('Proportion of buyers by capital vintage'
                                   ', steady state'), 0, 'vintage'))
        vintage_prop = vintage_final(durShift_input, ['agebin', 'vintage'])
        all_stats.append(cleaning(vintage_prop['id']/vintage_prop['id_ss'],
                                  'Proportion of buyers by capital vintage',
                                  0, 'vintage'))
        # reentered in later merged data
        view3.drop('agebin', axis=1, inplace=True)
    elif polType == 'First-time':
        # HOMEOWNERSHIP RATES OVER POLICY (FTHB only?)
        view7.appended = view7.appended.iloc[0:38]
        view7.appended['model'] = 'Steady-state'
        own_ss_tot = (view7.genMoments(99, ['fracOwn'])
                      .reset_index(-1, drop=True))  # Entire dataset
        own_ss = view7.genMoments(10, ['fracOwn'])  # 10-yr age bins
        own_ss.index.names = ['model', 'agebin']  # "ageb/Bin" name discrepancy
        all_stats.append(cleaning(own_ss_tot, 'Per-period homeownership rates',
                                  'fracOwn', 'model'))
        all_stats.append(cleaning(own_ss, 'Per-period homeownership rates',
                                  'fracOwn', 'model'))

        genbins(view6)
        pol_exposed = view6.loc[view6['age'] > view6['period']]
        collapsed = (pol_exposed.groupby('period')
                     [['renters', 'alive']].sum())
        collapsed['hrate'] = 1 - collapsed['renters']/collapsed['alive']
        collapsed2 = (pol_exposed.groupby(['period', 'agebin'])
                      [['renters', 'alive']].sum())
        collapsed2['hrate'] = 1 - collapsed2['renters']/collapsed2['alive']

        all_stats.append(cleaning(collapsed['hrate'],
                                  'Per-period homeownership rate',
                                  'hrate', 'period'))
        all_stats.append(cleaning(collapsed2['hrate'],
                                  'Per-period homeownership rate',
                                  'hrate', 'period'))

# %>


# %< PART 5: TYPES OF MARGINAL BUYERS
def stats5(polType, durSS, polEnd):
    print("Part 5 \n \n")
    if polType == 'Repeat':
        cap = (view4.loc[view4['PolTaken'] == polType, ['nextDurables']]
               .describe().loc['max', 'nextDurables'])
        durBool = ((view4['durablesOrig'] <= cap) & (view4['policy'] < polEnd))
    else:
        durBool = (view4['policy'] == 0)
    durStats = view4.loc[durBool & (view4['PolTaken'] == polType),
                         ['id', 'ageBought', 'nextDurables', 'nextAssets',
                          'consumption']]
    merged = merge(durStats, durSS, on=['id', 'ageBought'], how='left',
                   suffixes=('', '_ss'))
    merged = merge(merged, view3, suffixes=('', '_ss'),
                   left_on=['id', 'ageBought'], right_on=['id', 'age'])
    merged['marginal'] = 0
    merged.loc[merged['durables'].isnull(), 'marginal'] = 1
    merged['ageSD'] = merged['ageBought']
    genbins(merged, col='ageBought')

    aggdict = {'id': len, 'consumption': np.mean, 'consumption_ss': np.mean,
               'nextDurables': np.mean, 'nextDurables_ss': np.mean,
               'nextAssets': np.mean, 'nextAssets_ss': np.mean,
               'income_val': np.mean, 'assets': np.mean,
               'ageBought': np.mean, 'age': np.median, 'ageSD': np.std}
    collapsed = merged.groupby(['marginal']).agg(aggdict)
    collapsed2 = merged.groupby(['agebin', 'marginal']).agg(aggdict)
    all_stats.append(cleaning(collapsed.loc[idx[0], :],
                              'Buyer characteristics', 0, 'index'))
    all_stats.append(cleaning(melt(collapsed2.loc[idx[:, 0], :].reset_index(),
                                   'agebin'),
                              'Buyer characteristics', 'value', 'variable'))

    all_stats.append(cleaning(collapsed.loc[idx[1], :],
                              'Buyer characteristics, marginal buyers',
                              1, 'index'))
    all_stats.append(cleaning(melt(collapsed2.loc[idx[:, 1], :].reset_index(),
                                   'agebin'),
                              'Buyer characteristics, marginal buyers',
                              'value', 'variable'))
# %>


# %< PART 6: METRICS COMPARING DIFFERNT FTHB AGE DISTRIBUTIONS
def stats6(mdir, polType):
    print("Part 6 \n \n")
    if polType == 'First-time':
        ageDist = read_table(mdir + 'FTHB_age_trans.txt', sep='\t',
                             header=None, names=['age', 'distPol', 'distSta'])
        ageDistAgg = ageDist.loc[3:19, :]
        ageDist = ageDist[(ageDist['age'] > 3) & (ageDist['age'] < 38)]
        print ageDistAgg
        genbins(ageDist)
        all_stats.append(cleaning([entropy(ageDistAgg['distPol'],
                                           ageDistAgg['distSta'])],
                                  'KL divergence of FTHB distributions', 0))
        all_stats.append(cleaning(ageDist.groupby('agebin').apply(lambda x:
                                  entropy(x['distPol'], x['distSta'])),
                                  'KL divergence of FTHB distributions',
                                  0, 'index'))

        all_stats.append(cleaning([hellinger2(ageDistAgg['distPol'],
                                              ageDistAgg['distSta'])],
                                  'Hellinger dist. of FTHB distributions', 0))
        all_stats.append(cleaning(ageDist.groupby('agebin').apply(lambda x:
                                  hellinger2(x['distPol'], x['distSta'])),
                                  'Hellinger dist. of FTHB distributions',
                                  0, 'index'))
# %>


if __name__ == '__main__':
    print('Generating aggregate statistics: \n \n')
    print(argv)
    stats1()
    try:
        print('num periods {}'.format(argv[6]))
        stats2(end=int(float(argv[6])))
    except:
        stats2()

    durShift, durSS = stats3(argv[2], float(argv[3]),
                             float(argv[4]), float(argv[5]))
    try:
        stats4(argv[2], durShift, durSS)
        stats5(argv[2], durSS, float(argv[5]))
    except:
        stats4(argv[2], durShift, durSS)
    stats6('output/%s/' % argv[1], argv[2])
    final = concat(all_stats).drop(['index'], axis=1)
    try:
        final = final.drop(['adjustment_ss'], axis=1)
    except:
        pass
    final = final.swaplevel(-2, -1).sort_index()
    final.to_csv('stats_final_%s.csv' % argv[1])
