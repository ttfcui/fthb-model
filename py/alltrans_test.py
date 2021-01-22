from model.model_iterate import lifecycle_iterate
from numpy import nan
from pandas import merge, DataFrame, isnull

lc = lifecycle_iterate.appendModels('fthb', model=['experiment_monetary_nodown']).appended
lc2 = lifecycle_iterate.appendModels('transition_fthb', model=['experiment_monetary_nodown']).appended
life_adjust = (lc.loc[(lc['adjust']==2) & (lc['age'] <= 39),
               ['id', 'age', 'nextDurables']].set_index(['age', 'id']).unstack(-2))
life_adjust.columns = life_adjust.columns.droplevel()
lc2['policy'] = lc2['ageBought'] - lc2['age']
life_adjust2 = (lc2.loc[lc2['ageBought'] <= 39, ['id', 'age', 'policy', 'nextDurables']]
                .set_index(['policy', 'age', 'id']).unstack(-3).loc[2:39])
life_adjust2.columns = life_adjust2.columns.droplevel()

life_adjust = life_adjust.reset_index()
life_adjust2 = life_adjust2.reset_index()

matched = DataFrame()
for cohort in xrange(2, 40):
    sample = life_adjust2[life_adjust2['age'] == cohort]
    cohort_drop = min(-cohort + 2, -1)
    sample.iloc[:,cohort_drop:-1] = nan
    cf_sample = life_adjust.drop(range(1, cohort), axis=1)
    cf_sample.columns = ['id'] + range(0, cf_sample.shape[1] - 1)
    merged = merge(sample, cf_sample, on=['id'], suffixes=('', '_cf'))
    matched = matched.append(merged)
    
print(matched.iloc[:,37:].describe())
print(matched.loc[~isnull(matched['0'])].iloc[:,37:].describe())
