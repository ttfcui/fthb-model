import pandas as pd
from model_iterate import lifecycle_iterate


baseFile = 'lifecycle_TC04032016.f90'

for model, value in [('experiment_full', 1), ('experiment_little', 0.02),
                     ('experiment_theta', 0.15), ('experiment_zero', 0)]:
#for model, value in [('experiment_full_delay', 1)]:
    T = 1
    modelCal = lifecycle_iterate({'r_rentalO': 0.065, 'elasticity2O': 0.85,
				  'hpmin': -0.04, 'numhouseholds': 15000,
                                  'rentelasticity': 0,
				  'Probhp(1,1)': 1.0 - 1.0 / T, 
				  'Probhp(1,2)': 1.0 / T, 
				  'PolYrs': T, 'downtransfer': value},
				  baseFile)
    modelCal.execSh(model=model)


model = lifecycle_iterate({}, baseFile)

orig = model.appendModels('fthb', model=['experiment_little']).appended

orig = orig.loc[orig['age'] <= 60].set_index(['id', 'age']).sort_index()
orig['durablesl'] = (orig.groupby(level=0)['durables']
		     .transform(lambda x: x.shift()))
orig['durablesl'] = orig['durablesl'].fillna(method='bfill')
trans = (orig.loc[abs(orig['durables'] - orig['durablesl']) > 2e-2]
	 .groupby(level=0)['durables'].count())
pd.value_counts(trans).sort_index()

transR = (orig.loc[(orig['durables'] == 0) & (orig['durablesl'] != 0)]
	  .groupby(level=0)['durables'].count())
pd.value_counts(transR).sort_index()


orig = model.appendModels('fthb', model=['experiment_little']).appended

B = orig.loc[orig['durables'] == 0].groupby('id')['durables'].count()
orig_c = orig.loc[orig['id'].isin(B[B > 39].index)]
orig_c = orig_c.set_index(['id', 'age']).sort_index()
orig_c.index.names = ['id', 'ageBought']

orig_p = model.appendModels('dist_fthb', model=['experiment_little']).appended
orig_p.set_index(['id'], inplace=True)

model.appendModels('transition_fthb', model=['experiment_little'])

merged = pd.merge(orig_c, model.appended, left_index=True, 
                  right_on=['id', 'ageBought'])
len(merged.loc[(merged['durables_x'] == 0) & (merged['durables_y'] != 0),
               'id'].unique())

merged = pd.merge(orig_p, model.appended, left_index=True, 
                  right_on=['id'])
merged.groupby(['ageBought_x', 'ageBought_y'])['id'].count()
