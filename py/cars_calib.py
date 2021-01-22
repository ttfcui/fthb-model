from model.model_iterate import lifecycle_iterate
fthb = lifecycle_iterate.appendModels('fthb', model=['experiment_CARS']).appended
fthb['vintage'] = fthb['age'] - fthb['durTime']

def descriptive(data): 
    print(data['vintage'].loc[(data['vintage'] >= 0)].describe())
    print(data.loc[fthb['adjust'] > 1].shape[0]/
          float(data.loc[(fthb['rent']==0)].shape[0]))
    print(data.loc[:,'nextDurables'].describe())
    print(data.loc[(data['vintage']==2) ,'nextDurables'].describe())

descriptive(fthb.loc[fthb['age'] <= 39])
median = fthb.loc[fthb['age'] <= 39, 'income_val'].describe()['50%']
descriptive(fthb.loc[(fthb['age'] <= 39) & (fthb['income_val'] <= median)])
descriptive(fthb.loc[(fthb['age'] <= 39) & (fthb['income_val'] >= median)])
