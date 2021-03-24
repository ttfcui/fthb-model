FILENAME = 'fthb_model.csv'


def get_truncated_fthb(save=True):
    from model.model_iterate import lifecycle_iterate
    fthb = lifecycle_iterate.readModel('fthb').appended
    df = fthb[['age', 'id', 'income_val']]
    if save:
        df.to_csv(FILENAME, index=False)
    return df
