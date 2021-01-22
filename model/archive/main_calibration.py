#!/usr/bin/python

from scipy.optimize import minimize, brute
from model.model_iterate import lifecycle_iterate
from pandas import get_dummies, merge, read_csv
from subprocess import call
import numpy as np

""" %< Wrapper functions for model calibration over multiple parameters """


def momOut_allHH(data, workmax):

    assert data['age'].max() == workmax  # max age in model
    print('test2.1')
    data = data.loc[data['age'] >= 2]
    data.loc[data['age'] <= 39, 'NetWToIncome'] = data['netWorth']/data['income_val']
    eqPrice = ((data['netWorth'] - data['nextAssets'])/data['nextDurables']
                ).max()
    print('test2.2')
    netwMed = data['NetWToIncome'].median()
    incomeMean = data['income_val'].mean()
    print('test2.3')
    mediansOwn = data.loc[data['adjust'] == 2, 'nextDurables'].median()

    adjcounts = data['adjust'].value_counts().sort_index()

    # Moments 1, 2, 3: Median net worth to income ratio, the
    # ratio of average owned house size to average income, and
    # two statistics about turnover - HH who adjust and all HH. 
    moments1 = np.zeros([1, 4])
    moments1[0, :] = [netwMed/incomeMean, (eqPrice*mediansOwn)/incomeMean,
                      adjcounts.iloc[2].sum(), float(adjcounts.sum())]
    del data['NetWToIncome']

    return moments1


def momOut_FTHBs(data):

    assert 'ageBought' in data.columns  # ensures this is dist_fthb file
    dataMom = data.loc[data['ageBought'] <= 39]

    # Moments 1 : Median FTHB Age
    mediansF = dataMom.median()['ageBought']

    # 28-30 pooled net 23-25 pooled
    ageagg = dataMom.groupby('ageBought')['id'].count()
    print(ageagg)
    calibinfo_mom1 = (np.log(ageagg.loc[6:9].mean()) -
                      np.log(ageagg.loc[1:4].mean()))

    # FTHB income distribution stats
    incdist = dataMom['income_val'].describe()
    calibinfo_mom2 = incdist.loc[['50%', '75%']].values.flatten()

    return np.concatenate(([calibinfo_mom1, mediansF],
			    calibinfo_mom2, [float(data.shape[0])]))


def momOut_lifeCycle(model):

    # Moment 1: The average homeownership rate for workers.
    momentsAll = np.transpose(model.genMoments(40, ['fracOwn'])
                            .as_matrix())[:,0:1]
     
    # Moments 5, 6: The average homeownership rate for people beneath 30
    # and people ages above 65. This step involves some data manipulation.
    # People < 30 and > 65 are disjoint, so slice the data to only include
    # those two groups and then use an appropriate bin size in genMoments
    # so each bin only includes one of the two disjoint cohorts.
    model.appended = model.appended.loc[
        (10 >= model.appended['age']) | (45 <= model.appended['age'])]
    momentsPartial = np.transpose(model.genMoments(35, ['fracOwn']).as_matrix())

    return np.concatenate((momentsAll, momentsPartial), axis=1)


def momOut_marginal(data):

    assert 'marginal' in data.columns  # ensures this is propensity file
    dataMom = data.loc[data['marginal'] == 1]

    # Semielas - ratio of marginal buyers over inframarg buyers
    datainfra = data.loc[data['infraMarg'] == 1]
    hinfo_Elas = float(dataMom.shape[0])/float(datainfra.shape[0])
    
    # Median FTHB Age of Marginal buyers
    mediansF = dataMom.median()['age']

    # IQR distribution for rental size
    hinfo = dataMom[['durables', 'nextDurables']].describe()
    hinfo_mom1 = hinfo.loc[['25%', '75%'], 'durables'].values.flatten()

    # Proportion of marginal buyers bunching
    # TODO: Include inframarg too?
    view_bunch = dataMom.loc[(abs(dataMom['nextDurables'] -
                            hinfo.loc['min', 'nextDurables']) < 1e-2)]
    hinfo_mom2 = np.transpose([view_bunch.shape[0]/float(dataMom.shape[0])])

    # Proportion of marginal buyers *not* scaling durable level greatly
    dataMom['durdiff'] = np.log(dataMom['nextDurables']) - np.log(dataMom['durables'])
    view_nochg = dataMom.loc[(abs(dataMom['durdiff']) < 2.5e-1)]
    hinfo_mom3 = np.transpose([view_nochg.shape[0]/float(dataMom.shape[0])])

    return np.concatenate(([hinfo_Elas, mediansF], hinfo_mom1, hinfo_mom2),
                           axis=0)


def momOut_9mom(csvDir, **kwargs):
    """

    """
    from csv import reader

    print('test1')
    if 'model' in kwargs:
        view = (lifecycle_iterate.appendModels(
                'fthb', model=[kwargs['model']]).appended)
    else:
        view = lifecycle_iterate.readModel('fthb').appended
    print('test2')
    mom1 = momOut_allHH(view, 54)

    if 'model' in kwargs:
        view2 = (lifecycle_iterate.appendModels(
                'dist_fthb', model=[kwargs['model']]).appended)
    else:
        view2 = lifecycle_iterate.readModel('dist_fthb').appended
    mom2 = np.reshape(momOut_FTHBs(view2), (-1, 5))
    # A moment is constructed here, share of FTHBs among all transactions
    mom1h = mom2[0,-1]/mom1[0,-2]
    mom1 = np.delete(mom1, [2, 3], axis=1)
    mom2 = np.delete(mom2, [4], axis=1)
    print(mom2)
    print(mom1)
    mom1 = np.reshape(np.append(mom1, mom1h), (-1, 3))

    lifeMom = lifecycle_iterate.readModel('lifecycleprofiles')
    if 'ageScale' in kwargs and kwargs['ageScale'] == True:
        ageDist = read_csv('%sinput_data/ageDist.csv' % csvDir)
        ageDist['ageDist'] /= 100
        lifeMom.appended = merge(lifeMom.appended, ageDist, on='age')
        lifeData = lifeMom.appended
        lifeData['N'] = lifeData['N']*(lifeData['ageDist']/
				       (lifeData['N']/lifeData['N'].sum()))
    mom3 = momOut_lifeCycle(lifeMom)
   
    # Rescaling of moments
    mom1[0, -1] = 5.0*mom1[0, -1]  # FTHB share
    mom2[0, 0] = 10.0*mom2[0, 0]  # Age gradient
    mom2[0, -2:] = 5.0*mom2[0, -2:]  # FTHB Income percentiles

    momFinal = np.concatenate((mom1, mom2, mom3*10), axis=1)
    print(momFinal)
    momFinal = np.delete(momFinal, 1, 1)  # The "weird" house value/income ratio
    readf = reader(open('%sinput_data/moments2019.csv' % csvDir, 'r'))
    target = np.array(readf.next()).astype(float)
    print((momFinal - target))
    return momFinal, target


def momOut_6mom(csvDir, ageScale, **kwargs):
    """

    """
    from csv import reader
    from pandas import merge, read_csv

    # Moments 1, 2, 3: Median net worth to income ratio, FTHB age, + the
    # ratio of average owned house size to average rented house size.
    if 'model' in kwargs:
        view = (lifecycle_iterate.appendModels(
                'fthb', model=[kwargs['model']]).appended)
    else:
        view = lifecycle_iterate.readModel('fthb').appended
    eqPrice = ((view['netWorth'] - view['nextAssets'])/view['nextDurables']
                ).max()
    view['NetWToIncome'] = view['netWorth']/view['income_val']
    medians = view.median()
    means = view.mean()
    if 'model' in kwargs:
        view2 = (lifecycle_iterate.appendModels(
                'dist_fthb', model=[kwargs['model']]).appended)
    else:
        view2 = lifecycle_iterate.readModel('dist_fthb').appended
    medians2 = view2.loc[view2['ageBought'] <= 40].median()
    medians3 = view.loc[view['adjust'] == 2, 'nextDurables'].median()

    lifeMom = lifecycle_iterate.readModel('lifecycleprofiles')
    if ageScale:
        ageDist = read_csv('%sinput_data/ageDist.csv' % csvDir)
        ageDist['ageDist'] /= 100
        lifeMom.appended = merge(lifeMom.appended, ageDist, on='age')
        lifeData = lifeMom.appended
        lifeData['N'] = lifeData['N']*(lifeData['ageDist']/
				       (lifeData['N']/lifeData['N'].sum()))

    moments1 = np.zeros([1, 3])
    moments1[0, :] = [medians['NetWToIncome']/means['income_val'],
                      medians2['ageBought'],
                      (eqPrice*medians3)/medians['income_val']]  # the ratios

    # Moment 4: The average homeownership rate for workers.
    moments2 = np.transpose(lifeMom.genMoments(40, ['fracOwn'])
                            .as_matrix())[:,0:1]
     
    # Moments 5, 6: The average homeownership rate for people beneath 30
    # and people ages above 65. This step involves some data manipulation.
    # People < 30 and > 65 are disjoint, so slice the data to only include
    # those two groups and then use an appropriate bin size in genMoments
    # so each bin only includes one of the two disjoint cohorts.
    lifeMom.appended = lifeMom.appended.loc[
        (10 >= lifeMom.appended['age']) | (45 <= lifeMom.appended['age'])]
    moments3 = np.transpose(lifeMom.genMoments(35, ['fracOwn']).as_matrix())

    # Moment 7: The turnover rate in stationary equilibrium - agents
    # who adjust owned housing as a ratio of all agents.
    view3 = view  # [view['age'] <= 40]
    moments4 = np.transpose([[view3['adjust'].value_counts()[2]/
                            float(view3.shape[0])]])
    print([view3['adjust'].value_counts()[2], float(view3.shape[0])])

     
    # Aggregate, match to target moments, output
    moments5 = np.concatenate((moments1, moments2*10, moments3*10, moments4*10
                               ), axis=1)
    readf = reader(open('%sinput_data/moments.csv' % csvDir, 'r'))
    target = np.array(readf.next()).astype(float)
    print(moments5)
    print((moments5 - target))
    return moments5, target


def momOut_extension(model, **kwargs):
    """

    """
    from csv import reader
    from pandas import merge, read_csv

    # Model *must* be specified in function call
    try:
        view = read_csv('propensity_%s.csv' % model)
    except IOError:
        raise IOError('Propensity file not found. Are you sure you are in'
                      ' the right directory?')
        
    view['fthbs'] = (view['infraMarg'] == 1) | (view['marginal'] == 1)
    view2 = view.loc[view['marginal'] == 1]

    # FTHB age distribution gradient
    view_ageagg = view.loc[view['infraMarg'] == 1].groupby('age')['id'].count()
    # 27-29 pooled net 22-24 pooled
    print(view_ageagg)
    calibinfo_mom1 = np.transpose([
        np.log(view_ageagg.loc[7:10].mean()) - np.log(view_ageagg.loc[2:5].mean())])

    # FTHB income distribution stats
    view_incdist = view2['income_val'].describe()
    calibinfo_mom2 = view_incdist.loc[['50%', '75%']].values.flatten()

    # IQR distribution for rental size
    hinfo = view2[['durables', 'nextDurables']].describe()
    hinfo_mom1 = hinfo.loc[['25%', '75%'], 'durables'].values.flatten()

    # Proportion of marginal buyers bunching
    # TODO: Include inframarg too?
    view_bunch = view2.loc[(abs(view2['nextDurables'] -
                            hinfo.loc['min', 'nextDurables']) < 1e-2)]
    hinfo_mom2 = np.transpose([view_bunch.shape[0]/float(view2.shape[0])])

    # Proportion of marginal buyers *not* scaling durable level greatly
    view2['durdiff'] = (view2['nextDurables'] - view2['durables'])/view2['durables']
    view_nochg = view2.loc[(abs(view2['durdiff']) < 1e-1)]
    hinfo_mom3 = np.transpose([view_nochg.shape[0]/float(view2.shape[0])])

    # Aggregate, match to target moments, output
    print(calibinfo_mom1)
    print(calibinfo_mom2)
    moments3 = np.concatenate((calibinfo_mom1*10, calibinfo_mom2), axis=0)
    moments3Add = np.concatenate((hinfo_mom1, hinfo_mom2, hinfo_mom3), axis=0)
    print(moments3)
    return moments3, moments3Add


def calibration_9mom(params, *args):
    """ Calibrates the FTHB lifecycle model allowing five parameters
        to vary: EIS (elasticity), rent premium above the user
        cost of housing (rentPrem), initial income poisson parameter
        (poisMean), Minimum size of house (Dmin) and exogenous moving
        probability for renters (movprobR).
     
        The parameters are overidentified using nine moments
        observed in the model and data.
    Args:
        params: The values of the five parameters to be calibrated.
        args[0]: The "steady convergence threshold:" once guesses for the
            market clearing price changes less than this threshold in the
            Brent minimization algorithm, the algorithm stops. Should
            increase in magnitude when iterating this function.
        args[1]: The guess for the market clearing price.
        args[2], args[3]: The permitted maximum deviation, positive and
            negative, the actual market clearing price will be from the
            guess. A tighter bound means fewer iterations in finding the
            actual price.
        args[4]: If True, calibrates each age cohort's size to its proportion
            in the 2010 Census age distribution. This inflates the value
            of younger cohorts at the expense of near-retirement cohorts.
    """
    # Manually set upper bounds for the parameters (SUBJECT TO CHANGE).
    # If the minimization algorithm is unconstrained, still allows you
    # to input algorithm starting values using values from 0 to 1.
    scale = [1.0, 1.0e-2, 1.0, 1.0e-1, 1.0]
    ageScale = args[4]
    params = np.multiply(params, scale)
    print(params)
     
    # Run the model and get output
    iterDict = {'agridsize': 100, 'Dgridsize': 55,
                'zgridsize': 19, 'steady_conv_thres': args[0],
                'hpnodes(1)': args[1], 'EligYrsF': 3,
                'adjTransfer': 0.075,
                'startprice(1,1)': args[2],
                'startprice(2,1)': args[3],
                'elasticity': params[4], 'rentPrem': params[1],
                'beta2': params[0], 'rentUtil': params[2],
                'F2': params[3],
	        'ge_start': '.FALSE.', 'pe_start': '.TRUE.'}
    if (args[4]):
        iterDict['ss_only'] = '.TRUE.'
    else:
        iterDict['ss_only'] = '.FALSE.'

    modelCal = lifecycle_iterate(iterDict)
    modelCal.execSh(model='calib')
    modelOut, target = momOut_9mom(modelCal.dir, ageScale=ageScale)
    propsList = ['0.20', '%f' % args[1]]
    call(['./gen_propensity.py', 'calib', 'First-time'] + propsList)
    try:
        view = read_csv('propensity_calib.csv')
        momExtra = momOut_marginal(view)
    except IOError:
        raise IOError('Propensity file not found. Are you sure you are in'
                      ' the right directory?')
 
    try:
        np.savetxt(momFile, params, '%16.6f', ' ', '\n')
        np.savetxt(momFile, modelOut, fmt='%16.6f')
        np.savetxt(momFile, (modelOut - target), fmt='%16.6f')
        np.savetxt(momFile, np.array([momExtra]), fmt='%16.6f')
        np.savetxt(momFile, np.sum((modelOut - target)**2))
    except:
        pass
    return(np.sum((modelOut - target)**2))






def calibration_6mom(params, *args):
    """ Calibrates the FTHB lifecycle model allowing five parameters
        to vary: EIS (elasticity), rent premium above the user
        cost of housing (rentPrem), lump-sum transfer in retirement
        (ret_wealth), Minimum size of house (Dmin) and weight of
        bequests in lifetime utility (psi).
     
        The parameters are exactly identified using five simple moments
        observed in the model and data.
    Args:
        params: The values of the five parameters to be calibrated.
        args[0]: The "steady convergence threshold:" once guesses for the
            market clearing price changes less than this threshold in the
            Brent minimization algorithm, the algorithm stops. Should
            increase in magnitude when iterating this function.
        args[1]: The guess for the market clearing price.
        args[2], args[3]: The permitted maximum deviation, positive and
            negative, the actual market clearing price will be from the
            guess. A tighter bound means fewer iterations in finding the
            actual price.
        args[4]: If True, calibrates each age cohort's size to its proportion
            in the 2010 Census age distribution. This inflates the value
            of younger cohorts at the expense of near-retirement cohorts.
    """
    # Manually set upper bounds for the parameters (SUBJECT TO CHANGE).
    # If the minimization algorithm is unconstrained, still allows you
    # to input algorithm starting values using values from 0 to 1.
    scale = [3.0, 1.0e-2, 1.0, 1.0, 2.0, 1.0]
    ageScale = args[4]
    params = np.multiply(params, scale)
    print(params)
     
    # Run the model and get output
    modelCal = lifecycle_iterate({'agridsize': 120, 'Dgridsize': 55,
                                  'zgridsize': 13, 'steady_conv_thres': args[0],
                                  'hpnodes(1)': args[1], 'EligYrsF': 3,
                                  'startprice(1,1)': args[2],
                                  'startprice(2,1)': args[3],
                                  'elasticity': params[5], 'rentPrem': params[1],
                                  'Dmin': params[2], 'rentUtil': params[3],
                                  'ret_wealth': params[4], 'psi': params[0]})
    modelCal.execSh()
    modelOut, target = momOut_6mom(modelCal.dir, ageScale)
    propsList = ['0.20', '%f' % args[1]]
    call(['./gen_propensity.py', 'calib', 'First-time'] + propsList)
    try:
        view = read_csv('propensity_calib.csv')
        momExtra = momOut_marginal(view)
    except IOError:
        raise IOError('Propensity file not found. Are you sure you are in'
                      ' the right directory?')
 
    try:
        np.savetxt(momFile, params, '%16.6f', ' ', '\n')
        np.savetxt(momFile, modelOut, fmt='%16.6f')
        np.savetxt(momFile, (modelOut - target), fmt='%16.6f')
        np.savetxt(momFile, np.array([momExtra]), fmt='%16.6f')
        np.savetxt(momFile, np.sum((modelOut - target)**2))
    except:
        pass
    return(np.sum((modelOut - target)**2))


def calibration_2mom(params, *args):
    """ Calibrates the FTHB lifecycle model allowing two parameters
       to vary: rent premium above the user cost of housing (rentPrem)
       and disutility from housing (rentUtil).

       The parameters are overidentified using three simple moments
       observed in the model and data.
    Args:
       params: The values of the two parameters to be calibrated.
       args[0]: The "steady convergence threshold:" once guesses for the
           market clearing price changes less than this threshold in the
           Brent minimization algorithm, the algorithm stops. Should
           increase in magnitude when iterating this function.
       args[1]: The guess for the market clearing price.
       args[2], args[3]: The permitted maximum deviation, positive and
           negative, the actual market clearing price will be from the
           guess. A tighter bound means fewer iterations in finding the
           actual price.
    """
    from csv import reader
    # Manually set upper bounds for the parameters (SUBJECT TO CHANGE).
    # If the minimization algorithm is unconstrained, still allows you
    # to input algorithm starting values using values from 0 to 1.
    scale = [1.0e-2, 1.0e-1]
    params = np.multiply(params, scale)
    print(params)

    # Run the model and get output
    modelCal = lifecycle_iterate({'agridsize': 120, 'Dgridsize': 55,
                                  'zgridsize': 19, 'steady_conv_thres': args[0],
                                  'hpnodes(1)': args[1], 'EligYrsF': 3,
                                  'startprice(1,1)': args[2],
                                  'startprice(2,1)': args[3],
                                  'F2': params[1], 'rentPrem': params[0],
                                  'beta2': args[4], 'rentUtil': args[5],
                                  'elasticity': args[6], 'ge_start': '.FALSE.',
				  'ss_only': '.FALSE.', 'pe_start': '.TRUE.'})
    modelCal.execSh()
    modelOut, target = momOut_9mom(modelCal.dir, ageScale=True)
    propsList = ['0.20', '%f' % args[1]]
    call(['./gen_propensity.py', 'calib', 'First-time'] + propsList)
    try:
        view = read_csv('propensity_calib.csv')
        momExtra = momOut_marginal(view)
    except IOError:
        raise IOError('Propensity file not found. Are you sure you are in'
                      ' the right directory?')
    print(momExtra)
 
    modelOut = modelOut[0,np.array([2, 6, 7])]
    target = target[np.array([2, 6, 7])]
    print((modelOut - target))
    """
    lifeMom = lifecycle_iterate.readModel('lifecycleprofiles')

    # Moment 1: The average homeownership rate.
    moments1 = np.transpose(lifeMom.genMoments(40, ['fracOwn'])
                            .as_matrix())[:,0:1]
    moments1 = np.nan_to_num(moments1)  # safety measure

    # Moments 2, 3: The average homeownership rate for people beneath 30
    # and people over 65. This step involves some data manipulation.
    # People < 30 and > 65 are disjoint, so slice the data to only include
    # those two groups and then use an appropriate bin size in genMoments
    # so each bin only includes one of the two disjoint cohorts.
    lifeMom.appended = lifeMom.appended.loc[
        (10 >= lifeMom.appended['age']) | (50 < lifeMom.appended['age'])]
    moments2 = np.transpose(lifeMom.genMoments(35, ['fracOwn']).as_matrix())
    moments2 = moments2[0:1,0:1]
    moments3 = np.concatenate((moments1*10, moments2*10), axis=1)
    
    readf = reader(open('%sinput_data/moments.csv' % modelCal.dir, 'r'))
    target = np.array(readf.next()).astype(float)[3:-2]
    print((moments3 - target))
    # return((moments3 - target).flatten())
    return(np.sum((moments3 - target)**2))
    """
    return(np.sum(np.abs(modelOut - target)))


def calibrate_grid(steadyprice, prices):
    """ This is the big minimization call. Notice args[2], args[3], i.e. we
        think the market clearing price is no more than +/-n% of the test price
        we computed.

        These bounds on market clearing prices are also an implicit way to
        constrain the minimization problem. The Fortran model will quit before
        outputting simulation stats if a price is found not because it clears
        the market, but because the price hit the boundaries. Hence stats will
        not change relative to the last iteration, and the surface we minimize
        over is a bump function of sorts.
    """
    # Set up uneven grid for calibration: Variables are phi, rentPrem,
    # Dmin, rentUtil, ret_wealth and elasticity, in that order.
    # TODO: Instead of these grids generate a series of random numbers
    # (so like basinhopping)
    pRanges = (slice(0.865, 0.95, 0.025), slice(7.75, 9.00, 0.45),
               slice(0.90, 0.93, 0.04), slice(1.30, 2.15, 0.25),
               slice(2.50, 4.4, 1.0))

    res = brute(calibration_9mom, pRanges, args=(8e-2, steadyprice, prices[0],
                prices[1], False), full_output=True, finish=None)
    calibration_9mom(res[0], 5e-3, steadyprice, prices[0], prices[1], False)
    return res[0]


def calibrate_alg(params, steadyprice, prices):
    """ For the moment we use Nelder-Mead to find the minimizing parameter
        vector, because it spends less time evaluating the function than other
        gradient-based methods.
    """
    simplex = np.zeros([3, 2])
    simplex[0, :] = [8.00, 2.00]
    simplex[1, :] = [8.10, 3.00]
    simplex[2, :] = [7.50, 2.40]

    res = minimize(calibration_2mom, simplex[0, :],
                   args=(1e-2, steadyprice, prices[0], prices[1],
                         params[0], params[1], params[2]),
                   method='Nelder-Mead', options={'fatol': 1e-4,
                   'xatol': 1e-4, 'initial_simplex': simplex})
    print(res)
    return res.x


def calibrate_check(params, steadyprice, prices, recalib=True):
    """

    """
    # ModelMin feeds in the minimizing parameter found by earlier functions
    modelMin = lifecycle_iterate({'agridsize': 120, 'Dgridsize': 55,
                                  'zgridsize':19, 'steady_conv_thres': 1e-2,
                                  'hpnodes(1)': steadyprice, 'EligYrsF': 3,
                                  'numhouseholds': 50000,
                                  'adjTransfer': 0.075,
                                  'startprice(1,1)': prices[0],
                                  'startprice(2,1)': prices[1],
                                  'elasticity': params[4], 'rentPrem': params[1],
                                  'beta2': params[0], 'rentUtil': params[2],
                                  'F2': params[3], 'ge_start': '.FALSE.',
                                  'ss_only': '.FALSE.', 'pe_start': '.TRUE.'})
    if recalib:
        modelMin.execSh(model='calib')
        momOut_9mom(modelMin.dir, ageScale=False)

    # Get out lifecycle moments also generated from SCF calibration, split
    # over 5-year bins. Feed the CSVs into Stata file moment_graphs.do
    modelSS = lifecycle_iterate.readModel('lifecycleprofiles')
    data = modelSS.appended
    modelSS.appended['avgAssetsExDebt'] = (data['avgAssetsExDebtOwners'] *
        data['fracOwn'] + data['avgAssetsRental'] * (1-data['fracOwn']))

    momH = modelSS.genMoments(1, ['avgHouseOwners', 'avgAssetsExDebtOwners'],
                               'NOwn')
    mom = modelSS.genMoments(1, ['avgAssetsExDebt', 'fracOwn'])
    momH.to_csv('momentsH_1yearbin.csv')
    mom.to_csv('moments_1yearbin.csv')

    # Read from file recording every agent, not just FTHBs. Get raw data
    # for plotting wealth among renters and owners.
    modelSS = lifecycle_iterate.readModel('fthb')
    modelSS.wealthLorenz()
    view = modelSS.appended
    view[['renter', 'homeown', 'adjuster']] = get_dummies(view['adjust'])
    view = view.loc[view['age'] <= 39]
    view['adjcount'] = view.groupby('id')['adjuster'].transform(lambda x: x.cumsum())
    holdtimes = view.loc[view['rent'] == 0.0].groupby(['id', 'adjcount'])['adjust'].count().reset_index()
    print 'Average durable holding time: %6.2f' % holdtimes.groupby('id')['adjust'].mean().mean()
    print 'Median durable holding time: %6.2f' % holdtimes['adjust'].median()

    indivAssets = view.loc[(view['age'] < 54), ['id', 'age', 'rent',
                           'nextAssets', 'nextDurables', 'netWorth',
                           'income', 'income_val']]
    indivAssets.to_csv('indivAssets.csv', index=False)

    modelMin.matlabPlot(model=['calib'])

    # Output semielasticity estimates as in calibration
    propsList = ['0.20', '%f' % steadyprice]
    call(['./gen_propensity.py', 'calib', 'First-time'] + propsList)
    try:
        view = read_csv('propensity_calib.csv')
        momExtra = momOut_marginal(view)
    except IOError:
        raise IOError('Propensity file not found. Are you sure you are in'
                      ' the right directory?')
    print(momExtra)
 
# %>

if __name__ == "__main__":
    # Get an idea of what the steady state price should be.
    # init_params = [0.80, 3.000, 0.92, 1.00, 3.50]
    # intro = calibration_9mom(init_params, 1e-3, 0.00, 0.25, 0.00, True)
    steadyprice = 0.1000 # From checking model_log.txt

    """
    momFile = open('moments_102019v2.txt', 'a+')
    res = calibrate_grid(steadyprice, (0.08, -0.08))
    momFile.close()
    print(res)
    quit() 

    BEST PARAMETERS:

    res = [0.94, 8.20, 0.93, 2.05, 2.50]

    params_v1 = np.multiply([res[0], res[2], res[-1]], [1.0, 1.0, 1.0])
    res2 = calibrate_alg(params_v1, steadyprice, (0.05, -0.05))
    res2 = np.multiply([res2[0], res2[1]], [1.0e-2, 1.0])
    print res2
    """
    res_final = np.multiply([0.94, 7.70, 0.93, 2.38, 2.50],
                            [1.0, 1.0e-2, 1.0, 1.0e-1, 1.0])
    print res_final
    calibrate_check(res_final, steadyprice, (0.05, -0.60), False)

"""
res = minimize(calibration_6mom, np.zeros([1, 5]), method='Nelder-Mead',
               options={'fatol': 1e-4, 'initial_simplex': simplex},
               args=(1e-2, steadyprice, 0.03, -0.03))
simplex[0,:] = [0.05, 0.8]
simplex[1,:] = [0.08, 0.9]
simplex[2,:] = [0.2, 0.95]

res = minimize(calibration_2mom, np.zeros([1,2]), method='Nelder-Mead',
               options={'fatol':1e-5, 'initial_simplex': simplex},
               args=(steadyprice, steadydeviate[0], steadydeviate[1]))

res = least_squares(calibration_4mom, [0.98, 0.5, 0.75, 0.6], diff_step=1e-3,
                    bounds=([0.95, 0, 0, 0], [1, 1, 1, 1]), ftol=1e-4,
                    args=(steadyprice, steadydeviate[0], steadydeviate[1]))
"""
