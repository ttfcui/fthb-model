#!/usr/bin/python

from scipy.optimize import minimize, brute
from model.model_iterate import lifecycle_iterate
from pandas import get_dummies, merge, read_csv, read_table, concat
from subprocess import call
import numpy as np


""" Wrapper functions for moment calculations """
# %<

def momOut_allHH(data, workmax):
    """ Generates moments computed over all simulated agents (+ renters
        and retirees), hence "allHH."

    Args:
        data: A preloaded output file from the last model run, namely
        the "fthb.txt" file.
    """

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
    """ Generates moments computed over only FTHBs that are not retired,
        hence "FTHBs."

    Args:
        data: A preloaded output file from the last model run, namely
        the "dist_fthb.txt" file.
    """

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
    """ Generates homeownership moments binned over subpopulations,
        hence "lifeCycle."

    Args:
        model: *Not* an output file, but a call to a lifecycle_iterate
        moment. That class has a method that computes the moments.
    """

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
    """ Generates homeownership moments generated for just marginal
        FTHBs incentivized by the policy, hence "marginal."

    Args:
        data: A preloaded output file from the last model run, namely
        the latest "propensity_*" file.
    """

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
    """ Generates homeownership moments generated for just marginal
        FTHBs incentivized by the policy, hence "marginal."

    Args:
        csvDir: Directory with all the model f90 files and the CSV files.
        **model: This function allows for arbitrary models outputted
            after a model execution to be called.
        **ageScale: If True, calibrates each age cohort's size to its proportion
            in the 2010 Census age distribution. This inflates the value
            of younger cohorts at the expense of near-retirement cohorts.
    """
    from csv import reader

    if 'model' in kwargs:
        view = (lifecycle_iterate.appendModels(
                'fthb', model=[kwargs['model']]).appended)
    else:
        view = lifecycle_iterate.readModel('fthb').appended
    mom1 = momOut_allHH(view, 54)
    print('Processed all HH moments')

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
    print('Processed all FTHB moments')

    lifeMom = lifecycle_iterate.readModel('lifecycleprofiles')
    if 'ageScale' in kwargs and kwargs['ageScale'] == True:
        ageDist = read_csv('%sinput_data/ageDist.csv' % csvDir)
        ageDist['ageDist'] /= 100
        lifeMom.appended = merge(lifeMom.appended, ageDist, on='age')
        lifeData = lifeMom.appended
        lifeData['N'] = lifeData['N']*(lifeData['ageDist']/
				       (lifeData['N']/lifeData['N'].sum()))
    mom3 = momOut_lifeCycle(lifeMom)
    print('Processed all HO Rate moments')
   
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

# %>


""" Functions for individual calls in calibration process """
# %<
def calibration_execModel(paramsDict, *args):
    """ Given a dict of calibrated params specified earlier, runs the model
        with this specification and outputs simulated moments.

    Args:
        paramsDict: Dictionary of parameter values to be calibrated,
            with the key being the name of the parameter.
        args[0]: The "steady convergence threshold:" once guesses for the
            market clearing price changes less than this threshold in the
            Brent minimization algorithm, the algorithm stops. Should
            increase in magnitude when iterating this function.
        args[1]: The guess for the market clearing price.
        args[2], args[3]: The permitted maximum deviation, positive and
            negative, the actual market clearing price will be from the
            guess. A tighter bound means fewer iterations in finding the
            actual price.
        args[4]: If True, skips over the policy simulation half of the model
            execution. This means policy period moments outputted afterwards
            are not appropriate for calibration.
    """
    algoDict = {'agridsize': 100, 'Dgridsize': 55,
                'zgridsize': 19, 'steady_conv_thres': args[0],
                'hpnodes(1)': args[1], 'EligYrsF': 3,
                'adjTransfer': 0.075,
                'startprice(1,1)': args[2],
                'startprice(2,1)': args[3],
	        'ge_start': '.FALSE.', 'pe_start': '.TRUE.'}
    if (args[4]):
        paramsDict['ss_only'] = '.TRUE.'
    else:
        paramsDict['ss_only'] = '.FALSE.'

    # With both calibration parameters and algorithm params chosen...
    paramsDict.update(algoDict)
    # Model runs after this line
    modelCal = lifecycle_iterate(paramsDict)
    # Extracting all relevant moments
    modelCal.execSh(model='calib')
    modelOut, target = momOut_9mom(modelCal.dir, ageScale=False)
    propsList = ['0.20', '%f' % args[1]]
    call(['./gen_propensity.py', 'calib', 'First-time'] + propsList)
    try:
        view = read_csv('propensity_calib.csv')
        momExtra = momOut_marginal(view)
    except IOError:
        raise IOError('Propensity file not found. Are you sure you are in'
                      ' the right directory?')
    return modelOut, momExtra, target
 

def calibration_9mom(params, *args):
    """ Calibrates the FTHB lifecycle model allowing five parameters
        to vary: EIS (elasticity), rent premium above the user
        cost of housing (rentPrem), quadratic transaction cost term (F2)
        (poisMean), Minimum size of house (Dmin) and disutility of
        rental housing services (rentUtil).
     
        The parameters are overidentified using nine moments
        observed in the model and data.

    Args:
        params: The values of the five parameters to be calibrated.
        args: See calibration_execModel for details.
    """
    # Manually set upper bounds for the parameters (SUBJECT TO CHANGE).
    # If the minimization algorithm is unconstrained, still allows you
    # to input algorithm starting values using values from 0 to 1.
    scale = [1.0, 1.0e-2, 1.0, 1.0e-1, 1.0]
    params = np.multiply(params, scale)
    print(params)
     
    # Run the model and get output
    iterDict = {'elasticity': params[4], 'rentPrem': params[1],
                'beta2': params[0], 'rentUtil': params[2],
                'F2': params[3]}
    modelOut, momExtra, target = calibration_execModel(
            iterDict, *args)

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

       The parameters are overidentified using four simple moments
       observed in the model and data.

    Args:
       params: The values of the two parameters to be calibrated.
       args[5:]: These aren't algorithmic parameters, but calibrated
           parameters from an earlier loop step that we hold fixed.
           Hence, even though they're passed through to
           calibration_execModel, they're left unused.
    """
    # Manually set upper bounds for the parameters (SUBJECT TO CHANGE).
    # If the minimization algorithm is unconstrained, still allows you
    # to input algorithm starting values using values from 0 to 1.
    scale = [1.0e-2, 1.0e-1]
    params = np.multiply(params, scale)
    print(params)

    # Run the model and get output
    iterDict = {'F2': params[1], 'rentPrem': params[0],
                'beta2': args[5], 'rentUtil': args[6],
                'elasticity': args[7]}
    modelOut, momExtra, target = calibration_execModel(
            iterDict, *args)

    # Pick only the FTHB/HO rates moments we set as targets
    modelOut = modelOut[0,np.array([2, 6, 7])]
    target = target[np.array([2, 6, 7])]
    print((modelOut - target))
    return(np.sum(np.abs(modelOut - target)))

# %>


""" Wrapper functions for calibration iterations """
# %<
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

    Args:
        steadyprice: The approximately market-clearing house price level
            calculated at start of the file.
        prices: If GE simulations are run, this governs the interval
            of prices over which the market clearing price is searched.
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


def clean_calibgrid(calibfile):
    """ Small data cleaning function taking iterated Fortran outputs
        of model moments into more readable spreadsheets.

    Args:
        calibfile: Name of the file with iterated moments.
    """

    # We deliberate don't parse the space delimiters in orig file
    infile = read_table(calibfile, sep=',', header=None)
    # Why "8"? 5 lines of parameters + 3 lines of moments
    params = infile.loc[infile.index % 8 < 5].reset_index()
    params['iter'] = np.floor(params.index / 5.0)
    params = (params.groupby('iter')[0].apply(' '.join)
              .str.split('\s+', expand=True))

    # Parsing of moment rows
    mom1 = infile.loc[infile.index % 8 == 5].reset_index()
    mom1 = mom1.iloc[:,-1].str.split('\s+', expand=True)
    mom2 = infile.loc[infile.index % 8 == 6].reset_index()
    mom2 = mom2.iloc[:,-1].str.split('\s+', expand=True)
    mom3 = infile.loc[infile.index % 8 == 7].reset_index()
    mom3 = mom3.iloc[:,-1].str.split('\s+', expand=True)
    
    # Output cleaned spreadsheets
    concat([params, mom1, mom3], axis=1, ignore_index=True).to_csv(
            'calib_moments.csv')
    concat([params, mom2], axis=1, ignore_index=True).to_csv(
            'calib_diffs.csv')


def calibrate_alg(params, steadyprice, prices):
    """ For the moment we use Nelder-Mead to find the minimizing parameter
        vector, because it spends less time evaluating the function than other
        gradient-based methods.

    Args:
        params: Calibrated parameters from an earlier loop step
            that we hold fixed.
        steadyprice: The approximately market-clearing house price level
            calculated at start of the file.
        prices: If GE simulations are run, this governs the interval
            of prices over which the market clearing price is searched.
    """
    simplex = np.zeros([3, 2])
    simplex[0, :] = [8.00, 2.00]
    simplex[1, :] = [8.10, 3.00]
    simplex[2, :] = [7.50, 2.40]

    res = minimize(calibration_2mom, simplex[0, :],
                   args=(1e-2, steadyprice, prices[0], prices[1], False,
                         params[0], params[1], params[2]),
                   method='Nelder-Mead', options={'fatol': 1e-4,
                   'xatol': 1e-4, 'initial_simplex': simplex})
    print(res)
    return res.x


def calibrate_check(params, steadyprice, prices, recalib=True):
    """ Runs the model with final calibrated parameters one more time,
        generating moments and distributional stats to be pushed into
        steady_programs.do.

    Args:
        params: Calibrated parameters from an earlier loop step
            that we hold fixed.
        steadyprice: The approximately market-clearing house price level
            calculated at start of the file.
        prices: If GE simulations are run, this governs the interval
            of prices over which the market clearing price is searched.
        recalib: If True, the model is rerun beforehand. Mostly used
            to debug this function.

    """
    # ModelMin feeds in the minimizing parameter found by earlier functions
    modelMin = lifecycle_iterate({'agridsize': 120, 'Dgridsize': 55,
                                  'zgridsize': 19, 'steady_conv_thres': 1e-2,
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
        modelMin.execSh()
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
    clean_calibgrid(momFile)   
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
    calibrate_check(res_final, steadyprice, (0.05, -0.60))
