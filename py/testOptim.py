import sys, traceback
import pandas as pd
import csv
import numpy as np
import scipy
from scipy.optimize import minimize
import random
from model.model_iterate import lifecycle_iterate
from main_calibration import *
import time

output_file = 'myOptim_out.csv'
#col = ['elasticity2', 'beta2',  'delta', 'dtau', 'tax1', 'tax2', 'F', 'F2', 'Dmin', 'elasticity', 'rentPrem', 'rentUtil']
#col = ['elasticity2', 'F', 'F2', 'Dmin', 'rentPrem', 'rentUtil', 'beta2', 'tax1', 'thetamatlab']
col = ['r', 'rborrow', 'rentPrem', 'rentUtil']
moments2019 = pd.read_csv("model/input_data/moments.csv", header=None).loc[0,:]

def scale_factor():
    return 10

def scale():
    return [1,1,1,1,0.2,1,1,1,1]

def myCol():
    #return ['rborrow', 'rentPrem', 'rentUtil']
    return ['r', 'rborrow', 'rentPrem', 'rentUtil']

def myStart():
    #return np.array([0.018, 7.70e-2, 0.93])
    return np.array([0.01, 0.018, 7.70e-2, 0.93])

def myMoments2019():
    return moments2019

def get_dict(guess_array):
    d = {'elasticity': 4, 'Dmin': 0.25}
    col = myCol()
    real_array = myStart()
    for i in range(len(guess_array)):
        if abs(guess_array[i] -real_array[i]) < 0.000001:
            continue
        else:
            d[col[i]] = guess_array[i]
    return d

def myMoments(myParams):
    """
    @param <list> myParams: input parameters
    Returns:
    - <list> modelOut: moments calculated from the model
    """
    d = get_dict(myParams)
    print(d)
    obj = lifecycle_iterate(d)
    time.sleep(3)
    obj.execSh(model='calib')
    modelOut, target, weights = momOut_9mom(obj.dir, ageScale=False)
    modelOut = list(modelOut.reshape([9]))
    return modelOut

def my_mse(x):
#    x /= scale_factor()
    mine = myMoments(x)
    #Rescale moments
    mine = np.array(mine)
    moments = myMoments2019() * scale()
    mine *= scale()
    print('my moments, actual moments, difference')
    print(mine)
    print(moments.values.tolist())
    print(mine-moments)
    op = (mine-moments)**2
    print("Current SSE")
    op = sum(op.values.tolist())
    print(op)
    return op


def main():
    col = myCol()
    random_guess = np.array(myStart())
#    random_guess *= scale_factor()
    try:
        res = scipy.optimize.minimize(my_mse, random_guess, method = 'COBYLA', options={'maxiter':100000})
        print(res)
        myLoss = [res.fun]
        myLoss.extend(res.x)
        with open(output_file, 'a') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow(myLoss)
    except Exception as e:
        print(e)

def main_old():
    #done = pd.read_csv(output_file, header=None).iloc[:,0].values.tolist()
    #toDo = list(filter(lambda x: x not in done, mySeeds))
    col = myCol()
    #bnds = ((0, 1), (0, None), (0, None), (0, None), (None, None), (None, None))
    #random_guess = np.array([.759, 0.06, 0.238, 0, 7.70e-2, 0.93, 0.94, 0.175, 0.20])
    random_guess = np.array(myStart())
    random_guess *= 10
    for i in toDo:
        try:
            print(i)
            def my_smm(x):
                print('Raw params:')
                print(x)
                print('Scaled params:')
                x = x/10
                print(x)
                return my_mse(x, i)
            res = scipy.optimize.minimize(my_smm, random_guess, method = 'COBYLA', options={'maxiter':100000})
            print(res)
            myLoss = [i, res.fun]
            myLoss.extend(res.x)
            with open(output_file, 'a') as f:
                csv_writer = csv.writer(f)
                csv_writer.writerow(myLoss)
        except Exception as e:
            print(e)
            #traceback.print_exc()
            print('Seed '+ str(i)+ ' failed')


main()
