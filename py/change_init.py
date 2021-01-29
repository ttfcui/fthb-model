import pandas as pd
import os
import csv
import time
import numpy as np
import itertools
from model.model_iterate import lifecycle_iterate
from main_calibration import *

initial_grid = pd.read_csv("model/input_data/initial_grid_calibration_original.txt", sep='\t', lineterminator='\n', header=None)
shape = list(np.arange(0.25, 4, 0.5))
scale = list(np.arange(0.25, 4, 0.5))
combinations = list(itertools.product(shape, scale))
getdone = pd.read_csv("initAssets_moments.csv", header=None)
done = getdone.iloc[:,[0,1]].values.tolist()
toDo = list(filter(lambda x: list(x) not in done, combinations))
toDo = [[1,1]]


for i in toDo:
    print(i)
    try:
        start_time = time.time()
        initial_grid.iloc[:,2] *= i[0]
        initial_grid.iloc[:,3] *= i[1]
        initial_grid = initial_grid.round(10)
        initial_grid.to_csv("model/input_data/initial_grid_calibration.txt", sep='\t', header=False, index=False)
        obj = lifecycle_iterate({'Dmin': 0.2})
        obj.execSh(model='calib')
        modelOut, target, weights = momOut_9mom(obj.dir, ageScale=False)
        modelOut = list(modelOut.reshape([9]))
        target = list(target.reshape([9]))
        print('My moments, target moments')
        print(modelOut)
        print(target)
        obj.matlabPlot(model=['calib'])
        name = 'Dmin_plots/Dmin_' + str(i[0]) + '_' + str(i[1]) + '_fracFthb.pdf'
        os.system('cp output/calib/fracFthb.pdf '+ name)
        with open('initAssets_moments.csv', 'a') as f:
            csv_writer = csv.writer(f)
            myOutput = list(i)
            myOutput.extend(modelOut)
            csv_writer.writerow(myOutput)
        end_time = time.time()
        print(start_time-end_time)
    except Exception as e:
        print(e)
    #break

