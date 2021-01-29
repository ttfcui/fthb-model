import os
import csv
import numpy as np
from model.model_iterate import lifecycle_iterate
from main_calibration import *
import pandas as pd

param = 'F2'
param_values = list(np.arange(0.025, 0.15, 0.025))
myName = param + '_moments.csv'
#done = pd.read_csv(myName, header=None).iloc[:,0].values.tolist()
#toDo = list(filter(lambda x: x not in done, param_values))
#print(done)
#print(toDo)
toDo = param_values

def scale():
    return [1,1,1,1,0.1,1,1,1,1]

for i in toDo:
    print(i)
    try:
        obj = lifecycle_iterate({param: i})
        obj.execSh(model='calib')
        modelOut, target, weights = momOut_9mom(obj.dir, ageScale=False)
        modelOut = list(modelOut.reshape([9]))
        mine = np.array(modelOut)
        target = list(target.reshape([9]))
        target = np.array(target)
        print('My moments, target moments')
        print(mine)
        print(target)
        mine *= scale()
        target *= scale()
        op = (mine-target)**2
        op = sum(op)
        print("Current SSE:")
        print(op)
        obj.matlabPlot(model=['calib'])
        name = 'F2_plots/' + param + '_' + str(i) + '_fracFthb.pdf'
        os.system('cp output/calib/fracFthb.pdf '+ name)
        with open(myName, 'a') as f:
            csv_writer = csv.writer(f)
            myOutput = [i]
            myOutput.extend(modelOut)
            csv_writer.writerow(myOutput)
    except Exception as e:
        print(e)
    break
