import csv
#from subprocess import call
from model.model_iterate import lifecycle_iterate
from main_calibration import *

vals = 0.01
#d = {param: vals}
#params = [0.64670902600932323, 1.3159970003349724, 0.64184312623670015, 0.018115997502330707, -0.035177830159456261, 0.84610654323794254]
#d= dict(zip(['elasticity2', 'F', 'F2', 'Dmin', 'rentPrem', 'rentUtil'], params))
d={'Dmin': 0.25}
#d['seedvalue(:)']='{}'.format(6774)
obj = lifecycle_iterate(d)
obj.execSh(model='calib')
modelOut, target, weights = momOut_9mom(obj.dir, ageScale=False)
modelOut = list(modelOut.reshape([9]))
obj.matlabPlot(model=['calib'])
#call(['matlab', '-r', "T=38;POL=-1;transition_plots;exit;"])
print('My moments, target moments')
print(modelOut)
print(target)

param='new'
with open(param + '_moments.csv', 'a') as f:
	csv_writer = csv.writer(f)
	myOutput = [vals]
	myOutput.extend(modelOut)
	csv_writer.writerow(myOutput)

