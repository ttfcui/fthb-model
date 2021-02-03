#!/usr/bin/python

from types import *
from os import remove
import pandas as pd
import re
import pdb
import glob


class lifecycle_iterate:

    """ Class for testing the housing model first used in Berger, Guerrieri,
        Lorenzoni and Vavra, now used in Berger, Turner and Zwick.

    Initialization requires a dictionary of which parameters are to be
    adjusted. There is a lot to play with: calibrated and uncalibrated,
    Substantive and algorithmic. From there, the code can be compiled
    and relevant output files are categorized. After a multitude of models
    have been solved, class methods can be used to aggregate moments
    and statistics across models.

    Attrs:
        self.dir: the directory with the original model code.
        self.mainD: the main repository directory for code and project
            outputs, accessible through subdirectories.
        self.Keys: the parameters adjusted from the baseline model
            in this call.
        self.Vals: new values for each parameter stated in self.Keys.
        self.OKFiles: Output files for which aggregation methods are
            currently supported. Only these files get copied when
            looping over different model specifications.

    Args:
        paramSub: Dictionary of parameters to be changed. The key is the
            exact name of the parameter, most likely to be found in module
            share. The value is the new value, up to six
            significant digits.
        file_path: The name of the latest original model code. By "latest"
            it is meant code that does not just have changed parameters.
        new: Boolean. If True, generates all files into which the Fortran model
            exports outputs. Intended when the model is run for the first
            time on a machine.
        sub: Boolean. If True, swaps out new parameters from the
            initialization file share.f90, reparametrizing the model.

    """

    # Regex search. First group is the old parameter value.
    # Second group are inline comments, if any.
    def paramReplace(self, paramSpec, newParam, line):
        paramFlag = paramSpec[1]
        if ('%s=' % paramSpec[0] in line):
            if isinstance(newParam, float):
                newParam = '%.6f' % newParam
            line = re.sub('\s*=([0-9\-\.\+\s]*|\.[TF][A-Z]\+.*)[^!]*(?P<comm>!|\n)',
                          '= %s \g<comm>' % newParam, line)
            paramFlag = True
        else:
            pass
        return line, paramFlag

    def __init__(self, paramSub, file_path='share.f90', new=False, sub=True):
        from tempfile import mkstemp
        from shutil import move
        from os import getcwd, close, chmod

        assert type(paramSub) is DictType
        assert all(type(key) is StringType for key, val in paramSub.items())
        self.Keys = paramSub.keys()
        self.Vals = paramSub.values()
        self.dir = '%s/model/' % getcwd()
        self.mainD = '%s/' % getcwd()

        if sub:
            fh, abs_path = mkstemp()
            own_path = self.dir + file_path
            with open(abs_path, 'w+') as newFile:
                with open(own_path, 'r+') as oldFile:
                    # This is where the format of Fortran code is exploited.
                    # Only interested in the first parameter call, so once
                    # that is found and replaced the value is changed to
                    # indicate search is unnecessary.
                    for line in oldFile:
                        # If all values changed to True just parse the rest
                        # of the file (few hundred lines)
                        if [phrase for phrase, value in paramSub.items()
                                if value is not True] != []:
                            lineChk = line
                            for phrase, value in paramSub.items():
                                if (value is not True):
                                    lineChk, paramSub[phrase] = (
                                        self.paramReplace((phrase, value),
                                                          value, line))
                                    if lineChk != line:
                                        print '#   Substituted ' + phrase
                                        line = lineChk
                                        break
                        newFile.write(line)
                    oldFile.close()
            close(fh)
            chmod(abs_path, 0o775)
            move(abs_path, 'shareIter.f90')
            move('shareIter.f90', self.dir + 'shareIter.f90')
            

        self.tBegin = 1  # -1
        self.OKFiles = ['fthb', 'income_val', 'dist_fthb',
                        'transition_fthb', 'lifecycleprofiles',
                        'lifecycle_adjust', 'transition_adjust',
                        'transition_lives', 'housingstock', 'housing_transit']
        # Create all of the model output files in advance
        access = 'w+' if new else 'r+'
        for file in (['model_log', 'policy_constantwealth',
                      'transition_lagsleads', 'transition_debug'] + self.OKFiles):
            with open('%s%s.txt' % (self.dir, file), access) as newfile:
                newfile.close()

    def batchSubCheck(self, proc):
        from time import sleep
        # In this setup, the proc submits a batch job to a node.
        # How do we tell when the batch is done? The script ends
        # by echoing an output to a file.
        out, err = proc.communicate()
        self.processing = True
        while self.processing:
            try:
                with open(self.own_path, 'r+') as outFile:
                    outLines = outFile.readlines()
                    if re.match('Done!', outLines[0]):
                        self.processing = False
                        remove(self.own_path)
            except Exception as e:
                print(str(e))
            sleep(10.0)


    def execSh(self, thread=40, **kwargs):
        """ Executes a shell script that compiles the Fortran code
            and executes it. Of the 10+ output files, those in
            self.OKFiles get moved to a directory of its own, if
            multiple models are requested.

        Args:
            thread: OpenMP threads called on this operation.
                Default is 40 because the execution is on Z, but
                could be downgraded to 8 on a laptop.
            **kwargs['model']: Optional name ascribed to the model, used
                to distinguish output files.
            **kwargs['stop_model']: If specified and is True, the model is
                never run (but files in model dir could still be moved)
        """
        from subprocess import Popen, PIPE
        from os import makedirs, chdir, remove
        from time import sleep
        from os.path import exists
        from shutil import copyfile
        from subprocess import call

        self.own_path = self.dir + 'sbatch_out.txt'
        try:
            remove(self.own_path)
        except:
            pass

        if 'stop_model' not in kwargs:
            print('################################################# \n'
                   '#   FORTRAN CODE RUNNING... \n'
                   '#################################################')
            proc = Popen(['sbatch', self.mainD + 'lifecycle_build.sh',
                          '%d' % thread], stdout=PIPE, stderr=PIPE)
            self.batchSubCheck(proc)
            copyfile(self.dir + 'model_log.txt', self.dir + 'model_log_prev.txt')
        else:
            print('Model was not executed')

        # If crunching several models in a loop
        try:
            print('#   Moving files...')
            directory = self.mainD + 'output/%s/' % kwargs['model']
            if not exists(directory):
                makedirs(directory)
            chdir(self.dir)
            for files in glob.glob('*.txt'):
                if any(fname in files for fname in self.OKFiles):
                    print(files[:-4] + '_%s' % kwargs['model'] + '.txt')
                    copyfile(files, directory + files[:-4] + '_%s'
                             % kwargs['model'] + '.txt')
            chdir(self.mainD)
            print('#   Done')
        except KeyError:
            print('#   Error with moving files into output directory')
            pass

    def readParams(self, name):
        """ Pulls out column names to distinguish the Fortran output,
            which is otherwise just a space delimited matrix.
            May be deprecated.

        Args:
            name: Output file name before any model suffixes.
        """
        # Get column headers
        self.headers = {'fthb':
                        ['id', 'age', 'nextDurables', 'nextAssets',
                         'durables', 'netWorth', 'rent', 'adjust', 'moved', 'assets',
                         'consumption', 'consumptionmpc', 'durablesmpc', 'rentalmpc',
                         'welfare', 'durTime', 'income',
                         'income_val', 'income_l', 'income_val_l', 'income_l2',
                         'income_val_l2', 'consumption_l', 'assets_l'],
                        'income_val':
                        ['age', 'income', 'income_val'],
                        'dist_fthb':
                        ['id', 'ageBought', 'nextDurables', 'nextAssets',
                         'durables', 'durTime', 'moved', 'assets', 'consumption', 'welfare',
                         'income', 'income_val', 'income_l', 'income_val_l',
                         'consumption_l', 'assets_l', 'PolTaken'],
                        'transition_fthb':
                        ['id', 'age', 'ageBought', 'nextDurables',
                         'nextAssets', 'durables', 'durTime', 'subVal',
                         'assets', 'consumption', 'welfare', 'income', 'income_val',
                         'income_l', 'income_val_l', 'consumption_l',
                         'assets_l', 'PolTaken'],
                        'transition_lives':
                        ['id', 'cohort', 'policy', 'income', 'h_ss', 'h_pol',
                         'cons_ss', 'cons_pol', 'assets_ss', 'assets_pol',
                         'rent', 'rent_pol', 'adjust_ss', 'adjust_pol'],
                        'lifecycleprofiles':
                        ['age', 'N', 'NRent', 'NOwn', 'avgConsumption',
                         'avgMPC', 'avgConsumptionOwners', 'avgHouseOwners',
                         'avgRental', 'avgAssets', 'avgAssetsOwners',
                         'avgAssetsRental', 'avgAssetsExDebtOwners',
                         'avgHouseToNetWOwners', 'aggNetW',
                         'aggIncome', 'fracOwn'],
                        'lifecycle_adjust':
                        ['id', 'age', 'adjustment'],
                        'transition_adjust':
                        ['id', 'age', 'adjustment'],
                        'transition_lagsleads':
                        ['id', 'age', 'variable', 'var_l4', 'var_l3', 'var_l2',
		         'var_l1', 'var_ss', 'var_pol0', 'var_pol1', 'var_pol2',
                         'var_pol3', 'var_pol4', 'var_pol5', 'var_pol9'],
                        'housingstock':
                        ['age', 'houseWealth', 'housInv', 'houseCount'],
                        'housing_transit':
                        ['period', 'age', 'houseWealth', 'houseInv',
                         'houseTrans', 'renters', 'alive']}

        if name not in self.OKFiles + ['transition_lagsleads']:
            raise NameError('File requested currently not an output file.')
        columns = self.headers[name]
        return columns

    def treatParams(self, outputName):
        """

        """
        view = self.appended
        if outputName in ['dist_fthb', 'transition_fthb']:
            view.loc[view['PolTaken'] == -1, 'PolTaken'] = 'First-time'
            view.loc[view['PolTaken'] == -2, 'PolTaken'] = 'Repeat'
        # Default is that the first age cohort simulated in model is age 20,
        # but can adjust this if tBegin is nonzero
        for column in ['age', 'ageBought']:
            try:
                view.loc[:, column] += self.tBegin
            except:
                pass
        return view

    @classmethod
    def readModel(cls, outputName):
        """ Read one output file of one model into memory.

        """
        temp = cls({}, sub=False)
        columns = temp.readParams(outputName)
        temp.appended = pd.read_table('%s%s.txt' % (temp.dir, outputName),
                                      sep='\s+', header=None, names=columns)
        temp.appended['model'] = 'unique'
        temp.appended = temp.treatParams(outputName)
        return temp

    @classmethod
    def appendModels(cls, outputName, **kwargs):
        """ Read one output file of multiple models into memory, appending
            them into a long data table. Should be called on the last object
            of this class created in a parameter loop.

        Args:
            outputName:
            kwargs['model']:
        """
        temp = cls({}, sub=False)
        if kwargs['model'] is StringType:
            kwargs['model'] = list(kwargs['model'])
        temp.appended = pd.DataFrame()
        columns = temp.readParams(outputName)
        try:
            for name in kwargs['model']:
                modelDir = temp.mainD + 'output/%s/' % name
                f = pd.read_table(modelDir + '%s_%s.txt' % (outputName, name),
                                  sep='\s+', header=None, names=columns)
                f['model'] = '%s' % name
                temp.appended = temp.appended.append(f)
            temp.appended = temp.treatParams(outputName)
        except:
            raise IOError('Output data could not be merged together.')
        return temp

    def addLog(self):
        """

        """
        from numpy import log
        try:
            self.appended['logAge'] = log(self.appended['ageBought'])
        except:
            raise IOError('Column ageBought not found.'
                          'Are you looking at dist_fthb.txt, where age'
                          'is interpreted as a stopping time?')

    def summaryTable(self, colName):
        """ Takes one column of a one/many model data file and produces
            some summary statistics for later storage.

        Args:
            colName: list of column names, with values taken from the
                attribute associated with outputName.
        """
        try:
            self.appended = self.appended.set_index(['model', 'id'])
        except KeyError:
            pass
        # This is where Pandas shines!
        return self.appended[colName].unstack('model').describe()

    def genMoments(self, bin, colName, wgt='N'):
        """ Reads in a panel output file and create aggregate statistics
            over 1 or several age profiles.

        Args:
            bin: integer describing how many years are in each age bin.
            colName: list of column names, with values taken from the
                attribute associated with outputName.
            wgt: column name for weighing variable between observations,
                e.g. unit count if observation is already one year.
        """
        assert type(colName) is ListType
        try:
            from numpy import sum  # , average
            self.appended['ageBin'] = 20 + bin * (self.appended['age']
                                                  // bin)
            self.weighted = self.appended.copy()
            g = self.weighted.groupby(['model', 'ageBin'])
            self.weighted['Ndenom'] = g[wgt].transform(sum)
            self.weighted[colName] = self.weighted[colName].apply(
                lambda x: x*(self.weighted[wgt] / self.weighted['Ndenom']))
            return g[colName].sum()
        except:
            raise IOError('Are you looking at one of the panel files outputted'
                          'by the model, such as lifecycleprofiles.txt?')

    def percentiles(self, dataset, colName, newName):
        """ Generates a single-column dataset that contains a
            normalized running sum of a variable from the model
            simulation. Along with an index of the number of agents
            with that variable, the dataset can generate a Lorenz
            curve.

        Args:
            dataset: The dataset to be processed. Intended in case
                a slice of the dataset needs to be considered.
            colName: A *string* that indicates the column to be
                processed into a running sum. One at a time.
            newName: The new name for the column post conversion
                into a running sum.

        """
        assert isinstance(colName, str)
        cum = dataset.sort([colName]).loc[:, [colName]]
        cum.columns = [newName]
        cum = cum.cumsum()/cum.sum()*100
        cum.loc[:, 'Pct'] = cum.rank(method='first')
        cum = cum.set_index('Pct')
        return cum

    def wealthLorenz(self):
        """ Generates a collection of normalized running sums of variables,
            which can be used to construct Lorenz curves illustrating
            wealth inequality in the simulation.

        """
        from numpy import mod
        try:
            view = self.appended
            view = view.loc[view['age'] < 40]
            view['houses'] = view['nextDurables']*(1.0-view['rent'])
            self.Lorenz = self.percentiles(view, 'houses', 'houses').join(
                self.percentiles(view, 'income_val', 'income')).join(
                self.percentiles(view, 'netWorth', 'netWorth')).join(
                self.percentiles(view.loc[view['rent'] == 0],
                                 'houses', 'houses_h')).join(
                self.percentiles(view.loc[view['rent'] == 0],
                                 'income', 'income_h')).join(
                self.percentiles(view.loc[view['rent'] == 0],
                                 'netWorth', 'netWorth_h'))
            # Clears up graph, fewer points
            self.Lorenz = self.Lorenz[mod(self.Lorenz.index, 200) == 0]
            self.Lorenz.to_csv('lorenzTest.csv')
        except:
            raise IOError('Dataset generation failed.'
                          'Are you looking at fthb.txt, the only dataset'
                          'that contains all relevant variables?')

    def matlabPlot(self, life=40, end=40, pol='FTHB', **kwargs):
        """ Assemble model outputs with a MATLAB script in the same directory
            and runs it, which builds several plots of the steady-state
            and of dynamics after the policy is enacted.

        """
        from subprocess import call
        from shutil import copyfile
        from os import remove
        try:
            for name in kwargs['model']:
                directory = self.mainD + 'output/%s/' % name
                copyfile(self.mainD + 'matlab/transition_plots.m',
                         directory + 'transition_plots.m')
                copyfile(self.dir + 'input_data/FTHB_dist_data.csv',
                         directory + 'FTHB_dist_data.csv')
                for plotFile in ['dist_fthb', 'housing_transit',
                                 'housingstock', 'transition_fthb']:
                    copyfile(directory + '%s_%s.txt' % (plotFile, name),
                             directory + '%s.txt' % plotFile)
                time = 'Tretire=%d;' % (life - self.tBegin - 1)
                polEnd = 'End=%d;' % end
                policy = 'POL=%d;' % -1 if pol == 'FTHB' else 'POL=%d;' % -2
                suffixes = "suffix='%s';" % name
                # Execute the Matlab code in the model subdirectory
                call(['matlab', '-nodisplay', '-nodesktop', '-nosplash',
                      '-r', policy + time + polEnd + suffixes +
                      'cd(\'%s\');' % directory + 'transition_plots;exit'])
                for plotFile in ['dist_fthb', 'housing_transit',
                                 'housingstock', 'transition_fthb']:
                    remove(directory + '%s.txt' % plotFile)
                for plotOut in ['FthbShock', 'FthbShockCum', 'FthbShockAge',
                                'HouseInvShockAge', 'FthbShockAxisAdj',
                                'FthbShockCumAxisAdj']:
                    copyfile(directory + '%s.pdf' % plotOut,
                             directory + '../%s_%s.pdf' % (plotOut, name))
                copyfile(directory + 'PolSeries.txt',
                         directory + '../PolSeries_%s.txt' % name)

        except IOError:
            print('Directory not found, make sure you specified the model '
                  'argument correctly.')
            pass
        except:
            pass


if __name__ == "__main__":
    base = lifecycle_iterate({'thetamatlab': 0.20})
    base.execSh(model='base')
    base_theta = lifecycle_iterate({'thetamatlab': 0.05})
    base_theta.execSh(model='base_theta')
    base.appendModels('lifecycleprofiles', model=['base', 'base_theta'])
    momH = base.genMoments(5, ['avgHouseOwners', 'avgAssetsOwners',
                           'avgAssetsExDebtOwners'], 'NOwn')
    mom = base.genMoments(5, ['avgAssets', 'avgAssetsExDebt', 'fracOwn'])
    momH.to_csv('momentsH_5yearbin.csv')
    mom.to_csv('moments_5yearbin.csv')
