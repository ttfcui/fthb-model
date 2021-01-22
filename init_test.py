#!/usr/bin/python

""" Test script for successful execution of Fortran model

    To see if the model is running properly, pull the latest commit
    to the repo to your computer. Change the directory to the folder
    containing this file. Then, just execute this script from command
    line.

    If any error appears, that shows the model files (in
    subdirectory "model") are incomplete. Let me (TC) know about this.
"""

from model.model_iterate import lifecycle_iterate
test = lifecycle_iterate({}, new=True)
test.execSh(model='test')
test.appendModels('lifecycleprofiles', model=['test'])
