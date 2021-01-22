/* 
  Title: master.do
  Author: EZ, DB, TC
  Created: 2016-12-05
  Purpose: Collected analysis tasks for FTHB model.
  
*/

*******************************************************************************
** Environment variables (work in progress, change as necessary)
*******************************************************************************
global userprofile "C:/Users/ttfcui"
global repodir "$userprofile/Documents/student-debt-fthb/code/model"
global librarydir "$userprofile/Documents/stata-library"
global dboxdir "$userprofile/Dropbox/fthb_model/data"

* These should be system invariant
global dboxdatadir "$dboxdir/scf/raw_data"
global dbox_sda_datadir "$dboxdir/scf/sda_output"
global initdir "$repodir/model/input_data"
global outdir "$repodir/output"

* Shortcut for running interactive even if code breaks
global root $localcodedir


set more 1
*----------------------------------------
clear all
set mem 100m
set logtype text
global graphconfig graphregion(c(white)) bgcolor(white)
*----------------------------------------

**************************************************************************** %<
** Calibration prep
** Cleans datasets to produce files the model needs as inputs.
** Outputs moments that are used to calibrate steady-state model,
** and eventually the model in policy transition.
*******************************************************************************

do $repodir/do/IRS_clean.do $dboxdir/IRS $initdir
do $repodir/do/SCF_calibration.do $dboxdir/SCF $initdir
/* %> */

**************************************************************************** %<
** Calibration process
** This step accesses an external Python interface for the model, written in Fortran.
** You should probably not call these scripts from this file; more about
** knowing they're there.
*******************************************************************************

* Python script for calibration of the steady-state
!python ./main_calibration.py
* Policy/value functions checking
do $repodir/do/fthb_checkvfuncs.do
* Check other calibration (non-moment) targets
do $repodir/do/steady_programs.do

**************************************************************************** %<
** Policy simulation process
*******************************************************************************

* Python script for simulating the temporary policies
* Note that this step will start producing figures through calls to Matlab.
!python ./experiment.py
/* %> */

**************************************************************************** %<
** Slides/Draft 1 programs
** Code for visualization or tabulation of data (except Matlab visualizations
** directly called from Python scripts)
*******************************************************************************

do $repodir/do/policy_programs_042019.do // Statistics on aggregate data

/* %> */
