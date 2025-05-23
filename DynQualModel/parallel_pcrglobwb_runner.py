#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

#
# PCR-GLOBWB2 (PCRaster Global Water Balance) Global Hydrological Model
#
# Copyright (C) 2016, Edwin H. Sutanudjaja, Rens van Beek, Niko Wanders, Yoshihide Wada, 
# Joyce H. C. Bosmans, Niels Drost, Ruud J. van der Ent, Inge E. M. de Graaf, Jannis M. Hoch, 
# Kor de Jong, Derek Karssenberg, Patricia López López, Stefanie Peßenteiner, Oliver Schmitz, 
# Menno W. Straatsma, Ekkamol Vannametee, Dominik Wisser, and Marc F. P. Bierkens
# Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands
#
# DynQual (Dynamic Quality) Global Water Quality Model v1.0
# Edward R. Jones, Michelle T.H. van Vliet, Niko Wanders, Edwin H. Sutanudjaja, Rens van Beek, and Marc F. P. Bierkens
# Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands
#
# This script is designed for running online DynQual runs in parallel (self.offline = False)
# In this configuration, hydrology is simulated with the PCR-GLOBWB2 GHM.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import shutil
import glob
import datetime

import pcraster as pcr

import configuration
import virtualOS as vos

import logging
logger = logging.getLogger(__name__)

import disclaimer

# print disclaimer
disclaimer.print_disclaimer()

# get the configuration/ini file given in the system argument
iniFileName = os.path.abspath(sys.argv[1])

# option for debuging and paralelization 
debug_option = "parallel"
if len(sys.argv) > 2: debug_option = str(sys.argv[2]) 

# object to handle configuration/ini file
generalConfiguration = configuration.Configuration(iniFileName = iniFileName, debug_mode = False, no_modification = False)

clean_previous_output = True
if generalConfiguration.globalOptions['cloneAreas'] == "part_one" or\
   generalConfiguration.globalOptions['cloneAreas'] == "part_two":
    clean_previous_output = False

# clean any files exists on the ouput directory (this can be done for global runs)
if clean_previous_output and os.path.exists(generalConfiguration.globalOptions['outputDir']): shutil.rmtree(generalConfiguration.globalOptions['outputDir'])

# make log folder and initialize logging
generalOutputFolder = generalConfiguration.globalOptions['outputDir']
logFileFolder = generalOutputFolder + "/global/log/"
if os.path.exists(logFileFolder) and clean_previous_output:  
    shutil.rmtree(logFileFolder)
    os.makedirs(logFileFolder)
if os.path.exists(logFileFolder) == False: 
    os.makedirs(logFileFolder)
generalConfiguration.initialize_logging(logFileFolder)

# copy ini file to the log folder:
timestamp = datetime.datetime.now()
logger.info('Copying ini file to the folder %s', logFileFolder)
shutil.copy(iniFileName, logFileFolder + \
                         os.path.basename(iniFileName) + '_' +  str(timestamp.isoformat()).replace(":",".") + '.ini')

# make global output maps directory (it will contain a "maps" directory that will contain merged pcraster maps)
global_maps_folder = generalConfiguration.globalOptions['outputDir']+"/global/maps/"
if os.path.exists(global_maps_folder):
    shutil.rmtree(global_maps_folder)
    os.makedirs(global_maps_folder)
else:
    os.makedirs(global_maps_folder)
    
# make global netcdf directory (it will contain a "netcdf" directory that will contain merged netcdfs after post-processing)
global_netcdf_folder = generalConfiguration.globalOptions['outputDir']+"/global/netcdf/"
if os.path.exists(global_netcdf_folder):
    shutil.rmtree(global_netcdf_folder)
    os.makedirs(global_netcdf_folder)
else:
    os.makedirs(global_netcdf_folder)

# make the backup of these python scripts to a specific backup folder and go to the backup folder
scriptDir = generalConfiguration.globalOptions['outputDir'] + "/global/scripts/"
if generalConfiguration.globalOptions['cloneAreas'] == "part_one" or\
   generalConfiguration.globalOptions['cloneAreas'] == "part_two":
    scriptDir = scriptDir + "/" + generalConfiguration.globalOptions['cloneAreas'] + "/"
logger.info('Making the backup of the python scripts used to the folder %s', scriptDir)
if os.path.exists(scriptDir): shutil.rmtree(scriptDir)
os.makedirs(scriptDir)
path_of_this_module = os.path.abspath(os.path.dirname(__file__))
# working/starting directory where all script files are used
for filename in glob.glob(os.path.join(path_of_this_module, '*.py')):
    shutil.copy(filename, scriptDir)
os.chdir(scriptDir)    


# option to include merging process:
with_merging_or_modflow = True
# - in the ini file, we can also skip merging (e.g. for runs with spin-ups): 


# Note that for parallel runs with spin-up, we cannot do any merging and/or combine them with modflow 
if float(generalConfiguration.globalOptions['maxSpinUpsInYears']) > 0:
    msg = "The run is set with some spin-ups. We can NOT combine this with merging processes and modflow calculation. "
    logger.warning(msg)
    logger.warning(msg)
    logger.warning(msg)
    logger.warning(msg)
    logger.warning(msg)
    if "with_merging" in list(generalConfiguration.globalOptions.keys()) and generalConfiguration.globalOptions["with_merging"] == "False":
        with_merging_or_modflow = False
    else:
        msg = "You set this run (with spin-ups) either with modflow or merging processes. That is not possible."
        logger.warning(msg)
        sys.exit()


# pcr-globwb clone areas (for pcr-globwb multiple runs)
clone_codes = list(set(generalConfiguration.globalOptions['cloneAreas'].split(",")))
#
# - for one global run (that should be using a fat node):
if clone_codes[0] == "Global": 
    clone_codes = ['M%02d'%i for i in range(1,54,1)]

# command line(s) for PCR-GLOBWB 
logger.info('Running transient PCR-GLOBWB with/without MODFLOW ')
i_clone = 0
cmd = ''
for clone_code in clone_codes:

   cmd += "python deterministic_runner_glue_with_parallel_and_modflow_options.py " + iniFileName  + " " + debug_option + " " + clone_code + " "
   cmd = cmd + " & "
   i_clone += 1


# Note that for runs with spin-up, we should not combine it with modflow 


# command line(s) for merging and MODFLOW processes:       
if with_merging_or_modflow:

   logger.info('Also with merging and/or MODFLOW processes ')
   
   cmd += "python deterministic_runner_for_monthly_modflow_and_merging.py " + iniFileName +" "+debug_option +" transient"

   cmd = cmd + " & "       


# don't foget to add the following line
cmd = cmd + "wait"       

print(cmd)
msg = "Call: "+str(cmd)
logger.debug(msg)


# execute PCR-GLOBWB and MODFLOW
vos.cmd_line(cmd, using_subprocess = False)      
