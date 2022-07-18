#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
# PCR-GLOBWB2 (PCRaster Global Water Balance) Global Hydrological Model
#
# Copyright (C) 2016, Edwin H. Sutanudjaja, Rens van Beek, Niko Wanders, Yoshihide Wada, 
# Joyce H. C. Bosmans, Niels Drost, Ruud J. van der Ent, Inge E. M. de Graaf, Jannis M. Hoch, 
# Kor de Jong, Derek Karssenberg, Patricia López López, Stefanie Peßenteiner, Oliver Schmitz, 
# Menno W. Straatsma, Ekkamol Vannametee, Dominik Wisser, and Marc F. P. Bierkens
# Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands
#
# DynQual (Dynamic Quality) Global Water Quality Model
# Edward R. Jones, Michelle T.H. van Vliet, Niko Wanders, Edwin H. Sutanudjaja, Rens van Beek, and Marc F. P. Bierkens
# Faculty of Geosciences, Utrecht University, Utrecht, The Netherlands
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
import datetime

import pcraster as pcr
from pcraster.framework import DynamicModel
from pcraster.framework import DynamicFramework

from configuration import Configuration
from currTimeStep import ModelTime
from reporting import Reporting

from spinUp_offline import SpinUp
from pcrglobwb_offline import PCRGlobWB

import logging
logger = logging.getLogger(__name__)

import disclaimer

class DeterministicRunner(DynamicModel):

    def __init__(self, configuration, modelTime, initialState = None, system_argument = None):
        DynamicModel.__init__(self)

        self.modelTime = modelTime        
        self.model = PCRGlobWB(configuration, modelTime, initialState)
        self.reporting = Reporting(configuration, self.model, modelTime)
        
        # the model will set paramaters based on global pre-multipliers given in the argument:
        if system_argument != None: self.adusting_parameters(configuration, system_argument)

        # option to include merging processes for pcraster maps and netcdf files:
        self.with_merging = True
        if ('with_merging' in configuration.globalOptions.keys()) and (configuration.globalOptions['with_merging'] == "False"):
            self.with_merging = False

        # make the configuration available for the other method/function
        self.configuration = configuration
    
    def adusting_parameters(self, configuration, system_argument):
        pass
        
    def initial(self): 
        pass

    def dynamic(self):

        # re-calculate current model time using current pcraster timestep value
        self.modelTime.update(self.currentTimeStep())

        # read model forcing (will pick up current model time from model time object)
        self.model.read_forcings()
        
		# adjust the reference potential ET according to the given pre-multiplier
        try:
            self.model.meteo.referencePotET = self.model.meteo.referencePotET * self.multiplier_for_refPotET
        except:
            pass
                    
        # update model (will pick up current model time from model time object)
        # - for a run coupled to MODFLOW, water balance checks are not valid due to lateral flow. 
        if self.configuration.online_coupling_between_pcrglobwb_and_modflow:
            self.model.update(report_water_balance = False)
        else:
            self.model.update(report_water_balance = True)
		
        # do any needed reporting for this time step        
        self.reporting.report()

        # at the last day of the month, stop calculation until modflow and related merging process are ready (only for a run with modflow) 
        if self.modelTime.isLastDayOfMonth() and (self.configuration.online_coupling_between_pcrglobwb_and_modflow or\
                                                  self.with_merging):
            
            # wait until modflow files are ready
            if self.configuration.online_coupling_between_pcrglobwb_and_modflow:
                modflow_is_ready = False
                self.count_check = 0
                while modflow_is_ready == False:
                    if datetime.datetime.now().second == 14 or\
                       datetime.datetime.now().second == 29 or\
                       datetime.datetime.now().second == 34 or\
                       datetime.datetime.now().second == 59:\
                       modflow_is_ready = self.check_modflow_status()
                
            # wait until merged files are ready
            merged_files_are_ready = False
            while merged_files_are_ready == False:
                self.count_check = 0
                if datetime.datetime.now().second == 14 or\
                   datetime.datetime.now().second == 29 or\
                   datetime.datetime.now().second == 34 or\
                   datetime.datetime.now().second == 59:\
                   merged_files_are_ready = self.check_merging_status()

    def check_modflow_status(self):

        status_file = str(self.configuration.main_output_directory) + "/modflow/transient/maps/modflow_files_for_" + str(self.modelTime.fulldate) + "_are_ready.txt"
        msg = 'Waiting for the file: ' + status_file
        if self.count_check == 1: logger.info(msg)
        if self.count_check < 7:
            #~ logger.debug(msg)			# INACTIVATE THIS AS THIS MAKE A HUGE DEBUG (dbg) FILE
            self.count_check += 1
        status = os.path.exists(status_file)
        if status == False: return status	
        if status: self.count_check = 0            
        return status

    def check_merging_status(self):

        status_file = str(self.configuration.main_output_directory) + "/global/maps/merged_files_for_"    + str(self.modelTime.fulldate) + "_are_ready.txt"
        msg = 'Waiting for the file: ' + status_file
        if self.count_check == 1: logger.info(msg)
        if self.count_check < 7:
            #~ logger.debug(msg)			# INACTIVATE THIS AS THIS MAKE A HUGE DEBUG (dbg) FILE
            self.count_check += 1
        status = os.path.exists(status_file)
        if status == False: return status	
        if status: self.count_check = 0            
        return status
 

def main():
    
    # get the full path of configuration/ini file given in the system argument
    iniFileName   = os.path.abspath(sys.argv[1])
    
    # debug option
    debug_mode = False
    if len(sys.argv) > 2: 
        if sys.argv[2] == "debug" or sys.argv[2] == "debug_parallel": debug_mode = True
    
    # object to handle configuration/ini file
    configuration = Configuration(iniFileName = iniFileName, \
                                  debug_mode = debug_mode, \
                                  no_modification = False)      

    # parallel option
    this_run_is_part_of_a_set_of_parallel_run = False    
    if len(sys.argv) > 2: 
        if sys.argv[2] == "parallel" or sys.argv[2] == "debug_parallel": this_run_is_part_of_a_set_of_parallel_run = True
    
    # for a non parallel run (usually 30min), a specific directory given in the system argument (sys.argv[3]) will be assigned for a given parameter combination:
    if this_run_is_part_of_a_set_of_parallel_run == False:
        # modfiying 'outputDir' (based on the given system argument)
        configuration.globalOptions['outputDir'] += "/"+str(sys.argv[3])+"/" 

    # for a parallel run (usually 5min), we assign a specific directory based on the clone number/code:
    if this_run_is_part_of_a_set_of_parallel_run:
        # modfiying outputDir, clone-map and landmask (based on the given system arguments)
        clone_code = str(sys.argv[3])
        configuration.globalOptions['outputDir'] += "/"+clone_code+"/" 
        configuration.globalOptions['cloneMap']   = configuration.globalOptions['cloneMap'] %(clone_code)
        configuration.globalOptions['landmask']   = configuration.globalOptions['landmask'] %(clone_code)

    # set configuration
    configuration.set_configuration(system_arguments = sys.argv)

    # timeStep info: year, month, day, doy, hour, etc
    currTimeStep = ModelTime() 

    # object for spin_up
    spin_up = SpinUp(configuration)            
    
    # spinning-up 
    noSpinUps = int(configuration.globalOptions['maxSpinUpsInYears'])
    initial_state = None
    if noSpinUps > 0:
        
        logger.info('Spin-Up #Total Years: '+str(noSpinUps))

        spinUpRun = 0 ; has_converged = False
        while spinUpRun < noSpinUps and has_converged == False:
            spinUpRun += 1
            currTimeStep.getStartEndTimeStepsForSpinUp(
                    configuration.globalOptions['startTime'],
                    spinUpRun, noSpinUps)
            logger.info('Spin-Up Run No. '+str(spinUpRun))
            deterministic_runner = DeterministicRunner(configuration, currTimeStep, initial_state, sys.argv)
            
            all_state_begin = deterministic_runner.model.getAllState() 
            
            dynamic_framework = DynamicFramework(deterministic_runner,currTimeStep.nrOfTimeSteps)
            dynamic_framework.setQuiet(True)
            dynamic_framework.run()
            
            all_state_end = deterministic_runner.model.getAllState() 
            
            has_converged = spin_up.checkConvergence(all_state_begin, all_state_end, spinUpRun, deterministic_runner.model.routing.cellArea)
            
            initial_state = deterministic_runner.model.getState()
    #
    # Running the deterministic_runner (excluding DA scheme)
    currTimeStep.getStartEndTimeSteps(configuration.globalOptions['startTime'],
                                      configuration.globalOptions['endTime'])
    
    logger.info('Transient simulation run started.')
    deterministic_runner = DeterministicRunner(configuration, currTimeStep, initial_state, sys.argv)
    
    dynamic_framework = DynamicFramework(deterministic_runner,currTimeStep.nrOfTimeSteps)
    dynamic_framework.setQuiet(True)
    dynamic_framework.run()

if __name__ == '__main__':
    # print disclaimer
    disclaimer.print_disclaimer(with_logger = True)
    sys.exit(main())
