#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
# This script is designed for running online DynQual runs (self.offline = False)
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

from pcraster.framework import DynamicModel
from pcraster.framework import DynamicFramework

from configuration import Configuration
from currTimeStep import ModelTime
from reporting import Reporting

from spinUp import SpinUp
from pcrglobwb import PCRGlobWB

    
####----------------------------------------####

import logging
logger = logging.getLogger(__name__)

import oldcalc_framework
import disclaimer  

class DeterministicRunner(DynamicModel):

    def __init__(self, configuration, modelTime, initialState = None):
        DynamicModel.__init__(self)

        self.modelTime = modelTime                
        self.model = PCRGlobWB(configuration, modelTime, initialState)
        self.reporting = Reporting(configuration, self.model, modelTime)       
        
    def initial(self): 
        pass

    def dynamic(self):

        # re-calculate current model time using current pcraster timestep value
        self.modelTime.update(self.currentTimeStep())

        # update model (will pick up current model time from model time object)
        self.model.read_forcings()
        self.model.update(report_water_balance=True)
        
        #do any needed reporting for this time step        
        self.reporting.report()

def main():

    # print disclaimer
    disclaimer.print_disclaimer()
    
    # get the full path of configuration/ini file given in the system argument
    iniFileName   = os.path.abspath(sys.argv[1])
    
    # debug option
    debug_mode = False
    if len(sys.argv) > 2: 
        if sys.argv[2] == "debug": debug_mode = True
    
    # no modification in the given ini file, use it as it is
    no_modification = True
    
    # use the output directory as given in the system argument
    if len(sys.argv) > 3 and sys.argv[3] == "--output_dir": 
        no_modification = False
        output_directory = sys.argv[4]

    # object to handle configuration/ini file
    configuration = Configuration(iniFileName = iniFileName, \
                                  debug_mode = debug_mode, \
                                  no_modification = no_modification)      
    if no_modification == False:
        configuration.main_output_directory = output_directory
        configuration.globalOptions['outputDir'] = output_directory
        configuration.set_configuration()

    # timeStep info: year, month, day, doy, hour, etc
    currTimeStep = ModelTime() 
    
    # object for spin_up
    spin_up = SpinUp(configuration)            

    # spinningUp
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
            deterministic_runner = DeterministicRunner(configuration, currTimeStep, initial_state)
            
            all_state_begin = deterministic_runner.model.getAllState() 
            
            dynamic_framework = DynamicFramework(deterministic_runner,currTimeStep.nrOfTimeSteps)
            dynamic_framework.setQuiet(True)
            dynamic_framework.run()
            
            all_state_end = deterministic_runner.model.getAllState() 
            
            has_converged = spin_up.checkConvergence(all_state_begin, all_state_end, spinUpRun, deterministic_runner.model.routing.cellArea)
            
            initial_state = deterministic_runner.model.getState()
    
    # Running the deterministic_runner (excluding DA scheme)
    currTimeStep.getStartEndTimeSteps(configuration.globalOptions['startTime'],
                                      configuration.globalOptions['endTime'])
    logger.info('Transient simulation run started.')
    deterministic_runner = DeterministicRunner(configuration, currTimeStep, initial_state)
    dynamic_framework = DynamicFramework(deterministic_runner,currTimeStep.nrOfTimeSteps)
    dynamic_framework.setQuiet(True)
    dynamic_framework.run()

    if configuration.routingOptions["offlineRun"] == "False":
        # for debugging to PCR-GLOBWB version one
        if configuration.debug_to_version_one:
            logger.info('\n\n\n\n\n'+'Executing PCR-GLOBWB version 1.'+'\n\n\n\n\n')
    
            # reset modelTime object
            currTimeStep = None; currTimeStep = ModelTime() 
            currTimeStep.getStartEndTimeSteps(configuration.globalOptions['startTime'],
                                              configuration.globalOptions['endTime'])
            
            # execute PCR-GLOBWB version 1
            # - including comparing model outputs (from versions one and two)
            pcrglobwb_one = oldcalc_framework.PCRGlobWBVersionOne(configuration, \
                                                                  currTimeStep, \
                                                                  deterministic_runner.model.routing.landmask, \
                                                                  deterministic_runner.model.routing.cellArea)
            dynamic_framework = DynamicFramework(pcrglobwb_one, currTimeStep.nrOfTimeSteps)
            dynamic_framework.setQuiet(True)
            dynamic_framework.run()
        
if __name__ == '__main__':
    # print disclaimer
    disclaimer.print_disclaimer(with_logger = True)
    sys.exit(main())
