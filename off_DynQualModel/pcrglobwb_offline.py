import os
import sys
import math
import gc
import logging

import pcraster as pcr

import virtualOS as vos
import meteo
import routing


logger = logging.getLogger(__name__)

'''
Created on Jan 18, 2017

@author: Niko Wanders
'''
class PCRGlobWB(object):
    
    def __init__(self, configuration, currTimeStep, initialState = None):
        self._configuration = configuration
        self._modelTime = currTimeStep
        
        pcr.setclone(configuration.cloneMap)

        # Read the ldd map.
        self.lddMap = vos.readPCRmapClone(\
                  configuration.routingOptions['lddMap'],
                  configuration.cloneMap,configuration.tmpDir,configuration.globalOptions['inputDir'],True)
        #ensure ldd map is correct, and actually of type "ldd"
        self.lddMap = pcr.lddrepair(pcr.ldd(self.lddMap))
 
        if configuration.globalOptions['landmask'] != "None":
            self.landmask = vos.readPCRmapClone(\
            configuration.globalOptions['landmask'],
            configuration.cloneMap,configuration.tmpDir,configuration.globalOptions['inputDir'])
        else:
            self.landmask = pcr.defined(self.lddMap)
       
        # defining catchment areas
        self.catchment_class = 1.0
        
        self.createSubmodels(initialState)
         
    @property
    def configuration(self):
        return self._configuration
         
    def createSubmodels(self, initialState):

        # initializing sub modules
        self.meteo = meteo.Meteo(self._configuration,self.landmask,initialState)
        self.routing = routing.Routing(self._configuration, initialState, self.lddMap)
 
    def dumpState(self, outputDirectory):
        #write all state to disk to facilitate restarting

        if outputDirectory == None:
            return
        
        state = self.getState()
        
        routingState = state['routing']
        for variable, map in routingState.items():
            vos.writePCRmapToDir(\
             map,\
             str(variable)+"_"+
             str(self._modelTime.fulldate)+".map",\
             outputDirectory)
        
    def resume(self):
        #restore state from disk. used when restarting
        pass


    #FIXME: implement
    def setState(self, state):
        logger.info("cannot set state")

        
    def report(self):
        #report the state. which states are written when is based on the configuration

        #set total to 0 on first day of the year                             
        if self._modelTime.doy == 1 or self._modelTime.isFirstTimestep():

            # set all accumulated variables to zero
            self.precipitationAcc  = pcr.ifthen(self.landmask, pcr.scalar(0.0)) 

            self.baseflowAcc                  = pcr.ifthen(self.landmask, pcr.scalar(0.0))

            self.runoffAcc                    = pcr.ifthen(self.landmask, pcr.scalar(0.0))

        # accumulating until the last day of the year:
        self.precipitationAcc   += self.meteo.precipitation

        self.baseflowAcc         += self.routing.baseflow

        self.runoffAcc           += self.routing.runoff

        if self._modelTime.isLastDayOfMonth():
            self.dumpState(self._configuration.endStateDir)
            
            totalCellArea = vos.getMapTotal(pcr.ifthen(self.landmask,self.routing.cellArea))
            msg = 'Total area = %e km2'\
                    % (totalCellArea/1e6)
            logging.getLogger("model").info(msg)

            # reporting the endStates at the end of the Year:
            variableList = ['precipitation',
                            'baseflow',
                            'runoff']

            for var in variableList:
                volume = vos.getMapVolume(\
                            self.__getattribute__(var + 'Acc'),\
                            self.routing.cellArea)
                msg = 'Accumulated %s days 1 to %i in %i = %e km3 = %e mm'\
                    % (var,int(self._modelTime.doy),\
                           int(self._modelTime.year),volume/1e9,volume*1000/totalCellArea)
                logging.getLogger("model").info(msg)
        
    def getState(self):
        result = {}
        
        result['routing'] = self.routing.getState()
        
        return result
        
    def getPseudoState(self):
        result = {}
        
        result['routing'] = self.routing.getPseudoState()
        
        return result
    
    def getAllState(self):
        result = {}
        
        result['routing'] = self.routing.getState()
        result['routing'].update(self.routing.getPseudoState())
        
        return result
        
    
    def checkWaterBalance(self, storesAtBeginning, storesAtEnd):
		# for the entire modules: snow + interception + soil + groundwater + waterDemand
		# except: river/routing 

        precipitation   = pcr.ifthen(self.landmask,\
                                     self.meteo.precipitation)          # unit: m

        runoff           = pcr.ifthen(self.landmask,self.routing.runoff)
        
        vos.waterBalanceCheck([precipitation],\
                              [runoff],\
                              [storesAtBeginning],\
                              [storesAtEnd],\
                              'all modules (including water demand), but except river/routing',\
                               True,\
                               self._modelTime.fulldate,threshold=1e-3)

        
    def read_forcings(self):
        logger.info("reading forcings for time %s", self._modelTime)
        self.meteo.read_forcings(self._modelTime)
    
    def update(self, report_water_balance=False):
        logger.info("updating model to time %s", self._modelTime)
        
        #if (report_water_balance):
        #    storesAtBeginning = self.totalLandStores()

        self.meteo.update(self._modelTime)                                         
        self.routing.update_routing_only(self._modelTime,self.meteo)

        #if (report_water_balance):
        #    storesAtEnd = self.totalLandStores()
        #    self.checkWaterBalance(storesAtBeginning, storesAtEnd)
        
        if (report_water_balance):    
            self.report()
