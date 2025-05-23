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

import os
from pcraster.framework import *
import pcraster as pcr

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos
from ncConverter import *
import ETPFunctions as refPotET

class Meteo(object):

    def __init__(self,iniItems,landmask,spinUp):
        object.__init__(self)

        self.cloneMap = iniItems.cloneMap
        self.tmpDir = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
        
        # landmask/area of interest
        self.landmask = landmask
        if iniItems.globalOptions['landmask'] != "None":
           self.landmask = vos.readPCRmapClone(\
           iniItems.globalOptions['landmask'],
           self.cloneMap,self.tmpDir,self.inputDir)	

        self.preFileNC = iniItems.meteoOptions['precipitationNC']        # starting from 19 Feb 2014, we only support netcdf input files
        self.tmpFileNC = iniItems.meteoOptions['temperatureNC']

        self.refETPotMethod = iniItems.meteoOptions['referenceETPotMethod']
        if self.refETPotMethod == 'Hamon': self.latitudes = \
                                           pcr.ycoordinate(self.cloneMap) # needed to calculate 'referenceETPot'
        if self.refETPotMethod == 'Input': self.etpFileNC = \
                             iniItems.meteoOptions['refETPotFileNC']              

        try:
            self.precipitationCorrectionFactor = float(iniItems.meteoOptions['precipitationCorrectionFactor'])
        except:
            self.precipitationCorrectionFactor = float(1.0)
        try:
            self.precipitationVarName = iniItems.meteoOptions['precipitationVarName']
        except:
            self.precipitationVarName = "precipitation"

        try:
            self.temperatureVarName = iniItems.meteoOptions['temperatureVarName']
        except:
            self.temperatureVarName = "temperature"
        
        try:
            self.referenceEvaporationVarName = iniItems.meteoOptions['referenceETPotVarName']
        except:
            self.referenceEvaporationVarName = "evapotranspiration"


        # daily time step
        self.usingDailyTimeStepForcingData = False
        if iniItems.timeStep == 1.0 and iniItems.timeStepUnit == "day":
            self.usingDailyTimeStepForcingData = True
        
        # forcing downscaling options:
        self.forcingDownscalingOptions(iniItems)

        self.report = True
        try:
            self.outDailyTotNC = iniItems.meteoOptions['outDailyTotNC'].split(",")
            self.outMonthTotNC = iniItems.meteoOptions['outMonthTotNC'].split(",")
            self.outMonthAvgNC = iniItems.meteoOptions['outMonthAvgNC'].split(",")
            self.outMonthEndNC = iniItems.meteoOptions['outMonthEndNC'].split(",")
            self.outAnnuaTotNC = iniItems.meteoOptions['outAnnuaTotNC'].split(",")
            self.outAnnuaAvgNC = iniItems.meteoOptions['outAnnuaAvgNC'].split(",")
            self.outAnnuaEndNC = iniItems.meteoOptions['outAnnuaEndNC'].split(",")
        except:
            self.report = False
        if self.report == True:
            # daily output in netCDF files:
            self.outNCDir  = iniItems.outNCDir
            self.netcdfObj = PCR2netCDF(iniItems)
            #
            if self.outDailyTotNC[0] != "None":
                for var in self.outDailyTotNC:
                    # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_dailyTot.nc",\
                                                    var,"undefined")
            # MONTHly output in netCDF files:
            # - cummulative
            if self.outMonthTotNC[0] != "None":
                for var in self.outMonthTotNC:
                    # initiating monthlyVarTot (accumulator variable):
                    vars(self)[var+'MonthTot'] = None
                    # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_monthTot.nc",\
                                                    var,"undefined")
            # - average
            if self.outMonthAvgNC[0] != "None":
                for var in self.outMonthAvgNC:
                    # initiating monthlyTotAvg (accumulator variable)
                    vars(self)[var+'MonthTot'] = None
                    # initiating monthlyVarAvg:
                    vars(self)[var+'MonthAvg'] = None
                     # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_monthAvg.nc",\
                                                    var,"undefined")
            # - last day of the month
            if self.outMonthEndNC[0] != "None":
                for var in self.outMonthEndNC:
                     # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_monthEnd.nc",\
                                                    var,"undefined")
            # YEARly output in netCDF files:
            # - cummulative
            if self.outAnnuaTotNC[0] != "None":
                for var in self.outAnnuaTotNC:
                    # initiating yearly accumulator variable:
                    vars(self)[var+'AnnuaTot'] = None
                    # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_annuaTot.nc",\
                                                    var,"undefined")
            # - average
            if self.outAnnuaAvgNC[0] != "None":
                for var in self.outAnnuaAvgNC:
                    # initiating annualyVarAvg:
                    vars(self)[var+'AnnuaAvg'] = None
                    # initiating annualyTotAvg (accumulator variable)
                    vars(self)[var+'AnnuaTot'] = None
                     # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_annuaAvg.nc",\
                                                    var,"undefined")
            # - last day of the year
            if self.outAnnuaEndNC[0] != "None":
                for var in self.outAnnuaEndNC:
                     # creating the netCDF files:
                    self.netcdfObj.createNetCDF(str(self.outNCDir)+"/"+ \
                                                str(var)+"_annuaEnd.nc",\
                                                    var,"undefined")


    def forcingDownscalingOptions(self, iniItems):

        self.downscalePrecipitationOption  = False
        self.downscaleTemperatureOption    = False
        self.downscaleReferenceETPotOption = False

        if 'meteoDownscalingOptions' in iniItems.allSections:

            # downscaling options
            if iniItems.meteoDownscalingOptions['downscalePrecipitation']  == "True":
                self.downscalePrecipitationOption  = True  
                logger.info("Precipitation forcing will be downscaled to the cloneMap resolution.")

            if iniItems.meteoDownscalingOptions['downscaleTemperature']    == "True":
                self.downscaleTemperatureOption    = True  
                logger.info("Temperature forcing will be downscaled to the cloneMap resolution.")

            if iniItems.meteoDownscalingOptions['downscaleReferenceETPot'] == "True" and self.refETPotMethod != 'Hamon':
                self.downscaleReferenceETPotOption = True 
                logger.info("Reference potential evaporation will be downscaled to the cloneMap resolution.")

            # Note that for the Hamon method: referencePotET will be calculated based on temperature,  
            # therefore, we do not have to downscale it (particularly if temperature is already provided at high resolution). 

        if self.downscalePrecipitationOption or\
           self.downscaleTemperatureOption   or\
           self.downscaleReferenceETPotOption:

            # creating anomaly DEM
            highResolutionDEM = vos.readPCRmapClone(\
               iniItems.meteoDownscalingOptions['highResolutionDEM'],
               self.cloneMap,self.tmpDir,self.inputDir)
            highResolutionDEM = pcr.cover(highResolutionDEM, 0.0)
            highResolutionDEM = pcr.max(highResolutionDEM, 0.0)
            self.meteoDownscaleIds = vos.readPCRmapClone(\
               iniItems.meteoDownscalingOptions['meteoDownscaleIds'],
               self.cloneMap,self.tmpDir,self.inputDir,isLddMap=False,cover=None,isNomMap=True)
            self.cellArea = vos.readPCRmapClone(\
               iniItems.routingOptions['cellAreaMap'],
               self.cloneMap,self.tmpDir,self.inputDir)
            loweResolutionDEM = pcr.areatotal(pcr.cover(highResolutionDEM*self.cellArea, 0.0),\
                                              self.meteoDownscaleIds)/\
                                pcr.areatotal(pcr.cover(self.cellArea, 0.0),\
                                              self.meteoDownscaleIds)                  
            self.anomalyDEM = highResolutionDEM - loweResolutionDEM    # unit: meter  

            # temperature lapse rate (netCDF) file 
            self.temperLapseRateNC = vos.getFullPath(iniItems.meteoDownscalingOptions[\
                                        'temperLapseRateNC'],self.inputDir)                         
            self.temperatCorrelNC  = vos.getFullPath(iniItems.meteoDownscalingOptions[\
                                        'temperatCorrelNC'],self.inputDir)                    # TODO: Remove this criteria.                         

            # precipitation lapse rate (netCDF) file 
            self.precipLapseRateNC = vos.getFullPath(iniItems.meteoDownscalingOptions[\
                                        'precipLapseRateNC'],self.inputDir)
            self.precipitCorrelNC  = vos.getFullPath(iniItems.meteoDownscalingOptions[\
                                        'precipitCorrelNC'],self.inputDir)                    # TODO: Remove this criteria.                           

        else:
            logger.info("No forcing downscaling is implemented.")

        # forcing smoothing options: - THIS is still experimental. PS: MUST BE TESTED.
        self.forcingSmoothing = False
        if 'meteoDownscalingOptions' in iniItems.allSections and \
           'smoothingWindowsLength' in iniItems.meteoDownscalingOptions.keys():

            if float(iniItems.meteoDownscalingOptions['smoothingWindowsLength']) > 0.0:
                self.forcingSmoothing = True
                self.smoothingWindowsLength = vos.readPCRmapClone(\
                   iniItems.meteoDownscalingOptions['smoothingWindowsLength'],
                   self.cloneMap,self.tmpDir,self.inputDir)
                msg = "Forcing data are smoothed with 'windowaverage' using the window length:"+str(iniItems.meteoDownscalingOptions['smoothingWindowsLength'])
                logger.info(msg)   
 
    def perturb(self, name, **parameters):

        if name == "precipitation":

            # perturb the precipitation
            self.precipitation = self.precipitation * \
            pcr.min(pcr.max((1 + mapnormal() * parameters['standard_deviation']),0.01),2.0)
            #TODO: Please also make sure that precipitation >= 0
            #TODO: Add minimum and maximum 

        else:
            print("Error: only precipitation may be updated at this time")
            return -1


    def update(self,currTimeStep):
        #TODO: calculate  referencePotET
        pass

    def downscalePrecipitation(self, currTimeStep, useFactor = True, minCorrelationCriteria = 0.85):
        
        preSlope = 0.001 * vos.netcdf2PCRobjClone(\
                           self.precipLapseRateNC,self.precipitationVarName,\
                           currTimeStep.month, useDoy = "Yes",\
                           cloneMapFileName=self.cloneMap,\
                           LatitudeLongitude = True)
        preSlope = pcr.cover(preSlope, 0.0)
        preSlope = pcr.max(0.,preSlope)
        
        preCriteria = vos.netcdf2PCRobjClone(\
                     self.precipitCorrelNC,self.precipitationVarName,\
                     currTimeStep.month, useDoy = "Yes",\
                     cloneMapFileName=self.cloneMap,\
                     LatitudeLongitude = True)
        preSlope = pcr.ifthenelse(preCriteria > minCorrelationCriteria,\
                   preSlope, 0.0)             
        preSlope = pcr.cover(preSlope, 0.0)
    
        if useFactor == True:
            factor = pcr.max(0.,self.precipitation + preSlope*self.anomalyDEM)
            factor = factor / \
                     pcr.areaaverage(factor, self.meteoDownscaleIds)
            factor = pcr.cover(factor, 1.0)
            self.precipitation = factor * self.precipitation
        else:
            self.precipitation = self.precipitation + preSlope*self.anomalyDEM

        self.precipitation = pcr.max(0.0, self.precipitation)

    def downscaleTemperature(self, currTimeStep, useFactor = False, maxCorrelationCriteria = -0.75, zeroCelciusInKelvin = 273.15):
        
        tmpSlope = 1.000 * vos.netcdf2PCRobjClone(\
                           self.temperLapseRateNC,self.temperatureVarName,\
                           currTimeStep.month, useDoy = "Yes",\
                           cloneMapFileName=self.cloneMap,\
                           LatitudeLongitude = True)
        tmpSlope = pcr.min(0.,tmpSlope)  # must be negative
        tmpCriteria = vos.netcdf2PCRobjClone(\
                      self.temperatCorrelNC,self.temperatureVarName,\
                      currTimeStep.month, useDoy = "Yes",\
                      cloneMapFileName=self.cloneMap,\
                      LatitudeLongitude = True)
        tmpSlope = pcr.ifthenelse(tmpCriteria < maxCorrelationCriteria,\
                   tmpSlope, 0.0)             
        tmpSlope = pcr.cover(tmpSlope, 0.0)
    
        if useFactor == True:
            temperatureInKelvin = self.temperature + zeroCelciusInKelvin
            factor = pcr.max(0.0, temperatureInKelvin + tmpSlope * self.anomalyDEM)
            factor = factor / \
                     pcr.areaaverage(factor, self.meteoDownscaleIds)
            factor = pcr.cover(factor, 1.0)
            self.temperature = factor * temperatureInKelvin - zeroCelciusInKelvin
        else:
            self.temperature = self.temperature + tmpSlope*self.anomalyDEM

    def downscaleReferenceETPot(self, zeroCelciusInKelvin = 273.15):
        
        temperatureInKelvin = self.temperature + zeroCelciusInKelvin
        factor = pcr.max(0.0, temperatureInKelvin)
        factor = factor / \
                 pcr.areaaverage(factor, self.meteoDownscaleIds)
        factor = pcr.cover(factor, 1.0)
        self.referencePotET = pcr.max(0.0, factor * self.referencePotET)

    def read_forcings(self,currTimeStep):
        # reading precipitation:
        self.precipitation = vos.netcdf2PCRobjClone(\
                                  self.preFileNC,self.precipitationVarName,\
                                  str(currTimeStep.fulldate), 
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)

        self.precipitation = pcr.max(0.,self.precipitation*\
                self.precipitationCorrectionFactor)
        self.precipitation = pcr.cover( self.precipitation, 0.0)
        
        # ignore very small values of precipitation (less than 0.00001 m/day or less than 0.01 kg.m-2.day-1 )
        if self.usingDailyTimeStepForcingData:
            self.precipitation = pcr.rounddown(self.precipitation*100000.)/100000.

        # reading temperature
        self.temperature = vos.netcdf2PCRobjClone(\
                                 self.tmpFileNC,self.temperatureVarName,\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        self.temperature = self.temperature

        # Downscaling precipitation and temperature
        if self.downscalePrecipitationOption: self.downscalePrecipitation(currTimeStep)
        if self.downscaleTemperatureOption: self.downscaleTemperature(currTimeStep)

        # calculate or obtain referencePotET
        if self.refETPotMethod == 'Hamon': self.referencePotET = \
                                  refPotET.HamonPotET(self.temperature,\
                                                      currTimeStep.doy,\
                                                      self.latitudes)
        if self.refETPotMethod == 'Input': 
            self.referencePotET = vos.netcdf2PCRobjClone(\
                                     self.etpFileNC,self.referenceEvaporationVarName,\
                                     str(currTimeStep.fulldate), 
                                     useDoy = None,
                                      cloneMapFileName=self.cloneMap,\
                                      LatitudeLongitude = True)

        # Downscaling referenceETPot (based on temperature)
        if self.downscaleReferenceETPotOption: self.downscaleReferenceETPot()
        
        # smoothing:
        if self.forcingSmoothing == True:
            logger.info("Forcing data are smoothed.")   
            self.precipitation  = pcr.windowaverage(self.precipitation , self.smoothingWindowsLength)
            self.temperature    = pcr.windowaverage(self.temperature   , self.smoothingWindowsLength)
            self.referencePotET = pcr.windowaverage(self.referencePotET, self.smoothingWindowsLength)
        
        # define precipitation, temperature and referencePotET ONLY at landmask area (for reporting):
        self.precipitation  = pcr.ifthen(self.landmask, self.precipitation)
        self.temperature    = pcr.ifthen(self.landmask, self.temperature)
        self.referencePotET = pcr.ifthen(self.landmask, self.referencePotET)
 
        if self.report == True:
            timeStamp = datetime.datetime(currTimeStep.year,\
                                          currTimeStep.month,\
                                          currTimeStep.day,\
                                          0)
            # writing daily output to netcdf files
            timestepPCR = currTimeStep.timeStepPCR
            if self.outDailyTotNC[0] != "None":
                for var in self.outDailyTotNC:
                    self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_dailyTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,timestepPCR-1)

            # writing monthly output to netcdf files
            # -cummulative
            if self.outMonthTotNC[0] != "None":
                for var in self.outMonthTotNC:

                    # introduce variables at the beginning of simulation or
                    #     reset variables at the beginning of the month
                    if currTimeStep.timeStepPCR == 1 or \
                       currTimeStep.day == 1:\
                       vars(self)[var+'MonthTot'] = pcr.scalar(0.0)

                    # accumulating
                    vars(self)[var+'MonthTot'] += vars(self)[var]

                    # reporting at the end of the month:
                    if currTimeStep.endMonth == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'MonthTot'),\
                           vos.MV),timeStamp,currTimeStep.monthIdx-1)
            # -average
            if self.outMonthAvgNC[0] != "None":
                for var in self.outMonthAvgNC:
                    # only if a accumulator variable has not been defined: 
                    if var not in self.outMonthTotNC: 

                        # introduce accumulator at the beginning of simulation or
                        #     reset accumulator at the beginning of the month
                        if currTimeStep.timeStepPCR == 1 or \
                           currTimeStep.day == 1:\
                           vars(self)[var+'MonthTot'] = pcr.scalar(0.0)
                        # accumulating
                        vars(self)[var+'MonthTot'] += vars(self)[var]

                    # calculating average & reporting at the end of the month:
                    if currTimeStep.endMonth == True:
                        vars(self)[var+'MonthAvg'] = vars(self)[var+'MonthTot']/\
                                                     currTimeStep.day  
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthAvg.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'MonthAvg'),\
                           vos.MV),timeStamp,currTimeStep.monthIdx-1)
            #
            # -last day of the month
            if self.outMonthEndNC[0] != "None":
                for var in self.outMonthEndNC:
                    # reporting at the end of the month:
                    if currTimeStep.endMonth == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthEnd.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,currTimeStep.monthIdx-1)

            # writing yearly output to netcdf files
            # -cummulative
            if self.outAnnuaTotNC[0] != "None":
                for var in self.outAnnuaTotNC:

                    # introduce variables at the beginning of simulation or
                    #     reset variables at the beginning of the month
                    if currTimeStep.timeStepPCR == 1 or \
                       currTimeStep.doy == 1:\
                       vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)

                    # accumulating
                    vars(self)[var+'AnnuaTot'] += vars(self)[var]

                    # reporting at the end of the year:
                    if currTimeStep.endYear == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'AnnuaTot'),\
                           vos.MV),timeStamp,currTimeStep.annuaIdx-1)
            # -average
            if self.outAnnuaAvgNC[0] != "None":
                for var in self.outAnnuaAvgNC:
                    # only if a accumulator variable has not been defined: 
                    if var not in self.outAnnuaTotNC: 
                        # introduce accumulator at the beginning of simulation or
                        #     reset accumulator at the beginning of the year
                        if currTimeStep.timeStepPCR == 1 or \
                           currTimeStep.doy == 1:\
                           vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)
                        # accumulating
                        vars(self)[var+'AnnuaTot'] += vars(self)[var]
                    #
                    # calculating average & reporting at the end of the year:
                    if currTimeStep.endYear == True:
                        vars(self)[var+'AnnuaAvg'] = vars(self)[var+'AnnuaTot']/\
                                                     currTimeStep.doy  
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaAvg.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'AnnuaAvg'),\
                           vos.MV),timeStamp,currTimeStep.annuaIdx-1)
            #
            # -last day of the year
            if self.outAnnuaEndNC[0] != "None":
                for var in self.outAnnuaEndNC:
                    # reporting at the end of the year:
                    if currTimeStep.endYear == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaEnd.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,currTimeStep.annuaIdx-1)

    def read_forcings_forecasting(self,currTimeStep,sampleNumber):
        # reading precipitation:
        #try:
            #precipitationCorrectionFactor = float(iniItems.meteoOptions['precipitationCorrectionFactor'])
            #print "Nice"
        #except:
            #precipitationCorrectionFactor = float(1.0)
        #print iniItems.meteoOptions['precipitationCorrectionFactor']
        self.precipitation = vos.netcdf2PCRobjClone(\
                                  self.generate_Ensemble_Forcing_Name(self.preFileNC,sampleNumber),self.precipitationVarName,\
                                  str(currTimeStep.fulldate), 
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)

        self.precipitation = pcr.max(0.,self.precipitation*\
                self.precipitationCorrectionFactor)
        self.precipitation = pcr.cover( self.precipitation, 0.0)
        
        # ignore very small values of precipitation (less than 0.00001 m/day or less than 0.01 kg.m-2.day-1 )
        if self.usingDailyTimeStepForcingData:
            self.precipitation = pcr.rounddown(self.precipitation*100000.)/100000.

        # reading temperature
        self.temperature = vos.netcdf2PCRobjClone(\
                                 self.generate_Ensemble_Forcing_Name(self.tmpFileNC,sampleNumber),self.temperatureVarName,\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        self.temperature = self.temperature

        # Downscaling precipitation and temperature
        if self.downscalePrecipitationOption: self.downscalePrecipitation(currTimeStep)
        if self.downscaleTemperatureOption: self.downscaleTemperature(currTimeStep)

        # calculate or obtain referencePotET
        if self.refETPotMethod == 'Hamon': self.referencePotET = \
                                  refPotET.HamonPotET(self.temperature,\
                                                      currTimeStep.doy,\
                                                      self.latitudes)
        if self.refETPotMethod == 'Input': 
            self.referencePotET = vos.netcdf2PCRobjClone(\
                                     self.etpFileNC,'evapotranspiration',\
                                     str(currTimeStep.fulldate), 
                                     useDoy = None,
                                      cloneMapFileName=self.cloneMap,\
                                      LatitudeLongitude = True)

        # Downscaling referenceETPot (based on temperature)
        if self.downscaleReferenceETPotOption: self.downscaleReferenceETPot()
        
        # smoothing:
        if self.forcingSmoothing == True:
            logger.info("Forcing data are smoothed.")   
            self.precipitation  = pcr.windowaverage(self.precipitation , self.smoothingWindowsLength)
            self.temperature    = pcr.windowaverage(self.temperature   , self.smoothingWindowsLength)
            self.referencePotET = pcr.windowaverage(self.referencePotET, self.smoothingWindowsLength)
        
        # define precipitation, temperature and referencePotET ONLY at landmask area (for reporting):
        self.precipitation  = pcr.ifthen(self.landmask, self.precipitation)
        self.temperature    = pcr.ifthen(self.landmask, self.temperature)
        self.referencePotET = pcr.ifthen(self.landmask, self.referencePotET)
 
        if self.report == True:
            timeStamp = datetime.datetime(currTimeStep.year,\
                                          currTimeStep.month,\
                                          currTimeStep.day,\
                                          0)
            # writing daily output to netcdf files
            timestepPCR = currTimeStep.timeStepPCR
            if self.outDailyTotNC[0] != "None":
                for var in self.outDailyTotNC:
                    self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_dailyTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,timestepPCR-1)

            # writing monthly output to netcdf files
            # -cummulative
            if self.outMonthTotNC[0] != "None":
                for var in self.outMonthTotNC:

                    # introduce variables at the beginning of simulation or
                    #     reset variables at the beginning of the month
                    if currTimeStep.timeStepPCR == 1 or \
                       currTimeStep.day == 1:\
                       vars(self)[var+'MonthTot'] = pcr.scalar(0.0)

                    # accumulating
                    vars(self)[var+'MonthTot'] += vars(self)[var]

                    # reporting at the end of the month:
                    if currTimeStep.endMonth == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'MonthTot'),\
                           vos.MV),timeStamp,currTimeStep.monthIdx-1)
            # -average
            if self.outMonthAvgNC[0] != "None":
                for var in self.outMonthAvgNC:
                    # only if a accumulator variable has not been defined: 
                    if var not in self.outMonthTotNC: 

                        # introduce accumulator at the beginning of simulation or
                        #     reset accumulator at the beginning of the month
                        if currTimeStep.timeStepPCR == 1 or \
                           currTimeStep.day == 1:\
                           vars(self)[var+'MonthTot'] = pcr.scalar(0.0)
                        # accumulating
                        vars(self)[var+'MonthTot'] += vars(self)[var]

                    # calculating average & reporting at the end of the month:
                    if currTimeStep.endMonth == True:
                        vars(self)[var+'MonthAvg'] = vars(self)[var+'MonthTot']/\
                                                     currTimeStep.day  
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthAvg.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'MonthAvg'),\
                           vos.MV),timeStamp,currTimeStep.monthIdx-1)
            #
            # -last day of the month
            if self.outMonthEndNC[0] != "None":
                for var in self.outMonthEndNC:
                    # reporting at the end of the month:
                    if currTimeStep.endMonth == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_monthEnd.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,currTimeStep.monthIdx-1)

            # writing yearly output to netcdf files
            # -cummulative
            if self.outAnnuaTotNC[0] != "None":
                for var in self.outAnnuaTotNC:

                    # introduce variables at the beginning of simulation or
                    #     reset variables at the beginning of the month
                    if currTimeStep.timeStepPCR == 1 or \
                       currTimeStep.doy == 1:\
                       vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)

                    # accumulating
                    vars(self)[var+'AnnuaTot'] += vars(self)[var]

                    # reporting at the end of the year:
                    if currTimeStep.endYear == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaTot.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'AnnuaTot'),\
                           vos.MV),timeStamp,currTimeStep.annuaIdx-1)
            # -average
            if self.outAnnuaAvgNC[0] != "None":
                for var in self.outAnnuaAvgNC:
                    # only if a accumulator variable has not been defined: 
                    if var not in self.outAnnuaTotNC: 
                        # introduce accumulator at the beginning of simulation or
                        #     reset accumulator at the beginning of the year
                        if currTimeStep.timeStepPCR == 1 or \
                           currTimeStep.doy == 1:\
                           vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)
                        # accumulating
                        vars(self)[var+'AnnuaTot'] += vars(self)[var]
                    #
                    # calculating average & reporting at the end of the year:
                    if currTimeStep.endYear == True:
                        vars(self)[var+'AnnuaAvg'] = vars(self)[var+'AnnuaTot']/\
                                                     currTimeStep.doy  
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaAvg.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var+'AnnuaAvg'),\
                           vos.MV),timeStamp,currTimeStep.annuaIdx-1)
            #
            # -last day of the year
            if self.outAnnuaEndNC[0] != "None":
                for var in self.outAnnuaEndNC:
                    # reporting at the end of the year:
                    if currTimeStep.endYear == True: 
                        self.netcdfObj.data2NetCDF(str(self.outNCDir)+"/"+ \
                                         str(var)+"_annuaEnd.nc",\
                                         var,\
                          pcr2numpy(self.__getattribute__(var),vos.MV),\
                                         timeStamp,currTimeStep.annuaIdx-1)

    def generate_Ensemble_Forcing_Name(self,orgFile,sampleNumber):
        splitName = str.split(orgFile, "_")
        splitName[4] = "r"+str(sampleNumber)+"i1p1"
        newFile = "_".join(splitName)
        return newFile


    def matchCDF(data, orgDataCDF, refDataCDF, var="pr"):
        nx, ny, orgMax = orgDataCDF.shape
        refMax = refDataCDF.shape[2]
    
        if var == "pr":
            nonZero = orgDataCDF > 0.0
            nonZeroOrg = nonZero.argmax(2)-1
            nonZeroOrg[nonZeroOrg == 0.0] = 0.0001
            
            nonZero = refDataCDF > 0.0
            nonZeroRef = nonZero.argmax(2)-1
            
        optDif = np.zeros((nx,ny))+9e9
        lowDif = np.zeros((nx,ny))
        highDif = np.zeros((nx,ny))
        optVal = np.zeros((nx,ny))+9e10
        out = np.zeros((nx,ny))
        
        for p in range(orgMax):
            absDif = np.abs(orgDataCDF[:,:,p] - data)
            improveVal = optDif > absDif
            optVal[improveVal] = p
            optDif[improveVal] = absDif[improveVal]
            lowDif[improveVal] = np.abs(orgDataCDF[:,:,np.max([p-1, 0])] - data)[improveVal]
            highDif[improveVal] = np.abs(orgDataCDF[:,:,np.min([p+1, orgMax-1])] - data)[improveVal]
        
        upInt = highDif < lowDif
        lowInt = highDif > lowDif
        out = optVal
        out[upInt] = optVal[upInt] + optDif[upInt]/(highDif[upInt] + optDif[upInt])
        out[lowInt] = optVal[lowInt] - optDif[lowInt]/(lowDif[lowInt] + optDif[lowInt])
        out = out/(orgMax-1.)
        out[out < 0.0] = 0.0
        out[out > 100.0] = 100.0

        transData= np.zeros((nx,ny))
        out = out * (refMax-1.)
        for p in range(refMax):
            selVal = np.floor(out) == p
            maxVal = refDataCDF[:,:,np.min([p+1,refMax-1])][selVal]
            minVal = refDataCDF[:,:,np.max([p,0])][selVal]
            lowData = minVal + (out[selVal] - p) * (maxVal - minVal)
            highData = minVal + ((p+1) - out[selVal]) * (maxVal - minVal)
            transData[selVal] = (lowData + highData)/2
            if p == 0 and var == "pr":
                randomRainChance = 1-(nonZeroOrg - nonZeroRef)/nonZeroOrg
                randomRainChance[randomRainChance >= 0.99999] = 2
                randomRain = (np.random.random((nx,ny)) - randomRainChance) / (1.-randomRainChance)
                randomRain[randomRain < 0] = 0.0
                rainPercentile = (nonZeroRef-0.001) + np.ceil((nonZeroOrg - nonZeroRef) * randomRain)
                out[selVal] = rainPercentile[selVal]
        
        return(transData)
        
