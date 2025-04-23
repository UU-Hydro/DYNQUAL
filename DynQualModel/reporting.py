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

'''
Reporting (writing) of PCR-GLOBWB2 and DynQual output variables to netcdf files. Aggregates totals and averages for various time periods.
@authors (PCR-GLOBWB2): Edwin H. Sutanudjaja
@authors (DynQual)    : Edward R. Jones
'''

import os

import logging
logger = logging.getLogger(__name__)

import pcraster as pcr
import virtualOS as vos

from ncConverter import *

import variable_list as varDicts

class Reporting(object):

    #list of all output variables
    
    def __init__(self, configuration, model, modelTime, sampleNumber = None):

        self._model = model
        self._modelTime = modelTime

        self.offlineRun = configuration.routingOptions['offlineRun'] #for offline DynQual runs
        self.quality = configuration.routingOptions['quality'] #for water quality modelling
        self.calculateLoads = configuration.routingOptions['calculateLoads'] #for calculating pollutant loads in model run
        self.loadsPerSector = configuration.routingOptions['loadsPerSector'] #for calculating pollutant loads per sector

        # output directory storing netcdf files:
        self.outNCDir  = str(configuration.outNCDir)
        if sampleNumber != None:
            self.outNCDir = self.outNCDir+str(sampleNumber)+"/"
            os.mkdir(self.outNCDir)
            logger.info("Creating folder "+str(sampleNumber)+" in netcdf output directory")
        print("Output directory for netcdfs is located in: ", self.outNCDir)        

        # object for reporting:
        self.netcdfObj = PCR2netCDF(configuration)

        # initiating netcdf files for reporting
        #
        # - daily output in netCDF files:
        self.outDailyTotNC = ["None"]
        try:
            self.outDailyTotNC = configuration.reportingOptions['outDailyTotNC'].split(",")
        except:
            pass
        #
        if self.outDailyTotNC[0] != "None":
            for var in self.outDailyTotNC:
                
                logger.info("Creating the netcdf file for daily reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  
                
                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_dailyTot_output.nc",\
                                            short_name,unit,long_name)

        # -- weekly average
        self.outWeekTotNC = ["None"]
        try:
            self.outWeekTotNC = configuration.reportingOptions['outWeekTotNC'].split(",")
        except:
            pass
        if self.outWeekTotNC[0] != "None":
            for var in self.outWeekTotNC:

                # initiating monthlyVarTot (accumulator variable):
                vars(self)[var+'WeekTot'] = None

                logger.info("Creating the netcdf file for weekly accumulation reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_monthly_total_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_weekTot_output.nc",\
                                            short_name,unit,long_name)
        
        self.outWeekAvgNC = ["None"]
        try:
            self.outWeekAvgNC = configuration.reportingOptions['outWeekAvgNC'].split(",")
        except:
            pass
        if self.outWeekAvgNC[0] != "None":

            for var in self.outWeekAvgNC:

                # initiating monthlyTotAvg (accumulator variable)
                vars(self)[var+'WeekTot'] = None

                # initiating monthlyVarAvg:
                vars(self)[var+'WeekAvg'] = None

                logger.info("Creating the netcdf file for weekly average reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_weekAvg_output.nc",\
                                            short_name,unit,long_name)
        #
        # - MONTHly output in netCDF files:
        # -- cummulative
        self.outMonthTotNC = ["None"]
        try:
            self.outMonthTotNC = configuration.reportingOptions['outMonthTotNC'].split(",")
        except:
            pass
        if self.outMonthTotNC[0] != "None":
            for var in self.outMonthTotNC:

                # initiating monthlyVarTot (accumulator variable):
                vars(self)[var+'MonthTot'] = None

                logger.info("Creating the netcdf file for monthly accumulation reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_monthly_total_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_monthTot_output.nc",\
                                            short_name,unit,long_name)
        #
        # -- average
        self.outMonthAvgNC = ["None"]
        try:
            self.outMonthAvgNC = configuration.reportingOptions['outMonthAvgNC'].split(",")
        except:
            pass
        if self.outMonthAvgNC[0] != "None":

            for var in self.outMonthAvgNC:

                # initiating monthlyTotAvg (accumulator variable)
                vars(self)[var+'MonthTot'] = None

                # initiating monthlyVarAvg:
                vars(self)[var+'MonthAvg'] = None

                logger.info("Creating the netcdf file for monthly average reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_monthAvg_output.nc",\
                                            short_name,unit,long_name)
        #
        # -- last day of the month
        self.outMonthEndNC = ["None"]
        try:
            self.outMonthEndNC = configuration.reportingOptions['outMonthEndNC'].split(",")
        except:
            pass
        if self.outMonthEndNC[0] != "None":

            for var in self.outMonthEndNC:

                logger.info("Creating the netcdf file for monthly end reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_monthEnd_output.nc",\
                                            short_name,unit,long_name)
        
        # -- maximum of the month
        self.outMonthMaxNC = ["None"]
        try:
            self.outMonthMaxNC = configuration.reportingOptions['outMonthMaxNC'].split(",")
        except:
            pass
        if self.outMonthMaxNC[0] != "None":
            for var in self.outMonthMaxNC:
                logger.info("Creating the netcdf file for monthly maximum reporting for variable %s.", str(var))
                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name
                
                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_monthMax_output.nc",\
                                            short_name,unit,long_name)
        #
        # - Yearly output in netCDF files:
        # -- cumulative
        self.outAnnuaTotNC = ["None"]
        try:
            self.outAnnuaTotNC = configuration.reportingOptions['outAnnuaTotNC'].split(",")
        except:
            pass
        if self.outAnnuaTotNC[0] != "None":

            for var in self.outAnnuaTotNC:

                # initiating yearly accumulator variable:
                vars(self)[var+'AnnuaTot'] = None

                logger.info("Creating the netcdf file for annual accumulation reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_yearly_total_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_annuaTot_output.nc",\
                                            short_name,unit,long_name)
        #
        # -- average
        self.outAnnuaAvgNC = ["None"]
        try:
            self.outAnnuaAvgNC = configuration.reportingOptions['outAnnuaAvgNC'].split(",")
        except:
            pass
        if self.outAnnuaAvgNC[0] != "None":

            for var in self.outAnnuaAvgNC:

                # initiating annualyVarAvg:
                vars(self)[var+'AnnuaAvg'] = None

                # initiating annualyTotAvg (accumulator variable)
                vars(self)[var+'AnnuaTot'] = None

                logger.info("Creating the netcdf file for annual average reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_annuaAvg_output.nc",\
                                            short_name,unit,long_name)
        #
        # -- last day of the year
        self.outAnnuaEndNC = ["None"]
        try:
            self.outAnnuaEndNC = configuration.reportingOptions['outAnnuaEndNC'].split(",")
        except:
            pass
        if self.outAnnuaEndNC[0] != "None":

            for var in self.outAnnuaEndNC:

                logger.info("Creating the netcdf file for annual end reporting for variable %s.", str(var))

                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name  

                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_annuaEnd_output.nc",\
                                            short_name,unit,long_name)
        
        # -- maximum of the year
        self.outAnnuaMaxNC = ["None"]
        try:
            self.outAnnuaMaxNC = configuration.reportingOptions['outAnnuaMaxNC'].split(",")
        except:
            pass
        if self.outAnnuaMaxNC[0] != "None":
            for var in self.outAnnuaMaxNC:
                logger.info("Creating the netcdf file for annual maximum reporting for variable %s.", str(var))
                short_name = varDicts.netcdf_short_name[var]
                unit       = varDicts.netcdf_unit[var]      
                long_name  = varDicts.netcdf_long_name[var]
                if long_name == None: long_name = short_name
                
                # creating netCDF files:
                self.netcdfObj.createNetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_annuaMax_output.nc",\
                                            short_name,unit,long_name)
        
        # list of variables that will be reported:
        self.variables_for_report = self.outDailyTotNC +\
                                    self.outWeekTotNC +\
                                    self.outWeekAvgNC +\
                                    self.outMonthTotNC +\
                                    self.outMonthAvgNC +\
                                    self.outMonthEndNC +\
                                    self.outMonthMaxNC +\
                                    self.outAnnuaTotNC +\
                                    self.outAnnuaAvgNC +\
                                    self.outAnnuaEndNC +\
                                    self.outAnnuaMaxNC

    def post_processing(self):

        self.basic_post_processing() 
        self.additional_post_processing() 

    def basic_post_processing(self):

        self.precipitation  = self._model.meteo.precipitation 
        self.temperature    = self._model.meteo.temperature
        self.referencePotET = self._model.meteo.referencePotET 

        if self.offlineRun == "False":
        # routing options only required when coupled with PCR-GLOBWB
            
            self.totalLandSurfacePotET = self._model.landSurface.totalPotET 
            self.totLandSurfaceActuaET = self._model.landSurface.actualET
            
            self.fractionLandSurfaceET = vos.getValDivZero(self.totLandSurfaceActuaET,\
                                                           self.totalLandSurfacePotET,\
                                                           vos.smallNumber)
            
            self.interceptStor = self._model.landSurface.interceptStor
    
            self.snowCoverSWE  = self._model.landSurface.snowCoverSWE
            self.snowFreeWater = self._model.landSurface.snowFreeWater
    
            self.topWaterLayer = self._model.landSurface.topWaterLayer
            self.storUppTotal  = self._model.landSurface.storUppTotal
            self.storLowTotal  = self._model.landSurface.storLowTotal
            
            self.interceptEvap        = self._model.landSurface.interceptEvap
            self.actSnowFreeWaterEvap = self._model.landSurface.actSnowFreeWaterEvap
            self.topWaterLayerEvap    = self._model.landSurface.openWaterEvap
            self.actBareSoilEvap      = self._model.landSurface.actBareSoilEvap
            
            self.actTranspiTotal      = self._model.landSurface.actTranspiTotal
            self.actTranspiUppTotal   = self._model.landSurface.actTranspiUppTotal
            self.actTranspiLowTotal   = self._model.landSurface.actTranspiLowTotal
                                      
            self.directRunoff         = self._model.landSurface.directRunoff
            self.interflowTotal       = self._model.landSurface.interflowTotal
            
            self.infiltration         = self._model.landSurface.infiltration
            self.gwRecharge           = self._model.landSurface.gwRecharge
            self.gwNetCapRise         = pcr.ifthen(self._model.landSurface.gwRecharge < 0.0, self.gwRecharge*(-1.0))
            
            self.irrGrossDemand       = self._model.landSurface.irrGrossDemand    
            self.nonIrrGrossDemand    = self._model.landSurface.nonIrrGrossDemand
            self.totalGrossDemand     = self._model.landSurface.totalPotentialGrossDemand
            
            self.satDegUpp            = self._model.landSurface.satDegUppTotal
            self.satDegLow            = self._model.landSurface.satDegLowTotal
            
            self.storGroundwater      = self._model.groundwater.storGroundwater
            
            self.baseflow             = self._model.groundwater.baseflow
    
            self.surfaceWaterAbstraction         = self._model.landSurface.actSurfaceWaterAbstract
            self.nonFossilGroundWaterAbstraction = self._model.groundwater.nonFossilGroundwaterAbs
            self.otherWaterSourceAbstraction     = self._model.groundwater.unmetDemand
            self.totalAbstraction                = self.surfaceWaterAbstraction +\
                                                   self.nonFossilGroundWaterAbstraction +\
                                                   self.otherWaterSourceAbstraction
            
            # water body evaporation (m) - from surface water fractions only
            self.waterBodyActEvaporation = self._model.routing.waterBodyEvaporation
            self.waterBodyPotEvaporation = self._model.routing.waterBodyPotEvap
            #
            self.fractionWaterBodyEvaporation = vos.getValDivZero(self.waterBodyActEvaporation,\
                                                                  self.waterBodyPotEvaporation,\
                                                                  vos.smallNumber)
            # total evaporation (m), from land and water fractions
            self.totalEvaporation = self._model.landSurface.actualET + \
                                    self._model.routing.waterBodyEvaporation
        
        if self.offlineRun == "True":
        # for reporting out hydrology when DynQual is running in offline configuration
            self.directRunoff         = self._model.routing.directRunoff
            self.interflowTotal       = self._model.routing.interflowTotal
            self.baseflow             = self._model.routing.baseflow        
        
        # runoff (m) from land surface - not including local changes in water bodies
        self.runoff = self._model.routing.runoff
        self.frac_surfaceRunoff = self._model.routing.runoff
        
        # discharge (unit: m3/s)
        self.discharge = self._model.routing.disChanWaterBody
        self.dynamicFracWat = self._model.routing.dynamicFracWat
        
        # dynamic flooded fraction (-)
        self.dynamicFracWat = self._model.routing.dynamicFracWat
        
        # channel storage (m3)
        self.channelStorage = self._model.routing.channelStorage
        
        # channel storage with discharge constraint (for calculating in-stream concentrations)
        channelStorage_QThres = pcr.cover(0.1) #in m3/s
        self.channelStorage_Qthres = pcr.ifthenelse(pcr.cover(self.discharge,vos.MV) > channelStorage_QThres, self.channelStorage, vos.MV)
        
        # water temperature (K)
        self.waterTemp = self._model.routing.waterTemp
        
        # water height (m)
        self.waterHeight = self._model.routing.water_height
       
        # ice thickness (m)
        self.iceThickness = self._model.routing.iceThickness
        
        # Aspects related to salinity pollution
        self.TDSload = self._model.routing.TDSload #in grams
        self.routedTDS = self._model.routing.routedTDS #in grams
        self.salinity = pcr.ifthenelse(self.channelStorage_Qthres != vos.MV, self._model.routing.routedTDS / self.channelStorage_Qthres, 0.) #non-natural salinity in mg/L
        self.salinity = pcr.ifthenelse(self.salinity != vos.MV, self.salinity + self._model.routing.backgroundSalinity,self._model.routing.backgroundSalinity) # +background salinity  
        
        if self.loadsPerSector == "True":
            #-TDS
            self.Dom_TDSload = self._model.routing.Dom_TDSload
            self.Man_TDSload = self._model.routing.Man_TDSload
            self.USR_TDSload = self._model.routing.USR_TDSload
            self.Irr_TDSload = self._model.routing.Irr_TDSload
            
            self.routedDomTDS = self._model.routing.routedDomTDS
            self.routedManTDS = self._model.routing.routedManTDS
            self.routedUSRTDS = self._model.routing.routedUSRTDS
            self.routedIrrTDS = self._model.routing.routedIrrTDS
        
        # Aspects related to organic pollution
        self.BODload = self._model.routing.BODload #in grams
        self.routedBOD = self._model.routing.routedBOD #in grams
        self.organic = pcr.ifthenelse(self.channelStorage_Qthres != vos.MV, self._model.routing.routedBOD / self.channelStorage_Qthres, vos.MV) #in mg/l

        if self.loadsPerSector == "True":
            #-BOD
            self.Dom_BODload = self._model.routing.Dom_BODload
            self.Man_BODload = self._model.routing.Man_BODload
            self.USR_BODload = self._model.routing.USR_BODload
            self.intLiv_BODload = self._model.routing.intLiv_BODload
            self.extLiv_BODload = self._model.routing.extLiv_BODload
            
            self.routedDomBOD = self._model.routing.routedDomBOD
            self.routedManBOD = self._model.routing.routedManBOD
            self.routedUSRBOD = self._model.routing.routedUSRBOD
            self.routedintLivBOD = self._model.routing.routedintLivBOD
            self.routedextLivBOD = self._model.routing.routedextLivBOD

        # Aspects related to pathogen pollution
        self.FCload = self._model.routing.FCload #in million cfu
        self.routedFC = self._model.routing.routedFC #in million cfu
        self.pathogen = pcr.ifthenelse(self.channelStorage_Qthres != vos.MV, self._model.routing.routedFC * 100. / self.channelStorage_Qthres, vos.MV) # in cfu/100ml

        if self.loadsPerSector  == "True":            
            #-FC
            self.Dom_FCload = self._model.routing.Dom_FCload
            self.Man_FCload = self._model.routing.Man_FCload
            self.USR_FCload = self._model.routing.USR_FCload
            self.intLiv_FCload = self._model.routing.intLiv_FCload
            self.extLiv_FCload = self._model.routing.extLiv_FCload
            
            self.routedDomFC = self._model.routing.routedDomFC
            self.routedManFC = self._model.routing.routedManFC
            self.routedUSRFC = self._model.routing.routedUSRFC
            self.routedintLivFC = self._model.routing.routedintLivFC
            self.routedextLivFC = self._model.routing.routedextLivFC

    def additional_post_processing(self):
        # In this method/function, users can add their own post-processing.
        
        # consumption for and return flow from non irrigation water demand (unit: m/day)  
        ## TODO self.nonIrrWaterConsumption = self._model.routing.nonIrrWaterConsumption
        ## TODO self.nonIrrReturnFlow       = self._model.routing.nonIrrReturnFlow
        
        # accumulated runoff (m3/s) along the drainage network - not including local changes in water bodies
        if "accuRunoff" in self.variables_for_report:
            self.accuRunoff = pcr.catchmenttotal(self.runoff * self._model.routing.cellArea, self._model.routing.lddMap) / vos.secondsPerDay()
        
        # accumulated baseflow (m3) along the drainage network
        if "accuBaseflow" in self.variables_for_report:
            self.accuBaseflow = pcr.catchmenttotal(self.baseflow * self._model.routing.cellArea, self._model.routing.lddMap)

        # local changes in water bodies (i.e. abstraction, return flow, evaporation, bed exchange), excluding runoff
        self.local_water_body_flux = self._model.routing.local_input_to_surface_water / self._model.routing.cellArea - self.runoff
        
        # total runoff (m) from local land surface runoff and local changes in water bodies 
        self.totalRunoff = self.runoff + self.local_water_body_flux     # actually this is equal to self._model.routing.local_input_to_surface_water / self._model.routing.cellArea
        
        # accumulated total runoff (m3) along the drainage network - not including local changes in water bodies
        if "accuTotalRunoff" in self.variables_for_report:
            self.accuTotalRunoff = pcr.catchmenttotal(self.totalRunoff * self._model.routing.cellArea, self._model.routing.lddMap) / vos.secondsPerDay()

        # surfaceWaterStorage (unit: m) - negative values may be reported
        self.surfaceWaterStorage = self._model.routing.channelStorage / self._model.routing.cellArea

        # Stefanie's post processing: reporting lake and reservoir storage (unit: m3)
        self.waterBodyStorage = pcr.ifthen(self._model.routing.landmask, \
                                pcr.ifthen(\
                                pcr.scalar(self._model.routing.WaterBodies.waterBodyIds) > 0.,\
                                           self._model.routing.WaterBodies.waterBodyStorage))     # Note: This value is after lake/reservoir outflow.

        if self.offlineRun == "False":
        # routing options only required when coupled with PCR-GLOBWB
        
            # fossil groundwater storage
            self.storGroundwaterFossil = self._model.groundwater.storGroundwaterFossil
            
            # total groundwater storage: (non fossil and fossil)
            self.storGroundwaterTotal  = self._model.groundwater.storGroundwater + \
                                         self._model.groundwater.storGroundwaterFossil
            
            # total active storage thickness (m) for the entire water column - not including fossil groundwater (unmetDemand) 
            # - including: interception, snow, soil and non fossil groundwater 
            self.totalActiveStorageThickness = pcr.ifthen(\
                                               self._model.routing.landmask, \
                                               self._model.routing.channelStorage / self._model.routing.cellArea + \
                                               self._model.landSurface.totalSto + \
                                               self._model.groundwater.storGroundwater)
    
            # total water storage thickness (m) for the entire water column: 
            # - including: interception, snow, soil, non fossil groundwater and fossil groundwater (unmetDemand)
            # - this is usually used for GRACE comparison  
            self.totalWaterStorageThickness  = self.totalActiveStorageThickness + \
                                               self._model.groundwater.storGroundwaterFossil
    
            # Menno's post proccessing: fractions of water sources (allocated for) satisfying water demand in each cell
            self.fracSurfaceWaterAllocation = pcr.ifthen(self._model.routing.landmask, \
                                              vos.getValDivZero(\
                                              self._model.landSurface.allocSurfaceWaterAbstract, self.totalGrossDemand, vos.smallNumber))
            self.fracSurfaceWaterAllocation = pcr.ifthenelse(self.totalGrossDemand < vos.smallNumber, 1.0, self.fracSurfaceWaterAllocation)
            #
            self.fracNonFossilGroundwaterAllocation = pcr.ifthen(self._model.routing.landmask, \
                                                      vos.getValDivZero(\
                                                      self._model.groundwater.allocNonFossilGroundwater, self.totalGrossDemand, vos.smallNumber))
            self.fracNonFossilGroundwaterAllocation = pcr.ifthenelse(self.totalGrossDemand < vos.smallNumber, 0.0, self.fracNonFossilGroundwaterAllocation)
            #
            self.fracOtherWaterSourceAllocation = pcr.ifthen(self._model.routing.landmask, \
                                                  vos.getValDivZero(\
                                                  self._model.groundwater.unmetDemand, self.totalGrossDemand, vos.smallNumber))
            self.totalFracWaterSourceAllocation = self.fracSurfaceWaterAllocation + \
                                                  self.fracNonFossilGroundwaterAllocation + \
                                                  self.fracOtherWaterSourceAllocation 
    
            #
            # snowMelt (m/day)
            self.snowMelt = self._model.landSurface.snowMelt
    
            # soil moisture state from (approximately) the first 5 cm soil  
            if self._model.landSurface.numberOfSoilLayers == 3:
                self.storUppSurface   = self._model.landSurface.storUpp000005    # unit: m
                self.satDegUppSurface = self._model.landSurface.satDegUpp000005  # unit: percentage
            
            # reporting water balance from the land surface part (excluding surface water bodies)
            self.land_surface_water_balance = self._model.waterBalance
            
            # evaporation from irrigation areas (m/day) - values are average over the entire cell area
            if self._model.landSurface.includeIrrigation:\
               self.evaporation_from_irrigation = self._model.landSurface.landCoverObj['irrPaddy'].actualET * \
                                                  self._model.landSurface.landCoverObj['irrPaddy'].fracVegCover + \
                                                  self._model.landSurface.landCoverObj['irrNonPaddy'].actualET * \
                                                  self._model.landSurface.landCoverObj['irrNonPaddy'].fracVegCover 

    def report(self):

        self.post_processing()

        # time stamp for reporting
        timeStamp = datetime.datetime(self._modelTime.year,\
                                      self._modelTime.month,\
                                      self._modelTime.day,\
                                      0)

        # writing daily output to netcdf files
        if self.outDailyTotNC[0] != "None":
            for var in self.outDailyTotNC:
                
                short_name = varDicts.netcdf_short_name[var]
                self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                            "_dailyTot_output.nc",\
                                            short_name,\
                  pcr2numpy(self.__getattribute__(var),vos.MV),\
                                            timeStamp)

        # - weekly total
        if self.outWeekTotNC[0] != "None":
            for var in self.outWeekTotNC:

                # introduce accumulator at the beginning of simulation or
                if self._modelTime.timeStepPCR == 1 or self._modelTime.doy == 1:
                    vars(self)[var+'WeekTot'] = pcr.scalar(0.0)

                # accumulating
                valid = pcr.ifthen(pcr.defined(vars(self)[var]), vars(self)[var] != vos.MV)
                vars(self)[var+'WeekTot'] += pcr.ifthenelse(valid, vars(self)[var], pcr.scalar(0.))
                vars(self)[var+'_ndays_week'] += pcr.ifthenelse(valid, pcr.scalar(1.0), pcr.scalar(0.))

                # calculating total & reporting (53 weeks per year):
                if self._modelTime.doy % 7 == 0 or self._modelTime.endYear == True: 

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                               str(var)+\
                                               "_weekTot_output.nc",\
                                               short_name,\
                      pcr2numpy(self.__getattribute__(var+'WeekTot'),\
                       vos.MV),timeStamp)
                    vars(self)[var+'WeekTot'] = pcr.scalar(0.0)


        # - weekly average
        if self.outWeekAvgNC[0] != "None":
            for var in self.outWeekAvgNC:

                # introduce accumulator at the beginning of simulation or
                if self._modelTime.timeStepPCR == 1 or self._modelTime.doy == 1:
                    vars(self)[var+'WeekTot'] = pcr.scalar(0.0)
                    vars(self)[var+'_ndays_week'] = pcr.scalar(0.0)

                # accumulating
                valid = pcr.ifthen(pcr.defined(vars(self)[var]), vars(self)[var] != vos.MV)
                vars(self)[var+'WeekTot'] += pcr.ifthenelse(valid, vars(self)[var], pcr.scalar(0.))
                vars(self)[var+'_ndays_week'] += pcr.ifthenelse(valid, pcr.scalar(1.0), pcr.scalar(0.))

                # calculating average & reporting (53 weeks per year):
                if self._modelTime.doy % 7 == 0 or self._modelTime.endYear == True: 

                    vars(self)[var+'WeekAvg'] = pcr.ifthenelse(vars(self)[var+'_ndays_week'] > 0.0,\
                                                  vars(self)[var+'WeekTot'] / vars(self)[var+'_ndays_week'], pcr.scalar(vos.MV))

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                               str(var)+\
                                               "_weekAvg_output.nc",\
                                               short_name,\
                      pcr2numpy(self.__getattribute__(var+'WeekAvg'),\
                       vos.MV),timeStamp)
                    
                    vars(self)[var+'WeekTot'] = pcr.scalar(0.0)
                    vars(self)[var+'_ndays_week'] = pcr.scalar(0.0)

        # writing monthly output to netcdf files
        # - cummulative
        if self.outMonthTotNC[0] != "None":
            for var in self.outMonthTotNC:

                # introduce variables at the beginning of simulation or
                #     reset variables at the beginning of the month
                if self._modelTime.timeStepPCR == 1 or self._modelTime.day == 1:
                   vars(self)[var+'MonthTot'] = pcr.scalar(0.0)
                   vars(self)[var+'MonthDayCount'] = pcr.scalar(0.0)

                # accumulating
                valid = pcr.ifthen(pcr.defined(vars(self)[var]), vars(self)[var] != vos.MV)
                vars(self)[var+'MonthTot'] += pcr.ifthenelse(valid, vars(self)[var], pcr.scalar(0.))
                vars(self)[var+'MonthDayCount'] += pcr.ifthenelse(valid, pcr.scalar(1.0), pcr.scalar(0.))

                # reporting at the end of the month:
                if self._modelTime.endMonth == True: 

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                               "_monthTot_output.nc",\
                                               short_name,\
                      pcr2numpy(self.__getattribute__(var+'MonthTot'),\
                       vos.MV),timeStamp)
        #
        # - average
        if self.outMonthAvgNC[0] != "None":
            for var in self.outMonthAvgNC:

                # only if a accumulator variable has not been defined: 
                if var not in self.outMonthTotNC: 

                    # introduce accumulator at the beginning of simulation or
                    #     reset accumulator at the beginning of the month
                    if self._modelTime.timeStepPCR == 1 or self._modelTime.day == 1:
                       vars(self)[var+'MonthTot'] = pcr.scalar(0.0)
                       vars(self)[var+'MonthDayCount'] = pcr.scalar(0.0)

                    # accumulating
                    valid = pcr.ifthen(pcr.defined(vars(self)[var]), vars(self)[var] != vos.MV)
                    vars(self)[var+'MonthTot'] += pcr.ifthenelse(valid, vars(self)[var], pcr.scalar(0.))
                    vars(self)[var+'MonthDayCount'] += pcr.ifthenelse(valid, pcr.scalar(1.0), pcr.scalar(0.))

                # calculating average & reporting at the end of the month:
                if self._modelTime.endMonth == True:

                    vars(self)[var+'MonthAvg'] = pcr.ifthenelse(vars(self)[var+'MonthDayCount'] > 0.0,\
                                                  vars(self)[var+'MonthTot'] / vars(self)[var+'MonthDayCount'], pcr.scalar(vos.MV))

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                               str(var)+\
                                               "_monthAvg_output.nc",\
                                               short_name,\
                      pcr2numpy(self.__getattribute__(var+'MonthAvg'),\
                       vos.MV),timeStamp)
        #
        # - last day of the month
        if self.outMonthEndNC[0] != "None":
            for var in self.outMonthEndNC:

                # reporting at the end of the month:
                if self._modelTime.endMonth == True: 

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                               str(var)+\
                                               "_monthEnd_output.nc",\
                                               short_name,\
                      pcr2numpy(self.__getattribute__(var),\
                       vos.MV),timeStamp)
        
        # - maximum
        if self.outMonthMaxNC[0] != "None":
            for var in self.outMonthMaxNC:

                # introduce variables at the beginning of simulation or
                #     reset variables at the beginning of the month
                if self._modelTime.timeStepPCR == 1 or \
                   self._modelTime.day == 1:
                    vars(self)[var+'MonthMax'] = pcr.scalar(0.)

                # find the maximum
                valid = pcr.ifthen(pcr.defined(vars(self)[var]), vars(self)[var] != vos.MV)
                vars(self)[var+'MonthMax'] = pcr.ifthenelse(valid,pcr.max(vars(self)[var], vars(self)[var+'MonthMax']),vars(self)[var+'MonthMax'])

                # reporting at the end of the month:
                if self._modelTime.endMonth == True: 

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                               "_monthMax_output.nc",\
                                               short_name,\
                      pcr.pcr2numpy(self.__getattribute__(var+'MonthMax'),\
                       vos.MV),timeStamp)
        
        # writing yearly output to netcdf files
        # - cumulative
        if self.outAnnuaTotNC[0] != "None":
            for var in self.outAnnuaTotNC:

                # introduce variables at the beginning of simulation or
                #     reset variables at the beginning of the year
                if self._modelTime.timeStepPCR == 1 or self._modelTime.doy == 1:
                   vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)
                   vars(self)[var+'AnnualDayCount'] = pcr.scalar(0.0)

                # accumulating
                valid = pcr.ifthen(pcr.defined(vars(self)[var]), vars(self)[var] != vos.MV)
                vars(self)[var+'AnnuaTot'] += pcr.ifthenelse(valid, vars(self)[var], pcr.scalar(0.))
                vars(self)[var+'AnnualDayCount'] += pcr.ifthenelse(valid, pcr.scalar(1.0), pcr.scalar(0.))

                # reporting at the end of the year:
                if self._modelTime.endYear == True: 

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                               str(var)+\
                                               "_annuaTot_output.nc",\
                                               short_name,\
                      pcr2numpy(self.__getattribute__(var+'AnnuaTot'),\
                       vos.MV),timeStamp)

        # - average
        if self.outAnnuaAvgNC[0] != "None":
            for var in self.outAnnuaAvgNC:

                # only if a accumulator variable has not been defined: 
                if var not in self.outAnnuaTotNC: 

                    # introduce accumulator at the beginning of simulation or
                    #     reset accumulator at the beginning of the year
                    if self._modelTime.timeStepPCR == 1 or self._modelTime.doy == 1:
                       vars(self)[var+'AnnuaTot'] = pcr.scalar(0.0)
                       vars(self)[var+'AnnualDayCount'] = pcr.scalar(0.0)

                    # accumulating
                    valid = pcr.ifthen(pcr.defined(vars(self)[var]), vars(self)[var] != vos.MV)
                    vars(self)[var+'AnnuaTot'] += pcr.ifthenelse(valid, vars(self)[var], pcr.scalar(0.))
                    vars(self)[var+'AnnualDayCount'] += pcr.ifthenelse(valid, pcr.scalar(1.0), pcr.scalar(0.))

                # calculating average & reporting at the end of the year:
                if self._modelTime.endYear == True:

                    vars(self)[var+'AnnuaAvg'] = pcr.ifthenelse(vars(self)[var+'AnnualDayCount'] > 0.0,\
                                                  vars(self)[var+'AnnuaTot'] / vars(self)[var+'AnnualDayCount'], pcr.scalar(vos.MV))

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                               str(var)+\
                                               "_annuaAvg_output.nc",\
                                               short_name,\
                      pcr2numpy(self.__getattribute__(var+'AnnuaAvg'),\
                       vos.MV),timeStamp)
        #
        # -last day of the year
        if self.outAnnuaEndNC[0] != "None":
            for var in self.outAnnuaEndNC:

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                               str(var)+\
                                               "_annuaEnd_output.nc",\
                                               short_name,\
                      pcr2numpy(self.__getattribute__(var),\
                       vos.MV),timeStamp)
        #
        # - maximum
        if self.outAnnuaMaxNC[0] != "None":
            for var in self.outAnnuaMaxNC:

                # introduce variables at the beginning of simulation or
                #     reset variables at the beginning of the year
                if self._modelTime.timeStepPCR == 1 or \
                   self._modelTime.doy == 1:
                    vars(self)[var+'AnnuaMax'] = pcr.scalar(0.)

                # find the maximum
                valid = pcr.ifthen(pcr.defined(vars(self)[var]), vars(self)[var] != vos.MV)
                vars(self)[var+'AnnuaMax'] = pcr.ifthenelse(valid,pcr.max(vars(self)[var], vars(self)[var+'AnnuaMax']),vars(self)[var+'AnnuaMax'])

                # reporting at the end of the year:
                if self._modelTime.endYear == True: 

                    short_name = varDicts.netcdf_short_name[var]
                    self.netcdfObj.data2NetCDF(self.outNCDir+"/"+ \
                                            str(var)+\
                                               "_annuaMax_output.nc",\
                                               short_name,\
                      pcr.pcr2numpy(self.__getattribute__(var+'AnnuaMax'),\
                       vos.MV),timeStamp)

        logger.info("reporting for time %s", self._modelTime.currTime)