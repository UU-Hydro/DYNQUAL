#!/usr/bin/ python
# -*- coding: utf-8 -*-

#PCR-GLOBWB2 and DynQual routing.
#@authors (PCR-GLOBWB2): Edwin H. Sutanudjaja
#@authors (DynQual)    : Edward R. Jones, Niko Wanders

import os
import types
import math

from pcraster.framework import *
import pcraster as pcr

import logging
logger = logging.getLogger(__name__)

import virtualOS as vos
from ncConverter import *

import waterBodies

class Routing(object):
    
    #TODO: remove
    def getPseudoState(self):
        result = {}
        return result

    #TODO: remove
    def getVariables(self, names):
        result = {}
        return result

    def getState(self):
        result = {}
        
        # Hydrology elements
        result['timestepsToAvgDischarge']  = self.timestepsToAvgDischarge    #  day 
        result['channelStorage']           = self.channelStorage             #  m3     ; channel storage, including lake and reservoir storage 
        result['readAvlChannelStorage']    = self.readAvlChannelStorage      #  m3     ; readily available channel storage that can be extracted to satisfy water demand
        result['avgDischargeLong']         = self.avgDischarge               #  m3/s   ;  long term average discharge
        result['m2tDischargeLong']         = self.m2tDischarge               #  (m3/s)^2
        result['avgBaseflowLong']          = self.avgBaseflow                #  m3/s   ;  long term average baseflow
        result['riverbedExchange']         = self.riverbedExchange           #  m3/day : river bed infiltration (from surface water bdoies to groundwater)
        result['waterBodyStorage']            = self.waterBodyStorage        #  m3     ; storages of lakes and reservoirs            # values given are per water body id (not per cell)
        result['avgLakeReservoirOutflowLong'] = self.avgOutflow              #  m3/s   ; long term average lake & reservoir outflow  # values given are per water body id (not per cell)
        result['avgLakeReservoirInflowShort'] = self.avgInflow               #  m3/s   ; short term average lake & reservoir inflow  # values given are per water body id (not per cell)
        result['avgDischargeShort']        = self.avgDischargeShort          #  m3/s   ; short term average discharge 
        result['subDischarge']             = self.subDischarge               #  m3/s   ; sub-time step discharge (needed for kinematic wave methods/approaches)
        
        #for irrigation return flows
        try:
            result['avg_irrGrossDemand']       = self.avg_irrGrossDemand         #  m/day  ; average irrigation gross demand    
            result['avg_netLqWaterToSoil']     = self.avg_netLqWaterToSoil       #  m/day  ; average net liquid transferred to the soil
        except:
            logger.info("Irrigation return flows not estimated")        
        
        # Water quality elements
        try:
            result['waterTemperature']        = self.waterTemp                   #  K      ; water temperature
            result['iceThickness']            = self.iceThickness                #  m      ; ice thickness
            result['routedTDS']               = self.routedTDS                   #  g TDS  ; routed TDS load (for conversion to salinity pollution in mg/L)
            result['routedBOD']               = self.routedBOD                   #  g BOD  ; routed BOD load (for conversion to organic pollution mg/L)
            result['routedFC']                = self.routedFC                    #  cfu    ; routed FC load (for conversion to pathogen pollution in cfu/100mL)
            
            try:
                #- Route pollutants individually per sector (for analysis of sectoral contributions)
                #Domestic sector
                result['routedDomTDS']               = self.routedDomTDS                   #  g TDS
                result['routedDomBOD']               = self.routedDomBOD                   #  g BOD 
                result['routedDomFC']                = self.routedDomFC                    #  cfu   
                #Manufacturing sector
                result['routedManTDS']               = self.routedManTDS                   #  g TDS
                result['routedManBOD']               = self.routedManBOD                   #  g BOD 
                result['routedManFC']                = self.routedManFC                    #  cfu
                #Urban surface runoff
                result['routedUSRTDS']               = self.routedUSRTDS                   #  g TDS
                result['routedUSRBOD']               = self.routedUSRBOD                   #  g BOD 
                result['routedUSRFC']                = self.routedUSRFC                    #  cfu
                #Intensive livestock
                result['routedintLivBOD']            = self.routedintLivBOD                #  g BOD 
                result['routedintLivFC']             = self.routedintLivFC                 #  cfu
                #Extensive livestock
                result['routedextLivBOD']            = self.routedextLivBOD                #  g BOD 
                result['routedextLivFC']             = self.routedextLivFC                 #  cfu
                #Irrigation
                result['routedIrrTDS']               = self.routedIrrTDS                   #  g TDS
            except:
                logger.info("Water quality elements per sector not simulated")
            
        except:
            logger.info("Water quality elements not simulated")
        
        return result

    def __init__(self,iniItems,initialConditions,lddMap):
        object.__init__(self)

        self.lddMap = lddMap

        self.cloneMap = iniItems.cloneMap
        self.tmpDir = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']
       
        self.elevation_path = vos.getFullPath(iniItems.landSurfaceOptions['topographyNC'],self.inputDir)
        self.elevation = vos.netcdf2PCRobjCloneWithoutTime(self.elevation_path,'dem_average', self.cloneMap, True, None, self.inputDir)  
                     
        # option to activate water balance check
        self.debugWaterBalance = True
        if iniItems.routingOptions['debugWaterBalance'] == "False":
            self.debugWaterBalance = False

        self.method = iniItems.routingOptions['routingMethod']
        
        # TODO: 26 Feb 2014, EdeltaIceThicknessn found that reasonable runs are only found 
        # if all of these options = True.                    
        self.includeWaterBodies = True
        self.includeLakes = True
        self.includeReservoirs =  True

        # local drainage direction:
        self.lddMap = vos.readPCRmapClone(\
                  iniItems.routingOptions['lddMap'],
                  self.cloneMap,self.tmpDir,self.inputDir,True)
        self.lddMap = pcr.lddrepair(pcr.ldd(self.lddMap))
        self.lddMap = pcr.lddrepair(self.lddMap)

        # landmask:
        if iniItems.globalOptions['landmask'] != "None":
           self.landmask = vos.readPCRmapClone(\
           iniItems.globalOptions['landmask'],
           self.cloneMap,self.tmpDir,self.inputDir)
        else:   	
           self.landmask = pcr.defined(self.lddMap)
        self.landmask = pcr.ifthen(pcr.defined(self.lddMap), self.landmask)
        self.landmask = pcr.cover(self.landmask, pcr.boolean(0))   

        # ldd mask 
        self.lddMap = pcr.lddmask(self.lddMap, self.landmask)

        self.cellArea = vos.readPCRmapClone(\
                  iniItems.routingOptions['cellAreaMap'],
                  self.cloneMap,self.tmpDir,self.inputDir)

        self.cellSizeInArcDeg = vos.getMapAttributes(self.cloneMap,"cellsize")  

        # maximum number of days (timesteps) to calculate long term average flow values (default: 5 years = 5 * 365 days = 1825)
        self.maxTimestepsToAvgDischargeLong  = 1825.

        # maximum number of days (timesteps) to calculate short term average values (default: 1 month = 1 * 30 days = 30)
        self.maxTimestepsToAvgDischargeShort = 30.                            

        routingParameters = ['gradient','manningsN']
        for var in routingParameters:
            input = iniItems.routingOptions[str(var)]
            vars(self)[var] = vos.readPCRmapClone(input,\
                            self.cloneMap,self.tmpDir,self.inputDir)

        # parameters needed to estimate channel dimensions/parameters   
        self.eta = 0.25
        self.nu  = 0.40
        self.tau = 8.00
        self.phi = 0.58

        # option to use minimum channel width (m)
        self.minChannelWidth = pcr.scalar(0.0)
        if "minimumChannelWidth" in iniItems.routingOptions.keys():
            if iniItems.routingOptions['minimumChannelWidth'] != "None":\
               self.minChannelWidth = vos.readPCRmapClone(\
                                      iniItems.routingOptions['minimumChannelWidth'],
                                      self.cloneMap,self.tmpDir,self.inputDir)
        
        # option to use constant channel width (m)
        self.constantChannelWidth = None
        if "constantChannelWidth" in iniItems.routingOptions.keys():
            if iniItems.routingOptions['constantChannelWidth'] != "None":\
               self.constantChannelWidth = vos.readPCRmapClone(\
                                           iniItems.routingOptions['constantChannelWidth'],
                                           self.cloneMap,self.tmpDir,self.inputDir)

        # an assumption for broad sheet flow in kinematic wave methods/approaches        
        self.beta = 0.6 

        # cellLength (m) is approximated cell diagonal   
        #
        cellSizeInArcMin    =  self.cellSizeInArcDeg*60.
        verticalSizeInMeter =  cellSizeInArcMin*1852.                            
        #
        self.cellLengthFD = ((self.cellArea/verticalSizeInMeter)**(2)+\
                                           (verticalSizeInMeter)**(2))\
                                                                **(0.5) 
        nrCellsDownstream  = pcr.ldddist(self.lddMap,\
                                         pcr.nominal(self.lddMap) == 5,1.)
        distanceDownstream = pcr.ldddist(self.lddMap,\
                                         pcr.nominal(self.lddMap) == 5,\
                                         self.cellLengthFD)
        channelLengthDownstream = \
                (self.cellLengthFD + distanceDownstream)/\
                (nrCellsDownstream + 1)                 # unit: m
        self.dist2celllength  = channelLengthDownstream /\
                                  self.cellSizeInArcDeg # unit: m/arcDegree  

        # the channel gradient must be >= minGradient 
        minGradient   = 0.00001
        self.gradient = pcr.max(minGradient,\
                        pcr.cover(self.gradient, minGradient))

        # initiate/create WaterBody class
        self.WaterBodies = waterBodies.WaterBodies(iniItems,self.landmask)

        self.fileCropKC = vos.getFullPath(\
                     iniItems.routingOptions['cropCoefficientWaterNC'],\
                     self.inputDir)

        # courantNumber criteria for numerical stability in kinematic wave methods/approaches
        self.courantNumber = 0.50

        # empirical values for minimum number of sub-time steps:
        design_flood_speed = 5.00 # m/s
        design_length_of_sub_time_step   = pcr.cellvalue(
                                           pcr.mapminimum(
                                           self.courantNumber * self.cellLengthFD / design_flood_speed),1)[0]
        self.limit_num_of_sub_time_steps = np.ceil(
                                           vos.secondsPerDay() / design_length_of_sub_time_step)
        #
        # minimum number of sub-time steps: 24 ; hourly resolution as used in Van Beek et al. (2011) 
        self.limit_num_of_sub_time_steps = max(24.0, self.limit_num_of_sub_time_steps) 
                
        # minimum number of a sub time step based on the configuration/ini file:  
        if 'maxiumLengthOfSubTimeStep' in iniItems.routingOptions.keys():
            maxiumLengthOfSubTimeStep = float(iniItems.routingOptions['maxiumLengthOfSubTimeStep'])
            minimum_number_of_sub_time_step  = np.ceil(
                                               vos.secondsPerDay() / maxiumLengthOfSubTimeStep )
            self.limit_num_of_sub_time_steps = max(\
                                               minimum_number_of_sub_time_step, \
                                               self.limit_num_of_sub_time_steps)                                 
        # 
        self.limit_num_of_sub_time_steps = np.int(self.limit_num_of_sub_time_steps)
        
        # critical water height (m) used to select stable length of sub time step in kinematic wave methods/approaches
        self.critical_water_height = 0.25;					                                                          # used in Van Beek et al. (2011)

        # assumption for the minimum fracwat value used for calculating water height
        self.min_fracwat_for_water_height = 0.01
        self.max_water_height = 50000000
        
        # assumption for minimum crop coefficient for surface water bodies 
        self.minCropWaterKC = 0.00
        if 'minCropWaterKC' in iniItems.routingOptions.keys():
            self.minCropWaterKC = float(iniItems.routingOptions['minCropWaterKC'])
            
        self.floodPlain = iniItems.routingOptions['dynamicFloodPlain'] == "True"
        
        if self.floodPlain:
            # get the elevation profile per grid cell
            self.relZFileName= iniItems.routingOptions['relativeElevationFiles']
            areaFractions= iniItems.routingOptions['relativeElevationLevels']
            #self.areaFractions = self.areaFractions.split(',')
            self.areaFractions = []
            for level in areaFractions.split(','):
                self.areaFractions.append(float(level))
            self.nrZLevels = len(self.areaFractions)
            # reduction parameter of smoothing interval and error threshold
            self.reductionKK= float(iniItems.routingOptions['reductionKK'])
            self.criterionKK= float(iniItems.routingOptions['criterionKK'])
            self.floodplainManN= float(iniItems.routingOptions['floodplainManningsN'])
            self.channelStorageCapacity = vos.readPCRmapClone(iniItems.routingOptions['maxChannelCapacity'],self.cloneMap,self.tmpDir,self.inputDir)
            try:
               self.channelLength = vos.readPCRmapClone(iniItems.routingOptions['channelLength'],self.cloneMap,self.tmpDir,self.inputDir)
            except:
               self.channelLength = self.cellLengthFD
            self.channelDepth = vos.readPCRmapClone(iniItems.routingOptions['channelDepth'],self.cloneMap,self.tmpDir,self.inputDir)
            self.channelGradient = vos.readPCRmapClone(iniItems.routingOptions['channelGradient'],self.cloneMap,self.tmpDir,self.inputDir)
            self.lddMap = vos.readPCRmapClone(iniItems.routingOptions['channelLDD'],self.cloneMap,self.tmpDir,self.inputDir, True)
            self.lddMap = pcr.lddrepair(pcr.ldd(self.lddMap))
            self.lddMap = pcr.lddrepair(self.lddMap)
            self.lddMap = pcr.lddmask(self.lddMap, self.landmask)

            self.channelStorageCapacity = ifthen(self.landmask, pcr.cover(self.channelStorageCapacity, 0.0))
            self.channelLength = ifthen(self.landmask, pcr.cover(self.channelLength, 0.0))
            self.channelDepth = ifthen(self.landmask, pcr.cover(self.channelDepth, 0.0))
            self.channelGradient = ifthen(self.landmask, pcr.cover(self.channelGradient, 0.0))
            self.channelGradient = pcr.max(self.channelGradient, 0.00005)
            pcr.report(self.channelStorageCapacity, "channelStorageCapacity.map")
            pcr.report(self.channelLength, "channelLength.map")
            pcr.report(self.channelDepth, "channelDepth.map")
            pcr.report(self.channelGradient, "channelGradient.map")

            #-patch elevations: those that are part of sills are updated on the basis of the floodplain gradient
            # using local distances deltaX per increment upto z[N] and the sum over sills
            #-fill all lists including smoothing interval and slopes
            self.relZ= [0.]*self.nrZLevels
            for iCnt in range(1,self.nrZLevels):
                inputName = self.relZFileName %(self.areaFractions[iCnt]*100)
                self.relZ[iCnt] = vos.readPCRmapClone(inputName,self.cloneMap,self.tmpDir,self.inputDir)
                pcr.report(self.relZ[iCnt], "elev"+str(iCnt)+".map")
                self.relZ[iCnt] = ifthen(self.landmask, pcr.cover(self.relZ[iCnt], 0.0))
                pcr.report(self.relZ[iCnt], "Covelev"+str(iCnt)+".map")
            #-minimum slope of floodplain, being defined as the longest sill, first used to retrieve
            # longest cumulative distance 
            deltaX= [self.cellArea**0.5]*self.nrZLevels
            deltaX[0]= 0.
            sumX= deltaX[:]
            minSlope= 0.
            for iCnt in range(self.nrZLevels):
                if iCnt < self.nrZLevels-1:
                    deltaX[iCnt]= (self.areaFractions[iCnt+1]**0.5-self.areaFractions[iCnt]**0.5)*deltaX[iCnt]
                else:
                    deltaX[iCnt]= (1.-self.areaFractions[iCnt-1]**0.5)*deltaX[iCnt]
                if iCnt > 0:
                    sumX[iCnt]= pcr.ifthenelse(self.relZ[iCnt] == self.relZ[iCnt-1],sumX[iCnt-1]+deltaX[iCnt],0.)
                    minSlope= pcr.ifthenelse(self.relZ[iCnt] == self.relZ[iCnt-1],\
                        pcr.max(sumX[iCnt],minSlope),minSlope)
            minSlope= pcr.min(self.gradient,0.5*pcr.max(deltaX[1],minSlope)**-1.)
            #-add small increment to elevations to each sill except in the case of lakes
            for iCnt in range(self.nrZLevels):
                self.relZ[iCnt]= self.relZ[iCnt]+sumX[iCnt]*pcr.ifthenelse(self.relZ[self.nrZLevels-1] > 0.,\
                    minSlope,0.)
            #-set slope and smoothing interval between dy= y(i+1)-y(i) and dx= x(i+1)-x(i)
            # on the basis of volume
            #-slope and smoothing interval
            self.kSlope=  [0.]*(self.nrZLevels)
            self.mInterval= [0.]*(self.nrZLevels)
            self.floodVolume= [0.]*(self.nrZLevels)
            for iCnt in range(1,self.nrZLevels):
                self.floodVolume[iCnt]= self.floodVolume[iCnt-1]+\
                    0.5*(self.areaFractions[iCnt]+self.areaFractions[iCnt-1])*\
                    (self.relZ[iCnt]-self.relZ[iCnt-1])*self.cellArea
                self.kSlope[iCnt-1]= (self.areaFractions[iCnt]-self.areaFractions[iCnt-1])/\
                    pcr.max(0.001,self.floodVolume[iCnt]-self.floodVolume[iCnt-1])
            for iCnt in range(1,self.nrZLevels):
                if iCnt < (self.nrZLevels-1):
                    self.mInterval[iCnt]= 0.5*self.reductionKK*pcr.min(self.floodVolume[iCnt+1]-self.floodVolume[iCnt],\
                        self.floodVolume[iCnt]-self.floodVolume[iCnt-1])
                else:
                    self.mInterval[iCnt]= 0.5*self.reductionKK*(self.floodVolume[iCnt]-self.floodVolume[iCnt-1])        
        
        try:
          self.quality = iniItems.routingOptions['quality'] == "True"
          logger.info("Water quality modelling initiated.")
        
        except:
          self.quality = False
          logger.info("Water quality modelling not initiated.")

        print("waterTemperature =",self.quality)
        print("Salinity = ", self.quality)
        print("Organic = ", self.quality)
        print("Pathogen = ", self.quality)
        
        if self.quality:
            
            #-Define constants for the energy balance 
            self.iceThresTemp= pcr.scalar(273.15) # threshold temperature for snowmelt (degK)
            self.densityWater= 1000.0 # density of water [kg/m3]
            self.latentHeatVapor= pcr.scalar(2.5e6) # latent heat of vaporization [J/kg]
            self.latentHeatFusion= pcr.scalar(3.34e5) # latent heat of fusion [J/kg]
            self.specificHeatWater= pcr.scalar(4190.0) # specific heat of water [J/kg/degC]
            self.heatTransferWater= pcr.scalar(20.0) # heat transfer coefficient for water [W/m2/degC]
            self.heatTransferIce= pcr.scalar(8.0) # heat transfer coefficient for ice [W/m2/degC]
            self.albedoWater= pcr.scalar(0.15) # albedo of water [-]
            self.albedoSnow= pcr.scalar(0.50) # albedo of snow and ice [-]         
            self.deltaTPrec= pcr.scalar(1.5) #-energy balance, proxy for temperature of groundwater store: mean annual temperature and reduction in the temperature for falling rain 

            # - sunshine fraction table
            if iniItems.meteoOptions['sunhoursTable'] != "Default":
                self.sunFracTBL = vos.getFullPath(iniItems.meteoOptions['sunhoursTable'], self.inputDir) #convert cloud cover to sunshine hours (Doornkamp & Pruitt)
            else:               
                self.sunFracTBL = vos.getFullPath("sunhoursfrac.tbl", os.path.abspath(os.path.dirname( __file__ )))
                msg = "Using the default sunhoursfrac.tbl stored on " + self.sunFracTBL
                logger.info(msg)

            self.radCon= 0.25
            self.radSlope= 0.50
            self.stefanBoltzman= 5.67e-8 # [W/m2/K]
            self.maxThresTemp = pcr.scalar(322.15) # max river temperature set to 322.15 K (or 50C)

            #- Paths to cloudcover file
            self.cloudFileNC = vos.getFullPath(iniItems.meteoOptions['cloudcoverNC'], self.inputDir)
            self.radFileNC = vos.getFullPath(iniItems.meteoOptions['radiationNC'], self.inputDir)
            self.vapFileNC = vos.getFullPath(iniItems.meteoOptions['vaporNC'], self.inputDir)
            self.annualTFileNC = vos.getFullPath(iniItems.meteoOptions['annualAvgTNC'], self.inputDir)
            self.maxIceThickness= 3.0
            self.deltaIceThickness = 0.0
            
            #- Path to powerplant data
            self.PowRFNC = vos.getFullPath(iniItems.routingOptions["powerplantNC"], self.inputDir) #Power return flows (m3 day) at 5arc-min (Lohrmann et al., 2019)             
            
            #- Background TDS concentration (mg l-1)
            self.backgroundSalinityNC = vos.getFullPath(iniItems.routingOptions['backgroundSalinity'], self.inputDir)
            self.backgroundSalinity   = vos.netcdf2PCRobjCloneWithoutTime(self.backgroundSalinityNC,"bgTDS",self.cloneMap) # mg l-1
            
            #- Path to TSS loadings (from Beusen et al., 2005)
            self.tss = vos.readPCRmapClone(iniItems.routingOptions['TSSmap'],self.cloneMap,self.tmpDir,self.inputDir)
            
            #- Options required for offline runs
            if iniItems.routingOptions['offlineRun'] == "True":
                self.offlineRun = True
                logger.info("DynQual running in offline configuration")
                
                #Baseflow, interflow and direct runoff required for offline DynQual runs.
                self.baseflowNC = vos.getFullPath(iniItems.routingOptions['baseflowNC'], self.inputDir)
                self.interflowNC = vos.getFullPath(iniItems.routingOptions['interflowNC'], self.inputDir)
                self.directRunoffNC = vos.getFullPath(iniItems.routingOptions['directRunoffNC'], self.inputDir)
                    
                if iniItems.routingOptions['calculateLoads'] == "True":
                    logger.info("WARNING: Cannot calculate pollutant loadings in offline configuration.")
                    logger.info("Switch to online configuration or prescribe loadings directly.")
            
            else:
                self.offlineRun = False
                logger.info("DynQual running online")
             
            # For calculating loadings within model runs   
            if iniItems.routingOptions['calculateLoads'] == "True" and self.offlineRun == False:
                self.calculateLoads = True
                logger.info("Loadings calculated within model run.")
                
                if iniItems.routingOptions['loadsPerSector'] == "True":
                    self.loadsPerSector = True
                    logger.info("Option to report loads per sector enabled")
                
            else:
                self.calculateLoads = False
                self.loadsPerSector = False
                logger.info("Loadings are prescribed (i.e. akin to a forcing) to the model")
            
            if self.calculateLoads:
                
                ###File pathways and constant pollutant loading input data
                
                #Domestic
                self.PopulationNC = vos.getFullPath(iniItems.routingOptions["PopulationNC"], self.inputDir) #gridded population, annual, 5 arc-min (Lange & Geiger, 2020)
                self.Dom_ExcrLoadNC = vos.getFullPath(iniItems.routingOptions["Dom_ExcrLoadNC"], self.inputDir) #average (regional) excretion rates:
                self.DomTDS_ExcrLoad = vos.netcdf2PCRobjCloneWithoutTime(self.Dom_ExcrLoadNC,"Dom_Fixed_TDSload",self.cloneMap) # g/capita/day
                self.DomBOD_ExcrLoad = vos.netcdf2PCRobjCloneWithoutTime(self.Dom_ExcrLoadNC,"Dom_Fixed_BODload",self.cloneMap) # g/capita/day
                self.DomFC_ExcrLoad = vos.netcdf2PCRobjCloneWithoutTime(self.Dom_ExcrLoadNC,"Dom_Fixed_FCload",self.cloneMap)   # cfu/capita/day
                
                #Manufacturing
                self.factorInd_ManNC = vos.getFullPath(iniItems.routingOptions["factorInd_ManNC"], self.inputDir) #factor to converting industrial rf to man wastewater (-)
                self.Man_EfflConcNC = vos.getFullPath(iniItems.routingOptions["Man_EfflConcNC"], self.inputDir) #average (regional) man effluent concentrations:
                self.ManTDS_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.Man_EfflConcNC,"Man_Fixed_TDSload",self.cloneMap) # mg/L [i.e. g/m3]
                self.ManBOD_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.Man_EfflConcNC,"Man_Fixed_BODload",self.cloneMap) # mg/L [i.e. g/m3] 
                self.ManFC_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.Man_EfflConcNC,"Man_Fixed_FCload",self.cloneMap)   # cfu/100ml                
                
                #Urban suface runoff          
                self.UrbanFractionNC = vos.getFullPath(iniItems.routingOptions["UrbanFractionNC"], self.inputDir) #ratio 0 (no urban) - 1 (all urban)
                self.USR_EfflConcNC = vos.getFullPath(iniItems.routingOptions["USR_EfflConcNC"], self.inputDir) #average (regional) USR effluent concentration:
                self.USRTDS_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.USR_EfflConcNC,"USR_Fixed_TDSload",self.cloneMap) # mg/L [i.e. g/m3]
                self.USRBOD_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.USR_EfflConcNC,"USR_Fixed_BODload",self.cloneMap) # mg/L [i.e. g/m3]
                self.USRFC_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.USR_EfflConcNC,"USR_Fixed_FCload",self.cloneMap)   # cfu/100ml 
                
                #Livestock
                self.LivPopulationNC = vos.getFullPath(iniItems.routingOptions["LivPopulationNC"], self.inputDir) #Gridded livestock populations, 2010, 5 arc-min (Gilbert et al., 2018)
                self.Liv_ExcrLoadNC = vos.getFullPath(iniItems.routingOptions["Liv_ExcrLoadNC"], self.inputDir)  #average (regional) excretion rates:
                self.Bufallo_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"bufallo_BODload",self.cloneMap) # g/stock/day
                self.Chicken_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"chicken_BODload",self.cloneMap) # g/stock/day
                self.Cow_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"cow_BODload",self.cloneMap) # g/stock/day
                self.Duck_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"duck_BODload",self.cloneMap) # g/stock/day
                self.Goat_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"goat_BODload",self.cloneMap) # g/stock/day
                self.Horse_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"horse_BODload",self.cloneMap) # g/stock/day
                self.Pig_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"pig_BODload",self.cloneMap) # g/stock/day
                self.Sheep_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"sheep_BODload",self.cloneMap) # g/stock/day
                self.Bufallo_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"bufallo_FCload",self.cloneMap) # cfu/stock/day
                self.Chicken_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"chicken_FCload",self.cloneMap) # cfu/stock/day
                self.Cow_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"cow_FCload",self.cloneMap) # cfu/stock/day
                self.Duck_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"duck_FCload",self.cloneMap) # cfu/stock/day
                self.Goat_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"goat_FCload",self.cloneMap) # cfu/stock/day
                self.Horse_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"horse_FCload",self.cloneMap) # cfu/stock/day
                self.Pig_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"pig_FCload",self.cloneMap) # cfu/stock/day
                self.Sheep_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"sheep_FCload",self.cloneMap) # cfu/stock/day    
                
                #Irrigation
                self.Irr_EfflConcNC = vos.getFullPath(iniItems.routingOptions["Irr_EfflConcNC"], self.inputDir) #average soil concentration averaged over the topsoil and subsoil 
                self.IrrTDS_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.Irr_EfflConcNC,"soil_TDS",self.cloneMap) # mg/L
                                                    
                #Wastewater pathways and removal efficiencies (treatment [tertiary, secondary, primary], collected but untreated, basic sanitation, open defecation, direct)
                self.WWtPathwaysNC = vos.getFullPath(iniItems.routingOptions["WWtPathwaysNC"], self.inputDir)
                self.TDS_Ter_RemEff = pcr.scalar(0.)
                self.TDS_Sec_RemEff = pcr.scalar(0.)
                self.TDS_Pri_RemEff = pcr.scalar(0.)
                self.TDS_WWcut_RemEff = pcr.scalar(0.)
                self.TDS_dom_WWbs_RemEff = pcr.scalar(0.)
                self.TDS_dom_WWod_RemEff = pcr.scalar(0.)
                self.TDS_man_WWdirect_RemEff = pcr.scalar(0.)
                  
                self.BOD_Ter_RemEff = pcr.scalar(0.99)
                self.BOD_Sec_RemEff = pcr.scalar(0.85)
                self.BOD_Pri_RemEff = pcr.scalar(0.25)
                self.BOD_WWcut_RemEff = pcr.scalar(0.)
                self.BOD_dom_WWbs_RemEff = pcr.scalar(0.)
                self.BOD_dom_WWod_RemEff = pcr.scalar(0.)
                self.BOD_man_WWdirect_RemEff = pcr.scalar(0.)
                  
                self.FC_Ter_RemEff = pcr.scalar(0.9999)
                self.FC_Sec_RemEff = pcr.scalar(0.9745)
                self.FC_Pri_RemEff = pcr.scalar(0.4279)
                self.FC_dom_WWbs_RemEff = pcr.scalar(0.)
                self.FC_WWcut_RemEff = pcr.scalar(0.)
                self.FC_dom_WWod_RemEff = pcr.scalar(0.)
                self.FC_man_WWdirect_RemEff = pcr.scalar(0.)
                
            else:
            #- Path to (non-natural) TDS, BOD, FC loading inputs
                self.TDSloadNC = vos.getFullPath(iniItems.routingOptions["TDSloadNC"], self.inputDir)
                self.BODloadNC = vos.getFullPath(iniItems.routingOptions["BODloadNC"], self.inputDir)
                self.FCloadNC = vos.getFullPath(iniItems.routingOptions["FCloadNC"], self.inputDir)
                               
        # get the initialConditions
        self.getICs(iniItems, initialConditions)
        
        # initiate old style reporting
        self.initiate_old_style_routing_reporting(iniItems)
        

    def getICs(self,iniItems,iniConditions = None):

        if iniConditions == None:
            logger.info("Reading initial conditions from pcraster maps listed in the .ini file")
            # read initial conditions from pcraster maps listed in the ini file (for the first time step of the model; when the model just starts)
            self.timestepsToAvgDischarge = vos.readPCRmapClone(iniItems.routingOptions['timestepsToAvgDischargeIni'] ,self.cloneMap,self.tmpDir,self.inputDir)  
            self.channelStorage          = vos.readPCRmapClone(iniItems.routingOptions['channelStorageIni']          ,self.cloneMap,self.tmpDir,self.inputDir)
            self.readAvlChannelStorage   = vos.readPCRmapClone(iniItems.routingOptions['readAvlChannelStorageIni']   ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.avgDischarge            = vos.readPCRmapClone(iniItems.routingOptions['avgDischargeLongIni']        ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.m2tDischarge            = vos.readPCRmapClone(iniItems.routingOptions['m2tDischargeLongIni']        ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.avgBaseflow             = vos.readPCRmapClone(iniItems.routingOptions['avgBaseflowLongIni']         ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.riverbedExchange        = vos.readPCRmapClone(iniItems.routingOptions['riverbedExchangeIni']        ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.avgDischargeShort       = vos.readPCRmapClone(iniItems.routingOptions['avgDischargeShortIni']       ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.subDischarge            = vos.readPCRmapClone(iniItems.routingOptions['subDischargeIni'],self.cloneMap,self.tmpDir,self.inputDir)

            # Initial conditions needed for water quality module
            if self.quality:
                self.waterTemp    = vos.readPCRmapClone(iniItems.routingOptions['waterTemperatureIni'],self.cloneMap,self.tmpDir,self.inputDir)
                self.DO = (1-0.0001148*self.elevation)*exp(-139.34411+(157570.1)/(self.waterTemp)-(66423080)/(self.waterTemp**2)+(12438000000)/(self.waterTemp**3)-(862194900000)/(self.waterTemp**4))
                self.iceThickness   = vos.readPCRmapClone(iniItems.routingOptions['iceThicknessIni'],self.cloneMap,self.tmpDir,self.inputDir)
                self.routedTDS = vos.readPCRmapClone(iniItems.routingOptions['routedTDSIni'],self.cloneMap,self.tmpDir,self.inputDir) #initial conditions for salinity pollution
                self.routedBOD = vos.readPCRmapClone(iniItems.routingOptions['routedBODIni'],self.cloneMap,self.tmpDir,self.inputDir) #initial conditions for organic pollution
                self.routedFC = vos.readPCRmapClone(iniItems.routingOptions['routedFCIni'],self.cloneMap,self.tmpDir,self.inputDir) #initial conditions for pathogen pollution
                
                # Initial conditions for calculating average irrigation demand and net liquid transferred to the soil for irrigation return flow calculations
                if self.calculateLoads and self.offlineRun == False:
                    self.avg_irrGrossDemand   = vos.readPCRmapClone(iniItems.routingOptions['avg_irrGrossDemandIni'],self.cloneMap,self.tmpDir,self.inputDir)      
                    self.avg_netLqWaterToSoil = vos.readPCRmapClone(iniItems.routingOptions['avg_netLqWaterToSoilIni'], self.cloneMap,self.tmpDir,self.inputDir)
                    
                    # Per water quality sector
                    if self.loadsPerSector:
                        self.routedDomTDS = vos.readPCRmapClone(iniItems.routingOptions['routedDomTDSIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedDomBOD = vos.readPCRmapClone(iniItems.routingOptions['routedDomBODIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedDomFC = vos.readPCRmapClone(iniItems.routingOptions['routedDomFCIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedManTDS = vos.readPCRmapClone(iniItems.routingOptions['routedManTDSIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedManBOD = vos.readPCRmapClone(iniItems.routingOptions['routedManBODIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedManFC = vos.readPCRmapClone(iniItems.routingOptions['routedManFCIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedUSRTDS = vos.readPCRmapClone(iniItems.routingOptions['routedUSRTDSIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedUSRBOD = vos.readPCRmapClone(iniItems.routingOptions['routedUSRBODIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedUSRFC = vos.readPCRmapClone(iniItems.routingOptions['routedUSRFCIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedintLivBOD = vos.readPCRmapClone(iniItems.routingOptions['routedintLivBODIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedintLivFC = vos.readPCRmapClone(iniItems.routingOptions['routedintLivFCIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedextLivBOD = vos.readPCRmapClone(iniItems.routingOptions['routedextLivBODIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedextLivFC = vos.readPCRmapClone(iniItems.routingOptions['routedextLivFCIni'],self.cloneMap,self.tmpDir,self.inputDir)
                        self.routedIrrTDS = vos.readPCRmapClone(iniItems.routingOptions['routedIrrTDSIni'],self.cloneMap,self.tmpDir,self.inputDir)

        else:              
            logger.info("Reading initial conditions from memory.")
            # read initial conditions from the memory
            self.timestepsToAvgDischarge = iniConditions['routing']['timestepsToAvgDischarge']                                                      
            self.channelStorage          = iniConditions['routing']['channelStorage']
            self.readAvlChannelStorage   = iniConditions['routing']['readAvlChannelStorage']
            self.avgDischarge            = iniConditions['routing']['avgDischargeLong']
            self.m2tDischarge            = iniConditions['routing']['m2tDischargeLong']
            self.avgBaseflow             = iniConditions['routing']['avgBaseflowLong']
            self.riverbedExchange        = iniConditions['routing']['riverbedExchange']
            self.avgDischargeShort       = iniConditions['routing']['avgDischargeShort']
            self.subDischarge            = iniConditions['routing']['subDischarge']

            # Initial conditions needed for water quality module
            if self.quality:
                self.waterTemp               = iniConditions['routing']['waterTemperature']
                self.iceThickness            = iniConditions['routing']['iceThickness']
                self.routedTDS               = iniConditions['routing']['routedTDS']
                self.routedBOD               = iniConditions['routing']['routedBOD']
                self.routedFC                = iniConditions['routing']['routedFC']
                
                # Initial conditions for calculating average irrigation demand and net liquid transferred to the soil for irrigation return flow calculations                
                if self.calculateLoads and self.offlineRun == False:    
                    self.avg_irrGrossDemand   = iniConditions['routing']['avg_irrGrossDemand']      
                    self.avg_netLqWaterToSoil = iniConditions['routing']['avg_netLqWaterToSoil']
                    
                    #Per water quality sector
                    if self.loadsPerSector:
                        self.routedDomTDS = iniConditions['routing']['routedDomTDS']
                        self.routedDomBOD = iniConditions['routing']['routedDomBOD']
                        self.routedDomFC = iniConditions['routing']['routedDomFC']
                        self.routedManTDS = iniConditions['routing']['routedManTDS']
                        self.routedManBOD = iniConditions['routing']['routedManBOD'] 
                        self.routedManFC = iniConditions['routing']['routedManFC']
                        self.routedUSRTDS = iniConditions['routing']['routedUSRTDS']
                        self.routedUSRBOD = iniConditions['routing']['routedUSRBOD']
                        self.routedUSRFC =  iniConditions['routing']['routedUSRFC']
                        self.routedintLivBOD = iniConditions['routing']['routedintLivBOD']
                        self.routedintLivFC = iniConditions['routing']['routedintLivFC']
                        self.routedextLivBOD = iniConditions['routing']['routedextLivBOD']
                        self.routedextLivFC = iniConditions['routing']['routedextLivFC']
                        self.routedIrrTDS = iniConditions['routing']['routedIrrTDS']


        # Get values associated with landmask
        self.channelStorage        = pcr.ifthen(self.landmask, pcr.cover(self.channelStorage, 0.0))
        self.readAvlChannelStorage = pcr.ifthen(self.landmask, pcr.cover(self.readAvlChannelStorage, 0.0))
        self.avgDischarge          = pcr.ifthen(self.landmask, pcr.cover(self.avgDischarge, 0.0))
        self.m2tDischarge          = pcr.ifthen(self.landmask, pcr.cover(self.m2tDischarge, 0.0))
        self.avgDischargeShort     = pcr.ifthen(self.landmask, pcr.cover(self.avgDischargeShort, 0.0))
        self.avgBaseflow           = pcr.ifthen(self.landmask, pcr.cover(self.avgBaseflow, 0.0))
        self.riverbedExchange      = pcr.ifthen(self.landmask, pcr.cover(self.riverbedExchange, 0.0))
        self.subDischarge          = pcr.ifthen(self.landmask, pcr.cover(self.subDischarge , 0.0))        
        self.readAvlChannelStorage = pcr.max(pcr.min(self.readAvlChannelStorage, self.channelStorage),0.0)

        if self.floodPlain:
            self.subDischarge          = pcr.ifthen(self.landmask, pcr.cover(self.subDischarge , 0.0))

        if self.quality:
            self.waterTemp             = pcr.ifthen(self.landmask, pcr.cover(self.waterTemp, 0.0))
            self.iceThickness          = pcr.ifthen(self.landmask, pcr.cover(self.iceThickness , 0.0))
            self.channelStorageTimeBefore = self.channelStorage
            self.totEW = self.channelStorage * self.waterTemp*self.specificHeatWater * self.densityWater
            self.temp_water_height = yMean = self.eta * pow (self.avgDischarge, self.nu)
            self.routedTDS = pcr.ifthen(self.landmask, pcr.cover(self.routedTDS, 0.0))
            self.routedBOD = pcr.ifthen(self.landmask, pcr.cover(self.routedBOD, 0.0))
            self.routedFC = pcr.ifthen(self.landmask, pcr.cover(self.routedFC,  0.0))
            
            #Per water quality sector
            if self.calculateLoads and self.offlineRun == False: 
                self.avg_irrGrossDemand   = pcr.ifthen(self.landmask, pcr.cover(self.avg_irrGrossDemand, 0.0))
                self.avg_netLqWaterToSoil = pcr.ifthen(self.landmask, pcr.cover(self.avg_netLqWaterToSoil, 0.0))
                
                if self.loadsPerSector:          
                    self.routedDomTDS = pcr.ifthen(self.landmask, pcr.cover(self.routedDomTDS, 0.0))
                    self.routedDomBOD = pcr.ifthen(self.landmask, pcr.cover(self.routedDomBOD, 0.0))
                    self.routedDomFC = pcr.ifthen(self.landmask, pcr.cover(self.routedDomFC, 0.0))
                    self.routedManTDS = pcr.ifthen(self.landmask, pcr.cover(self.routedManTDS, 0.0))
                    self.routedManBOD = pcr.ifthen(self.landmask, pcr.cover(self.routedManBOD, 0.0))
                    self.routedManFC = pcr.ifthen(self.landmask, pcr.cover(self.routedManFC, 0.0))
                    self.routedUSRTDS = pcr.ifthen(self.landmask, pcr.cover(self.routedUSRTDS, 0.0))
                    self.routedUSRBOD = pcr.ifthen(self.landmask, pcr.cover(self.routedUSRBOD, 0.0))
                    self.routedUSRFC = pcr.ifthen(self.landmask, pcr.cover(self.routedUSRFC, 0.0))
                    self.routedintLivBOD = pcr.ifthen(self.landmask, pcr.cover(self.routedintLivBOD, 0.0))
                    self.routedintLivFC = pcr.ifthen(self.landmask, pcr.cover(self.routedintLivFC, 0.0))
                    self.routedextLivBOD = pcr.ifthen(self.landmask, pcr.cover(self.routedextLivBOD, 0.0))
                    self.routedextLivFC = pcr.ifthen(self.landmask, pcr.cover(self.routedextLivFC, 0.0))
                    self.routedIrrTDS = pcr.ifthen(self.landmask, pcr.cover(self.routedIrrTDS, 0.0))
        
        # make sure that timestepsToAvgDischarge is consistent (or the same) for the entire map:
        try:
            self.timestepsToAvgDischarge = pcr.mapmaximum(self.timestepsToAvgDischarge)
        except:    
            pass # We have to use 'try/except' because 'pcr.mapmaximum' cannot handle scalar value

        # for netcdf reporting, we have to make sure that timestepsToAvgDischarge is spatial and scalar (especially while performing pcr2numpy operations)
        self.timestepsToAvgDischarge = pcr.spatial(pcr.scalar(self.timestepsToAvgDischarge))
        self.timestepsToAvgDischarge = pcr.ifthen(self.landmask, self.timestepsToAvgDischarge)

        # Initial conditions needed for water bodies:
        # - initial short term average inflow (m3/s) and 
        #           long term average outflow (m3/s)
        if iniConditions == None:
            # read initial conditions from pcraster maps listed in the ini file (for the first time step of the model; when the model just starts)
            self.avgInflow  = vos.readPCRmapClone(iniItems.routingOptions['avgLakeReservoirInflowShortIni'],self.cloneMap,self.tmpDir,self.inputDir)
            self.avgOutflow = vos.readPCRmapClone(iniItems.routingOptions['avgLakeReservoirOutflowLongIni'],self.cloneMap,self.tmpDir,self.inputDir)
            if iniItems.routingOptions['waterBodyStorageIni'] is not None:
                self.waterBodyStorage = vos.readPCRmapClone(iniItems.routingOptions['waterBodyStorageIni'],self.cloneMap,self.tmpDir,self.inputDir)
                self.waterBodyStorage = pcr.ifthen(self.landmask, self.waterBodyStorage)
            else:
                self.waterBodyStorage = None
        else:
            # read initial conditions from the memory
            self.avgInflow        = iniConditions['routing']['avgLakeReservoirInflowShort']
            self.avgOutflow       = iniConditions['routing']['avgLakeReservoirOutflowLong']
            self.waterBodyStorage = iniConditions['routing']['waterBodyStorage']


    def getRoutingParamAvgDischarge(self, avgDischarge, dist2celllength):
        # obtain routing parameters based on average (longterm) discharge
        # output: channel dimensions and 
        #         characteristicDistance (for accuTravelTime input)
        
        yMean = self.eta * pow (avgDischarge, self.nu ) # avgDischarge in m3/s
        wMean = self.tau * pow (avgDischarge, self.phi)
 
        # option to use constant channel width (m)
        if self.constantChannelWidth is not None:\
           wMean = pcr.cover(self.constantChannelWidth, wMean)

        # minimum channel width (m)
        wMean = pcr.max(self.minChannelWidth, wMean)

        yMean =   pcr.max(yMean,0.01) # channel depth (m)
        wMean =   pcr.max(wMean,0.01) # channel width (m)
        yMean = pcr.cover(yMean,0.01)
        wMean = pcr.cover(wMean,0.01)

        # characteristicDistance (dimensionless)
        # - This will be used for accutraveltimeflux & accutraveltimestate
        # - discharge & storage = accutraveltimeflux & accutraveltimestate
        # - discharge = the total amount of material flowing through the cell (m3/s)
        # - storage   = the amount of material which is deposited in the cell (m3)
        #
        characteristicDistance = \
             ( (yMean *   wMean)/ \
               (wMean + 2*yMean) )**(2./3.) * \
              ((self.gradient)**(0.5))/ \
                self.manningsN * \
                vos.secondsPerDay()                         #  meter/day

        characteristicDistance = \
         pcr.max((self.cellSizeInArcDeg)*0.000000001,\
                 characteristicDistance/dist2celllength)    # arcDeg/day
        
        # charateristicDistance for each lake/reservoir:
        lakeReservoirCharacteristicDistance = pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
                                              pcr.areaaverage(characteristicDistance, self.WaterBodies.waterBodyIds))
        #
        # - make sure that all outflow will be released outside lakes and reservoirs
        outlets = pcr.cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyOut) > 0, pcr.boolean(1)), pcr.boolean(0))
        distance_to_outlets = pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
                              pcr.ldddist(self.lddMap, outlets, pcr.scalar(1.0)))
        lakeReservoirCharacteristicDistance = pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
                                              pcr.max(distance_to_outlets + pcr.downstreamdist(self.lddMap)*1.50, lakeReservoirCharacteristicDistance))
        #
        # TODO: calculate lakeReservoirCharacteristicDistance while obtaining lake & reservoir parameters
        
        characteristicDistance = pcr.cover(lakeReservoirCharacteristicDistance, characteristicDistance)                      
        
        # PS: In accutraveltime function: 
        #     If characteristicDistance (velocity) = 0 then:
        #     - accutraveltimestate will give zero 
        #     - accutraveltimeflux will be very high 
        
        # TODO: Consider to use downstreamdist function.
        
        # current solution: using the function "roundup" to ignore 
        #                   zero and very small values 
        characteristicDistance = \
         pcr.roundup(characteristicDistance*100.)/100.      # arcDeg/day
        
        # and set minimum value of characteristicDistance:
        characteristicDistance = pcr.cover(characteristicDistance, 0.1*self.cellSizeInArcDeg)
        characteristicDistance = pcr.max(0.100*self.cellSizeInArcDeg, characteristicDistance) # TODO: check what the minimum distance for accutraveltime function

        return (yMean, wMean, characteristicDistance)
        
    def estimate_length_of_sub_time_step(self): 

        # estimate the length of sub-time step (unit: s):
        # - the shorter is the better
        # - estimated based on the initial or latest sub-time step discharge (unit: m3/s)
        #
        length_of_sub_time_step = pcr.ifthenelse(self.subDischarge > 0.0, 
                                  self.water_height * self.dynamicFracWat * self.cellArea / self.subDischarge, vos.secondsPerDay())

        # determine the number of sub time steps (based on Rens van Beek's method - check this method with him)
        #
        critical_condition = (length_of_sub_time_step < vos.secondsPerDay())  & \
                             (self.water_height > self.critical_water_height) & \
                             (self.lddMap != pcr.ldd(5))
        #
        number_of_sub_time_steps = vos.secondsPerDay() /\
                                   pcr.cover(
                                   pcr.areaminimum(\
                                   pcr.ifthen(critical_condition, \
                                              length_of_sub_time_step),self.landmask),\
                                             vos.secondsPerDay()/self.limit_num_of_sub_time_steps)   
        number_of_sub_time_steps = 1.25 * number_of_sub_time_steps + 1
        number_of_sub_time_steps = pcr.roundup(number_of_sub_time_steps)
        #
        number_of_loops = max(1.0, pcr.cellvalue(pcr.mapmaximum(number_of_sub_time_steps),1)[1])     # minimum number of sub_time_steps = 1 
        number_of_loops = int(max(self.limit_num_of_sub_time_steps, number_of_loops))
        
        # actual length of sub-time step (s)
        length_of_sub_time_step = vos.secondsPerDay() / number_of_loops

        return (length_of_sub_time_step, number_of_loops)                               

    def simplifiedKinematicWave(self, meteo, landSurface, groundwater): 
        """
        The 'simplifiedKinematicWave':
        1. First, assume that all local fluxes has been added to 'channelStorage'. This is done outside of this function/method.
        2. Then, the 'channelStorage' is routed by using 'pcr.kinematic function' with 'lateral_inflow' = 0.0.
        """

        ##########################################################################################################################

        logger.info("Using the simplifiedKinematicWave method!")
        
        # route only non negative channelStorage (otherwise stay):
        self.channelStorage = self.channelStorage
        
        channelStorageThatWillNotMove = pcr.ifthenelse(self.channelStorage < 0.0, self.channelStorage, 0.0)
        
        # channelStorage that will be given to the ROUTING operation:
        channelStorageForRouting = pcr.max(0.0, self.channelStorage)                              # unit: m3
        
        # water height (m)
        self.water_height = pcr.min(self.max_water_height, channelStorageForRouting / (pcr.max(self.min_fracwat_for_water_height, self.dynamicFracWat) * self.cellArea))
        self.iniWater = self.water_height
        # estimate the length of sub-time step (unit: s):
        length_of_sub_time_step, number_of_loops = \
          self.estimate_length_of_sub_time_step()
              
        for i_loop in range(number_of_loops):
            
            #msg = "sub-daily time step "+str(i_loop+1)+" from "+str(number_of_loops)
            #logger.info(msg)
            
            if self.floodPlain:
                self.dynamicFracWat, self.water_height, alpha, dischargeInitial = self.kinAlphaComposite(channelStorageForRouting)
                self.dynamicFracWat = pcr.min(pcr.max(self.dynamicFracWat, self.WaterBodies.fracWat),1.0)
            else:    
                #alpha parameter and initial discharge variable needed for kinematic wave
                alpha, dischargeInitial = self.calculate_alpha_and_initial_discharge_for_kinematic_wave() 

            # at the lake/reservoir outlets, use the discharge of water bofy outflow
            waterBodyOutflowInM3PerSec = pcr.ifthen(\
                                         self.WaterBodies.waterBodyOut,
                                         self.WaterBodies.waterBodyOutflow) / vos.secondsPerDay()
            waterBodyOutflowInM3PerSec = pcr.ifthen(\
                                         pcr.scalar(self.WaterBodies.waterBodyIds) > 0.0,
                                         pcr.cover(waterBodyOutflowInM3PerSec,0.0))                                         
            dischargeInitial = pcr.cover(waterBodyOutflowInM3PerSec, dischargeInitial)                             

            # discharge (m3/s) based on kinematic wave approximation
            #logger.info('start pcr.kinematic')
            self.subDischarge = pcr.kinematic(self.lddMap, dischargeInitial, 0.0, 
                                              alpha, self.beta, \
                                              1, length_of_sub_time_step, self.cellLengthFD)
            #logger.info('done')
            
            # update channelStorage (m3)
            storage_change_in_volume  = pcr.upstream(self.lddMap, self.subDischarge * length_of_sub_time_step) - self.subDischarge * length_of_sub_time_step 
            channelStorageForRouting += storage_change_in_volume 
            #
            # route only non negative channelStorage (otherwise stay):
            channelStorageThatWillNotMove += pcr.ifthenelse(channelStorageForRouting < 0.0, channelStorageForRouting, 0.0)
            channelStorageForRouting       = pcr.max(0.000, channelStorageForRouting)
            #
            # update water_height (this will be passed to the next loop)
            if self.floodPlain != True:
                self.water_height = pcr.min(self.max_water_height, channelStorageForRouting / (pcr.max(self.min_fracwat_for_water_height, self.dynamicFracWat) * self.cellArea))

            # total discharge_volume (m3) until this present i_loop
            if i_loop == 0: discharge_volume = pcr.scalar(0.0)
            discharge_volume += self.subDischarge * length_of_sub_time_step
            
            if self.quality:
                self.channelStorageNow = pcr.max(0.0, channelStorageForRouting)
                self.qualityRouting(length_of_sub_time_step)
                self.channelStorageTimeBefore = pcr.max(0.0, self.channelStorageNow)

        # channel discharge (m3/day) = self.Q
        self.Q = discharge_volume

        # updating channelStorage (after routing)
        self.channelStorage = channelStorageForRouting

        # return channelStorageThatWillNotMove to channelStorage:
        self.channelStorage += channelStorageThatWillNotMove    
        

    def update(self,landSurface,groundwater,currTimeStep,meteo):

        logger.info("Routing in progress (one-way coupling with PCR-GLOBWB2)")

        # waterBodies: 
        # - get parameters at the beginning of each year or simulation
        # - note that the following function should be called first, specifically because  
        #   we have to define initial conditions at the beginning of simulaution, 
        #
        if currTimeStep.timeStepPCR == 1:
            initial_conditions_for_water_bodies = self.getState()
            self.WaterBodies.getParameterFiles(currTimeStep,\
                                               self.cellArea,\
                                               self.lddMap,\
                                               self.cellLengthFD,\
                                               self.cellSizeInArcDeg,\
                                               initial_conditions_for_water_bodies)               # the last line is for the initial conditions of lakes/reservoirs

        if (currTimeStep.doy == 1) and (currTimeStep.timeStepPCR > 1):
            self.WaterBodies.getParameterFiles(currTimeStep,\
                                               self.cellArea,\
                                               self.lddMap,\
                                               self.cellLengthFD,\
                                               self.cellSizeInArcDeg)
        #
        #self.WaterBodies.waterBodyIds = pcr.ifthen(self.landmask, pcr.nominal(-1)) #TODO
        # downstreamDemand (m3/s) for reservoirs 
        # - this one must be called before updating timestepsToAvgDischarge
        # - estimated based on environmental flow discharge 
        self.downstreamDemand = self.estimate_discharge_for_environmental_flow(self.channelStorage)
        
        # get routing/channel parameters/dimensions (based on avgDischarge)
        # and estimating water bodies fraction ; this is needed for calculating evaporation from water bodies
        # 
        self.yMean, self.wMean, self.characteristicDistance = \
                self.getRoutingParamAvgDischarge(self.avgDischarge,\
                self.dist2celllength)
        # 
        channelFraction = pcr.max(0.0, pcr.min(1.0,\
                          self.wMean * self.cellLengthFD / (self.cellArea)))
        if currTimeStep.timeStepPCR == 1:
            if self.floodPlain:
                self.dynamicFracWat, self.water_height = self.returnFloodedFraction(self.channelStorage)
                self.dynamicFracWat = pcr.min(pcr.max(self.dynamicFracWat, self.WaterBodies.fracWat),1.0)
            else: 
                self.dynamicFracWat = pcr.max(channelFraction, self.WaterBodies.fracWat)
            self.dynamicFracWat = pcr.ifthen(self.landmask, self.dynamicFracWat)
        
        # routing methods
        if self.method == "accuTravelTime" or self.method == "simplifiedKinematicWave": \
           self.simple_update(landSurface,groundwater,currTimeStep,meteo)
        #
        if self.method == "kinematicWave": \
           self.kinematic_wave_update(landSurface,groundwater,currTimeStep,meteo)                 
        # NOTE that this method require abstraction from fossil groundwater.
        
        # infiltration from surface water bodies (rivers/channels, as well as lakes and/or reservoirs) to groundwater bodies
        # - this exchange fluxes will be handed in the next time step
        # - in the future, this will be the interface between PCR-GLOBWB & MODFLOW (based on the difference between surface water levels & groundwater heads)
        #
        self.calculate_exchange_to_groundwater(groundwater,currTimeStep) 

        # volume water released in pits (losses: to the ocean / endorheic basin)
        self.outgoing_volume_at_pits = pcr.ifthen(self.landmask,
                                       pcr.cover(
                                       pcr.ifthen(self.lddMap == pcr.ldd(5), self.Q), 0.0))
        #
        # TODO: accumulate water in endorheic basins that are considered as lakes/reservoirs
        
        # estimate volume of water that can be extracted for abstraction in the next time step
        self.readAvlChannelStorage = self.estimate_available_volume_for_abstraction(self.channelStorage)
        
        # old-style reporting                             
        self.old_style_routing_reporting(currTimeStep)                 # TODO: remove this one

        # make sure that waterBodyStorage is only in the landmask region (as suggested by E.S. 24 okt 2018)
        self.waterBodyStorage = pcr.ifthen(self.landmask, self.waterBodyStorage)
        
    def update_routing_only(self,currTimeStep,meteo):
        
        logger.info("Routing in progress (offline configuration)")

        # waterBodies: 
        # - get parameters at the beginning of each year or simulation
        # - note that the following function should be called first, specifically because  
        #   we have to define initial conditions at the beginning of simulaution, 
        #
        if currTimeStep.timeStepPCR == 1:
            initial_conditions_for_water_bodies = self.getState()
            self.WaterBodies.getParameterFiles(currTimeStep,\
                                               self.cellArea,\
                                               self.lddMap,\
                                               self.cellLengthFD,\
                                               self.cellSizeInArcDeg,\
                                               initial_conditions_for_water_bodies)               # the last line is for the initial conditions of lakes/reservoirs

        if (currTimeStep.doy == 1) and (currTimeStep.timeStepPCR > 1):
            self.WaterBodies.getParameterFiles(currTimeStep,\
                                               self.cellArea,\
                                               self.lddMap,\
                                               self.cellLengthFD,\
                                               self.cellSizeInArcDeg)
        #
        #self.WaterBodies.waterBodyIds = pcr.ifthen(self.landmask, pcr.nominal(-1)) #TODO
        # downstreamDemand (m3/s) for reservoirs 
        # - this one must be called before updating timestepsToAvgDischarge
        # - estimated based on environmental flow discharge 
        self.downstreamDemand = self.estimate_discharge_for_environmental_flow(self.channelStorage)
        
        # get routing/channel parameters/dimensions (based on avgDischarge)
        # and estimating water bodies fraction ; this is needed for calculating evaporation from water bodies
        # 
        self.yMean, self.wMean, self.characteristicDistance = \
                self.getRoutingParamAvgDischarge(self.avgDischarge,\
                self.dist2celllength)
        # 
        channelFraction = pcr.max(0.0, pcr.min(1.0,\
                          self.wMean * self.cellLengthFD / (self.cellArea)))
        if currTimeStep.timeStepPCR == 1:
            if self.floodPlain:
                self.dynamicFracWat, self.water_height = self.returnFloodedFraction(self.channelStorage)
                self.dynamicFracWat = pcr.min(pcr.max(self.dynamicFracWat, self.WaterBodies.dynamicFracWat),1.0)
            else: 
                self.dynamicFracWat = pcr.max(channelFraction, self.WaterBodies.dynamicFracWat)
            self.dynamicFracWat = pcr.ifthen(self.landmask, self.dynamicFracWat)
        
        # routing methods
        if self.method == "accuTravelTime" or self.method == "simplifiedKinematicWave": \
           self.simple_update_routing_only(currTimeStep,meteo)
        #
        if self.method == "kinematicWave": \
           self.kinematic_wave_update(currTimeStep,meteo)                 

        # volume water released in pits (losses: to the ocean / endorheic basin)
        self.outgoing_volume_at_pits = pcr.ifthen(self.landmask,
                                       pcr.cover(
                                       pcr.ifthen(self.lddMap == pcr.ldd(5), self.Q), 0.0))
        
        # estimate volume of water that can be extracted for abstraction in the next time step
        self.readAvlChannelStorage = self.estimate_available_volume_for_abstraction(self.channelStorage)
                
        # old-style reporting                             
        self.old_style_routing_reporting(currTimeStep)

    def calculate_potential_evaporation(self,landSurface,currTimeStep,meteo):

        # potential evaporation from water bodies
        # current principle: 
        # - if landSurface.actualET < waterKC * meteo.referencePotET * self.fracWat
        #   then, we add more evaporation
        #
        if (currTimeStep.day == 1) or (currTimeStep.timeStepPCR == 1):
            waterKC = vos.netcdf2PCRobjClone(self.fileCropKC,'kc', \
                               currTimeStep.fulldate, useDoy = 'month',\
                                       cloneMapFileName = self.cloneMap)
            self.waterKC = pcr.ifthen(self.landmask,\
                           pcr.cover(waterKC, 0.0))
            self.waterKC = pcr.max(self.minCropWaterKC, self.waterKC)
            
        # potential evaporation from water bodies (m/day)) - reduced by evaporation that has been calculated in the landSurface module
        waterBodyPotEvapOvesSurfaceWaterArea = pcr.ifthen(self.landmask, \
                                               pcr.max(0.0,\
                                               self.waterKC * meteo.referencePotET -\
                                               landSurface.actualET ))              # These values are NOT over the entire cell area.
        
        # potential evaporation from water bodies over the entire cell area (m/day)
        waterBodyPotEvap = waterBodyPotEvapOvesSurfaceWaterArea * self.dynamicFracWat
        return waterBodyPotEvap
        
    def calculate_potential_evaporation_routing_only(self,currTimeStep,meteo):

        # potential evaporation from water bodies
        # current principle: 
        # - if landSurface.actualET < waterKC * meteo.referencePotET * self.fracWat
        #   then, we add more evaporation
        #
        if (currTimeStep.day == 1) or (currTimeStep.timeStepPCR == 1):
            waterKC = vos.netcdf2PCRobjClone(self.fileCropKC,'kc', \
                               currTimeStep.fulldate, useDoy = 'month',\
                                       cloneMapFileName = self.cloneMap)
            self.waterKC = pcr.ifthen(self.landmask,\
                           pcr.cover(waterKC, 0.0))
            self.waterKC = pcr.max(self.minCropWaterKC, self.waterKC)
            
        # potential evaporation from water bodies (m/day)) - reduced by evaporation that has been calculated in the landSurface module
        waterBodyPotEvapOvesSurfaceWaterArea = pcr.ifthen(self.landmask, \
                                               pcr.max(0.0,\
                                               self.waterKC * meteo.referencePotET))  # These values are NOT over the entire cell area.
        
        # potential evaporation from water bodies over the entire cell area (m/day)
        waterBodyPotEvap = waterBodyPotEvapOvesSurfaceWaterArea * self.dynamicFracWat
        return waterBodyPotEvap

    def calculate_evaporation(self,landSurface,groundwater,currTimeStep,meteo):

        # calculate potential evaporation from water bodies OVER THE ENTIRE CELL AREA (m/day) ; not only over surface water bodies
        self.waterBodyPotEvap = self.calculate_potential_evaporation(landSurface,currTimeStep,meteo)
        
        # evaporation volume from water bodies (m3)
        # - not limited to available channelStorage 
        volLocEvapWaterBody = self.waterBodyPotEvap * self.cellArea
        # - limited to available channelStorage
        volLocEvapWaterBody = pcr.min(\
                              pcr.max(0.0,self.channelStorage), volLocEvapWaterBody)

        # update channelStorage (m3) after evaporation from water bodies
        self.channelStorage = self.channelStorage -\
                              volLocEvapWaterBody
        self.local_input_to_surface_water -= volLocEvapWaterBody
        
        # evaporation (m) from water bodies                             
        self.waterBodyEvaporation = volLocEvapWaterBody / self.cellArea
        self.waterBodyEvaporation = pcr.ifthen(self.landmask, self.waterBodyEvaporation)

        # remaining potential evaporation (m) from water bodies
        self.remainWaterBodyPotEvap = pcr.max(0.0, self.waterBodyPotEvap - self.waterBodyEvaporation)
        
    def calculate_evaporation_routing_only(self,currTimeStep,meteo):

        # calculate potential evaporation from water bodies OVER THE ENTIRE CELL AREA (m/day) ; not only over surface water bodies
        self.waterBodyPotEvap = self.calculate_potential_evaporation_routing_only(currTimeStep,meteo)
        
        # evaporation volume from water bodies (m3)
        # - not limited to available channelStorage 
        volLocEvapWaterBody = self.waterBodyPotEvap * self.cellArea
        # - limited to available channelStorage
        volLocEvapWaterBody = pcr.min(\
                              pcr.max(0.0,self.channelStorage), volLocEvapWaterBody)

        # update channelStorage (m3) after evaporation from water bodies
        self.channelStorage = self.channelStorage -\
                              volLocEvapWaterBody
        self.local_input_to_surface_water -= volLocEvapWaterBody
        
        # evaporation (m) from water bodies                             
        self.waterBodyEvaporation = volLocEvapWaterBody / self.cellArea
        self.waterBodyEvaporation = pcr.ifthen(self.landmask, self.waterBodyEvaporation)

        # remaining potential evaporation (m) from water bodies
        self.remainWaterBodyPotEvap = pcr.max(0.0, self.waterBodyPotEvap - self.waterBodyEvaporation)
        
    def calculate_extra_evaporation(self):
		# limited to self.remainWaterBodyPotEvap: remaining potential evaporation (m) from water bodies

        # evaporation volume from water bodies (m3) - limited to available channelStorage
        volLocEvapWaterBody = pcr.min(\
                              pcr.max(0.0,self.channelStorage),
                              self.remainWaterBodyPotEvap * self.dynamicFracWat * self.cellArea)

        # update channelStorage (m3) after evaporation from water bodies
        self.channelStorage = self.channelStorage -\
                              volLocEvapWaterBody
        self.local_input_to_surface_water -= volLocEvapWaterBody
        
        # update evaporation (m) from water bodies                             
        self.waterBodyEvaporation += volLocEvapWaterBody / self.cellArea

        # remaining potential evaporation (m) from water bodies
        self.remainWaterBodyPotEvap = pcr.max(0.0, self.remainWaterBodyPotEvap - volLocEvapWaterBody / self.cellArea)

    def calculate_exchange_to_groundwater(self,groundwater,currTimeStep):

        if self.debugWaterBalance:\
           preStorage = self.channelStorage                            # unit: m3

        # riverbed infiltration (m3/day):
        #
        # - current implementation based on Inge's principle (later, will be based on groundater head (MODFLOW) and can be negative)
        # - happening only if 0.0 < baseflow < total_groundwater_abstraction
        # - total_groundwater_abstraction = groundwater.nonFossilGroundwaterAbs + groundwater.unmetDemand
        # - infiltration rate will be based on aquifer saturated conductivity
        # - limited to fracWat
        # - limited to available channelStorage
        # - this infiltration will be handed to groundwater in the next time step
        # - References: de Graaf et al. (2014); Wada et al. (2012); Wada et al. (2010)
        #
        riverbedConductivity  = groundwater.kSatAquifer # unit: m/day
        total_groundwater_abstraction = pcr.max(0.0, groundwater.nonFossilGroundwaterAbs + groundwater.unmetDemand)   # unit: m
        self.riverbedExchange = pcr.max(0.0,\
                                pcr.min(pcr.max(0.0,self.channelStorage),\
                                pcr.ifthenelse(groundwater.baseflow > 0.0, \
                                pcr.ifthenelse(total_groundwater_abstraction > groundwater.baseflow, \
                                riverbedConductivity * self.dynamicFracWat * self.cellArea, \
                                0.0), 0.0)))
        self.riverbedExchange = pcr.cover(self.riverbedExchange, 0.0)                         
        factor = 0.05 # to avoid flip flop
        self.riverbedExchange = pcr.min(self.riverbedExchange, (1.0-factor)*pcr.max(0.0,self.channelStorage))                                                             
        self.riverbedExchange = pcr.ifthenelse(self.channelStorage < 0.0, 0.0, self.riverbedExchange)
        self.riverbedExchange = pcr.cover(self.riverbedExchange, 0.0)
        self.riverbedExchange = pcr.ifthen(self.landmask, self.riverbedExchange)

        # update channelStorage (m3) after riverbedExchange (m3)
        self.channelStorage  -= self.riverbedExchange
        self.local_input_to_surface_water -= self.riverbedExchange

        if self.debugWaterBalance:\
           vos.waterBalanceCheck([pcr.scalar(0.0)],\
                                 [self.riverbedExchange/self.cellArea],\
                                 [           preStorage/self.cellArea],\
                                 [  self.channelStorage/self.cellArea],\
                                   'channelStorage after surface water infiltration',\
                                  True,\
                                  currTimeStep.fulldate,threshold=1e-4)


    def reduce_unmet_demand(self,landSurface,groundwater,currTimeStep):

        logger.info("Reducing unmetDemand by allowing extra surface water abstraction.")

        extra_surface_water_abstraction = pcr.scalar(0.0)
        reduction_for_unmetDemand       = pcr.scalar(0.0)
        
        # estimate channel storage that can be extracted (unit: m3)
        self.readAvlChannelStorage = self.estimate_available_volume_for_abstraction(self.channelStorage)
        
        # potential_unmet_demand (unit: m) 
        potential_unmet_demand = landSurface.totalPotentialGrossDemand -\
                                 landSurface.allocSurfaceWaterAbstract -\
                                 groundwater.allocNonFossilGroundwater

        if self.debugWaterBalance:
            test = pcr.ifthen(potential_unmet_demand < 0.0, potential_unmet_demand)
            a,b,c = vos.getMinMaxMean(pcr.scalar(test),True)
            threshold = 1e-3
            if abs(a) > threshold or abs(b) > threshold:
                logger.info("WARNING !!!!! Water Balance Error. There is negative unmetDemand ... Min %f Max %f Mean %f" %(a,b,c))

        if landSurface.usingAllocSegments == False and landSurface.limitAbstraction == False:
        
            logger.info("Surface water abstraction is only to satisfy local demand. No network.")
            
            # reduction_for_unmetDemand
            reduction_for_unmetDemand = pcr.min(self.readAvlChannelStorage / self.cellArea, \
                                                potential_unmet_demand)                           # unit: m

            # actual extra surface water abstraction in meter 
            extra_surface_water_abstraction = pcr.ifthen(self.landmask, reduction_for_unmetDemand)
                                                
            
        if landSurface.usingAllocSegments == True and landSurface.limitAbstraction == False:
        
            # TODO: Assuming that there is also network for distributing groundwater abstractions.
            # Notes: Incorporating distribution network of groundwater source is possible only if limitAbstraction = False.  

            logger.info("Using allocation to reduce unmetDemand.")

            # gross/potential demand volume in each cell (unit: m3) - ignore small values (less than 1 m3)
            cellVolGrossDemand = pcr.rounddown(
                                 potential_unmet_demand*self.cellArea)
            
            # demand in each segment/zone (unit: m3)
            segTtlGrossDemand  = pcr.areatotal(cellVolGrossDemand, landSurface.allocSegments)
            
            # total available water volume in each cell - ignore small values (less than 1 m3)
            cellAvlWater = pcr.rounddown(pcr.max(0.00, self.readAvlChannelStorage))
            
            # total available surface water volume in each segment/zone  (unit: m3)
            segAvlWater  = pcr.areatotal(cellAvlWater, landSurface.allocSegments)
            
            # total actual extra surface water abstraction volume in each segment/zone (unit: m3)
            # - limited to available water
            segActWaterAbs = pcr.min(segAvlWater, segTtlGrossDemand)
            
            # actual extra surface water abstraction volume in each cell (unit: m3)
            volActWaterAbstract = vos.getValDivZero(\
                                  cellAvlWater, segAvlWater, vos.smallNumber) * \
                                  segActWaterAbs                                                 
            volActWaterAbstract = pcr.min(cellAvlWater,volActWaterAbstract)                               # unit: m3
            
            # actual extra surface water abstraction in meter 
            extra_surface_water_abstraction    = pcr.ifthen(self.landmask, volActWaterAbstract) /\
                                                                             self.cellArea                # unit: m

            # allocation extra surface water abstraction volume to each cell (unit: m3)
            extraVolAllocSurfaceWaterAbstract  = vos.getValDivZero(\
                                                 cellVolGrossDemand, segTtlGrossDemand, vos.smallNumber) *\
                                                 segActWaterAbs                                           # unit: m3 
            # reduction for unmetDemand (unit: m)
            reduction_for_unmetDemand = pcr.ifthen(self.landmask, 
                                        extraVolAllocSurfaceWaterAbstract / self.cellArea)                # unit: m
            
            if self.debugWaterBalance:
    
                abstraction = pcr.cover(pcr.areatotal(volActWaterAbstract              , landSurface.allocSegments)/landSurface.segmentArea, 0.0)
                allocation  = pcr.cover(pcr.areatotal(extraVolAllocSurfaceWaterAbstract, landSurface.allocSegments)/landSurface.segmentArea, 0.0)
            
                vos.waterBalanceCheck([pcr.ifthen(self.landmask,abstraction)],\
                                      [pcr.ifthen(self.landmask, allocation)],\
                                      [pcr.scalar(0.0)],\
                                      [pcr.scalar(0.0)],\
                                      'extra surface water abstraction - allocation per zone/segment (PS: Error here may be caused by rounding error.)' ,\
                                       True,\
                                       "",threshold=5e-4)

        # correcting surface water abstraction 
        landSurface.actSurfaceWaterAbstract   += extra_surface_water_abstraction                # unit: m
            
        # update channelStorage (m3) after extra_surface_water_abstraction
        self.channelStorage = self.channelStorage -\
                              extra_surface_water_abstraction * self.cellArea
        self.local_input_to_surface_water -= extra_surface_water_abstraction * self.cellArea

        # correcting surface water allocation after reduction of unmetDemand
        landSurface.allocSurfaceWaterAbstract += reduction_for_unmetDemand                      # unit: m

        # recalculating unmetDemand (m)
        groundwater.unmetDemand =  landSurface.totalPotentialGrossDemand -\
                                   landSurface.allocSurfaceWaterAbstract -\
                                   groundwater.allocNonFossilGroundwater

        if self.debugWaterBalance:

            test = pcr.ifthen(groundwater.unmetDemand < 0.0, groundwater.unmetDemand)
            a,b,c = vos.getMinMaxMean(pcr.scalar(test),True)
            threshold = 1e-3
            if abs(a) > threshold or abs(b) > threshold:
                logger.info("WARNING !!!!! Water Balance Error. There is negative unmetDemand ... Min %f Max %f Mean %f" %(a,b,c))

    def simple_update(self,landSurface,groundwater,currTimeStep,meteo):

        # updating timesteps to calculate long and short term statistics values of avgDischarge, avgInflow, avgOutflow, etc.
        self.timestepsToAvgDischarge += 1.

        if self.debugWaterBalance:\
           preStorage = self.channelStorage                                                        # unit: m3

        # the following variable defines total local change (input) to surface water storage bodies # unit: m3 
        # - only local processes; therefore not considering any routing processes
        self.local_input_to_surface_water = pcr.scalar(0.0)          # initiate the variable, start from zero

        # runoff from landSurface cells (unit: m/day)
        self.runoff = landSurface.landSurfaceRunoff +\
                      groundwater.baseflow   
        
        # update channelStorage (unit: m3) after runoff
        self.channelStorage += self.runoff * self.cellArea
        self.local_input_to_surface_water += self.runoff * self.cellArea

        # update channelStorage (unit: m3) after actSurfaceWaterAbstraction 
        self.channelStorage -= landSurface.actSurfaceWaterAbstract * self.cellArea
        self.local_input_to_surface_water -= landSurface.actSurfaceWaterAbstract * self.cellArea

        # reporting channelStorage after surface water abstraction (unit: m3)
        self.channelStorageAfterAbstraction = pcr.ifthen(self.landmask, self.channelStorage) 
        
        # return flow from domestic water demand (m day-1)
        self.DomesticReturnFlow = landSurface.domesticGrossDemand * landSurface.DomReturnFlowFraction 
        self.DomesticReturnFlowVol = self.DomesticReturnFlow * self.cellArea #in m3 day (assuming all demands are met.)
        
        # return flow from industry water demands (m day-1)
        self.IndustryReturnFlow = landSurface.industryGrossDemand * landSurface.IndReturnFlowFraction
        self.IndustryReturnFlowVol = self.IndustryReturnFlow * self.cellArea #in m3 day (assuming all demands are met.)
        
        # return flow from (m) non irrigation water demand
        self.nonIrrReturnFlow = landSurface.nonIrrReturnFlowFraction * landSurface.potentialNonIrrGrossWaterDemand
        nonIrrReturnFlowVol   = self.nonIrrReturnFlow*self.cellArea
        self.channelStorage  += nonIrrReturnFlowVol
        self.local_input_to_surface_water += nonIrrReturnFlowVol
        
        # water consumption for non irrigation water demand (m) - this water is removed from the water balance
        self.nonIrrWaterConsumption = pcr.max(0.0,\
                                      landSurface.potentialNonIrrGrossWaterDemand - \
                                      self.nonIrrReturnFlow)

        # Note that in case of limitAbstraction = True ; landSurface.nonIrrGrossDemand has been reduced by available water                               
        
        # calculate evaporation from water bodies - this will return self.waterBodyEvaporation (unit: m)
        self.calculate_evaporation(landSurface,groundwater,currTimeStep,meteo)
        
        if self.debugWaterBalance:\
           vos.waterBalanceCheck([self.runoff,\
                                  self.nonIrrReturnFlow],\
                                 [landSurface.actSurfaceWaterAbstract,self.waterBodyEvaporation],\
                                 [           preStorage/self.cellArea],\
                                 [  self.channelStorage/self.cellArea],\
                                   'channelStorage (unit: m) before lake/reservoir outflow',\
                                  True,\
                                  currTimeStep.fulldate,threshold=5e-3)
        
        # LAKE AND RESERVOIR OPERATIONS
        ##########################################################################################################################
        if self.debugWaterBalance: \
           preStorage = self.channelStorage                                  # unit: m3

        # at cells where lakes and/or reservoirs defined, move channelStorage to waterBodyStorage
        #
        storageAtLakeAndReservoirs = \
         pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
                               self.channelStorage)
        storageAtLakeAndReservoirs = pcr.cover(storageAtLakeAndReservoirs,0.0)
        #
        # - move only non negative values and use rounddown values
        storageAtLakeAndReservoirs = pcr.max(0.00, pcr.rounddown(storageAtLakeAndReservoirs))
        self.channelStorage -= storageAtLakeAndReservoirs                    # unit: m3
        
        # update waterBodyStorage (inflow, storage and outflow)
        self.WaterBodies.update(storageAtLakeAndReservoirs,\
                                self.timestepsToAvgDischarge,\
                                self.maxTimestepsToAvgDischargeShort,\
                                self.maxTimestepsToAvgDischargeLong,\
                                currTimeStep,\
                                self.avgDischarge,\
                                vos.secondsPerDay(),\
                                self.downstreamDemand)

        # waterBodyStorage (m3) after outflow:                               # values given are per water body id (not per cell)
        self.waterBodyStorage = pcr.ifthen(self.landmask,self.WaterBodies.waterBodyStorage)
        
        if self.quality:    
            self.waterBodyStorageTimeBefore = self.waterBodyStorage + self.WaterBodies.waterBodyOutflow
            self.waterBodyOutFlowDay = pcr.cover(\
                           pcr.ifthen(\
                           self.WaterBodies.waterBodyOut,
                           self.WaterBodies.waterBodyOutflow), 0.0)  
        
        # transfer outflow from lakes and/or reservoirs to channelStorages
        waterBodyOutflow = pcr.cover(\
                           pcr.ifthen(\
                           self.WaterBodies.waterBodyOut,
                           self.WaterBodies.waterBodyOutflow), 0.0)          # unit: m3/day
        
        if self.method == "accuTravelTime":
            # distribute outflow to water body storage
            # - this is to avoid 'waterBodyOutflow' skipping cells 
            # - this is done by distributing waterBodyOutflow within lake/reservoir cells 
            #
            waterBodyOutflow = pcr.areaaverage(waterBodyOutflow, self.WaterBodies.waterBodyIds)
            waterBodyOutflow = pcr.ifthen(\
                               pcr.scalar(self.WaterBodies.waterBodyIds) > 0.0,
                               waterBodyOutflow)                                 
        self.waterBodyOutflow = pcr.cover(waterBodyOutflow, 0.0)             # unit: m3/day

        # update channelStorage (m3) after waterBodyOutflow (m3)
        self.channelStorage += self.waterBodyOutflow
        # Note that local_input_to_surface_water does not include waterBodyOutflow
        
        if self.debugWaterBalance:\
           vos.waterBalanceCheck([self.waterBodyOutflow/self.cellArea],\
                                 [storageAtLakeAndReservoirs/self.cellArea],\
                                 [           preStorage/self.cellArea],\
                                 [  self.channelStorage/self.cellArea],\
                                   'channelStorage (unit: m) after lake reservoir/outflow fluxes (errors here are most likely due to pcraster implementation in float_32)',\
                                  True,\
                                  currTimeStep.fulldate,threshold=1e-3)

        if self.quality:
            # Input data provided at daily resolution, unless otherwise adjusted in function
            self.readExtensiveMeteo(currTimeStep)
            
            # Input data for pollutant loadings typically provided at monthly resolution
            if currTimeStep.day == 1:
                self.readPowerplantData(currTimeStep)
            
                if self.calculateLoads:
                    self.readPollutantLoadingsInputData(currTimeStep)  
            
            # Pollutant loadings calculated (or read in) at daily resolution, unless otherwise adjusted in function
            if self.calculateLoads:
                self.calculatePollutantLoadings(currTimeStep,landSurface,groundwater)
            else:
                self.readPollutantLoadings(currTimeStep)  
                
            self.channelStorageTimeBefore = pcr.max(0.0, self.channelStorage)
            self.energyLocal(meteo, landSurface, groundwater)
            self.energyWaterBody()
            self.qualityLocal()
            self.qualityWaterBody()            
                                  
        # ROUTING OPERATION:
        ##########################################################################################################################
        # - this will return new self.channelStorage (but still without waterBodyStorage)
        # - also, this will return self.Q which is channel discharge in m3/day
        #
        if self.method == "accuTravelTime":          self.accuTravelTime() 		
        if self.method == "simplifiedKinematicWave": self.simplifiedKinematicWave(meteo, landSurface, groundwater) 		
        #
        #
        # channel discharge (m3/s): for current time step
        #
        self.discharge = self.Q / vos.secondsPerDay()
        self.discharge = pcr.max(0., self.discharge)                   # reported channel discharge cannot be negative
        self.discharge = pcr.ifthen(self.landmask, self.discharge)
        #
        self.disChanWaterBody = pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,\
                                pcr.areamaximum(self.discharge,self.WaterBodies.waterBodyIds))
        self.disChanWaterBody = pcr.cover(self.disChanWaterBody, self.discharge)
        self.disChanWaterBody = pcr.ifthen(self.landmask, self.disChanWaterBody)
        #
        self.disChanWaterBody = pcr.max(0.,self.disChanWaterBody)      # reported channel discharge cannot be negative
        #
        #
        ##########################################################################################################################

        # calculate the statistics of long and short term flow values
        self.calculate_statistics(groundwater, landSurface)
        
        self.allow_extra_evaporation_and_abstraction = False # This option is still EXPERIMENTAL (and not recommended)
        if self.allow_extra_evaporation_and_abstraction:\
           self.update_with_extra_evaporation_and_unmet_demand_reduction()

        if self.quality:
            self.energyWaterBodyAverage()
            self.qualityWaterBodyAverage()       
         
        self.channelStorage = self.return_water_body_storage_to_channel(self.channelStorage) # return waterBodyStorage to channelStorage 

    def simple_update_routing_only(self,currTimeStep,meteo):

        if self.offlineRun:
            logger.info("Reading extensive hydrology: total runoff = baseflow + interflow + surface runoff")
            self.readExtensiveHydro(currTimeStep) #daily timestep

        # updating timesteps to calculate long and short term statistics values of avgDischarge, avgInflow, avgOutflow, etc.
        self.timestepsToAvgDischarge += 1.

        if self.debugWaterBalance:\
           preStorage = self.channelStorage                                                        # unit: m3

        # the following variable defines total local change (input) to surface water storage bodies # unit: m3 
        # - only local processes; therefore not considering any routing processes
        self.local_input_to_surface_water = pcr.scalar(0.0)          # initiate the variable, start from zero
        
        # update channelStorage (unit: m3) after runoff
        self.channelStorage += self.runoff * self.cellArea
        self.local_input_to_surface_water += self.runoff * self.cellArea

        # reporting channelStorage after surface water abstraction (unit: m3)
        self.channelStorageAfterAbstraction = pcr.ifthen(self.landmask, self.channelStorage) 

        # calculate evaporation from water bodies - this will return self.waterBodyEvaporation (unit: m)
        self.calculate_evaporation_routing_only(currTimeStep,meteo)
        
        if self.debugWaterBalance:\
           vos.waterBalanceCheck([self.runoff],\
                                 [self.waterBodyEvaporation],\
                                 [           preStorage/self.cellArea],\
                                 [  self.channelStorage/self.cellArea],\
                                   'channelStorage (unit: m) before lake/reservoir outflow',\
                                  True,\
                                  currTimeStep.fulldate,threshold=1e-4)
        
        # LAKE AND RESERVOIR OPERATIONS
        ##########################################################################################################################
        if self.debugWaterBalance: \
           preStorage = self.channelStorage                                  # unit: m3

        # at cells where lakes and/or reservoirs defined, move channelStorage to waterBodyStorage
        #
        storageAtLakeAndReservoirs = \
         pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
                               self.channelStorage)
        storageAtLakeAndReservoirs = pcr.cover(storageAtLakeAndReservoirs,0.0)
        #
        # - move only non negative values and use rounddown values
        storageAtLakeAndReservoirs = pcr.max(0.00, pcr.rounddown(storageAtLakeAndReservoirs))
        self.channelStorage -= storageAtLakeAndReservoirs                    # unit: m3
        
        # update waterBodyStorage (inflow, storage and outflow)
        self.WaterBodies.update(storageAtLakeAndReservoirs,\
                                self.timestepsToAvgDischarge,\
                                self.maxTimestepsToAvgDischargeShort,\
                                self.maxTimestepsToAvgDischargeLong,\
                                currTimeStep,\
                                self.avgDischarge,\
                                vos.secondsPerDay(),\
                                self.downstreamDemand)

        # waterBodyStorage (m3) after outflow:                               # values given are per water body id (not per cell)
        self.waterBodyStorage = self.WaterBodies.waterBodyStorage
        
        if self.quality:    
            self.waterBodyStorageTimeBefore = self.waterBodyStorage + self.WaterBodies.waterBodyOutflow
            self.waterBodyOutFlowDay = pcr.cover(\
                           pcr.ifthen(\
                           self.WaterBodies.waterBodyOut,
                           self.WaterBodies.waterBodyOutflow), 0.0)  
        
        # transfer outflow from lakes and/or reservoirs to channelStorages
        waterBodyOutflow = pcr.cover(\
                           pcr.ifthen(\
                           self.WaterBodies.waterBodyOut,
                           self.WaterBodies.waterBodyOutflow), 0.0)          # unit: m3/day
        
        if self.method == "accuTravelTime":
            # distribute outflow to water body storage
            # - this is to avoid 'waterBodyOutflow' skipping cells 
            # - this is done by distributing waterBodyOutflow within lake/reservoir cells 
            #
            waterBodyOutflow = pcr.areaaverage(waterBodyOutflow, self.WaterBodies.waterBodyIds)
            waterBodyOutflow = pcr.ifthen(\
                               pcr.scalar(self.WaterBodies.waterBodyIds) > 0.0,
                               waterBodyOutflow)                                 
        self.waterBodyOutflow = pcr.cover(waterBodyOutflow, 0.0)             # unit: m3/day

        # update channelStorage (m3) after waterBodyOutflow (m3)
        self.channelStorage += self.waterBodyOutflow
        # Note that local_input_to_surface_water does not include waterBodyOutflow
        
        if self.debugWaterBalance:\
           vos.waterBalanceCheck([self.waterBodyOutflow/self.cellArea],\
                                 [storageAtLakeAndReservoirs/self.cellArea],\
                                 [           preStorage/self.cellArea],\
                                 [  self.channelStorage/self.cellArea],\
                                   'channelStorage (unit: m) after lake reservoir/outflow fluxes (errors here are most likely due to pcraster implementation in float_32)',\
                                  True,\
                                  currTimeStep.fulldate,threshold=1e-3)

        if self.quality:
            # Input data provided at daily resolution, unless otherwise adjusted in function
            self.readExtensiveMeteo(currTimeStep)
            self.readPowerplantData(currTimeStep)
            self.readPollutantLoadings(currTimeStep)

            self.channelStorageTimeBefore = pcr.max(0.0, self.channelStorage)
            self.energyLocal_routing_only(meteo)
            self.energyWaterBody()
            self.qualityLocal()
            self.qualityWaterBody()
        
        # ROUTING OPERATION:
        ##########################################################################################################################
        # - this will return new self.channelStorage (but still without waterBodyStorage)
        # - also, this will return self.Q which is channel discharge in m3/day
        #
        if self.method == "accuTravelTime":          self.accuTravelTime() 		
        if self.method == "simplifiedKinematicWave": self.simplifiedKinematicWave(meteo, 0.0,0.0) 		
        #
        #
        # channel discharge (m3/s): for current time step
        #
        self.discharge = self.Q / vos.secondsPerDay()
        self.discharge = pcr.max(0., self.discharge)                   # reported channel discharge cannot be negative
        self.discharge = pcr.ifthen(self.landmask, self.discharge)
        #
        self.disChanWaterBody = pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,\
                                pcr.areamaximum(self.discharge,self.WaterBodies.waterBodyIds))
        self.disChanWaterBody = pcr.cover(self.disChanWaterBody, self.discharge)
        self.disChanWaterBody = pcr.ifthen(self.landmask, self.disChanWaterBody)
        #
        self.disChanWaterBody = pcr.max(0.,self.disChanWaterBody)      # reported channel discharge cannot be negative
        #
        #
        ##########################################################################################################################

        self.allow_extra_evaporation_and_abstraction = False # This option is still EXPERIMENTAL (and not recommended)
        if self.allow_extra_evaporation_and_abstraction:\
           self.update_with_extra_evaporation_and_unmet_demand_reduction()          

        if self.quality:
            self.energyWaterBodyAverage() 
            self.qualityWaterBodyAverage()
           
        # return waterBodyStorage to channelStorage  
        self.channelStorage = self.return_water_body_storage_to_channel(self.channelStorage)

    def update_with_extra_evaporation_and_unmet_demand_reduction(self): 
        # This function is still EXPERIMENTAL (and not recommended)

        # add extra evaporation
        self.calculate_extra_evaporation()
        # reduce fossil groundwater storage abstraction (unmetDemand)
        if groundwater.limitAbstraction == False: self.reduce_unmet_demand(landSurface,groundwater,currTimeStep) 

    def calculate_alpha_and_initial_discharge_for_kinematic_wave(self): 

        # calculate alpha (dimensionless), which is the roughness coefficient 
        # - for kinewatic wave (see: http://pcraster.geo.uu.nl/pcraster/4.0.0/doc/manual/op_kinematic.html)
        # - based on wetted area (m2) and wetted perimeter (m), as well as self.beta (dimensionless)
        # - assuming rectangular channel with channel_width = self.wMean and channel_length = self.dist2celllength
        #
        if self.quality: 
            manIce= pcr.max(self.manningsN,\
              0.0493*pcr.max(0.01,self.channelStorage/(self.dynamicFracWat*self.cellArea))**\
              (-0.23)*self.iceThickness**0.57)
            manningsWithIce= (0.5*(self.manningsN**1.5+manIce**1.5))**(2./3.)
            wetA= self.channelStorage/self.cellLengthFD
            wetP= 2.*wetA/self.wMean+self.wMean
            alpha= (manningsWithIce*wetP**(2./3.)*self.gradient**-0.5)**self.beta
            #TODO alpha = (self.manningsN*wetP**(2./3.)*self.gradient**(-0.5))**self.beta            
        else:
            channel_wetted_area      =   self.water_height * self.wMean                                  # unit: m2
            channel_wetted_perimeter = 2*self.water_height + self.wMean                                  # unit: m  
        #
        
            alpha = (self.manningsN*channel_wetted_perimeter**(2./3.)*self.gradient**(-0.5))**self.beta  # dimensionless
        
        # estimate of channel discharge (m3/s) based on water height
        #
        dischargeInitial = pcr.ifthenelse(alpha > 0.0,\
                                         (self.water_height * self.wMean / alpha)**(1/self.beta),0.0)
        return (alpha, dischargeInitial)    


    def return_water_body_storage_to_channel(self, channelStorage):

        # return waterBodyStorage to channelStorage  
        #
        waterBodyStorageTotal = \
         pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areaaverage(\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyStorage),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)) + \
         pcr.areatotal(pcr.cover(\
         pcr.ifthen(self.landmask,channelStorage), 0.0),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)))
        waterBodyStoragePerCell = \
         waterBodyStorageTotal*\
                       self.cellArea/\
         pcr.areatotal(pcr.cover(\
         self.cellArea, 0.0),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))
        waterBodyStoragePerCell = \
         pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         waterBodyStoragePerCell)                                                      # unit: m3
        #
        channelStorage = pcr.cover(waterBodyStoragePerCell, channelStorage)            # unit: m3
        channelStorage = pcr.ifthen(self.landmask, channelStorage)
        return channelStorage

        
    def calculate_statistics(self, groundwater, landSurface):

        # short term average inflow (m3/s) and long term average outflow (m3/s) from lake and reservoirs
        self.avgInflow  = pcr.ifthen(self.landmask, pcr.cover(self.WaterBodies.avgInflow , 0.0)) 
        self.avgOutflow = pcr.ifthen(self.landmask, pcr.cover(self.WaterBodies.avgOutflow, 0.0))

        # short term and long term average discharge (m3/s)
        # - see: online algorithm on http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        #
        # - long term average diadvectedEnergyPreciprge
        #
        diadvectedEnergyPreciprgeUsed      = pcr.max(0.0, self.discharge)
        diadvectedEnergyPreciprgeUsed      = pcr.max(diadvectedEnergyPreciprgeUsed, self.disChanWaterBody)
        #
        deltaAnoDischarge = diadvectedEnergyPreciprgeUsed - self.avgDischarge  
        self.avgDischarge = self.avgDischarge +\
                            deltaAnoDischarge/\
                            pcr.min(self.maxTimestepsToAvgDischargeLong, self.timestepsToAvgDischarge)
        self.avgDischarge = pcr.max(0.0, self.avgDischarge)                                    
        self.m2tDischarge = self.m2tDischarge + pcr.abs(deltaAnoDischarge*(diadvectedEnergyPreciprgeUsed - self.avgDischarge))                             
        #
        # - short term average diadvectedEnergyPreciprge
        #
        deltaAnoDischargeShort = diadvectedEnergyPreciprgeUsed - self.avgDischargeShort  
        self.avgDischargeShort = self.avgDischargeShort +\
                                 deltaAnoDischargeShort/\
                                 pcr.min(self.maxTimestepsToAvgDischargeShort, self.timestepsToAvgDischarge)
        self.avgDischargeShort = pcr.max(0.0, self.avgDischargeShort)                         

        # long term average baseflow (m3/s) ; used as proxies for partitioning groundwater and surface water abstractions
        #
        baseflowM3PerSec = groundwater.baseflow * self.cellArea / vos.secondsPerDay()
        deltaAnoBaseflow = baseflowM3PerSec - self.avgBaseflow  
        self.avgBaseflow = self.avgBaseflow +\
                           deltaAnoBaseflow/\
                           pcr.min(self.maxTimestepsToAvgDischargeLong, self.timestepsToAvgDischarge)                
        self.avgBaseflow = pcr.max(0.0, self.avgBaseflow)

        # average irrigation water that is allocated from the last 30 days (needed for online runs where loadings are calculated in loop)
        if self.calculateLoads and self.offlineRun == False: 
            # calculate average irrigation gross demands from the last 30 days       
            deltaAno_irrGrossDemand = pcr.max(0.0, landSurface.irrGrossDemand) - self.avg_irrGrossDemand  
            self.avg_irrGrossDemand = self.avg_irrGrossDemand +\
                                      deltaAno_irrGrossDemand/\
                                      pcr.min(30.0, self.timestepsToAvgDischarge)
            self.avg_irrGrossDemand = pcr.max(0.0, self.avg_irrGrossDemand)                 
        
            # average netLqWaterToSoil from the last 30 days
            deltaAno_netLqWaterToSoil = pcr.max(0.0, landSurface.netLqWaterToSoil) - self.avg_netLqWaterToSoil  
            self.avg_netLqWaterToSoil = self.avg_netLqWaterToSoil +\
                                        deltaAno_netLqWaterToSoil/\
                                        pcr.min(30.0, self.timestepsToAvgDischarge)
            self.avg_netLqWaterToSoil = pcr.max(0.0, self.avg_netLqWaterToSoil)                         

    def estimate_discharge_for_environmental_flow(self, channelStorage):

        # long term variance and standard deviation of discharge values
        varDischarge = self.m2tDischarge / \
                       pcr.max(1.,\
                       pcr.min(self.maxTimestepsToAvgDischargeLong, self.timestepsToAvgDischarge)-1.)                             
                       # see: online algorithm on http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        stdDischarge = pcr.max(varDischarge**0.5, 0.0)
        
        # calculate minimum discharge for environmental flow (m3/s)
        minDischargeForEnvironmentalFlow = pcr.max(0.001, self.avgDischarge - 3.5*stdDischarge)
        factor = 0.10 # to avoid flip flop
        minDischargeForEnvironmentalFlow = pcr.max(factor*self.avgDischarge, minDischargeForEnvironmentalFlow)   # unit: m3/s
        return minDischargeForEnvironmentalFlow


    def estimate_available_volume_for_abstraction(self, channelStorage):
        # input: channelStorage    in m3

        # estimate minimum discharge for environmental flow (m3/s)
        minDischargeForEnvironmentalFlow = self.estimate_discharge_for_environmental_flow(channelStorage)
        
        # available channelStorage that can be extracted for surface water abstraction
        readAvlChannelStorage = pcr.max(0.0,channelStorage)                                                             
        
        # safety factor to reduce readAvlChannelStorage
        safety_factor = vos.getValDivZero(pcr.max(0.0, pcr.min(self.avgDischargeShort, self.avgDischarge)), \
                                          minDischargeForEnvironmentalFlow, vos.smallNumber)
        safety_factor = pcr.min(1.00, pcr.max(0.00, safety_factor))
        readAvlChannelStorage = safety_factor * pcr.max(0.0, readAvlChannelStorage)                                                             

        # ignore small volume values - less than 1 m3
        readAvlChannelStorage = pcr.rounddown(readAvlChannelStorage*1.)/1.
        readAvlChannelStorage = pcr.ifthen(self.landmask, readAvlChannelStorage)
        return readAvlChannelStorage       # unit: m3

    def initiate_old_style_routing_reporting(self,iniItems):

        self.report = True
        try:
            self.outDailyTotNC = iniItems.routingOptions['outDailyTotNC'].split(",")
            self.outMonthTotNC = iniItems.routingOptions['outMonthTotNC'].split(",")
            self.outMonthAvgNC = iniItems.routingOptions['outMonthAvgNC'].split(",")
            self.outMonthEndNC = iniItems.routingOptions['outMonthEndNC'].split(",")
            self.outAnnuaTotNC = iniItems.routingOptions['outAnnuaTotNC'].split(",")
            self.outAnnuaAvgNC = iniItems.routingOptions['outAnnuaAvgNC'].split(",")
            self.outAnnuaEndNC = iniItems.routingOptions['outAnnuaEndNC'].split(",")
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

    def old_style_routing_reporting(self,currTimeStep):

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

    def returnFloodedFraction(self,channelStorage):
        #-returns the flooded fraction given the flood volume and the associated water height
        # using a logistic smoother near intersections (K&K, 2007)
        #-find the match on the basis of the shortest distance to the available intersections or steps
        deltaXMin= self.floodVolume[self.nrZLevels-1]
        y_i= pcr.scalar(1.)
        k= [pcr.scalar(0.)]*2
        mInt= pcr.scalar(0.)
        for iCnt in range(self.nrZLevels-1,0,-1):
            #-find x_i for current volume and update match if applicable
            # also update slope and intercept
            deltaX= channelStorage-self.floodVolume[iCnt]
            mask= pcr.abs(deltaX) < pcr.abs(deltaXMin)
            deltaXMin= pcr.ifthenelse(mask,deltaX,deltaXMin)
            y_i= pcr.ifthenelse(mask,self.areaFractions[iCnt],y_i)
            k[0]= pcr.ifthenelse(mask,self.kSlope[iCnt-1],k[0])
            k[1]= pcr.ifthenelse(mask,self.kSlope[iCnt],k[1])
            mInt= pcr.ifthenelse(mask,self.mInterval[iCnt],mInt)
        #-all values returned, process data: calculate scaled deltaX and smoothed function
        # on the basis of the integrated logistic functions PHI(x) and 1-PHI(x)
        deltaX= deltaXMin
        deltaXScaled= pcr.ifthenelse(deltaX < 0.,pcr.scalar(-1.),1.)*\
            pcr.min(self.criterionKK,pcr.abs(deltaX/pcr.max(1.,mInt)))
        logInt= self.integralLogisticFunction(deltaXScaled)
        #-compute fractional flooded area and flooded depth
        floodedFraction= pcr.ifthenelse(channelStorage > 0.,\
            pcr.ifthenelse(pcr.abs(deltaXScaled) < self.criterionKK,\
            y_i-k[0]*mInt*logInt[0]+k[1]*mInt*logInt[1],\
            y_i+pcr.ifthenelse(deltaX < 0.,k[0],k[1])*deltaX),0.)
        floodedFraction= pcr.max(0.,pcr.min(1.,floodedFraction))
        floodDepth= pcr.ifthenelse(floodedFraction > 0.,channelStorage/(floodedFraction*self.cellArea),0.)
        floodDepth = pcr.min(self.max_water_height, floodDepth)
        #floodedFraction = pcr.ifthen(self.landmask, pcr.cover(floodedFraction , 0.0))
        #floodDepth = pcr.ifthen(self.landmask, pcr.cover(floodDepth , 0.0))
        return floodedFraction, floodDepth

    def integralLogisticFunction(self,x):
        #-returns a tupple of two values holding the integral of the logistic functions
        # of (x) and (-x)
        logInt=pcr.ln(pcr.exp(-x)+1)
        return logInt,x+logInt

    def kinAlphaStatic(self,channelStorage):
        #-given the total water storage in the cell, returns the Q-A relationiceHeatTransferp
        # for the kinematic wave and required parameters using a static floodplain extent
        if self.quality: 
            manIce= pcr.max(self.manningsN,\
              0.0493*pcr.max(0.01,self.channelStorage/(self.dynamicFracWat*self.cellArea))**\
              (-0.23)*self.iceThickness**0.57)
            manningsWithIce= (0.5*(self.manningsN**1.5+manIce**1.5))**(2./3.)
            wetA= self.channelStorage/self.channelLength
            wetP= 2.*wetA/self.wMean+self.wMean
            alphaQ = (manningsWithIce*wetP**(2./3.)*self.channelGradient**-0.5)**self.beta
        else:
            wetA= channelStorage/self.channelLength
            wetP= 2.*wetA/self.wMean+self.wMean
            alphaQ= (self.manningsN*wetP**(2./3.)*self.channelGradient**-0.5)**self.beta	  
        #-returning variable of interest: flooded fraction, cross-sectional area
        # and alphaQ
        dischargeInitial = pcr.ifthenelse(alphaQ > 0.0,(wetA / alphaQ)**(1/self.beta),0.0)    
        return alphaQ, dischargeInitial

    def kinAlphaDynamic(self,channelStorage):
        #-given the total water storage in the cell, returns the Q-A relationiceHeatTransferp
        # for the kinematic wave and required parameters
        floodVol= pcr.max(0,channelStorage-self.channelStorageCapacity)
        floodFrac, floodZ= self.returnFloodedFraction(floodVol)
        channelFraction = pcr.max(0.0, pcr.min(1.0,\
             self.wMean * self.cellLengthFD / (self.cellArea)))
        floodFrac += channelFraction
        #-wetted perimeter, cross-sectional area and
        # corresponding mannings' n
        wetA= channelStorage/self.channelLength
        #-wetted perimeter, alpha and composite manning's n
        wetPFld= pcr.max(0.,floodFrac*self.cellArea/self.channelLength-\
            self.wMean)+2.*floodZ
        wetPCh= self.wMean+\
            2.*pcr.min(self.channelDepth,channelStorage/(self.channelLength*self.wMean))
        wetP= wetPFld+wetPCh
        manQ= (wetPCh/wetP*self.manningsN**1.5+\
            wetPFld/wetP*self.floodplainManN**1.5)**(2./3.)
        alphaQ= (manQ*wetP**(2./3.)*self.channelGradient**-0.5)**self.beta        
        # estimate of channel discharge (m3/s) based on water height
        #
        dischargeInitial = pcr.ifthenelse(alphaQ > 0.0,(wetA / alphaQ)**(1/self.beta),0.0)
        #-returning variables of interest: flooded fraction, cross-sectional area
        # and alphaQ
        return floodFrac,floodZ,alphaQ, dischargeInitial

    def kinAlphaComposite(self,channelStorage):
        #-given the total water storage and the mask specifying the occurrence of
        # floodplain conditions, retrns the Q-A relationiceHeatTransferp for the kinematic
        # wave and the associated parameters
        mask = pcr.boolean(1)
        #floodplainStorage= channelStorage
        floodFrac, floodZ, dynamicAlphaQ, dynamicDischargeInitial= self.kinAlphaDynamic(channelStorage)
        staticAlphaQ, staticDischargeInitial = self.kinAlphaStatic(channelStorage)
        floodFrac= pcr.ifthenelse(mask,floodFrac,0.)
        floodZ= pcr.ifthenelse(mask,floodZ,0.)
        alphaQ= pcr.ifthenelse(mask,dynamicAlphaQ,staticAlphaQ)
        dischargeInitial= pcr.ifthenelse(mask,dynamicDischargeInitial,staticDischargeInitial)
        return floodFrac,floodZ,alphaQ, dischargeInitial

    def readExtensiveMeteo(self, currTimeStep):
        #Read meteorological input directly from netCDF files
        
        if currTimeStep.day == 1:
            self.cloudCover = vos.netcdf2PCRobjClone(\
                                     self.cloudFileNC,'cld',\
                                     str(currTimeStep.fulldate), 
                                     useDoy = "monthly",
                                      cloneMapFileName=self.cloneMap,\
                                      LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)

            self.vaporPressure = vos.netcdf2PCRobjClone(\
                                     self.vapFileNC,'vap',\
                                     str(currTimeStep.fulldate), 
                                     useDoy = "monthly",
                                      cloneMapFileName=self.cloneMap,\
                                      LatitudeLongitude = True,specificFillValue = -999.)
                                      
            self.annualT = vos.netcdf2PCRobjClone(\
                                             self.annualTFileNC,'tas',\
                                             str(currTimeStep.fulldate), 
                                             useDoy = "yearly",
                                              cloneMapFileName=self.cloneMap,\
                                              LatitudeLongitude = True,specificFillValue = -999.) + pcr.scalar(273.15)                                      
                                  
        self.radiation =  vos.netcdf2PCRobjClone(\
                                 self.radFileNC,'rsds',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = "daily",
                                 cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)

        #-vapour pressure, used to return atmospheric emissivity [-]                          
        self.atmosEmis= pcr.scalar(1.0) #pcr.min(1.,(0.53+0.0065*(self.vaporPressure)**0.5)*(1.+0.4*self.cloudCover))
        cld1= pcr.roundoff(10*self.cloudCover+0.5)
        cld0= cld1-1
        sun0= pcr.lookupscalar(self.sunFracTBL, cld0)
        deltaSun= (pcr.lookupscalar(self.sunFracTBL,cld1)-sun0)/(cld1-cld0)
        sunFrac= sun0+(10*self.cloudCover-cld0)*deltaSun
        radFrac= self.radCon+self.radSlope*sunFrac
        self.rsw= radFrac*self.radiation
        
    def readExtensiveHydro(self, currTimeStep):
          
        #Read hydrological input directly from netCDF files (daily timestep)
        self.baseflow = vos.netcdf2PCRobjClone(\
                                  self.baseflowNC,"baseflow",\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        self.baseflow = pcr.ifthen(self.landmask,self.baseflow)
        
        self.interflowTotal = vos.netcdf2PCRobjClone(\
                                        self.interflowNC,"interflow",\
                                        str(currTimeStep.fulldate),
                                        useDoy = None,
                                        cloneMapFileName=self.cloneMap,\
                                        LatitudeLongitude = True)
        self.interflowTotal = pcr.ifthen(self.landmask,self.interflowTotal)
                                          
        self.directRunoff = vos.netcdf2PCRobjClone(\
                                  self.directRunoffNC, "direct_runoff",\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        self.directRunoff = pcr.ifthen(self.landmask,self.directRunoff)
                                               
        self.runoff = self.directRunoff + self.interflowTotal + self.baseflow

    def readPowerplantData(self, currTimeStep):
        
        #read power return flows (from Lohrmann et al., 2019)
        if currTimeStep.day == 1:
            logger.info("Reading powerplant data")
            self.PowRF =  vos.netcdf2PCRobjClone(\
                                     self.PowRFNC,'PowRF',\
                                     str(currTimeStep.fulldate), 
                                     useDoy = None,
                                     cloneMapFileName=self.cloneMap,\
                                      LatitudeLongitude = True,specificFillValue = None)
        
            #Define delta T heat dump
            self.deltaT = pcr.scalar(7.) #difference between effluent temperature and river temperature (K) (van Vliet et al., 2012) 
            self.PowTwload = self.PowRF * self.specificHeatWater * self.densityWater * self.deltaT #W
    
    def readPollutantLoadingsInputData(self, currTimeStep):
        logger.info("Loading input data required to calculate pollutant loadings")    
        
        #Domestic
        self.Population = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'Population',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for gridded population file (5 arc-mins)    
        
        #Manufacturing
        self.factorInd_Man = vos.netcdf2PCRobjCloneWithoutTime(self.factorInd_ManNC,"factorInd_Man",self.cloneMap) #factor for converting ind return flow to man wastewater (-)
                                      
        #Urban surface runoff
        self.urban_area_fraction = vos.netcdf2PCRobjClone(\
                                             self.UrbanFractionNC,'urban_fraction',\
                                             str(currTimeStep.fulldate), 
                                             useDoy = None,
                                             cloneMapFileName=self.cloneMap,\
                                             LatitudeLongitude = True,specificFillValue = None) #fraction urban area (0-1)
        self.urban_area_fraction = pcr.cover(self.urban_area_fraction,0) #if urban fraction missing (e.g. for lakes), make 0
        
        #Livestock
        self.BufalloPopulation = vos.netcdf2PCRobjClone(\
                                 self.LivPopulationNC,'BufalloPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for bufallo population (5 arc-mins)
        self.BufalloPopulation = pcr.cover(self.BufalloPopulation,0.)
        
        self.ChickenPopulation = vos.netcdf2PCRobjClone(\
                                 self.LivPopulationNC,'ChickenPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for chicken population (5 arc-mins)
        self.ChickenPopulation = pcr.cover(self.ChickenPopulation,0.)
                                          
        self.CowPopulation = vos.netcdf2PCRobjClone(\
                                 self.LivPopulationNC,'CowPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for cow population (5 arc-mins)
        self.CowPopulation = pcr.cover(self.CowPopulation,0.)
        
        self.DuckPopulation = vos.netcdf2PCRobjClone(\
                                 self.LivPopulationNC,'DuckPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for duck population (5 arc-mins)                 
        self.DuckPopulation = pcr.cover(self.DuckPopulation,0.)
                  
        self.GoatPopulation = vos.netcdf2PCRobjClone(\
                                 self.LivPopulationNC,'GoatPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for goat population (5 arc-mins)  
        self.GoatPopulation = pcr.cover(self.GoatPopulation,0.)
        
        self.HorsePopulation = vos.netcdf2PCRobjClone(\
                                 self.LivPopulationNC,'HorsePop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for horse population (5 arc-mins)  
        self.HorsePopulation = pcr.cover(self.HorsePopulation,0.)
                                          
        self.PigPopulation = vos.netcdf2PCRobjClone(\
                                 self.LivPopulationNC,'PigPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for pig population (5 arc-mins)  
        self.PigPopulation = pcr.cover(self.PigPopulation,0.)
        
        self.SheepPopulation = vos.netcdf2PCRobjClone(\
                                 self.LivPopulationNC,'SheepPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #input pathway for sheep population (5 arc-mins)            
        self.SheepPopulation = pcr.cover(self.SheepPopulation,0.)
        
        #Calculate livestock densities accounting for livestock units (Wen et al., 2018)
        self.cellArea_km2 = self.cellArea/1000000. #convert m2 to km2
        self.LivDensityThres = pcr.scalar(25.)
        
        self.BufalloDensity = self.BufalloPopulation/ self.cellArea_km2
        self.ChickenDensity = (self.ChickenPopulation * 0.01) /self.cellArea_km2
        self.CowDensity = self.CowPopulation/ self.cellArea_km2
        self.DuckDensity = (self.DuckPopulation * 0.01)/ self.cellArea_km2
        self.GoatDensity = (self.GoatPopulation * 0.1)/ self.cellArea_km2
        self.HorseDensity = self.HorsePopulation/ self.cellArea_km2
        self.PigDensity = (self.PigPopulation * 0.3)/ self.cellArea_km2
        self.SheepDensity = (self.SheepPopulation * 0.1)/ self.cellArea_km2
                                  
        #Wastewater Pathways
        self.ratio_WWt_Tertiary = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'WWt_Tertiary',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True, specificFillValue = None) #ratio of tertiary treatment (5 arc-min)
        self.ratio_WWt_Tertiary = pcr.cover(self.ratio_WWt_Tertiary,0)
                                  
        self.ratio_WWt_Secondary = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'WWt_Secondary',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #ratio of secondary treatment (5 arc-min)
        self.ratio_WWt_Secondary = pcr.cover(self.ratio_WWt_Secondary,0)
                                          
        self.ratio_WWt_Primary = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'WWt_Primary',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #ratio of primary treatment (5 arc-min)
        self.ratio_WWt_Primary = pcr.cover(self.ratio_WWt_Primary,0)
                                          
        self.ratio_WWcut = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'WWcut',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #fraction of wastewater collected, but untreated (5 arc-min)
        self.ratio_WWcut = pcr.cover(self.ratio_WWcut,0) 
                                          
        self.ratio_dom_WWbs = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'dom_WWbs',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #fraction of domestic wastewater via basic sanitation (5 arc-min)
        self.ratio_dom_WWbs = pcr.cover(self.ratio_dom_WWbs,0)
                                          
        self.ratio_dom_WWod = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'dom_WWod',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #fraction of domestic wastewater via open defecation (5 arc-min)
        self.ratio_dom_WWod = pcr.cover(self.ratio_dom_WWod,0) 
                                  
        self.ratio_man_WWdirect = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'man_WWdirect',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = None) #fraction of manufacturing wastewater directly discharged (5 arc-min)
        self.ratio_man_WWdirect = pcr.cover(self.ratio_man_WWdirect,0)

    def calculatePollutantLoadings(self, currTimeStep, landSurface, groundwater):
        #calculate pollutant loadings directly
        logger.info("Calculating pollutant loadings")
        
        #Municipal wastewater treatment     
        TDS_RemEff = (self.ratio_WWt_Tertiary * self.TDS_Ter_RemEff) +\
                     (self.ratio_WWt_Secondary * self.TDS_Sec_RemEff) +\
                     (self.ratio_WWt_Primary * self.TDS_Pri_RemEff) +\
                     (self.ratio_WWcut * self.TDS_WWcut_RemEff)
        BOD_RemEff = (self.ratio_WWt_Tertiary * self.BOD_Ter_RemEff) +\
                     (self.ratio_WWt_Secondary * self.BOD_Sec_RemEff) +\
                     (self.ratio_WWt_Primary * self.BOD_Pri_RemEff) +\
                     (self.ratio_WWcut * self.BOD_WWcut_RemEff)
        FC_RemEff =  (self.ratio_WWt_Tertiary * self.FC_Ter_RemEff) +\
                     (self.ratio_WWt_Secondary * self.FC_Sec_RemEff) +\
                     (self.ratio_WWt_Primary * self.FC_Pri_RemEff) +\
                     (self.ratio_WWcut * self.FC_WWcut_RemEff)
        
        #- Fraction of direct runoff directly from PCR-GLOBWB (direct runoff / total runoff)           
        self.frac_surfaceRunoff = vos.getValDivZero(landSurface.directRunoff, (landSurface.landSurfaceRunoff + groundwater.baseflow))
        
        #Domestic sector-specific wastewater pathways        
        TDS_Rem_BasicSan = self.ratio_dom_WWbs * self.TDS_dom_WWbs_RemEff     #Removal of TDS due to basic sanitation practices
        BOD_Rem_BasicSan = self.ratio_dom_WWbs * self.BOD_dom_WWbs_RemEff     #Removal of BOD due to basic sanitation practices
        FC_Rem_BasicSan = self.ratio_dom_WWbs * self.FC_dom_WWbs_RemEff       #Removal of FC due to basic sanitation practices
        Rem_OpenDef = self.ratio_dom_WWod * self.frac_surfaceRunoff           #Pollutants transferred to surface water proportional to fraction surface runoff
        
        #Domestic loadings: Gridded population (capita) * per Capita excretion rate [g/capita/day; cfu/capita/day]
        self.Dom_TDSload = self.Population * self.DomTDS_ExcrLoad * (1 - (TDS_RemEff + TDS_Rem_BasicSan + Rem_OpenDef)) #g/day
        self.Dom_BODload = self.Population * self.DomBOD_ExcrLoad * (1 - (BOD_RemEff + BOD_Rem_BasicSan + Rem_OpenDef)) #g/day
        self.Dom_FCload = (self.Population * self.DomFC_ExcrLoad * (1 - (FC_RemEff + FC_Rem_BasicSan + Rem_OpenDef))) /1000000. #million cfu/day
        
        #Manufacturing loadings: Manufacturing wastewater [m3/day] * average manufacturing effluent concentration [mg/L; cfu/100ml]
        self.ManWWp = self.IndustryReturnFlowVol * self.factorInd_Man #multipled by factor for estimated manufacturing wastewater (based on Jones et al., 2021 data)
        self.Man_TDSload = self.ManWWp * self.ManTDS_EfflConc * (1 - TDS_RemEff) #g/day
        self.Man_BODload = self.ManWWp * self.ManBOD_EfflConc * (1 - BOD_RemEff) #g/day
        self.Man_FCload = (self.ManWWp * self.ManFC_EfflConc * (1 - FC_RemEff))/ 100. #million cfu/day
        
        #Urban surface runoff loadings: USR return flow [m3/day] * average USR effluent concentration [mg/L; cfu/100ml]  
        self.USR_RF = landSurface.directRunoff * self.urban_area_fraction * self.cellArea #m3 day
        
        self.USR_TDSload = self.USR_RF * self.USRTDS_EfflConc * (1 - TDS_RemEff) #g/day
        self.USR_BODload = self.USR_RF * self.USRBOD_EfflConc * (1 - BOD_RemEff) #g/day
        self.USR_FCload = (self.USR_RF * self.USRFC_EfflConc * (1 - FC_RemEff))/ 100. #million cfu/day
                
        #Livestock loadings
        BOD_Rem_IntLiv = (self.ratio_WWt_Tertiary + self.ratio_WWt_Secondary) * self.BOD_Sec_RemEff
        FC_Rem_IntLiv = (self.ratio_WWt_Tertiary + self.ratio_WWt_Secondary) * self.FC_Sec_RemEff
        
        self.intLiv_Bufallo_BODload = pcr.ifthenelse(self.BufalloDensity > self.LivDensityThres,\
                                      self.BufalloPopulation * self.Bufallo_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Chicken_BODload = pcr.ifthenelse(self.ChickenDensity > self.LivDensityThres,\
                                      self.ChickenPopulation * self.Chicken_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Cow_BODload     = pcr.ifthenelse(self.CowDensity > self.LivDensityThres,\
                                      self.CowPopulation * self.Cow_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Duck_BODload    = pcr.ifthenelse(self.DuckDensity > self.LivDensityThres,\
                                      self.DuckPopulation * self.Duck_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Goat_BODload    = pcr.ifthenelse(self.GoatDensity > self.LivDensityThres,\
                                      self.GoatPopulation * self.Goat_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Horse_BODload   = pcr.ifthenelse(self.HorseDensity > self.LivDensityThres,\
                                      self.HorsePopulation * self.Horse_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Pig_BODload     = pcr.ifthenelse(self.PigDensity > self.LivDensityThres,\
                                      self.PigPopulation * self.Pig_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Sheep_BODload   = pcr.ifthenelse(self.SheepDensity > self.LivDensityThres,\
                                      self.SheepPopulation * self.Sheep_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        
        self.extLiv_Bufallo_BODload = pcr.ifthenelse(self.BufalloDensity <= self.LivDensityThres, self.BufalloPopulation * self.Bufallo_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Chicken_BODload = pcr.ifthenelse(self.ChickenDensity <= self.LivDensityThres, self.ChickenPopulation * self.Chicken_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Cow_BODload = pcr.ifthenelse(self.CowDensity <= self.LivDensityThres, self.CowPopulation * self.Cow_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Duck_BODload = pcr.ifthenelse(self.DuckDensity <= self.LivDensityThres, self.DuckPopulation * self.Duck_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Goat_BODload = pcr.ifthenelse(self.GoatDensity <= self.LivDensityThres, self.GoatPopulation * self.Goat_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Horse_BODload = pcr.ifthenelse(self.HorseDensity <= self.LivDensityThres, self.HorsePopulation * self.Horse_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Pig_BODload = pcr.ifthenelse(self.PigDensity <= self.LivDensityThres, self.PigPopulation * self.Pig_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Sheep_BODload = pcr.ifthenelse(self.SheepDensity <= self.LivDensityThres, self.SheepPopulation * self.Sheep_BODload * self.frac_surfaceRunoff, 0.0)
        
        self.intLiv_BODload = self.intLiv_Bufallo_BODload + self.intLiv_Chicken_BODload + self.intLiv_Cow_BODload +\
                              self.intLiv_Duck_BODload + self.intLiv_Goat_BODload + self.intLiv_Horse_BODload + self.intLiv_Pig_BODload + self.intLiv_Sheep_BODload #g/day        
                              
        self.extLiv_BODload = self.extLiv_Bufallo_BODload + self.extLiv_Chicken_BODload + self.extLiv_Cow_BODload +\
                              self.extLiv_Duck_BODload + self.extLiv_Goat_BODload + self.extLiv_Horse_BODload + self.extLiv_Pig_BODload + self.extLiv_Sheep_BODload #g/day

        self.intLiv_Bufallo_FCload = pcr.ifthenelse(self.BufalloDensity > self.LivDensityThres,\
                                     self.BufalloPopulation * self.Bufallo_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Chicken_FCload = pcr.ifthenelse(self.ChickenDensity > self.LivDensityThres,\
                                     self.ChickenPopulation * self.Chicken_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Cow_FCload     = pcr.ifthenelse(self.CowDensity > self.LivDensityThres,\
                                     self.CowPopulation * self.Cow_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Duck_FCload    = pcr.ifthenelse(self.DuckDensity > self.LivDensityThres,\
                                     self.DuckPopulation * self.Duck_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Goat_FCload    = pcr.ifthenelse(self.GoatDensity > self.LivDensityThres,\
                                     self.GoatPopulation * self.Goat_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Horse_FCload   = pcr.ifthenelse(self.HorseDensity > self.LivDensityThres,\
                                     self.HorsePopulation * self.Horse_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Pig_FCload     = pcr.ifthenelse(self.PigDensity > self.LivDensityThres,\
                                     self.PigPopulation * self.Pig_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        self.intLiv_Sheep_FCload   = pcr.ifthenelse(self.SheepDensity > self.LivDensityThres,\
                                     self.SheepPopulation * self.Sheep_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
        
        self.extLiv_Bufallo_FCload = pcr.ifthenelse(self.BufalloDensity <= self.LivDensityThres, self.BufalloPopulation * self.Bufallo_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Chicken_FCload = pcr.ifthenelse(self.ChickenDensity <= self.LivDensityThres, self.ChickenPopulation * self.Chicken_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Cow_FCload = pcr.ifthenelse(self.CowDensity <= self.LivDensityThres, self.CowPopulation * self.Cow_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Duck_FCload = pcr.ifthenelse(self.DuckDensity <= self.LivDensityThres, self.DuckPopulation * self.Duck_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Goat_FCload = pcr.ifthenelse(self.GoatDensity <= self.LivDensityThres, self.GoatPopulation * self.Goat_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Horse_FCload = pcr.ifthenelse(self.HorseDensity <= self.LivDensityThres, self.HorsePopulation * self.Horse_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Pig_FCload = pcr.ifthenelse(self.PigDensity <= self.LivDensityThres, self.PigPopulation * self.Pig_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Sheep_FCload = pcr.ifthenelse(self.SheepDensity <= self.LivDensityThres, self.SheepPopulation * self.Sheep_FCload * self.frac_surfaceRunoff, 0.0)
        
        self.intLiv_FCload = (self.intLiv_Bufallo_FCload + self.intLiv_Chicken_FCload + self.intLiv_Cow_FCload +\
                              self.intLiv_Duck_FCload + self.intLiv_Goat_FCload + self.intLiv_Horse_FCload + self.intLiv_Pig_FCload + self.intLiv_Sheep_FCload)\
                              / 1000000. #million cfu/day
        self.extLiv_FCload = (self.extLiv_Bufallo_FCload + self.extLiv_Chicken_FCload + self.extLiv_Cow_FCload +\
                              self.extLiv_Duck_FCload + self.extLiv_Goat_FCload + self.extLiv_Horse_FCload + self.extLiv_Pig_FCload + self.extLiv_Sheep_FCload)\
                              / 1000000. #million cfu/day
        
        #Irrigation loadings: Irrigation return flow (m3/day) * soil TDS (mg/L)
        #- Irrigation return flows calculated directly from PCR-GLOBWB (from direct runoff, interflow and baseflow)
        self.irr_rf_from_direct_runoff = landSurface.directRunoff * vos.getValDivZero(self.avg_irrGrossDemand, (self.avg_irrGrossDemand + self.avg_netLqWaterToSoil))
        self.irr_rf_from_interflow     = landSurface.interflowTotal * vos.getValDivZero(self.avg_irrGrossDemand, (self.avg_irrGrossDemand + self.avg_netLqWaterToSoil))
        self.irr_rf_from_baseflow     =  groundwater.baseflow * vos.getValDivZero(self.avg_irrGrossDemand, (self.avg_irrGrossDemand + self.avg_netLqWaterToSoil))

        self.Irr_RF = (self.irr_rf_from_direct_runoff +\
                       self.irr_rf_from_interflow +\
                       self.irr_rf_from_baseflow) * self.cellArea  # - total irirgation return flow in volume units (m3 day-1) 
        self.Irr_TDSload = self.Irr_RF * self.IrrTDS_EfflConc #g/day
        
        #Combined loadings
        self.TDSload = (self.Dom_TDSload + self.Man_TDSload + self.USR_TDSload + self.Irr_TDSload) #g day-1
        self.BODload = (self.Dom_BODload + self.Man_BODload + self.USR_BODload + self.intLiv_BODload + self.extLiv_BODload) #g day-1
        self.FCload = (self.Dom_FCload + self.Man_FCload + self.USR_FCload + self.intLiv_FCload + self.extLiv_FCload) #million cfu day-1
        
    def readPollutantLoadings(self, currTimeStep):
        #Use pre-calculated pollutant loadings
        logger.info("Reading loadings directly")        
        
        # read TDS loadings (combined for all sectoral activities)
        self.TDSload = vos.netcdf2PCRobjClone(\
                                 self.TDSloadNC,'TDSload',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                 cloneMapFileName=self.cloneMap,\
                                 LatitudeLongitude = True)
        self.TDSload = pcr.ifthen(self.landmask, self.TDSload)

        #read BOD loadings (combined for all sectoral activities)
        self.BODload = vos.netcdf2PCRobjClone(\
                                 self.BODloadNC,'BODload',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                 cloneMapFileName=self.cloneMap,\
                                 LatitudeLongitude = True)							 
        self.BODload = pcr.ifthen(self.landmask, self.BODload)
        
        #read FC loadings (combined for all sectoral activities)
        self.FCload = vos.netcdf2PCRobjClone(\
                                 self.FCloadNC,'FCload',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                 cloneMapFileName=self.cloneMap,\
                                 LatitudeLongitude = True)                                   
        self.FCload = pcr.ifthen(self.landmask, self.FCload)
                
    def energyLocal(self, meteo, landSurface, groundwater, timeSec = vos.secondsPerDay()):
        #-surface water energy fluxes [W/m2]
        # within the current time step
        #-ice formation evaluated prior to routing to account for loss
        # in water height, vertical change in energy evaluated,
        # warming capped to increase to air temperature
        # noIce:    boolean variable indicating the absence of ice,
        #           false when ice is present or when the flux to the
        #           ice layer is negative, indicating growth
        # SHI:      surface energy flux (heat transfer phi) of ice (+: melt)
        # SHW:      heat transfer to surface water
        # SHR:      heat transfer due to short and longwave radiation
        # SHA:      advected energy due to rain or snow
        # SHQ:      advected energy due to lateral inflow
        # SHL:      latent heat flux, based on actual open water evaporation
        #           to be evaluated when actual evap is known
        
        self.temperatureKelvin = meteo.temperature + pcr.scalar(273.15)
        landRunoff = self.runoff
        self.correctPrecip = pcr.scalar(0.0)
        self.dynamicFracWatBeforeRouting = self.dynamicFracWat
        
        landT= pcr.cover(landSurface.directRunoff/landRunoff*pcr.max(self.iceThresTemp+0.1,self.temperatureKelvin-self.deltaTPrec)+\
           landSurface.interflowTotal/landRunoff*pcr.max(self.iceThresTemp+0.1,self.temperatureKelvin)+\
           groundwater.baseflow/landRunoff*pcr.max(self.iceThresTemp+5.0,self.annualT),self.temperatureKelvin)
  
        iceHeatTransfer = self.heatTransferWater * (self.temperatureKelvin - self.iceThresTemp)
        waterHeatTransfer = self.heatTransferIce * (self.iceThresTemp - self.waterTemp)
        noIce = pcr.ifthenelse(self.iceThickness > 0,pcr.boolean(0),\
          pcr.ifthenelse(((iceHeatTransfer-waterHeatTransfer) < 0) & (self.temperatureKelvin < self.iceThresTemp),\
          pcr.boolean(0),pcr.boolean(1)))
        waterHeatTransfer = pcr.ifthenelse(noIce, self.heatTransferWater * (self.temperatureKelvin - self.waterTemp), waterHeatTransfer)
        radiativHeatTransfer = (1 - pcr.ifthenelse(noIce, self.albedoWater, self.albedoSnow)) * self.rsw
        radiativHeatTransfer = radiativHeatTransfer - self.stefanBoltzman * (pcr.ifthenelse(noIce,\
          self.waterTemp, self.iceThresTemp)**4 - self.atmosEmis * self.temperatureKelvin**4)
        advectedEnergyPrecip = pcr.max(0,self.correctPrecip) *\
          pcr.max(self.iceThresTemp+0.1,self.temperatureKelvin-self.deltaTPrec)*self.specificHeatWater*self.densityWater/timeSec
        advectedEnergyInflow = (1-self.dynamicFracWat)/self.dynamicFracWat * landRunoff * landT * self.specificHeatWater * self.densityWater / timeSec
        advectedEnergyInflow -= landSurface.actSurfaceWaterAbstract/self.dynamicFracWat * self.waterTemp * self.specificHeatWater * self.densityWater / timeSec
        #-ice formation
        # DSHI:     net flux for ice layer [W/m2]
        # wi:       thickness of ice cover [m]
        # wh:       available water height
        # deltaIceThickness:      change in thickness per day, melt negative
        #diceHeatTransfer= pcr.ifthenelse(noIce,0,iceHeatTransfer-waterHeatTransfer+advectedEnergyPrecip+radiativHeatTransfer)
        diceHeatTransfer= pcr.ifthenelse(noIce,0,iceHeatTransfer-waterHeatTransfer+radiativHeatTransfer)
        self.deltaIceThickness = -diceHeatTransfer * timeSec / (self.densityWater * self.latentHeatFusion)
        self.deltaIceThickness = pcr.max(-self.iceThickness,self.deltaIceThickness)
        self.deltaIceThickness = pcr.min(self.deltaIceThickness, pcr.max(0,self.maxIceThickness-self.iceThickness))
        #-returning direct gain over water surface
        watQ= pcr.ifthenelse(self.temperatureKelvin >= self.iceThresTemp,pcr.max(0,self.correctPrecip) - self.waterBodyEvaporation/self.dynamicFracWat,0)
        #watQ= pcr.ifthenelse(self.temperatureKelvin >= self.iceThresTemp,self.waterBodyEvaporation/self.dynamicFracWat,0)
        self.waterBodyEvaporation= pcr.ifthenelse(self.temperatureKelvin >= self.iceThresTemp,self.waterBodyEvaporation,0) # TODO move

        #-returning vertical gains/losses [m3] over lakes and channels
        # given their corresponding area
        deltaIceThickness_melt= pcr.max(0,-self.deltaIceThickness)
        #landQCont= pcr.max(0,landFrac/watFrac*landQ)
        verticalGain = (watQ+(landRunoff-landSurface.actSurfaceWaterAbstract)/(self.dynamicFracWat)+deltaIceThickness_melt)
        #channeldStor = verticalGain ### only ice in river
        #lakedStor = self.dynamicFracWat * verticalGain

        #self.waterBodies.waterBodyStorage = pcr.max(0.0, self.waterBodies.waterBodyStorage - lakedStor)
        #self.channelStorage = pcr.max(0.0, self.channelStorage - channeldStor)

        #-net cumulative input for mass balance check [m3]
        #self.dtotStor= channeldStor
        #-change in water storage due to vertical change only
        # used to limit heating and cooling of surface water
        dtotStorLoc= verticalGain
        totStorLoc = self.channelStorageTimeBefore/(self.dynamicFracWat*self.cellArea)-dtotStorLoc
        totStorLoc = self.return_water_body_storage_to_channel(self.channelStorageTimeBefore) / (self.dynamicFracWat*self.cellArea)-dtotStorLoc
        self.totEW = totStorLoc*self.waterTemp*self.specificHeatWater*self.densityWater
        dtotStorLoc = pcr.max(-totStorLoc,dtotStorLoc)
        
        #-latent heat flux due to evapotranspiration [W/m2]
        # and advected energy due to ice melt included here
        # to account for correction in water storage
        latentHeat= -self.waterBodyEvaporation/self.dynamicFracWat*self.densityWater*self.latentHeatVapor/timeSec
        advectedEnergyPrecip= advectedEnergyPrecip+deltaIceThickness_melt*self.iceThresTemp*self.specificHeatWater*self.densityWater/timeSec

        totEWC= totStorLoc * self.specificHeatWater * self.densityWater
        dtotEWC= dtotStorLoc * self.specificHeatWater * self.densityWater
        dtotEWLoc= (waterHeatTransfer + pcr.scalar(noIce) * (radiativHeatTransfer+latentHeat))*timeSec
        dtotEWAdv= (advectedEnergyInflow + pcr.scalar(noIce) * advectedEnergyPrecip)*timeSec
        dtotEWLoc= pcr.min(dtotEWLoc,\
          pcr.max(0,totEWC*self.temperatureKelvin-self.totEW)+pcr.ifthenelse(dtotStorLoc > 0,\
          pcr.max(0,dtotEWC*self.temperatureKelvin-dtotEWAdv),0))
        dtotEWLoc= pcr.ifthenelse(self.waterTemp > self.temperatureKelvin,\
          pcr.min(0,dtotEWLoc),dtotEWLoc)
        dtotEWLoc= pcr.max(dtotEWLoc,\
          pcr.min(0,(totEWC+dtotEWC)*pcr.max(self.temperatureKelvin,self.iceThresTemp+.1)-(self.totEW+dtotEWAdv)))
        self.temp_water_height = pcr.max(1e-16,totStorLoc + dtotStorLoc)
        
        #Calculate heat dump from powerplants
        self.maxTwload = pcr.cover(totStorLoc*self.deltaT*self.specificHeatWater*self.densityWater,0.) #calculate theoretical EW from powerplants to warm gridcell by 7K
        self.Twload = pcr.min(self.PowTwload, self.maxTwload)
        
        ###original DynWat code: -change in energy storage and resulting temperature
        self.totEW = pcr.max(0,self.totEW+dtotEWLoc+dtotEWAdv+self.Twload)
        self.waterTemp= pcr.ifthenelse(self.temp_water_height > self.critical_water_height,\
          self.totEW/self.temp_water_height/(self.specificHeatWater*self.densityWater),self.temperatureKelvin)        
        #self.waterTemp= pcr.ifthenelse(self.waterTemp < self.iceThresTemp+0.1,self.iceThresTemp+0.1,self.waterTemp)
        self.waterTemp= min(pcr.ifthenelse(self.waterTemp < self.iceThresTemp+0.1,self.iceThresTemp+0.1,self.waterTemp), self.maxThresTemp)               

    def energyLocal_routing_only(self, meteo, timeSec = vos.secondsPerDay()):
        #-surface water energy fluxes [W/m2]
        # within the current time step
        #-ice formation evaluated prior to routing to account for loss
        # in water height, vertical change in energy evaluated,
        # warming capped to increase to air temperature
        # noIce:    boolean variable indicating the absence of ice,
        #           false when ice is present or when the flux to the
        #           ice layer is negative, indicating growth
        # SHI:      surface energy flux (heat transfer phi) of ice (+: melt)
        # SHW:      heat transfer to surface water
        # SHR:      heat transfer due to short and longwave radiation
        # SHA:      advected energy due to rain or snow
        # SHQ:      advected energy due to lateral inflow
        # SHL:      latent heat flux, based on actual open water evaporation
        #           to be evaluated when actual evap is known
        
        self.temperatureKelvin = meteo.temperature + pcr.scalar(273.15)
        landRunoff = self.runoff
        self.correctPrecip = pcr.scalar(0.0)
        self.dynamicFracWatBeforeRouting = self.dynamicFracWat
        self.dynamicFracWat = pcr.max(self.dynamicFracWat,0.001)
        landT= pcr.cover(self.directRunoff/landRunoff*pcr.max(self.iceThresTemp+0.1,self.temperatureKelvin-self.deltaTPrec)+\
           self.interflowTotal/landRunoff*pcr.max(self.iceThresTemp+0.1,self.temperatureKelvin)+\
           self.baseflow/landRunoff*pcr.max(self.iceThresTemp+5.0,self.annualT),self.temperatureKelvin)
        iceHeatTransfer = self.heatTransferWater * (self.temperatureKelvin - self.iceThresTemp)
        waterHeatTransfer = self.heatTransferIce * (self.iceThresTemp - self.waterTemp)
        noIce = pcr.ifthenelse(self.iceThickness > 0,pcr.boolean(0),\
          pcr.ifthenelse(((iceHeatTransfer-waterHeatTransfer) < 0) & (self.temperatureKelvin < self.iceThresTemp),\
          pcr.boolean(0),pcr.boolean(1)))
        waterHeatTransfer = pcr.ifthenelse(noIce, self.heatTransferWater * (self.temperatureKelvin - self.waterTemp), waterHeatTransfer)
        radiativHeatTransfer = (1 - pcr.ifthenelse(noIce, self.albedoWater, self.albedoSnow)) * self.rsw
        radiativHeatTransfer = radiativHeatTransfer - self.stefanBoltzman * (pcr.ifthenelse(noIce,\
          self.waterTemp, self.iceThresTemp)**4 - self.atmosEmis * self.temperatureKelvin**4)
        advectedEnergyPrecip = pcr.max(0,self.correctPrecip) *\
          pcr.max(self.iceThresTemp+0.1,self.temperatureKelvin-self.deltaTPrec)*self.specificHeatWater*self.densityWater/timeSec
        advectedEnergyInflow = (1-self.dynamicFracWat)/self.dynamicFracWat * landRunoff * landT * self.specificHeatWater * self.densityWater / timeSec
        #-ice formation
        # DSHI:     net flux for ice layer [W/m2]
        # wi:       thickness of ice cover [m]
        # wh:       available water height
        # deltaIceThickness:      change in thickness per day, melt negative
        #diceHeatTransfer= pcr.ifthenelse(noIce,0,iceHeatTransfer-waterHeatTransfer+advectedEnergyPrecip+radiativHeatTransfer)
        diceHeatTransfer= pcr.ifthenelse(noIce,0,iceHeatTransfer-waterHeatTransfer+radiativHeatTransfer)
        self.deltaIceThickness = -diceHeatTransfer * timeSec / (self.densityWater * self.latentHeatFusion)
        self.deltaIceThickness = pcr.max(-self.iceThickness,self.deltaIceThickness)
        self.deltaIceThickness = pcr.min(self.deltaIceThickness, pcr.max(0,self.maxIceThickness-self.iceThickness))
        #-returning direct gain over water surface
        watQ= pcr.ifthenelse(self.temperatureKelvin >= self.iceThresTemp,pcr.max(0,self.correctPrecip) - pcr.cover(self.waterBodyEvaporation/self.dynamicFracWat,0),0)
        #watQ= pcr.ifthenelse(self.temperatureKelvin >= self.iceThresTemp,self.waterBodyEvaporation/self.dynamicFracWat,0)
        self.waterBodyEvaporation= pcr.ifthenelse(self.temperatureKelvin >= self.iceThresTemp,self.waterBodyEvaporation,0) # TODO move
        #-returning vertical gains/losses [m3] over lakes and channels
        # given their corresponding area
        deltaIceThickness_melt= pcr.max(0,-self.deltaIceThickness)
        #landQCont= pcr.max(0,landFrac/watFrac*landQ)
        verticalGain = (watQ+(landRunoff)/(self.dynamicFracWat)+deltaIceThickness_melt)
        #channeldStor = verticalGain ### only ice in river
        #lakedStor = self.dynamicFracWat * verticalGain

        #self.waterBodies.waterBodyStorage = pcr.max(0.0, self.waterBodies.waterBodyStorage - lakedStor)
        #self.channelStorage = pcr.max(0.0, self.channelStorage - channeldStor)

        #-net cumulative input for mass balance check [m3]
        #self.dtotStor= channeldStor
        #-change in water storage due to vertical change only
        # used to limit heating and cooling of surface water
        dtotStorLoc= verticalGain
        #totStorLoc = self.channelStorageTimeBefore/(self.dynamicFracWat*self.cellArea)-dtotStorLoc
        totStorLoc = self.return_water_body_storage_to_channel(self.channelStorageTimeBefore) / (self.dynamicFracWat*self.cellArea)-dtotStorLoc
        self.totEW = (totStorLoc*self.waterTemp*self.specificHeatWater*self.densityWater)
        dtotStorLoc = pcr.max(-totStorLoc,dtotStorLoc)
        
        #-latent heat flux due to evapotranspiration [W/m2]
        # and advected energy due to ice melt included here
        # to account for correction in water storage
        latentHeat= -self.waterBodyEvaporation/self.dynamicFracWat*self.densityWater*self.latentHeatVapor/timeSec
        advectedEnergyPrecip= advectedEnergyPrecip+deltaIceThickness_melt*self.iceThresTemp*self.specificHeatWater*self.densityWater/timeSec

        totEWC= totStorLoc * self.specificHeatWater * self.densityWater
        dtotEWC= dtotStorLoc * self.specificHeatWater * self.densityWater
        dtotEWLoc= (waterHeatTransfer + pcr.scalar(noIce) * (radiativHeatTransfer+latentHeat))*timeSec
        dtotEWAdv= (advectedEnergyInflow + pcr.scalar(noIce) * advectedEnergyPrecip)*timeSec
        dtotEWLoc= pcr.min(dtotEWLoc,\
          pcr.max(0,totEWC*self.temperatureKelvin-self.totEW)+pcr.ifthenelse(dtotStorLoc > 0,\
          pcr.max(0,dtotEWC*self.temperatureKelvin-dtotEWAdv),0))
        dtotEWLoc= pcr.ifthenelse(self.waterTemp > self.temperatureKelvin,\
          pcr.min(0,dtotEWLoc),dtotEWLoc)
        dtotEWLoc= pcr.max(dtotEWLoc,\
          pcr.min(0,(totEWC+dtotEWC)*pcr.max(self.temperatureKelvin,self.iceThresTemp+.1)-(self.totEW+dtotEWAdv)))

        #Calculate heat dump from powerplants
        self.maxTwload = totStorLoc*self.deltaT*self.specificHeatWater*self.densityWater #calculate theoretical EW from powerplants to warm gridcell by 7K
        self.Twload = pcr.min(self.PowTwload, self.maxTwload) 
        
        ###original DynWat code: -change in energy storage and resulting temperature
        self.totEW = pcr.max(0,self.totEW+dtotEWLoc+dtotEWAdv+self.Twload)
        self.waterTemp= pcr.ifthenelse(self.temp_water_height > self.critical_water_height,\
          self.totEW/self.temp_water_height/(self.specificHeatWater*self.densityWater),self.temperatureKelvin)        
        #self.waterTemp= pcr.ifthenelse(self.waterTemp < self.iceThresTemp+0.1,self.iceThresTemp+0.1,self.waterTemp)
        self.waterTemp= min(pcr.ifthenelse(self.waterTemp < self.iceThresTemp+0.1,self.iceThresTemp+0.1,self.waterTemp), self.maxThresTemp)  


    def qualityLocal(self, timeSec = vos.secondsPerDay()): 
        
        self.travel_time = 1 #model setup for daily timestep
        
        ###Salinity load (conservative substances approach)
        #---TDS routing
        self.routedTDS = self.routedTDS + self.TDSload 
        
        if self.loadsPerSector:
            self.routedDomTDS = self.routedDomTDS + self.Dom_TDSload
            self.routedManTDS = self.routedManTDS + self.Man_TDSload
            self.routedUSRTDS = self.routedUSRTDS + self.USR_TDSload
            self.routedIrrTDS = self.routedIrrTDS + self.Irr_TDSload
        
        ###Organic load (non-conservative, temperature dependent decay)
        self.routedBOD = self.routedBOD + self.BODload #get BOD load before decay
        
        #---Temperature dependent decay
        self.k_BOD = pcr.scalar(0.35)     #first-order degradation coefficient at 20C (van Vliet et al., 2021)
        self.watertempcorrection_BOD = pcr.scalar(1.047)     #temperature correction (van Vliet et al., 2021; Wen et al., 2017)
        self.waterTemp_BOD = self.waterTemp - pcr.scalar(273.15)   #water temperature in gridcell in C
        self.BODdecay_temperature = cover(self.k_BOD*(self.watertempcorrection_BOD**(self.waterTemp_BOD - 20)),0.0)
        
        #---Streeter-Phelps (Dissolved oxygen)
        self.BOD_conc = pcr.ifthenelse(self.channelStorage > 0, self.routedBOD / self.channelStorage, 0)  # BOD concentration in mg/l
        self.k1_BOD = self.BODdecay_temperature * self.BOD_conc
        self.DOsat = (1-0.0001148*self.elevation)*exp(-139.34411+(157570.1)/(self.waterTemp)-(66423080)/(self.waterTemp**2)+(12438000000)/(self.waterTemp**3)-(862194900000)/(self.waterTemp**4)) # oxygen saturation in mg/l
        self.velocity = self.avgDischarge / (self.yMean * self.wMean) # velocity assuming rectangular channel (m/s)
        self.k2 = 3.93 * ((self.velocity**0.5) / (self.yMean**1.5)) # reaeration rate in /d (O'Connor and Dobbins, 1958)
        self.k2 = pcr.ifthenelse(self.k2 > 1.5, 1.5, self.k2)
        self.k2 = pcr.ifthenelse(self.k2 < 0.4, 0.4, self.k2)
        self.DO = self.DO - self.k1_BOD + self.k2*(self.DOsat - self.DO) # DO concentration in mg/l
        self.DO = pcr.ifthenelse(self.DO < 0.0, 0.0, self.DO) # non-negative DO
        
        #---BOD routing
        self.BODdecay = exp(-self.BODdecay_temperature * self.travel_time) # /d
        self.routedBOD = self.routedBOD * self.BODdecay # calculate BOD load after decay (in grams)
                
        if self.loadsPerSector:
            self.routedDomBOD = (self.routedDomBOD + self.Dom_BODload) * self.BODdecay
            self.routedManBOD = (self.routedManBOD + self.Man_BODload) * self.BODdecay
            self.routedUSRBOD = (self.routedUSRBOD + self.USR_BODload) * self.BODdecay
            self.routedintLivBOD = (self.routedintLivBOD + self.intLiv_BODload) * self.BODdecay
            self.routedextLivBOD = (self.routedextLivBOD + self.extLiv_BODload) * self.BODdecay
        
        ###Pathogen load (non-conservative, decay coefficient a function of temperature, solar radiation and sedimentation)
        self.routedFC = self.routedFC + self.FCload #get FC load before decay
              
        #---Temperature dependent decay
        self.waterTemp_FC = self.waterTemp - pcr.scalar(273.15)    #water temperature in gridcell in C
        self.darkinactivation_FC = pcr.scalar(0.82)     #days-1; Reder et al., (2015)
        self.watertempcorrection_FC = pcr.scalar(1.07)     #Reder et al., (2015)
        FCdecay_temperature = cover(self.darkinactivation_FC * (self.watertempcorrection_FC**(self.waterTemp_FC - 20)),0.0)
        #---Solar radiation dependent decay
        self.water_height_pathogen = pcr.max(self.water_height, 0.1) #set minimum water depth of 0.1m for decay coefficients
        self.sunlightinactivation_FC = pcr.scalar(0.0068)     #m2 w-1     #Reder et al., (2015) 
        self.attenuation_FC = 0.0931 * self.tss + 0.881       #m-1; Reder et al., (2015)
        self.solarradiation_FC = self.rsw                     #w m-2
        FCdecay_solarradiation = cover(self.sunlightinactivation_FC * (self.solarradiation_FC/ (self.attenuation_FC * self.water_height_pathogen))*(1-(exp(-(self.attenuation_FC * self.water_height_pathogen)))),0.0)
        #---Sedimentation
        self.threshold_FC_settlingdepth = pcr.scalar(0.5) #m ; stream depth must exceed 50cm in order for sedimentation to occur
        self.settlingvelocity_FC = pcr.scalar(1.656)    #m/day; Reder et al., (2015)
        FCdecay_sedimentation = pcr.ifthenelse(self.water_height_pathogen > self.threshold_FC_settlingdepth, cover(self.settlingvelocity_FC / self.water_height,0.0), 0)  #day -1
        #---Combine decay coefficients
        self.FCdecay = exp(-(FCdecay_temperature + FCdecay_solarradiation + FCdecay_sedimentation)* self.travel_time)        
        
        #---FC routing
        self.routedFC = self.routedFC * self.FCdecay #calculate FC load after decay
        
        if self.loadsPerSector:
            self.routedDomFC = (self.routedDomFC + self.Dom_FCload) * self.FCdecay
            self.routedManFC = (self.routedManFC + self.Man_FCload) * self.FCdecay
            self.routedUSRFC = (self.routedUSRFC + self.USR_FCload) * self.FCdecay
            self.routedintLivFC = (self.routedintLivFC + self.intLiv_FCload) * self.FCdecay
            self.routedextLivFC = (self.routedextLivFC + self.extLiv_FCload) * self.FCdecay 
            
    def qualityRouting(self, timeSec):
        
        channelTransFrac = cover(pcr.max(pcr.min((self.subDischarge * timeSec) / self.channelStorageTimeBefore, 1.0),0.0), 0.0)
        
        #Energy (for water temperature) routing 
        dtotEWLat= channelTransFrac*self.volumeEW
        self.volumeEW = (self.volumeEW +pcr.upstream(self.lddMap,dtotEWLat)-dtotEWLat)
        
        #Salinity (TDS) routing
        dTDSLat = channelTransFrac*self.routedTDS
        self.routedTDS = (self.routedTDS +pcr.upstream(self.lddMap,dTDSLat)-dTDSLat)        
        
        if self.loadsPerSector:
            dDomTDSLat = channelTransFrac*self.routedDomTDS
            self.routedDomTDS = (self.routedDomTDS +pcr.upstream(self.lddMap,dDomTDSLat)-dDomTDSLat) 
            dManTDSLat = channelTransFrac*self.routedManTDS
            self.routedManTDS = (self.routedManTDS +pcr.upstream(self.lddMap,dManTDSLat)-dManTDSLat)
            dUSRTDSLat = channelTransFrac*self.routedUSRTDS
            self.routedUSRTDS = (self.routedUSRTDS +pcr.upstream(self.lddMap,dUSRTDSLat)-dUSRTDSLat)
            dIrrTDSLat = channelTransFrac*self.routedIrrTDS
            self.routedIrrTDS = (self.routedIrrTDS +pcr.upstream(self.lddMap,dIrrTDSLat)-dIrrTDSLat)
        
        #Organic (BOD) routing
        dBODLat = channelTransFrac*self.routedBOD
        self.routedBOD = (self.routedBOD +pcr.upstream(self.lddMap,dBODLat)-dBODLat)
        
        if self.loadsPerSector:          
            dDomBODLat = channelTransFrac*self.routedDomBOD
            self.routedDomBOD = (self.routedDomBOD +pcr.upstream(self.lddMap,dDomBODLat)-dDomBODLat) 
            dManBODLat = channelTransFrac*self.routedManBOD
            self.routedManBOD = (self.routedManBOD +pcr.upstream(self.lddMap,dManBODLat)-dManBODLat)
            dUSRBODLat = channelTransFrac*self.routedUSRBOD
            self.routedUSRBOD = (self.routedUSRBOD +pcr.upstream(self.lddMap,dUSRBODLat)-dUSRBODLat)
            dintLivBODLat = channelTransFrac*self.routedintLivBOD
            self.routedintLivBOD = (self.routedintLivBOD +pcr.upstream(self.lddMap,dintLivBODLat)-dintLivBODLat)
            dextLivBODLat = channelTransFrac*self.routedextLivBOD
            self.routedextLivBOD = (self.routedextLivBOD +pcr.upstream(self.lddMap,dextLivBODLat)-dextLivBODLat)
        
        #Pathogen (FC) routing
        dFCLat = channelTransFrac*self.routedFC
        self.routedFC = (self.routedFC +pcr.upstream(self.lddMap,dFCLat)-dFCLat)

        if self.loadsPerSector:          
            dDomFCLat = channelTransFrac*self.routedDomFC
            self.routedDomFC = (self.routedDomFC +pcr.upstream(self.lddMap,dDomFCLat)-dDomFCLat) 
            dManFCLat = channelTransFrac*self.routedManFC
            self.routedManFC = (self.routedManFC +pcr.upstream(self.lddMap,dManFCLat)-dManFCLat)
            dUSRFCLat = channelTransFrac*self.routedUSRFC
            self.routedUSRFC = (self.routedUSRFC +pcr.upstream(self.lddMap,dUSRFCLat)-dUSRFCLat)
            dintLivFCLat = channelTransFrac*self.routedintLivFC
            self.routedintLivFC = (self.routedintLivFC +pcr.upstream(self.lddMap,dintLivFCLat)-dintLivFCLat)
            dextLivFCLat = channelTransFrac*self.routedextLivFC
            self.routedextLivFC = (self.routedextLivFC +pcr.upstream(self.lddMap,dextLivFCLat)-dextLivFCLat)

    def energyWaterBody(self):
        lakeTransFrac = pcr.max(pcr.min((self.WaterBodies.waterBodyOutflow) / (self.waterBodyStorageTimeBefore), 1.0),0.0)
        lakeTransFrac = cover(ifthen(self.WaterBodies.waterBodyOut, lakeTransFrac), 0.0)
            
        energyTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.totEW * self.dynamicFracWat * self.cellArea),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.totEW * self.dynamicFracWat * self.cellArea)
        self.volumeEW = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., energyTotal*lakeTransFrac),energyTotal)
        self.remainingVolumeEW = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * energyTotal), 0.0)
        
    def qualityWaterBody(self):
        lakeTransFrac = pcr.max(pcr.min((self.WaterBodies.waterBodyOutflow) / (self.waterBodyStorageTimeBefore), 1.0),0.0)
        lakeTransFrac = cover(ifthen(self.WaterBodies.waterBodyOut, lakeTransFrac), 0.0)
        
        #Salinity (amount in water body)
        wbTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.routedTDS),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedTDS)
        self.wbVolumeTDS = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbTDSTotal*lakeTransFrac),wbTDSTotal)
        self.wbRemainingVolumeTDS = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbTDSTotal), 0.0)
        
        if self.loadsPerSector:
            wbDomTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedDomTDS),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedDomTDS)
            self.wbVolumeDomTDS = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbDomTDSTotal*lakeTransFrac),wbDomTDSTotal)
            self.wbRemainingVolumeDomTDS = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbDomTDSTotal), 0.0)
            
            wbManTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedManTDS),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedManTDS)
            self.wbVolumeManTDS = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbManTDSTotal*lakeTransFrac),wbManTDSTotal)
            self.wbRemainingVolumeManTDS = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbManTDSTotal), 0.0)
            
            wbUSRTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedUSRTDS),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedUSRTDS)
            self.wbVolumeUSRTDS = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbUSRTDSTotal*lakeTransFrac),wbUSRTDSTotal)
            self.wbRemainingVolumeUSRTDS = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbUSRTDSTotal), 0.0)
            
            wbIrrTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedIrrTDS),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedIrrTDS)
            self.wbVolumeIrrTDS = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbIrrTDSTotal*lakeTransFrac),wbIrrTDSTotal)
            self.wbRemainingVolumeIrrTDS = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbIrrTDSTotal), 0.0)
            
        #Organic (amount in water body)
        wbBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.routedBOD),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedBOD)
        self.wbVolumeBOD = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbBODTotal*lakeTransFrac),wbBODTotal)
        self.wbRemainingVolumeBOD = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbBODTotal), 0.0)

        if self.loadsPerSector: 
            wbDomBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedDomBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedDomBOD)
            self.wbVolumeDomBOD = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbDomBODTotal*lakeTransFrac),wbDomBODTotal)
            self.wbRemainingVolumeDomBOD = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbDomBODTotal), 0.0)
            
            wbManBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedManBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedManBOD)
            self.wbVolumeManBOD = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbManBODTotal*lakeTransFrac),wbManBODTotal)
            self.wbRemainingVolumeManBOD = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbManBODTotal), 0.0)
            
            wbUSRBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedUSRBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedUSRBOD)
            self.wbVolumeUSRBOD = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbUSRBODTotal*lakeTransFrac),wbUSRBODTotal)
            self.wbRemainingVolumeUSRBOD = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbUSRBODTotal), 0.0)
            
            wbintLivBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedintLivBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedintLivBOD)
            self.wbVolumeintLivBOD = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbintLivBODTotal*lakeTransFrac),wbintLivBODTotal)
            self.wbRemainingVolumeintLivBOD = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbintLivBODTotal), 0.0)
            
            wbextLivBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedextLivBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedextLivBOD)
            self.wbVolumeextLivBOD = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbextLivBODTotal*lakeTransFrac),wbextLivBODTotal)
            self.wbRemainingVolumeextLivBOD = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbextLivBODTotal), 0.0)
        
        #Pathogen (amount in water body)
        wbFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.routedFC),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedFC)
        self.wbVolumeFC = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbFCTotal*lakeTransFrac),wbFCTotal)
        self.wbRemainingVolumeFC = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbFCTotal), 0.0)
        
        if self.loadsPerSector:
            wbDomFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedDomFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedDomFC)
            self.wbVolumeDomFC = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbDomFCTotal*lakeTransFrac),wbDomFCTotal)
            self.wbRemainingVolumeDomFC = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbDomFCTotal), 0.0)
            
            wbManFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedManFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedManFC)
            self.wbVolumeManFC = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbManFCTotal*lakeTransFrac),wbManFCTotal)
            self.wbRemainingVolumeManFC = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbManFCTotal), 0.0)
            
            wbUSRFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedUSRFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedUSRFC)
            self.wbVolumeUSRFC = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbUSRFCTotal*lakeTransFrac),wbUSRFCTotal)
            self.wbRemainingVolumeUSRFC = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbUSRFCTotal), 0.0)
            
            wbintLivFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedintLivFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedintLivFC)
            self.wbVolumeintLivFC = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbintLivFCTotal*lakeTransFrac),wbintLivFCTotal)
            self.wbRemainingVolumeintLivFC = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbintLivFCTotal), 0.0)
            
            wbextLivFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.routedextLivFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedextLivFC)
            self.wbVolumeextLivFC = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., wbextLivFCTotal*lakeTransFrac),wbextLivFCTotal)
            self.wbRemainingVolumeextLivFC = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * wbextLivFCTotal), 0.0)

    def energyWaterBodyAverage(self):
        
        #Energy (averaged over water body)
        self.totalVolumeEW = self.volumeEW + self.remainingVolumeEW

        energyTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.totalVolumeEW),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.totalVolumeEW)
             
        energyAverageLakeCell = cover(energyTotal * self.cellArea \
          /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), energyTotal)
        self.totEW = cover(energyAverageLakeCell /(self.dynamicFracWat * self.cellArea), 1e-16)
            
        self.temp_water_height = self.return_water_body_storage_to_channel(self.channelStorageNow)/(self.dynamicFracWat * self.cellArea)
      
        iceReductionFactor = ifthen(self.landmask, cover(self.dynamicFracWatBeforeRouting/self.dynamicFracWat,1.0))
        
        self.deltaIceThickness = iceReductionFactor * self.deltaIceThickness
        self.deltaIceThickness= pcr.min(self.deltaIceThickness,self.temp_water_height)
        
        self.iceThickness = iceReductionFactor * self.iceThickness

        self.iceThickness= pcr.max(0,self.iceThickness+(self.deltaIceThickness+pcr.ifthenelse(self.temperatureKelvin >= self.iceThresTemp,0,self.correctPrecip)))
        self.iceThickness= pcr.ifthenelse((self.iceThickness <= 0.001) & (self.deltaIceThickness < 0),0,self.iceThickness)
        self.channelStorageNow = self.channelStorageNow - self.deltaIceThickness * self.dynamicFracWat * self.cellArea        
        
        self.waterTemp= pcr.ifthenelse(self.temp_water_height > self.critical_water_height,\
          self.totEW/self.temp_water_height/(self.specificHeatWater*self.densityWater),self.temperatureKelvin)
        #self.waterTemp= pcr.ifthenelse(self.waterTemp < self.iceThresTemp+0.1,self.iceThresTemp+0.1,self.waterTemp)
        self.waterTemp= min(pcr.ifthenelse(self.waterTemp < self.iceThresTemp+0.1,self.iceThresTemp+0.1,self.waterTemp), self.maxThresTemp)

    def qualityWaterBodyAverage(self):
        
        #Salinity (averaged over water body)
        wbTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeTDS),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedTDS)
        TDSAverageLakeCell = cover(wbTDSTotal * self.cellArea \
          /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbTDSTotal)
        self.routedTDS = cover(TDSAverageLakeCell, vos.MV)
        
        if self.loadsPerSector:        
            wbDomTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeDomTDS),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedDomTDS)
            DomTDSAverageLakeCell = cover(wbDomTDSTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbDomTDSTotal)
            self.routedDomTDS = cover(DomTDSAverageLakeCell, vos.MV)
            
            wbManTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeManTDS),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedManTDS)
            ManTDSAverageLakeCell = cover(wbManTDSTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbManTDSTotal)
            self.routedManTDS = cover(ManTDSAverageLakeCell, vos.MV)
            
            wbUSRTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeUSRTDS),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedUSRTDS)
            USRTDSAverageLakeCell = cover(wbUSRTDSTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbUSRTDSTotal)
            self.routedUSRTDS = cover(USRTDSAverageLakeCell, vos.MV)
            
            wbIrrTDSTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeIrrTDS),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedIrrTDS)
            IrrTDSAverageLakeCell = cover(wbIrrTDSTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbIrrTDSTotal)
            self.routedIrrTDS = cover(IrrTDSAverageLakeCell, vos.MV)
        
        #Organic (averaged over water body)
        wbBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeBOD),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedBOD)
        BODAverageLakeCell = cover(wbBODTotal * self.cellArea \
          /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbBODTotal)
        self.routedBOD = cover(BODAverageLakeCell, vos.MV)

        if self.loadsPerSector:        
            wbDomBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeDomBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedDomBOD)
            DomBODAverageLakeCell = cover(wbDomBODTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbDomBODTotal)
            self.routedDomBOD = cover(DomBODAverageLakeCell, vos.MV)
            
            wbManBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeManBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedManBOD)
            ManBODAverageLakeCell = cover(wbManBODTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbManBODTotal)
            self.routedManBOD = cover(ManBODAverageLakeCell, vos.MV)
            
            wbUSRBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeUSRBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedUSRBOD)
            USRBODAverageLakeCell = cover(wbUSRBODTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbUSRBODTotal)
            self.routedUSRBOD = cover(USRBODAverageLakeCell, vos.MV)
            
            wbintLivBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeintLivBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedintLivBOD)
            intLivBODAverageLakeCell = cover(wbintLivBODTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbintLivBODTotal)
            self.routedintLivBOD = cover(intLivBODAverageLakeCell, vos.MV)
            
            wbextLivBODTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeextLivBOD),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedextLivBOD)
            extLivBODAverageLakeCell = cover(wbextLivBODTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbextLivBODTotal)
            self.routedextLivBOD = cover(extLivBODAverageLakeCell, vos.MV)
        
        #Pathogen (averaged over water body)
        wbFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeFC),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedFC)
        FCAverageLakeCell = cover(wbFCTotal * self.cellArea \
          /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbFCTotal)
        self.routedFC = cover(FCAverageLakeCell, vos.MV)
        
        if self.loadsPerSector:        
            wbDomFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeDomFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedDomFC)
            DomFCAverageLakeCell = cover(wbDomFCTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbDomFCTotal)
            self.routedDomFC = cover(DomFCAverageLakeCell, vos.MV)
            
            wbManFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeManFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedManFC)
            ManFCAverageLakeCell = cover(wbManFCTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbManFCTotal)
            self.routedManFC = cover(ManFCAverageLakeCell, vos.MV)
            
            wbUSRFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeUSRFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedUSRFC)
            USRFCAverageLakeCell = cover(wbUSRFCTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbUSRFCTotal)
            self.routedUSRFC = cover(USRFCAverageLakeCell, vos.MV)
            
            wbintLivFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeintLivFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedintLivFC)
            intLivFCAverageLakeCell = cover(wbintLivFCTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbintLivFCTotal)
            self.routedintLivFC = cover(intLivFCAverageLakeCell, vos.MV)
            
            wbextLivFCTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
             pcr.areatotal(pcr.ifthen(self.landmask,self.wbRemainingVolumeextLivFC),\
             pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.routedextLivFC)
            extLivFCAverageLakeCell = cover(wbextLivFCTotal * self.cellArea \
              /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), wbextLivFCTotal)
            self.routedextLivFC = cover(extLivFCAverageLakeCell, vos.MV)
