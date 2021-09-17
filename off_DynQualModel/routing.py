#!/usr/bin/ python
# -*- coding: utf-8 -*-

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
        result['m2tDischargeLong']         = self.m2tDischarge               # (m3/s)^2
        result['avgBaseflowLong']          = self.avgBaseflow                #  m3/s   ;  long term average baseflow
        result['riverbedExchange']         = self.riverbedExchange           #  m3/day : river bed infiltration (from surface water bdoies to groundwater)
        result['waterBodyStorage']            = self.waterBodyStorage        #  m3     ; storages of lakes and reservoirs            # values given are per water body id (not per cell)
        result['avgLakeReservoirOutflowLong'] = self.avgOutflow              #  m3/s   ; long term average lake & reservoir outflow  # values given are per water body id (not per cell)
        result['avgLakeReservoirInflowShort'] = self.avgInflow               #  m3/s   ; short term average lake & reservoir inflow  # values given are per water body id (not per cell)
        result['avgDischargeShort']        = self.avgDischargeShort          #  m3/s   ; short term average discharge 
        result['subDischarge']             = self.subDischarge               #  m3/s   ; sub-time step discharge (needed for kinematic wave methods/approaches)
        result['avg_irrGrossDemand']       = self.avg_irrGrossDemand         #  m/day  ; average irrigation gross demand    
        result['avg_netLqWaterToSoil']     = self.avg_netLqWaterToSoil       #  m/day  ; average net liquid transferred to the soil
        
        # Water quality elements
        result['waterTemperature']        = self.waterTemp                   #  K      ; water temperature
        result['iceThickness']            = self.iceThickness                #  m      ; ice thickness
        result['salinity']                = self.salinity                    #  g TDS  ; routed salinity load (for conversion to salinity pollution in mg/L)
        result['organic']                 = self.organic                     #  g BOD  ; routed organic load (for conversion to organic pollution mg/L)
        result['pathogen']                = self.pathogen                    #  cfu    ; routed pathogen load [10^6] (for conversion to pathogen pollution in cfu/100mL)
        
        return result

    def __init__(self,iniItems,initialConditions,lddMap):
        object.__init__(self)

        self.lddMap = lddMap

        self.cloneMap = iniItems.cloneMap
        self.tmpDir = iniItems.tmpDir
        self.inputDir = iniItems.globalOptions['inputDir']

        # option to activate water balance check
        self.debugWaterBalance = True
        if iniItems.routingOptions['debugWaterBalance'] == "False":
            self.debugWaterBalance = False

        self.method = iniItems.routingOptions['routingMethod']
        
        self.baseflowNC = vos.getFullPath(iniItems.routingOptions['baseflowNC'], self.inputDir)
        self.interflowNC = vos.getFullPath(iniItems.routingOptions['interflowNC'], self.inputDir)
        self.directRunoffNC = vos.getFullPath(iniItems.routingOptions['directRunoffNC'], self.inputDir)

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
                                         self.lddMap == 5,1.)
        distanceDownstream = pcr.ldddist(self.lddMap,\
                                         self.lddMap == 5,\
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
          logger.info("Quality started")
        except:
          self.quality = False
        if self.method == "accuTravelTime" or "kinematicWave": self.quality = False
        self.quality = True
        print("waterTemperature =",self.quality)
        print("Salinity = ",self.quality)
        print("Organic = ", self.quality)
        print("Pathogen= ", self.quality)
        
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
            self.sunFracTBL= vos.getFullPath(iniItems.meteoOptions['sunhoursTable'], self.inputDir) #convert cloud cover to sunshine hours (Doornkamp & Pruitt)           
            self.radCon= 0.25
            self.radSlope= 0.50
            self.stefanBoltzman= 5.67e-8 # [W/m2/K]
            self.maxThresTemp = pcr.scalar(318.15) # max river temperature set to 318.15 K (or 45C)

            #- Paths to cloudcover file
            self.cloudFileNC = vos.getFullPath(iniItems.meteoOptions['cloudcoverNC'], self.inputDir)
            self.radFileNC = vos.getFullPath(iniItems.meteoOptions['radiationNC'], self.inputDir)
            self.vapFileNC = vos.getFullPath(iniItems.meteoOptions['vaporNC'], self.inputDir)
            self.annualT = vos.netcdf2PCRobjCloneWithoutTime(vos.getFullPath(iniItems.meteoOptions['annualAvgTNC'], self.inputDir), "temperature",self.cloneMap) + pcr.scalar(273.15)
            self.maxIceThickness= 3.0
            self.deltaIceThickness = 0.0
            
            #- Path to powerplant data
            self.PowRFNC = vos.getFullPath(iniItems.routingOptions["powerplantNC"], self.inputDir) #Power return flows aggregated to 5arc-min (Lohrmann et al., 2019)             
            
            #- Path to TSS loadings (from Beusen et al., 2005)
            self.tss = vos.readPCRmapClone(iniItems.routingOptions['TSSmap'],self.cloneMap,self.tmpDir,self.inputDir)
            
            try:
              self.calculateLoads = iniItems.routingOptions['calculateLoads'] == "True"
            except:
              self.calculateLoads = False
            
            if self.calculateLoads == True:
                
                ###File pathways and constant pollutant loading input data
                
                #Domestic
                self.PopulationNC = vos.getFullPath(iniItems.routingOptions["PopulationNC"], self.inputDir) #gridded population, annual, 5 arc-min (Lange & Geiger, 2020)
                self.Dom_ExcrLoadNC = vos.getFullPath(iniItems.routingOptions["Dom_ExcrLoadNC"], self.inputDir) #average (regional) excretion (g or cfu/capita/day)
                self.DomTDS_ExcrLoad = vos.netcdf2PCRobjCloneWithoutTime(self.Dom_ExcrLoadNC,"Dom_Fixed_TDSload",self.cloneMap)
                self.DomBOD_ExcrLoad = vos.netcdf2PCRobjCloneWithoutTime(self.Dom_ExcrLoadNC,"Dom_Fixed_BODload",self.cloneMap)
                self.DomFC_ExcrLoad = vos.netcdf2PCRobjCloneWithoutTime(self.Dom_ExcrLoadNC,"Dom_Fixed_FCload",self.cloneMap)
                
                #Manufacturing
                self.ManWWpNC = vos.getFullPath(iniItems.routingOptions["ManWWpNC"], self.inputDir)
                self.Man_EfflConcNC = vos.getFullPath(iniItems.routingOptions["Man_EfflConcNC"], self.inputDir) #average (regional) man effluent concentration (mg/L or cfu/100ml)
                self.ManTDS_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.Man_EfflConcNC,"Man_Fixed_TDSload",self.cloneMap) 
                self.ManBOD_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.Man_EfflConcNC,"Man_Fixed_BODload",self.cloneMap)
                self.ManFC_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.Man_EfflConcNC,"Man_Fixed_FCload",self.cloneMap)                
                
                #Urban suface runoff          
                self.UrbanFractionNC = vos.getFullPath(iniItems.routingOptions["UrbanFractionNC"], self.inputDir)
                self.USR_EfflConcNC = vos.getFullPath(iniItems.routingOptions["USR_EfflConcNC"], self.inputDir) #average (regional) USR effluent concentration (mg/L or cfu/100ml)
                self.USRTDS_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.USR_EfflConcNC,"USR_Fixed_TDSload",self.cloneMap)
                self.USRBOD_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.USR_EfflConcNC,"USR_Fixed_BODload",self.cloneMap)
                self.USRFC_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.USR_EfflConcNC,"USR_Fixed_FCload",self.cloneMap)
                
                #Livestock
                self.LivPopulationNC = vos.getFullPath(iniItems.routingOptions["LivPopulationNC"], self.inputDir) #Gilbert et al., (2010)
                self.Liv_ExcrLoadNC = vos.getFullPath(iniItems.routingOptions["Liv_ExcrLoadNC"], self.inputDir)  #Gridded livestock populations, 2010, 5 arc-min (Gilbert et al., 2018)
                self.Buffalo_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Buffalo_BODload",self.cloneMap)
                self.Chicken_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Chicken_BODload",self.cloneMap)
                self.Cow_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Cow_BODload",self.cloneMap)
                self.Duck_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Duck_BODload",self.cloneMap)
                self.Goat_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Goat_BODload",self.cloneMap)
                self.Horse_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Horse_BODload",self.cloneMap)
                self.Pig_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Pig_BODload",self.cloneMap)
                self.Sheep_BODload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Sheep_BODload",self.cloneMap)
                self.Buffalo_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Buffalo_FCload",self.cloneMap)
                self.Chicken_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Chicken_FCload",self.cloneMap)
                self.Cow_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Cow_FCload",self.cloneMap)
                self.Duck_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Duck_FCload",self.cloneMap)
                self.Goat_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Goat_FCload",self.cloneMap)
                self.Horse_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Horse_FCload",self.cloneMap)
                self.Pig_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Pig_FCload",self.cloneMap)
                self.Sheep_FCload = vos.netcdf2PCRobjCloneWithoutTime(self.Liv_ExcrLoadNC,"Sheep_FCload",self.cloneMap)                
                
                #Irrigation
                self.Irr_EfflConcNC = vos.getFullPath(iniItems.routingOptions["Irr_EfflConcNC"], self.inputDir) #average soil salinity averaged over the topsoil and subsoil (mg/L) 
                self.IrrTDS_EfflConc = vos.netcdf2PCRobjCloneWithoutTime(self.Irr_EfflConcNC,"soil_TDS",self.cloneMap)                 
                                                    
                #Wastewater pathways and removal efficiencies (treatment [tertiary, secondary, primary], collected but untreated, basic sanitation, open defecation, direct)
                self.WWtPathwaysNC = vos.getFullPath(iniItems.routingOptions["WWtPathwaysNC"], self.inputDir)
                #self.WWtRemEffTBL = iniItems.routingOptions["WWtRemEffs"] #TODO: Change this to a flexible table
                self.TDS_Ter_RemEff = 0.
                self.TDS_Sec_RemEff = 0.
                self.TDS_Pri_RemEff = 0.
                self.TDS_WWcut_RemEff = 0.
                self.TDS_dom_WWbs_RemEff = 0.
                self.TDS_dom_WWod_RemEff = 0.
                self.TDS_man_WWdirect_RemEff = 0.
                  
                self.BOD_Ter_RemEff = 0.99
                self.BOD_Sec_RemEff = 0.85
                self.BOD_Pri_RemEff = 0.25
                self.BOD_WWcut_RemEff = 0.
                self.BOD_dom_WWbs_RemEff = 0.
                self.BOD_dom_WWod_RemEff = 0.
                self.BOD_man_WWdirect_RemEff = 0.
                  
                self.FC_Ter_RemEff = 0.9999
                self.FC_Sec_RemEff = 0.9745
                self.FC_Pri_RemEff = 0.4279
                self.FC_dom_WWbs_RemEff = 0.
                self.FC_dom_WWod_RemEff = 0.
                self.FC_man_WWdirect_RemEff = 0.
                
            else:
            #- Path to (non-natural) salinity, organic, pathogen loading inputs
                self.salinityNC = vos.getFullPath(iniItems.routingOptions["salinityNC"], self.inputDir)
                self.organicNC = vos.getFullPath(iniItems.routingOptions["organicNC"], self.inputDir)
                self.pathogenNC = vos.getFullPath(iniItems.routingOptions["pathogenNC"], self.inputDir)
                               
        # get the initialConditions
        self.getICs(iniItems, initialConditions)
        
        # initiate old style reporting                                  # TODO: remove this!
        self.initiate_old_style_routing_reporting(iniItems)
        

    def getICs(self,iniItems,iniConditions = None):

        if iniConditions == None:

            # read initial conditions from pcraster maps listed in the ini file (for the first time step of the model; when the model just starts)
            self.timestepsToAvgDischarge = vos.readPCRmapClone(iniItems.routingOptions['timestepsToAvgDischargeIni'] ,self.cloneMap,self.tmpDir,self.inputDir)  
            self.channelStorage          = vos.readPCRmapClone(iniItems.routingOptions['channelStorageIni']          ,self.cloneMap,self.tmpDir,self.inputDir)
            self.readAvlChannelStorage   = vos.readPCRmapClone(iniItems.routingOptions['readAvlChannelStorageIni']   ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.avgDischarge            = vos.readPCRmapClone(iniItems.routingOptions['avgDischargeLongIni']        ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.m2tDischarge            = vos.readPCRmapClone(iniItems.routingOptions['m2tDischargeLongIni']        ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.avgBaseflow             = vos.readPCRmapClone(iniItems.routingOptions['avgBaseflowLongIni']         ,self.cloneMap,self.tmpDir,self.inputDir) 
            self.riverbedExchange        = vos.readPCRmapClone(iniItems.routingOptions['riverbedExchangeIni']        ,self.cloneMap,self.tmpDir,self.inputDir) 
            
            # New initial condition variable introduced in the version 2.0.2: avgDischargeShort 
            self.avgDischargeShort       = vos.readPCRmapClone(iniItems.routingOptions['avgDischargeShortIni']       ,self.cloneMap,self.tmpDir,self.inputDir) 

            # Initial conditions needed for kinematic wave methods
            self.subDischarge            = vos.readPCRmapClone(iniItems.routingOptions['subDischargeIni'],self.cloneMap,self.tmpDir,self.inputDir)

            # Initial conditions needed for water quality module
            if self.quality:
                self.waterTemp               = vos.readPCRmapClone(iniItems.routingOptions['waterTemperatureIni'],self.cloneMap,self.tmpDir,self.inputDir)
                self.iceThickness               = vos.readPCRmapClone(iniItems.routingOptions['iceThicknessIni'],self.cloneMap,self.tmpDir,self.inputDir)
                self.salinity = vos.readPCRmapClone(iniItems.routingOptions['salinityIni'],self.cloneMap,self.tmpDir,self.inputDir) #initial conditions for salinity pollution
                self.organic = vos.readPCRmapClone(iniItems.routingOptions['organicIni'],self.cloneMap,self.tmpDir,self.inputDir) #initial conditions for organic pollution
                self.pathogen = vos.readPCRmapClone(iniItems.routingOptions['pathogenIni'],self.cloneMap,self.tmpDir,self.inputDir) #initial conditions for pathogen pollution
                
            # Initial conditions for calculating average irrigation demand and net liquid transferred to the soil for irrigation return flow calculations
            if self.calculateLoads:
                self.avg_irrGrossDemand   = vos.readPCRmapClone(iniItems.routingOptions['avg_irrGrossDemandIni'],self.cloneMap,self.tmpDir,self.inputDir)      
                self.avg_netLqWaterToSoil = vos.readPCRmapClone(iniItems.routingOptions['avg_netLqWaterToSoil'], self.cloneMap,self.tmpDir,self.inputDir) 

        else:              

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

            if self.quality:
                self.waterTemp               = iniConditions['routing']['waterTemperature']
                self.iceThickness            = iniConditions['routing']['iceThickness']
                self.salinity               = iniConditions['routing']['salinity']
                self.organic               = iniConditions['routing']['organic']
                self.pathogen              = iniConditions['routing']['pathogen']
                
            if self.calculateLoads:    
                self.avg_irrGrossDemand   = iniConditions['routing']['avg_irrGrossDemand']      
                self.avg_netLqWaterToSoil = iniConditions['routing']['avg_netLqWaterToSoil']

        self.channelStorage        = pcr.ifthen(self.landmask, pcr.cover(self.channelStorage, 0.0))
        self.readAvlChannelStorage = pcr.ifthen(self.landmask, pcr.cover(self.readAvlChannelStorage, 0.0))
        self.avgDischarge          = pcr.ifthen(self.landmask, pcr.cover(self.avgDischarge, 0.0))
        self.m2tDischarge          = pcr.ifthen(self.landmask, pcr.cover(self.m2tDischarge, 0.0))
        self.avgDischargeShort     = pcr.ifthen(self.landmask, pcr.cover(self.avgDischargeShort, 0.0))
        self.avgBaseflow           = pcr.ifthen(self.landmask, pcr.cover(self.avgBaseflow, 0.0))
        self.riverbedExchange      = pcr.ifthen(self.landmask, pcr.cover(self.riverbedExchange, 0.0))
        self.subDischarge          = pcr.ifthen(self.landmask, pcr.cover(self.subDischarge , 0.0))
        
        self.readAvlChannelStorage = pcr.min(self.readAvlChannelStorage, self.channelStorage)
        self.readAvlChannelStorage = pcr.max(self.readAvlChannelStorage, 0.0)
        
        if self.floodPlain:
            self.subDischarge          = pcr.ifthen(self.landmask, pcr.cover(self.subDischarge , 0.0))

        if self.quality:
            self.waterTemp             = pcr.ifthen(self.landmask, pcr.cover(self.waterTemp, 0.0))
            self.iceThickness          = pcr.ifthen(self.landmask, pcr.cover(self.iceThickness , 0.0))
            self.channelStorageTimeBefore = self.channelStorage
            self.totEW = self.channelStorage * self.waterTemp*self.specificHeatWater * self.densityWater
            self.temp_water_height = yMean = self.eta * pow (self.avgDischarge, self.nu)
            self.salinity = pcr.ifthen(self.landmask, pcr.cover(self.salinity, 0.0))
            self.organic = pcr.ifthen(self.landmask, pcr.cover(self.organic, 0.0))
            self.pathogen = pcr.ifthen(self.landmask, pcr.cover(self.pathogen,  0.0))
            
        if self.calculateLoads:  
            self.avg_irrGrossDemand   = pcr.ifthen(self.landmask, pcr.cover(self.avg_irrGrossDemand,   0.0))
            self.avg_netLqWaterToSoil = pcr.ifthen(self.landmask, pcr.cover(self.avg_netLqWaterToSoil, 0.0))
        
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
        

    def accuTravelTime(self):
        		
        # accuTravelTime ROUTING OPERATIONS
        ##############n############################################################################################################

        # route only non negative channelStorage (otherwise stay):
        channelStorageThatWillNotMove = pcr.ifthenelse(self.channelStorage < 0.0, self.channelStorage, 0.0)
        self.channelStorage           = pcr.max(0.0, self.channelStorage)
        
        # also at least 1.0 m3 of water will stay - this is to minimize numerical errors due to float_32 pcraster implementations
        channelStorageThatWillNotMove += self.channelStorage - pcr.rounddown(self.channelStorage)
        self.channelStorage            = pcr.rounddown(self.channelStorage) 
        
        # channelStorage that will be given to the ROUTING operation:
        channelStorageForAccuTravelTime = pcr.max(0.0, self.channelStorage)
        channelStorageForAccuTravelTime = pcr.cover(channelStorageForAccuTravelTime,0.0)       # TODO: check why do we have to use the "cover" operation?

        # estimating channel discharge (m3/day)
        self.Q = pcr.accutraveltimeflux(self.lddMap,\
                                        channelStorageForAccuTravelTime,\
                                        self.characteristicDistance)
        self.Q = pcr.cover(self.Q, 0.0)
        # for very small velocity (i.e. characteristicDistanceForAccuTravelTime), discharge can be missing value.
        # see: http://sourceforge.net/p/pcraster/bugs-and-feature-requests/543/
        #      http://karssenberg.geo.uu.nl/tt/TravelTimeSpecification.htm
        #
        # and make sure that no negative discharge
        self.Q = pcr.max(0.0, self.Q)                                   # unit: m3/day        

        # updating channelStorage (after routing)
        #
        # - alternative 1: using accutraveltimestate
        self.channelStorage = pcr.accutraveltimestate(self.lddMap,\
                              channelStorageForAccuTravelTime,\
                              self.characteristicDistance)              # unit: m3

        #~ # - alternative 2: using the calculated Q (Can we do this?)
        #~ storage_change_in_volume  = pcr.upstream(self.lddMap, self.Q) - self.Q
        #~ channelStorageForRouting += storage_change_in_volume 

        # return channelStorageThatWillNotMove to channelStorage:
        self.channelStorage += channelStorageThatWillNotMove            # unit: m3

        # for non kinematic wave approach, set subDishcarge to missing values
        self.subDischarge = pcr.scalar(vos.MV) 

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

        logger.info("Using the simplifiedKinematicWave method ! ")
        
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
        

    def update_routing_only(self,currTimeStep,meteo):
        
        logger.info("routing only, reading fluxes")
        
        #Read hydrological input directly from netCDF files
        self.baseflow = vos.netcdf2PCRobjClone(\
                                  self.baseflowNC,"baseflow",\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        
        self.interflow = vos.netcdf2PCRobjClone(\
                                  self.interflowNC,"interflow",\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
                                  
        self.directRunoff = vos.netcdf2PCRobjClone(\
                                  self.directRunoffNC, "direct_runoff",\
                                  str(currTimeStep.fulldate),
                                  useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True)
        self.runoff = self.directRunoff + self.interflow + self.baseflow

        logger.info("routing in progress")

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

        # volume water released in pits (losses: to the ocean / endorheic basin)
        self.outgoing_volume_at_pits = pcr.ifthen(self.landmask,
                                       pcr.cover(
                                       pcr.ifthen(self.lddMap == pcr.ldd(5), self.Q), 0.0))
        
        # estimate volume of water that can be extracted for abstraction in the next time step
        self.readAvlChannelStorage = self.estimate_available_volume_for_abstraction(self.channelStorage)
                
        # old-style reporting                             
        self.old_style_routing_reporting(currTimeStep)                 # TODO: remove this one
        
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
                                               self.waterKC * meteo.referencePotET)) ##TODO -\
        ## TODO check                                       landSurface.actualET ))              # These values are NOT over the entire cell area.
        
        # potential evaporation from water bodies over the entire cell area (m/day)
        waterBodyPotEvap = waterBodyPotEvapOvesSurfaceWaterArea * self.dynamicFracWat
        return waterBodyPotEvap
        
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
        # - TODO: This concept should be IMPROVED. 
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

    def simple_update_routing_only(self,currTimeStep,meteo):

        # updating timesteps to calculate long and short term statistics values of avgDischarge, avgInflow, avgOutflow, etc.
        self.timestepsToAvgDischarge += 1.

        if self.debugWaterBalance:\
           preStorage = self.channelStorage                                                        # unit: m3

        # the following variable defines total local change (input) to surface water storage bodies # unit: m3 
        # - only local processes; therefore not considering any routing processes
        self.local_input_to_surface_water = pcr.scalar(0.0)          # initiate the variable, start from zero

        # runoff from landSurface cells (unit: m/day)
        #self.runoff = self.directRunoff + self.interflow + self.baseflow   
        
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
            if currTimeStep.day == 1 or currTimeStep.timeStepPCR == 1:
                self.readExtensiveMeteo(currTimeStep)
                self.readPowerplantData(currTimeStep)
                
            if self.calculateLoads == True:
                self.readPollutantLoadingsInputData(currTimeStep, landSurface, groundwater)
                self.calculatePollutantLoadings(currTimeStep)            
            else:
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


        # average irrigation water that is allocated from the last 30 days
        if self.calculateLoads:
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
        self.cloudCover = vos.netcdf2PCRobjClone(\
                                 self.cloudFileNC,'cld',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)
                                  
        self.radiation =  vos.netcdf2PCRobjClone(\
                                 self.radFileNC,'RSW',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = "month",
                                 cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)

        self.vaporPressure = vos.netcdf2PCRobjClone(\
                                 self.vapFileNC,'vap',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
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

    def readPowerplantData(self, currTimeStep):
        logger.info("Reading powerplant data")
        #read Power return flows (from Lohrmann et al., 2019)
        self.PowRF =  vos.netcdf2PCRobjClone(\
                                 self.PowRFNC,'PowRF',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                 cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)
        
        #Define delta T heat dump
        self.deltaT = pcr.scalar(7.) #difference between effluent temperature and river temperature (K) (van Vliet et al., 2012) 
        self.PowTwload = self.PowRF * self.specificHeatWater * self.densityWater * self.deltaT #W
    
    def readPollutantLoadingsInputData(self, currTimeStep, landSurface, groundwater):
        logger.info("Loading input data required to calculate pollutant loadings")    
        
        #Domestic
        self.Population = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'Population',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100) 
        self.frac_surfaceRunoff = vos.getValDivZero(landSurface.directRunoff , self.runoff)        
                        
        #Manufacturing
        self.ManWWp = vos.netcdf2PCRobjClone(\
                                 self.ManWWpNC,'man_WWp',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)
                                  
        #Urban surface runoff
        #- Urban surface runoff directly from PCR-GLOBWB (gridcell surface runoff [m/day] * fraction urban area in gridcell * cell area) [TODO ED]
        self.urban_area_fraction = vos.netcdf2PCRobjClone(\
                                             self.UrbanFractionNC,'urban_fraction',\
                                             str(currTimeStep.fulldate), 
                                             useDoy = None,
                                             cloneMapFileName=self.cloneMap,\
                                             LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)
        self.USR_RF = self.directRunoff * self.urban_area_fraction * self.cellArea #in m3/day

        #Livestock
        self.BuffaloPopulation = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'BuffaloPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100) 

        self.ChickenPopulation = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'ChickenPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100) 
                                  
        self.CowPopulation = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'CowPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100) 

        self.DuckPopulation = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'DuckPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)                                                                      
          
        self.GoatPopulation = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'GoatPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100) 

        self.HorsePopulation = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'HorsePop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100) 
                                  
        self.PigPopulation = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'PigPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100) 

        self.SheepPopulation = vos.netcdf2PCRobjClone(\
                                 self.PopulationNC,'SheepPop',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)           

        #Calculate livestock densities accounting for livestock units (Wen et al., 2018)
        self.cellArea_km2 = self.cellArea/1000000.
        self.LivDensityThres = 25.
        
        self.BuffaloDensity = self.BuffaloPopulation/ self.cellArea_km2
        self.ChickenDensity = (self.ChickenPopulation * 0.01) /self.cellArea_km2
        self.CowDensity = self.CowPopulation/ self.cellArea_km2
        self.DuckDensity = (self.DuckPopulation * 0.01)/ self.cellArea_km2
        self.GoatDensity = (self.GoatPopulation * 0.1)/ self.cellArea_km2
        self.HorseDensity = self.HorsePopulation/ self.cellArea_km2
        self.PigDensity = (self.PigPopulation * 0.3)/ self.cellArea_km2
        self.SheepDensity = (self.BuffaloPopulation * 0.1)/ self.cellArea_km2
    
        #Irrigation
        #- Irrigation return flows calculated directly from PCR-GLOBWB (from direct runoff, interflow and baseflow) [TODO ED]
        self.irr_rf_from_direct_runoff = self.directRunoff * vos.getValDivZero(landSurface.irrGrossDemand, ( landSurface.irrGrossDemand + landSurface.netLqWaterToSoil))
        self.irr_rf_from_interflow     = self.interflowTotal * vos.getValDivZero(self.avg_irrGrossDemand, ( self.avg_irrGrossDemand + self.avg_netLqWaterToSoil))
        self.irr_rf_from_baseflow     =  self.baseflow       * vos.getValDivZero(self.avg_irrGrossDemand, ( self.avg_irrGrossDemand + self.avg_netLqWaterToSoil))
        
        self.Irr_RF = (self.irr_rf_from_direct_runoff +\
                       self.irr_rf_from_interflow +\
                       self.irr_rf_from_baseflow) * self.cellArea  # - total irirgation return flow in volume units (i.e. m3 day-1) 
                                  
        #Wastewater Pathways
        self.ratio_WWt_Tertiary = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'WWt_Tertiary',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)
                                  
        self.ratio_WWt_Secondary = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'WWt_Secondary',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)
                                  
        self.ratio_WWt_Primary = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'WWt_Primary',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)                                                                            
        
        self.ratio_WWcut = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'WWcut',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100) 
                                  
        self.ratio_dom_WWbs = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'dom_WWbs',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)
                                  
        self.ratio_dom_WWod = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'dom_WWod',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)
                                  
        self.ratio_man_WWdirect = vos.netcdf2PCRobjClone(\
                                 self.WWtPathwaysNC,'man_WWdirect',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                  cloneMapFileName=self.cloneMap,\
                                  LatitudeLongitude = True,specificFillValue = -999.)/pcr.scalar(100)                                     

    def calculatePollutantLoadings(self, currTimeStep):
        #calculate pollutant loadings directly
        logger.info("Calculating pollutant loadings")
        
        #Municipal wastewater treatment     
        TDS_RemEff = (self.ratio_WWt_Tertiary * TDS_Ter_RemEff) +\
                     (self.ratio_WWt_Secondary * TDS_Sec_RemEff) +\
                     (self.ratio_WWt_Primary * TDS_Pri_RemEff) +\
                     (self.ratio_WWcut * TDS_WWcut_RemEff)
        BOD_RemEff = (self.ratio_WWt_Tertiary * BOD_Ter_RemEff) +\
                     (self.ratio_WWt_Secondary * BOD_Sec_RemEff) +\
                     (self.ratio_WWt_Primary * BOD_Pri_RemEff) +\
                     (self.ratio_WWcut * BOD_WWcut_RemEff)
        FC_RemEff =  (self.ratio_WWt_Tertiary * FC_Ter_RemEff) +\
                     (self.ratio_WWt_Secondary * FC_Sec_RemEff) +\
                     (self.ratio_WWt_Primary * FC_Pri_RemEff) +\
                     (self.ratio_WWcut * FC_WWcut_RemEff)
                
        #Domestic loadings: Gridded population (capita) * per Capita excretion rate [g/capita/day; cfu/capita/day]
        TDS_Rem_BasicSan = self.ratio_dom_WWbs * TDS_dom_WWbs_RemEff #Removal of TDS due to basic sanitation practices
        BOD_Rem_BasicSan = self.ratio_dom_WWbs * BOD_dom_WWbs_RemEff #Removal of BOD due to basic sanitation practices
        FC_Rem_BasicSan = self.ratio_dom_WWbs * FC_dom_WWbs_RemEff #Removal of FC due to basic sanitation practices
        Rem_OpenDef = self.ratio_dom_WWod * self.frac_surfaceRunoff #Pollutants transferred to surface water proportional to fraction surface runoff
        
        self.Dom_TDSload = self.Population * self.DomTDS_ExcrLoad * (1 - (TDS_RemEff + TDS_Rem_BasicSan + Rem_OpenDef)) #g/day
        self.Dom_BODload = self.Population * self.DomBOD_ExcrLoad * (1 - (BOD_RemEff + BOD_Rem_BasicSan + Rem_OpenDef)) #g/day
        self.Dom_FCload = (self.Population * self.DomFC_ExcrLoad * (1 - (FC_RemEff + FC_Rem_BasicSan + Rem_OpenDef))) /1000000. #million cfu/day
        
        #Manufacturing loadings: Manufacturing wastewater [m3/day] * average manufacturing effluent concentration [mg/L; cfu/100ml]
        self.Man_TDSload = self.ManWWp * self.ManTDS_EfflConc * (1 - TDS_RemEff) #g/day
        self.Man_BODload = self.ManWWp * self.ManBOD_EfflConc * (1 - BOD_RemEff) #g/day
        self.Man_FCload = self.ManWWp * self.ManFC_EfflConc * (1 - FC_RemEff)/ 100. #million cfu/day
        
        #Urban surface runoff loadings: USR return flow [m3/day] * average USR effluent concentration [mg/L; cfu/100ml]
        self.USR_TDSload = self.USR_RF * self.USRTDS_EfflConc * (1 - TDS_RemEff) #g/day
        self.USR_BODload = self.USR_RF * self.USRBOD_EfflConc * (1 - BOD_RemEff) #g/day
        self.USR_FCload = (self.USR_RF * self.USRFC_EfflConc * (1 - FC_RemEff))/ 100. #million cfu/day
                
        #Livestock loadings
        BOD_Rem_IntLiv = (self.ratio_WWt_Tertiary + ratio_WWt_Secondary) * BOD_Sec_RemEff
        FC_Rem_IntLiv = (self.ratio_WWt_Tertiary + ratio_WWt_Secondary) * FC_Sec_RemEff
        
        self.intLiv_Buffalo_BODload = pcr.ifthenelse(self.BuffaloDensity > self.LivDensityThres,\
                                      self.BuffaloPopulation * self.Buffalo_BODload * (1-BOD_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
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

        self.extLiv_Buffalo_BODload = pcr.ifthenelse(self.BuffaloDensity <= self.LivDensityThres, self.BuffaloPopulation * self.Buffalo_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Chicken_BODload = pcr.ifthenelse(self.ChickenDensity <= self.LivDensityThres, self.ChickenPopulation * self.Chicken_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Cow_BODload = pcr.ifthenelse(self.CowDensity <= self.LivDensityThres, self.CowPopulation * self.Cow_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Duck_BODload = pcr.ifthenelse(self.DuckDensity <= self.LivDensityThres, self.DuckPopulation * self.Duck_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Goat_BODload = pcr.ifthenelse(self.GoatDensity <= self.LivDensityThres, self.GoatPopulation * self.Goat_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Horse_BODload = pcr.ifthenelse(self.HorseDensity <= self.LivDensityThres, self.HorsePopulation * self.Horse_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Pig_BODload = pcr.ifthenelse(self.PigDensity <= self.LivDensityThres, self.PigPopulation * self.Pig_BODload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Sheep_BODload = pcr.ifthenelse(self.SheepDensity <= self.LivDensityThres, self.SheepPopulation * self.Sheep_BODload * self.frac_surfaceRunoff, 0.0)
        
        self.intLiv_BODload = self.intLiv_Buffalo_BODload + self.intLiv_Chicken_BODload + self.intLiv_Cow_BODload +\
                              self.intLiv_Duck_BODload + self.intLiv_Goat_BODload + self.intLiv_Horse_BODload + self.intLiv_Pig_BODload + self.intLiv_Sheep_BODload #g/day
        self.extLiv_BODload = self.extLiv_Buffalo_BODload + self.extLiv_Chicken_BODload + self.extLiv_Cow_BODload +\
                              self.extLiv_Duck_BODload + self.extLiv_Goat_BODload + self.extLiv_Horse_BODload + self.extLiv_Pig_BODload + self.extLiv_Sheep_BODload #g/day

        self.intLiv_Buffalo_FCload = pcr.ifthenelse(self.BuffaloDensity > self.LivDensityThres,\
                                     self.BuffaloPopulation * self.Buffalo_FCload * (1-FC_Rem_IntLiv) * self.frac_surfaceRunoff, 0.0)
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
        
        self.extLiv_Buffalo_FCload = pcr.ifthenelse(self.BuffaloDensity <= self.LivDensityThres, self.BuffaloPopulation * self.Buffalo_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Chicken_FCload = pcr.ifthenelse(self.ChickenDensity <= self.LivDensityThres, self.ChickenPopulation * self.Chicken_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Cow_FCload = pcr.ifthenelse(self.CowDensity <= self.LivDensityThres, self.CowPopulation * self.Cow_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Duck_FCload = pcr.ifthenelse(self.DuckDensity <= self.LivDensityThres, self.DuckPopulation * self.Duck_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Goat_FCload = pcr.ifthenelse(self.GoatDensity <= self.LivDensityThres, self.GoatPopulation * self.Goat_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Horse_FCload = pcr.ifthenelse(self.HorseDensity <= self.LivDensityThres, self.HorsePopulation * self.Horse_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Pig_FCload = pcr.ifthenelse(self.PigDensity <= self.LivDensityThres, self.PigPopulation * self.Pig_FCload * self.frac_surfaceRunoff, 0.0)
        self.extLiv_Sheep_FCload = pcr.ifthenelse(self.SheepDensity <= self.LivDensityThres, self.SheepPopulation * self.Sheep_FCload * self.frac_surfaceRunoff, 0.0)
        
        self.intLiv_FCload = (self.intLiv_Buffalo_FCload + self.intLiv_Chicken_FCload + self.intLiv_Cow_FCload +\
                              self.intLiv_Duck_FCload + self.intLiv_Goat_FCload + self.intLiv_Horse_FCload + self.intLiv_Pig_FCload + self.intLiv_Sheep_FCload)\
                              / 1000000. #million cfu/day
        self.extLiv_FCload = (self.extLiv_Buffalo_FCload + self.extLiv_Chicken_FCload + self.extLiv_Cow_FCload +\
                              self.extLiv_Duck_FCload + self.extLiv_Goat_FCload + self.extLiv_Horse_FCload + self.extLiv_Pig_FCload + self.extLiv_Sheep_FCload)\
                              / 1000000. #million cfu/day
        
        #Irrigation loadings: Irrigation return flow (m3/day) * soil TDS (mg/L)
        self.Irr_TDSload = self.Irr_RF * self.IrrTDS_EfflConc #g/day
        
        #Combined loadings
        self.TDSload = self.Dom_TDSload + self.Man_TDSload + self.USR_TDSload + self.Irr_TDSload
        self.BODload = self.Dom_BODload + self.Man_BODload + self.USR_BODload + self.intLiv_BODload + self.extLiv_BODload
        self.FCload = self.Dom_FCload + self.Man_FCload + self.USR_FCload + self.intLiv_FCload + self.extLiv_FCload
        
    def readPollutantLoadings(self, currTimeStep):
        #Use pre-calculated pollutant loadings
        logger.info("Reading loadings directly")        
        
        # read TDS loadings (combined from domestic, manufacturing, urban surface runoff and irrigation sources)
        self.TDSload = vos.netcdf2PCRobjClone(\
                                 self.salinityNC,'TDS',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                 cloneMapFileName=self.cloneMap,\
                                 LatitudeLongitude = True)

        #read BOD loadings (combined from domestic, manufacturing, urban surface runoff and livestock sources)
        self.BODload = vos.netcdf2PCRobjClone(\
                                 self.organicNC,'BOD',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                 cloneMapFileName=self.cloneMap,\
                                 LatitudeLongitude = True)							 

        #read FC loadings (combined from domestic, manufacturing, urban surface runoff and livestock sources)
        self.FCload = vos.netcdf2PCRobjClone(\
                                 self.pathogenNC,'FC',\
                                 str(currTimeStep.fulldate), 
                                 useDoy = None,
                                 cloneMapFileName=self.cloneMap,\
                                 LatitudeLongitude = True)                                                  

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
           self.interflow/landRunoff*pcr.max(self.iceThresTemp+0.1,self.temperatureKelvin)+\
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
        self.salinity = self.salinity + self.TDSload 
        
        ###Organic load (non-conservative, temperature dependent decay)
        self.organic = self.organic + self.BODload #get BOD load before decay
        
        #---Temperature
        self.k_BOD = pcr.scalar(0.35)     #first-order degradation coefficient at 20C (van Vliet et al., 2021)
        self.watertempcorrection_BOD = pcr.scalar(1.047)     #temperature correction (van Vliet et al., 2021; Wen et al., 2017)
        self.waterTemp_BOD = self.waterTemp - pcr.scalar(273.15)   #water temperature in gridcell in C
        BODdecay_temperature = cover(self.k_BOD*(self.watertempcorrection_BOD**(self.waterTemp_BOD - 20)),0.0)
        self.BODdecay = exp(-(BODdecay_temperature)* self.travel_time)
        
        self.organic = self.organic * self.BODdecay # calculate BOD load after decay
        
        ###Pathogen load (non-conservative, decay coefficient a function of temperature, solar radiation and sedimentation)
        self.pathogen = self.pathogen + self.FCload #get FC load before decay
              
        #---Temperature
        self.waterTemp_FC = self.waterTemp - pcr.scalar(273.15)    #water temperature in gridcell in C
        self.darkinactivation_FC = pcr.scalar(0.82)     #days-1; Reder et al., (2015)
        self.watertempcorrection_FC = pcr.scalar(1.07)     #Reder et al., (2015)
        FCdecay_temperature = cover(self.darkinactivation_FC * (self.watertempcorrection_FC**(self.waterTemp_FC - 20)),0.0)
        #---Solar radiation
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
        
        self.pathogen = self.pathogen * self.FCdecay #calculate FC load after decay
        
    def qualityRouting(self, timeSec):
        
        channelTransFrac = cover(pcr.max(pcr.min((self.subDischarge * timeSec) / self.channelStorageTimeBefore, 1.0),0.0), 0.0)
        
        #Energy (for water temperature) routing 
        dtotEWLat= channelTransFrac*self.volumeEW
        self.volumeEW = (self.volumeEW +pcr.upstream(self.lddMap,dtotEWLat)-dtotEWLat)
        
        #Salinity (TDS) routing
        dSalinityLat = channelTransFrac*self.salinity
        self.salinity = (self.salinity +pcr.upstream(self.lddMap,dSalinityLat)-dSalinityLat)        
        
        #Organic (BOD) routing
        dOrganicLat = channelTransFrac*self.organic
        self.organic = (self.organic +pcr.upstream(self.lddMap,dOrganicLat)-dOrganicLat)
                
        #Pathogen (FC) routing
        dPathogenLat = channelTransFrac*self.pathogen
        self.pathogen = (self.pathogen +pcr.upstream(self.lddMap,dPathogenLat)-dPathogenLat)

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
        salinityTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.salinity),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.salinity)
        self.volumeSalinity = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., salinityTotal*lakeTransFrac),salinityTotal)
        self.remainingVolumeSalinity = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * salinityTotal), 0.0)
        
        #Organic (amount in water body)
        organicTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.organic),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.organic)
        self.volumeOrganic = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., organicTotal*lakeTransFrac),organicTotal)
        self.remainingVolumeOrganic = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * organicTotal), 0.0)
        
        #Pathogen (amount in water body)
        lakeTransFrac = pcr.max(pcr.min((self.WaterBodies.waterBodyOutflow) / (self.waterBodyStorageTimeBefore), 1.0),0.0)
        lakeTransFrac = cover(ifthen(self.WaterBodies.waterBodyOut, lakeTransFrac), 0.0)

        pathogenTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.pathogen),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.pathogen)
        self.volumePathogen = cover(ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0., pathogenTotal*lakeTransFrac),pathogenTotal)
        self.remainingVolumePathogen = cover(ifthen(self.WaterBodies.waterBodyOut, (1-lakeTransFrac) * pathogenTotal), 0.0)

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
        salinityTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.remainingVolumeSalinity),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.salinity)
        salinityAverageLakeCell = cover(salinityTotal * self.cellArea \
          /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), salinityTotal)
        self.salinity = cover(salinityAverageLakeCell, vos.MV)
        
        #Organic (averaged over water body)
        organicTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.remainingVolumeOrganic),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.organic)
        organicAverageLakeCell = cover(organicTotal * self.cellArea \
          /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), organicTotal)
        self.organic = cover(organicAverageLakeCell, vos.MV)
        
        #Pathogen (averaged over water body)
        pathogenTotal = cover(pcr.ifthen(pcr.scalar(self.WaterBodies.waterBodyIds) > 0.,
         pcr.areatotal(pcr.ifthen(self.landmask,self.remainingVolumePathogen),\
         pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds))), self.pathogen)
        pathogenAverageLakeCell = cover(pathogenTotal * self.cellArea \
          /pcr.areatotal(pcr.cover(self.cellArea, 0.0),pcr.ifthen(self.landmask,self.WaterBodies.waterBodyIds)), pathogenTotal)
        self.pathogen = cover(pathogenAverageLakeCell, vos.MV)
