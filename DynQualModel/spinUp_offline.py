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
Initiate spin-up for DynQual in offline congifuration.
@authors (DynQual)    : Niko Wanders, Edward R. Jones
'''

import os
import math
import logging

import pcraster as pcr

import virtualOS as vos

class SpinUp(object):

    def __init__(self, iniItems):
        object.__init__(self)
        
        self.noSpinUps      = None

        self.setupConvergence(iniItems) 

    def setupConvergence(self,iniItems):

        self.noSpinUps         =   int(iniItems.globalOptions['maxSpinUpsInYears'])

        self.minConvForTotlSto = float(iniItems.globalOptions['minConvForTotlSto'])
        self.minConvForSoilSto = float(iniItems.globalOptions['minConvForSoilSto'])
        self.minConvForGwatSto = float(iniItems.globalOptions['minConvForGwatSto'])
        self.minConvForChanSto = float(iniItems.globalOptions['minConvForChanSto'])

        # TODO: including the convergence of ResvSto (reservoir storage)
        # self.minConvForResvSto = float(iniItems.globalOptions['minConvForResvSto'])
        
        self.endStateDir = iniItems.endStateDir     

    def getIniStates(self,model):

        self.iniChanSto = max(1E-20,\
                          vos.getMapVolume(\
                          model.routing.channelStorage,1))

    def channelStorageVolume(self, state, cellAreaMap):
        return vos.getMapVolume(state['routing']['channelStorage'], cellAreaMap) # unit: m3
    
    def checkConvergence(self,beginState, endState, spinUpRun, cellAreaMap):
        
        #calculate convergence of channel storage
        
        beginChanSto = max(1E-20,self.channelStorageVolume(beginState, cellAreaMap))
        endChanSto = self.channelStorageVolume(endState, cellAreaMap)
          
        convChanSto = math.fabs(100*(endChanSto-beginChanSto)/beginChanSto)
        
        logging.getLogger('spinup').info('Delta ChanStorage = %.2f percent' \
                    %(convChanSto))

        return convChanSto <= self.minConvForChanSto
