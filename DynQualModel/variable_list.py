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
List of variables for PCR-GLOBWB2 and DynQual
@authors (PCR-GLOBWB2): Edwin H. Sutanudjaja
@authors (DynQual)    : Edward R. Jones, Niko Wanders
'''

netcdf_short_name = {}
netcdf_unit       = {}
netcdf_monthly_total_unit = {} 
netcdf_yearly_total_unit  = {}
netcdf_long_name  = {}
description       = {}
comment           = {}
latex_symbol      = {}

# precipitation
pcrglobwb_variable_name = 'precipitation'
netcdf_short_name[pcrglobwb_variable_name] = 'precipitation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# temperature
pcrglobwb_variable_name = 'temperature'
netcdf_short_name[pcrglobwb_variable_name] = 'temperature'
netcdf_unit[pcrglobwb_variable_name]       = 'degrees Celcius'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'mean_air_temperature'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# referencePotET
pcrglobwb_variable_name = 'referencePotET'
netcdf_short_name[pcrglobwb_variable_name] = 'reference_potential_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# totalLandSurfacePotET
pcrglobwb_variable_name = 'totalLandSurfacePotET'
netcdf_short_name[pcrglobwb_variable_name] = 'land_surface_potential_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'total_potential_evaporation_and_transpiration_at_land_surface'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Not including water bodies. Values given are over the entire cell area.'
latex_symbol[pcrglobwb_variable_name]      = None

# totLandSurfaceActuaET
pcrglobwb_variable_name = 'totLandSurfaceActuaET'
netcdf_short_name[pcrglobwb_variable_name] = 'land_surface_actual_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'total_actual_evaporation_and_transpiration_at_land_surface'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Not including water bodies. Values given are over the entire cell area.'
latex_symbol[pcrglobwb_variable_name]      = None

# fractionLandSurfaceET
pcrglobwb_variable_name = 'fractionLandSurfaceET'
netcdf_short_name[pcrglobwb_variable_name] = 'land_surface_evaporation_fraction'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'ratio_between_actual_and_potential_values_of_evaporation_and_transpiration_at_land_surface'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Not including water bodies.'
latex_symbol[pcrglobwb_variable_name]      = None

# interceptStor
pcrglobwb_variable_name = 'interceptStor'
netcdf_short_name[pcrglobwb_variable_name] = 'interception_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# snowCoverSWE 
pcrglobwb_variable_name = 'snowCoverSWE'
netcdf_short_name[pcrglobwb_variable_name] = 'snow_water_equivalent'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'snow_cover_in_water_equivalent_amount'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# snowFreeWater
pcrglobwb_variable_name = 'snowFreeWater'
netcdf_short_name[pcrglobwb_variable_name] = 'snow_free_water'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'liquid_water_within_snowpack'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# topWaterLayer
pcrglobwb_variable_name = 'topWaterLayer'
netcdf_short_name[pcrglobwb_variable_name] = 'top_water_layer'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'water_layer_storage_above_soil'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# storUppTotal 
pcrglobwb_variable_name = 'storUppTotal'
netcdf_short_name[pcrglobwb_variable_name] = 'upper_soil_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'upper_soil_storage'       # first 30 cm of soil
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# storLowTotal 
pcrglobwb_variable_name = 'storLowTotal'
netcdf_short_name[pcrglobwb_variable_name] = 'lower_soil_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'lower_soil_storage'       # next 30-150 cm of soil
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# interceptEvap       
pcrglobwb_variable_name = 'interceptEvap'
netcdf_short_name[pcrglobwb_variable_name] = 'interception_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'evaporation_from_interception_storage'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# actSnowFreeWaterEvap
pcrglobwb_variable_name = 'actSnowFreeWaterEvap'
netcdf_short_name[pcrglobwb_variable_name] = 'snow_free_water_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'evaporation_from_liquid_water_within_snowpack'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# topWaterLayerEvap   
pcrglobwb_variable_name = 'topWaterLayerEvap'
netcdf_short_name[pcrglobwb_variable_name] = 'top_water_layer_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'evaporation_from_water_layer_storage_above_soil'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# actBareSoilEvap     
pcrglobwb_variable_name = 'actBareSoilEvap'
netcdf_short_name[pcrglobwb_variable_name] = 'bare_soil_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# actTranspiTotal     
pcrglobwb_variable_name = 'actTranspiTotal'
netcdf_short_name[pcrglobwb_variable_name] = 'total_transpiration'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'total_plant_transpiration_from_entire_soil_storages'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# actTranspiUppTotal
pcrglobwb_variable_name = 'actTranspiUppTotal'
netcdf_short_name[pcrglobwb_variable_name] = 'upper_soil_transpiration'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'total_plant_transpiration_from_upper_soil_storage(s)'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# actTranspiLowTotal
pcrglobwb_variable_name = 'actTranspiLowTotal'
netcdf_short_name[pcrglobwb_variable_name] = 'lower_soil_transpiration'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'total_plant_transpiration_from_lower_soil_storage'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# directRunoff                    
pcrglobwb_variable_name = 'directRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'direct_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# fraction direct runoff
pcrglobwb_variable_name = 'frac_surfaceRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'fraction_surface_runoff'
netcdf_unit[pcrglobwb_variable_name]       = '-'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "fraction of runoff that is surface runoff"
latex_symbol[pcrglobwb_variable_name]      = None

# interflowTotal                  
pcrglobwb_variable_name = 'interflowTotal'
netcdf_short_name[pcrglobwb_variable_name] = 'interflow'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# baseflow                  
pcrglobwb_variable_name = 'baseflow'
netcdf_short_name[pcrglobwb_variable_name] = 'baseflow'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# infiltration                    
pcrglobwb_variable_name = 'infiltration'
netcdf_short_name[pcrglobwb_variable_name] = 'infiltration'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# gwRecharge                      
pcrglobwb_variable_name = 'gwRecharge'
netcdf_short_name[pcrglobwb_variable_name] = 'groundwater_recharge'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "negative values indicating (net) capillary rise from groundater store"
latex_symbol[pcrglobwb_variable_name]      = None

# gwNetCapRise                      
pcrglobwb_variable_name = 'gwNetCapRise'
netcdf_short_name[pcrglobwb_variable_name] = 'groundwater_capillary_rise'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "values (positive) indicating (net) capillary rise from groundater store; only positive values given to the field."
latex_symbol[pcrglobwb_variable_name]      = None

# irrGrossDemand                  
pcrglobwb_variable_name = 'irrGrossDemand'
netcdf_short_name[pcrglobwb_variable_name] = 'irrigation_gross_demand'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# nonIrrGrossDemand                  
pcrglobwb_variable_name = 'nonIrrGrossDemand'
netcdf_short_name[pcrglobwb_variable_name] = 'non_irrigation_gross_demand'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# totalGrossDemand                  
pcrglobwb_variable_name = 'totalGrossDemand'
netcdf_short_name[pcrglobwb_variable_name] = 'total_gross_demand'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# satDegUpp                       
pcrglobwb_variable_name = 'satDegUpp'
netcdf_short_name[pcrglobwb_variable_name] = 'upper_soil_saturation_degree'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# satDegLow                       
pcrglobwb_variable_name = 'satDegLow'
netcdf_short_name[pcrglobwb_variable_name] = 'lower_soil_saturation_degree'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# storGroundwater                 
pcrglobwb_variable_name = 'storGroundwater'
netcdf_short_name[pcrglobwb_variable_name] = 'groundwater_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'non_fossil_groundwater_storage'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# storGroundwaterFossil                 
pcrglobwb_variable_name = 'storGroundwaterFossil'
netcdf_short_name[pcrglobwb_variable_name] = 'fossil_groundwater_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# storGroundwaterTotal                 
pcrglobwb_variable_name = 'storGroundwaterTotal'
netcdf_short_name[pcrglobwb_variable_name] = 'total_groundwater_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Non fossil and fossil groundwater storage.'
latex_symbol[pcrglobwb_variable_name]      = None

# surfaceWaterAbstraction         
pcrglobwb_variable_name = 'surfaceWaterAbstraction'
netcdf_short_name[pcrglobwb_variable_name] = 'surface_water_abstraction'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# nonFossilGroundWaterAbstraction 
pcrglobwb_variable_name = 'nonFossilGroundWaterAbstraction'
netcdf_short_name[pcrglobwb_variable_name] = 'non_fossil_groundwater_abstraction'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# otherWaterSourceAbstraction     
pcrglobwb_variable_name = 'otherWaterSourceAbstraction'
netcdf_short_name[pcrglobwb_variable_name] = 'other_water_source_abstraction'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# totalAbstraction
pcrglobwb_variable_name = 'totalAbstraction'
netcdf_short_name[pcrglobwb_variable_name] = 'total__abstraction'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = "Total abstraction from all water sources: surface water, non fossil groundwater and other water sources (e.g. fossil groundwater, desalinisation)."
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# fracSurfaceWaterAllocation
pcrglobwb_variable_name = 'fracSurfaceWaterAllocation'
netcdf_short_name[pcrglobwb_variable_name] = 'fraction_of_surface_water_allocation'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Values equal to 1 indicate either 100% allocation (from surface water) or zero water demand."
latex_symbol[pcrglobwb_variable_name]      = None

# fracNonFossilGroundwaterAllocation
pcrglobwb_variable_name = 'fracNonFossilGroundwaterAllocation'
netcdf_short_name[pcrglobwb_variable_name] = 'fraction_of_non_fossil_groundwater_allocation'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Values equal to 0 indicate either zero allocation or zero water demand."
latex_symbol[pcrglobwb_variable_name]      = None

# fracOtherWaterSourceAllocation
pcrglobwb_variable_name = 'fracOtherWaterSourceAllocation'
netcdf_short_name[pcrglobwb_variable_name] = 'fraction_of_other_water_source_allocation'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Values equal to 0 indicate either zero allocation or zero water demand."
latex_symbol[pcrglobwb_variable_name]      = None

# totalFracWaterSourceAllocation
pcrglobwb_variable_name = 'totalFracWaterSourceAllocation'
netcdf_short_name[pcrglobwb_variable_name] = 'total_fraction_water_allocation'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "All values must be equal to 1. Otherwise, water balance errors."
latex_symbol[pcrglobwb_variable_name]      = None

# waterBodyActEvaporation
pcrglobwb_variable_name = 'waterBodyActEvaporation'
netcdf_short_name[pcrglobwb_variable_name] = 'water_body_actual_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Flux values given are over the entire cell area (not only over surface water body fraction).'
latex_symbol[pcrglobwb_variable_name]      = None

# waterBodyPotEvaporation
pcrglobwb_variable_name = 'waterBodyPotEvaporation'
netcdf_short_name[pcrglobwb_variable_name] = 'water_body_potential_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Flux values given are over the entire cell area (not only over surface water body fraction).'
latex_symbol[pcrglobwb_variable_name]      = None

# fractionWaterBodyEvaporation
pcrglobwb_variable_name = 'fractionWaterBodyEvaporation'
netcdf_short_name[pcrglobwb_variable_name] = 'water_body_evaporation_fraction'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'ratio_between_actual_and_potential_values_of_evaporation_and_transpiration_at_surface_water_bodies'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# totalEvaporation
pcrglobwb_variable_name = 'totalEvaporation'
netcdf_short_name[pcrglobwb_variable_name] = 'total_evaporation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Including from water bodies.'
latex_symbol[pcrglobwb_variable_name]      = None

# runoff
pcrglobwb_variable_name = 'runoff'
netcdf_short_name[pcrglobwb_variable_name] = 'land_surface_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "direct_runoff + interflow + baseflow, but not including local runoff from water bodies."
latex_symbol[pcrglobwb_variable_name]      = None

# accuRunoff
pcrglobwb_variable_name = 'accuRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'accumulated_land_surface_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm3.s-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "direct_runoff + interflow + baseflow, but not including local runoff from water bodies."
latex_symbol[pcrglobwb_variable_name]      = None

# accuBaseflow
pcrglobwb_variable_name = 'accuBaseflow'
netcdf_short_name[pcrglobwb_variable_name] = 'accumulated_land_surface_baseflow'
netcdf_unit[pcrglobwb_variable_name]       = 'm3.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = "m3.month-1" 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = "m3.year-1"
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# discharge
pcrglobwb_variable_name = 'discharge'
netcdf_short_name[pcrglobwb_variable_name] = 'discharge'
netcdf_unit[pcrglobwb_variable_name]       = 'm3.s-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# totalRunoff
pcrglobwb_variable_name = 'totalRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'total_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Including local changes at water bodies."
latex_symbol[pcrglobwb_variable_name]      = None

# local_water_body_flux
pcrglobwb_variable_name = 'local_water_body_flux'
netcdf_short_name[pcrglobwb_variable_name] = 'local_water_body_flux'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# accuTotalRunoff
pcrglobwb_variable_name = 'accuTotalRunoff'
netcdf_short_name[pcrglobwb_variable_name] = 'accumulated_total_surface_runoff'
netcdf_unit[pcrglobwb_variable_name]       = 'm3.s-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Including runoff from water bodies."
latex_symbol[pcrglobwb_variable_name]      = None

# totalActiveStorageThickness
pcrglobwb_variable_name = 'totalActiveStorageThickness'
netcdf_short_name[pcrglobwb_variable_name] = 'total_thickness_of_active_water_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Not including fossil groundwater (unmetDemand)."
latex_symbol[pcrglobwb_variable_name]      = None

# totalWaterStorageThickness
pcrglobwb_variable_name = 'totalWaterStorageThickness'
netcdf_short_name[pcrglobwb_variable_name] = 'total_thickness_of_water_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = "Including fossil groundwater (unmetDemand)."
latex_symbol[pcrglobwb_variable_name]      = None

# surfaceWaterStorage
pcrglobwb_variable_name = 'surfaceWaterStorage'
netcdf_short_name[pcrglobwb_variable_name] = 'surface_water_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Negative values may be reported, due to excessive demands.'
latex_symbol[pcrglobwb_variable_name]      = None

# waterBodyStorage 
pcrglobwb_variable_name = 'waterBodyStorage'
netcdf_short_name[pcrglobwb_variable_name] = 'lake_and_reservoir_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm3'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'The values given are for every lake and reservoir ids (not per cells) and after lake/reservoir releases/outflows.'
latex_symbol[pcrglobwb_variable_name]      = None

# snowMelt
pcrglobwb_variable_name = 'snowMelt'
netcdf_short_name[pcrglobwb_variable_name] = 'snow_melt'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# satDegUppSurface                       
pcrglobwb_variable_name = 'satDegUppSurface'
netcdf_short_name[pcrglobwb_variable_name] = 'near_surface_soil_saturation_degree'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'This variable can only be reported if 3 layer soil model is used.'
latex_symbol[pcrglobwb_variable_name]      = None

# storUppSurface 
pcrglobwb_variable_name = 'storUppSurface'
netcdf_short_name[pcrglobwb_variable_name] = 'near_surface_soil_storage'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = None                       # first 5 cm of soil
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'This variable can only be reported if 3 layer soil model is used.'
latex_symbol[pcrglobwb_variable_name]      = None

# nonIrrWaterConsumption                  
pcrglobwb_variable_name = 'nonIrrWaterConsumption'
netcdf_short_name[pcrglobwb_variable_name] = 'consumption_for_non_irrigation_demand'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# nonIrrReturnFlow                  
pcrglobwb_variable_name = 'nonIrrReturnFlow'
netcdf_short_name[pcrglobwb_variable_name] = 'return_flow_from_non_irrigation_demand'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = 'return_flow_from_non_irrigation_demand_to_surface_water_bodies'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = None
latex_symbol[pcrglobwb_variable_name]      = None

# land_surface_water_balance                  
pcrglobwb_variable_name = 'land_surface_water_balance'
netcdf_short_name[pcrglobwb_variable_name] = 'land_surface_water_balance'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Excluding surface water bodies.'
latex_symbol[pcrglobwb_variable_name]      = None

# evaporation_from_irrigation
pcrglobwb_variable_name = 'evaporation_from_irrigation'
netcdf_short_name[pcrglobwb_variable_name] = 'evaporation_from_irrigation'
netcdf_unit[pcrglobwb_variable_name]       = 'm.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = 'm.month-1' 
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = 'm.year-1'
netcdf_long_name[pcrglobwb_variable_name]  = None
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Flux values given are over the entire cell area (not only irrigation fraction).'
latex_symbol[pcrglobwb_variable_name]      = None

# channel storage
pcrglobwb_variable_name = 'channelStorage'
netcdf_short_name[pcrglobwb_variable_name] = 'channelStorage'
netcdf_unit[pcrglobwb_variable_name]       = 'm3'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Total_in_channel_storage'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Riverine channel storage, values include reservoirs and lakes'
latex_symbol[pcrglobwb_variable_name]      = None

# water height
pcrglobwb_variable_name = 'waterHeight'
netcdf_short_name[pcrglobwb_variable_name] = 'waterHeight'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Total_in_channel_water_height'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Water height in channel, values include reservoirs and lakes'
latex_symbol[pcrglobwb_variable_name]      = None

# dynamic water fraction
pcrglobwb_variable_name = 'dynamicFracWat'
netcdf_short_name[pcrglobwb_variable_name] = 'dynamic_water_fraction'
netcdf_unit[pcrglobwb_variable_name]       = '1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Fraction_of_cell_flooded'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Flooded fraction is a combination of channel, waterbody and floodplain'
latex_symbol[pcrglobwb_variable_name]      = None

# Water temperature
pcrglobwb_variable_name = 'waterTemp'
netcdf_short_name[pcrglobwb_variable_name] = 'waterTemperature'
netcdf_unit[pcrglobwb_variable_name]       = 'K'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Temperature_of_surface_water'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Surface water temperature assuming fully mixed conditions'
latex_symbol[pcrglobwb_variable_name]      = None

# ice Thickness
pcrglobwb_variable_name = 'iceThickness'
netcdf_short_name[pcrglobwb_variable_name] = 'iceThickness'
netcdf_unit[pcrglobwb_variable_name]       = 'm'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'Thickness_of_ice_layer'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Thickness of ice layer on channel network'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted TDS loads
pcrglobwb_variable_name = 'TDSload'
netcdf_short_name[pcrglobwb_variable_name] = 'TDSload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_TDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted TDS loadings (for salinity pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Domestic TDS loads
pcrglobwb_variable_name = 'Dom_TDSload'
netcdf_short_name[pcrglobwb_variable_name] = 'DomTDSload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_DomTDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Domestic TDS loadings (for salinity pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Manufacturing TDS loads
pcrglobwb_variable_name = 'Man_TDSload'
netcdf_short_name[pcrglobwb_variable_name] = 'ManTDSload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_ManTDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Manufacturing TDS loadings (for salinity pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Urban Surface Runoff TDS loads
pcrglobwb_variable_name = 'USR_TDSload'
netcdf_short_name[pcrglobwb_variable_name] = 'USRTDSload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_USRTDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Urban Surface Runoff TDS loadings (for salinity pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# Irr RF
pcrglobwb_variable_name = 'Irr_RF'
netcdf_short_name[pcrglobwb_variable_name] = 'Irr_RF'
netcdf_unit[pcrglobwb_variable_name]       = 'm3 day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'irrigation_return_flow'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'irrigation return flows'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Irr TDS loads
pcrglobwb_variable_name = 'Irr_TDSload'
netcdf_short_name[pcrglobwb_variable_name] = 'IrrTDSload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_IrrTDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted irrigation TDS loadings (for salinity pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted BOD loads
pcrglobwb_variable_name = 'BODload'
netcdf_short_name[pcrglobwb_variable_name] = 'BODload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_BOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted BOD loadings (for organic pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Domestic BOD loads
pcrglobwb_variable_name = 'Dom_BODload'
netcdf_short_name[pcrglobwb_variable_name] = 'DomBODload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_DomBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Domestic BOD loadings (for organic pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Manufacturing BOD loads
pcrglobwb_variable_name = 'Man_BODload'
netcdf_short_name[pcrglobwb_variable_name] = 'ManBODload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_ManBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Manufacturing BOD loadings (for organic pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Urban Surface Runoff BOD loads
pcrglobwb_variable_name = 'USR_BODload'
netcdf_short_name[pcrglobwb_variable_name] = 'USRBODload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_USRBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Urban Surface Runoff BOD loadings (for organic pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted intensive livestock BOD loads
pcrglobwb_variable_name = 'intLiv_BODload'
netcdf_short_name[pcrglobwb_variable_name] = 'intLivBODload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_intLivBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Intensive Livestock BOD loadings (for organic pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted extensive livestock BOD loads
pcrglobwb_variable_name = 'extLiv_BODload'
netcdf_short_name[pcrglobwb_variable_name] = 'extLivBODload'
netcdf_unit[pcrglobwb_variable_name]       = 'g.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_extLivBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Extensive Livestock BOD loadings (for organic pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted FC loads
pcrglobwb_variable_name = 'FCload'
netcdf_short_name[pcrglobwb_variable_name] = 'FCload'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_FC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted FC loadings (for pathogen pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Domestic FC loads
pcrglobwb_variable_name = 'Dom_FCload'
netcdf_short_name[pcrglobwb_variable_name] = 'DomFCload'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_DomFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Domestic FC loadings (for pathogen pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Manufacturing FC loads
pcrglobwb_variable_name = 'Man_FCload'
netcdf_short_name[pcrglobwb_variable_name] = 'ManFCload'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_ManFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Manufacturing FC loadings (for pathogen pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Urban Surface Runoff FC loads
pcrglobwb_variable_name = 'USR_FCload'
netcdf_short_name[pcrglobwb_variable_name] = 'USRFCload'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_USRFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Urban Surface Runoff FC loadings (for pathogen pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Intensive Livestock FC loads
pcrglobwb_variable_name = 'intLiv_FCload'
netcdf_short_name[pcrglobwb_variable_name] = 'intLivFCload'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_intLivFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Intensive Livestock FC loadings (for pathogen pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted extensive Livestock FC loads
pcrglobwb_variable_name = 'extLiv_FCload'
netcdf_short_name[pcrglobwb_variable_name] = 'extLivFCload'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu.day-1'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_extLivFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'unrouted Extensive Livestock FC loadings (for pathogen pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# unrouted Tw loads
pcrglobwb_variable_name = 'Twload'
netcdf_short_name[pcrglobwb_variable_name] = 'Twload'
netcdf_unit[pcrglobwb_variable_name]       = 'W'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'unrouted_temperature_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'heat dumps from thermoelectric powerplants (for temperature pollution)'
latex_symbol[pcrglobwb_variable_name]      = None

# routed TDS loads
pcrglobwb_variable_name = 'routedTDS'
netcdf_short_name[pcrglobwb_variable_name] = 'routedTDS'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_TDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'TDS loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Domestic TDS loads
pcrglobwb_variable_name = 'routedDomTDS'
netcdf_short_name[pcrglobwb_variable_name] = 'routedDomTDS'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_DomTDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Domestic TDS loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Manufacturing TDS loads
pcrglobwb_variable_name = 'routedManTDS'
netcdf_short_name[pcrglobwb_variable_name] = 'routedManTDS'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_ManTDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Manufacturing TDS loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Urban Surface Runoff TDS loads
pcrglobwb_variable_name = 'routedUSRTDS'
netcdf_short_name[pcrglobwb_variable_name] = 'routedUSRTDS'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_USRTDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Urban Surface Runoff TDS loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Irrigation TDS loads
pcrglobwb_variable_name = 'routedIrrTDS'
netcdf_short_name[pcrglobwb_variable_name] = 'routedIrrTDS'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_IrrTDS_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Irrigation TDS loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed BOD loads
pcrglobwb_variable_name = 'routedBOD'
netcdf_short_name[pcrglobwb_variable_name] = 'routedBOD'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_BOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'BOD loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Domestic BOD loads
pcrglobwb_variable_name = 'routedDomBOD'
netcdf_short_name[pcrglobwb_variable_name] = 'routedDomBOD'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_DomBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Domestic BOD loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Manufacturing BOD loads
pcrglobwb_variable_name = 'routedManBOD'
netcdf_short_name[pcrglobwb_variable_name] = 'routedManBOD'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_ManBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Manufacturing BOD loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Urban Surface Runoff BOD loads
pcrglobwb_variable_name = 'routedUSRBOD'
netcdf_short_name[pcrglobwb_variable_name] = 'routedUSRBOD'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_USRBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Urban Surface Runoff BOD loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Intensive Livestock BOD loads
pcrglobwb_variable_name = 'routedintLivBOD'
netcdf_short_name[pcrglobwb_variable_name] = 'routedintLivBOD'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_intLivBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Intensive Livestock BOD loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Extensive Livestock BOD loads
pcrglobwb_variable_name = 'routedextLivBOD'
netcdf_short_name[pcrglobwb_variable_name] = 'routedextLivBOD'
netcdf_unit[pcrglobwb_variable_name]       = 'g'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_extLivBOD_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Extensive Livestock BOD loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed FC loads
pcrglobwb_variable_name = 'routedFC'
netcdf_short_name[pcrglobwb_variable_name] = 'routedFC'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_FC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'FC loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Domestic FC loads
pcrglobwb_variable_name = 'routedDomFC'
netcdf_short_name[pcrglobwb_variable_name] = 'routedDomFC'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_DomFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Domestic FC loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Manufacturing FC loads
pcrglobwb_variable_name = 'routedManFC'
netcdf_short_name[pcrglobwb_variable_name] = 'routedManFC'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_ManFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Manufacturing FC loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Urban Surface Runoff FC loads
pcrglobwb_variable_name = 'routedUSRFC'
netcdf_short_name[pcrglobwb_variable_name] = 'routedUSRFC'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_USRFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Urban Surface Runoff FC loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Intensive Livestock FC loads
pcrglobwb_variable_name = 'routedintLivFC'
netcdf_short_name[pcrglobwb_variable_name] = 'routedintLivFC'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_intLivFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Intensive Livestock FC loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# routed Extensive Livestock FC loads
pcrglobwb_variable_name = 'routedextLivFC'
netcdf_short_name[pcrglobwb_variable_name] = 'routedextLivFC'
netcdf_unit[pcrglobwb_variable_name]       = 'million_cfu'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'routed_extLivFC_loadings'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'Extensive Livestock FC loadings routed through surface water network'
latex_symbol[pcrglobwb_variable_name]      = None

# Salinity pollution (concentration in TDS mg.l)
pcrglobwb_variable_name = 'salinity'
netcdf_short_name[pcrglobwb_variable_name] = 'salinity'
netcdf_unit[pcrglobwb_variable_name]       = 'mg.l'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'salinity_concentration'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'In-stream salinity (TDS) concentration in mg.l'
latex_symbol[pcrglobwb_variable_name]      = None

# Organic pollution (concentration in BOD mg.l)
pcrglobwb_variable_name = 'organic'
netcdf_short_name[pcrglobwb_variable_name] = 'organic'
netcdf_unit[pcrglobwb_variable_name]       = 'mg.l'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'organic_concentration'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'In-stream organic (BOD) concentration in mg.l'
latex_symbol[pcrglobwb_variable_name]      = None

# Pathogen pollution (concentration in FC cfu.100ml)
pcrglobwb_variable_name = 'pathogen'
netcdf_short_name[pcrglobwb_variable_name] = 'pathogen'
netcdf_unit[pcrglobwb_variable_name]       = 'cfu.100ml'
netcdf_monthly_total_unit[pcrglobwb_variable_name] = None
netcdf_yearly_total_unit[pcrglobwb_variable_name]  = None
netcdf_long_name[pcrglobwb_variable_name]  = 'pathogen_concentration'
description[pcrglobwb_variable_name]       = None
comment[pcrglobwb_variable_name]           = 'In-stream pathogen (FC) concentration in cfu.100ml'
latex_symbol[pcrglobwb_variable_name]      = None


#~ # remove/clear pcrglobwb_variable_name 
#~ pcrglobwb_variable_name = None
#~ del pcrglobwb_variable_name
