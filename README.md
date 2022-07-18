# DYNQUAL
Water quality model for simulating water temperature (Tw) and salinity (TDS), organic (BOD) and pathogen (FC) pollution.

Builds on the global hydrological model PCR-GLOBWB (https://github.com/UU-Hydro/PCR-GLOBWB_model) and the water temperature model DynWat (https://github.com/wande001/dynWat).

Version allows for DynQual to be run in both 'online' and 'offline' configurations: 
- Online version is one-way coupled with PCR-GLOBWB: in .ini file set 'offlineRun = False'
- Offline version can be forced with hydrological (baseflow, interflow, direct runoff; in m day-1) input data from any GHM: in .ini file set 'offlineRun = True'

Pollutant loadings can be calculated within DynQual (only in online configuration) or prescribed as a forcing (online & offline configurations):
- Pre-calculated pollutant loadings prescribed as a forcing: in .ini file set 'calculateLoads = False'
- Pollutant loadings calculated within DynQual: in .ini file set 'calculateLoads = True'. Here, additional socio-economic data must be provided to enable calculation of pollutant loadings.
