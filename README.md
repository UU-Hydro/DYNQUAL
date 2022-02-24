# DYNQUAL
Water quality model for simulating water temperature and salinity (TDS), organic (BOD) and pathogen (FC) pollution.

Builds on the global hydrological model PCR-GLOBWB (https://github.com/UU-Hydro/PCR-GLOBWB_model) and the water temperature model DynWat (https://github.com/wande001/dynWat).

Version allows for Dynqual runs both 1) online vs. offline and 2) with pre-calculated loadings vs. loadings calculated in model run.
- Online version is coupled with PCR-GLOBWB: in .ini file set 'offlineRun = False'
- Offline version can be forced with hydrology (baseflow, interflow, direct runoff): in .ini file set 'offlineRun = True'
- (Pre-calculated) loadings can be prescribed as a forcing: in .ini file set 'calculateLoads = False'
- Loadings can also be calculated in loop: in .ini file set 'calculateLoads = True'
