# DynQual v1.1

The dynamical surface water quality model (DynQual v1.1) is a large scale water quality model for simulating daily dissolved oxygen (DO), water temperature (Tw), total dissolved solids (TDS), biological oxygen demand (BOD) and fecal coliform (FC) at 5 arc-minute spatial resolution with a daily timestep. 

The model builds on the global hydrological model PCR-GLOBWB2 (Sutanudjaja et al., 2018; https://github.com/UU-Hydro/PCR-GLOBWB_model) and the water temperature model DynWat (Wanders et al., 2019; https://github.com/wande001/dynWat), allowing for both the quantification of pollutant emissions from different human water use activities and their subsequent routing through the surface water network. Water temperature is simulated by solving the surface water energy balance, dissolved oxygen is simulated by solving the Streeter-Phelps equation, while TDS, BOD and FC concentrations are simulated accounting for both the dilution capacity and in-stream decay processes.

We offer two options for running DynQual:
1)	One-way coupled with PCR-GLOBWB2; or
2)	In a stand-alone configuration with user-defined hydrological input from any land surface or hydrological model (i.e. surface runoff, interflow and baseflow).

In both model configurations, pollutant loadings can be prescribed directly (akin to a forcing). Alternatively, when running DynQual coupled with PCR-GLOBWB2, pollutant loadings can be simulated within the model runs by providing additional socio-economic input data.

This repository holds an installation guide and the model scripts for running DynQual (DynQual_Manual.md).

An example set-up that includes all necessary input data for the Rhine-Meuse basin (~ 6GB) is provided through Zenodo (https://doi.org/10.5281/zenodo.13895791). The related .ini files necessary for running DynQual for the Rhine basin are also provided through Zenodo. 

Additionally, a global set-up is also provided that links to input files available on the OPeNDAP server (https://opendap.4tu.nl/thredds/catalog/data2/pcrglobwb/catalog.html). This allows users to access input files from a remote server and thus perform DynQual runs without needing to download the input files locally (> 250 GB). The related .ini files necessary for running DynQual linked to te OPenDAP server are provided through this GitHub. 


### Main references/ papers:
**DynQual**: Jones, E.R., Bierkens M. F. P., Wanders, N., Sutanudjaja, E. H., van Beek, R., van Vliet, M. T. H. (2023). DynQual v1.0: A high-resolution global surface water quality model. Geosci. Model Dev., 16, 4481–4500, https://doi.org/10.5194/gmd-16-4481-2023.

**PCR-GLOBWB2**: Sutanudjaja, E. H., van Beek, R., Wanders, N., Wada, Y., Bosmans, J. H. C., Drost, N., van der Ent, R. J., de Graaf, I. E. M., Hoch, J. M., de Jong, K., Karssenberg, D., López López, P., Peßenteiner, S., Schmitz, O., Straatsma, M. W., Vannametee, E., Wisser, D., and Bierkens, M. F. P. (2018). PCR-GLOBWB 2: a 5 arcmin global hydrological and water resources model, Geosci. Model Dev., 11, 2429-2453, https://doi.org/10.5194/gmd-11-2429-2018.

**DynWat**: Wanders, N., van Vliet, M. T., Wada, Y., Bierkens, M. F., & van Beek, L. P. (2019). High‐resolution global water temperature modeling. Water Resources Research, 55(4), 2760-2778, https://doi.org/10.1029/2018WR023250
