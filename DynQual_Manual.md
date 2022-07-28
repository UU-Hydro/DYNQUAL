# DynQual Manual 

This file contains information for installing and running DynQual.


## Installation Guide

DynQual can be installed following similar steps as for PCR-GLOBWB2 (https://github.com/UU-Hydro/PCR-GLOBWB_model)

1. You will need a working Python environment, we recommend to install Miniconda, particularly for Python 3. Follow their instructions given at https://docs.conda.io/en/latest/miniconda.html. The user guide and short reference on conda can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html).

2. Get the requirement or environment file from this repository conda_env/pcrglobwb_py3.yml and use it to install all modules required (e.g. PCRaster, netCDF4) to run DynQual:

`conda env create --name pcrglobwb_python3 -f pcrglobwb_py3.yml`

3. Activate the environment in a command prompt:

`conda activate pcrglobwb_python3`

4. Clone or download the DynQual repository into the current working directory.

`git clone https://github.com/UU-Hydro/DYNQUAL.git`



## DynQual configuration .ini file

For running DynQual, a configuration .ini file is required. In this file, you provide all the necessary information related to your run (e.g. time period, study extent, etc) and paths to your input data (e.g. climatological forcing, land cover, etc).

Some example configuration .ini files are provided in the 'ini' directory, both for DynQual runs that are one-way coupled to PCR-GLOBWB2 ("Online") and for runs using hydrological input as a forcing ("Offline"). We provide all the required input data files for the Rhine-Meuse basin which can be downloaded to your local machine (~6GB), facilitating a self-contained example run (**LINK**). Input data for DynQual runs at the global extent are availiable on the OPeNDAP server (**LINK**) - allowing for users to run DynQual for any land area **without** needing to download the input files (which exceed 250 GB).     

Some adjustments must be made:
- `inputDir =`  : this must be set to the directory where the input data is stored.
- `outputDir =` : this must be set to a directory you can access.
- `cloneMap =`  : please make sure this is stored locally in your computer. The cloneMap file defines the spatial resolution and extent of your study area and must be in the pcraster format. Some examples are given in this [repository](https://github.com/UU-Hydro/PCR-GLOBWB_model/blob/master/clone_landmask_maps/clone_landmask_examples.zip) 

Other adjustments are optional (e.g.):
- `startTime =` and `endTime =` : denote your simulation period.
- `precipitationNC =` : denotes the location of your precipitation input data file, relative to `inputDir`.
- `outDailyTotNC =` : denote which variables you would like written to netcdfs at daily resolution to `outputDir/netcdf/`. 

Some DynQual specific adjustments can also be made (e.g.):
- `calculateLoads =` : Set to TRUE for simulating pollutant loadings, FALSE for prescribing (pre-calculated) loadings. Note: TRUE is only avaliable for **online** DynQual runs.
- `loadsPerSector =` : Set to TRUE for segregation of pollutant loadings per sector (i.e. for attribution), FALSE for combined loadings only. Note: TRUE is only avaliable for **online** DynQual runs.


When adjusting input data files, remember to **check your units!**


## Running DynQual

Ensure the correct conda environment in a command prompt: `conda activate pcrglobwb_python3`

Navigate to the DynQual model directory (*DynQualModel*). You can start a DynQual run using the following commands:

For **online** run: `python deterministic_runner.py <ini_configuration_file>`

For **offline** run: `python deterministic_runner_offline.py <ini_configuration_file>`

where <ini_configuration_file> is your configuration file of DynQual.
