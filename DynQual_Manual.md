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

For running DynQual, a configuration .ini file is required. Some example configuration .ini files are provided in the 'ini' directory, both for DynQual runs that are one-way coupled to PCR-GLOBWB2 ("Online") and for runs using hydrological input as a forcing ("Offline").




## Running DynQual

Ensure the correct conda environment in a command prompt: `conda activate pcrglobwb_python3`

Navigate to the DynQual model directory (*DynQualModel*). You can start a DynQual run using the following commands:

For **online** run: `python deterministic_runner.py <ini_configuration_file>`

For **offline** run: `python deterministic_runner_offline.py <ini_configuration_file>`

where <ini_configuration_file> is your configuration file of DynQual.
