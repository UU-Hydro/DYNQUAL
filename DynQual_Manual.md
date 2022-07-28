# DynQual Manual 

This file contains information for installing and running DynQual.


## Installation Guide

DynQual can be installed following the same steps as reuqired for PCR-GLOBWB2 (https://github.com/UU-Hydro/PCR-GLOBWB_model)

Please follow the following steps:

1. You will need a working Python environment, we recommend to install Miniconda, particularly for Python 3. Follow their instructions given at https://docs.conda.io/en/latest/miniconda.html. The user guide and short reference on conda can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html).

2. Get the requirement or environment file from this repository conda_env/pcrglobwb_py3.yml and use it to install all modules required (e.g. PCRaster, netCDF4) to run PCR-GLOBWB:

conda env create --name pcrglobwb_python3 -f pcrglobwb_py3.yml

The requirements file will create a environment named pcrglobwb_python3.

3. Activate the environment in a command prompt:

conda activate pcrglobwb_python3

4. Clone or download this repository. We suggest to use the latest version of the model, which should also be in the default branch.

git clone https://github.com/UU-Hydro/PCR-GLOBWB_model.git

This will clone PCR-GLOBWB into the current working directory.
