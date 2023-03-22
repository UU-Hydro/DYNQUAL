#!/bin/bash

# we reserve one node of eejit, one contains 96 cores
#SBATCH -N 1

# we use all cores
#SBATCH -n 96

# activate the following if we want to reserve the entire node (which is the case if -n 96)
#SBATCH --exclusive

# time limit
#SBATCH -t 480:00:00

# the partition name, always use "defq" for eejit (we do not need gpu)
#SBATCH -p defq

# job name
#SBATCH -J dynqual


#~ # mail alert at start, end and abortion of execution
#~ #SBATCH --mail-type=ALL
#~ # send mail to this address
#~ #SBATCH --mail-user=XXXX@gmail.com


# please set where you stored DYNQUAL scripts
SCRIPT_FOLDER=$HOME/DYNQUAL/DynQualModel/

# please 
INI_FILE=$HOME/DYNQUAL/ini/global_run/DynQual_05min_global.ini
# - Note that I noticed that Edward still use quite old parameter files from PCR-GLOBWB. Shall we update this? If yes, I'll allocate some hours around next week. 

# we activate the correct conda environment on eejit and many other settings
# - abandon any existing PYTHONPATH (recommended, if you want to use miniconda or anaconda)
unset PYTHONPATH
# - load anaconda from the module
module load opt/all
module load anaconda3/2021.11
# - activate the correct conda env
source activate /quanta1/home/sutan101/opt/miniconda3/envs/pcrglobwb_python3_2022-10-17
#~ # - use at least 8 workers
#~ export PCRASTER_NR_WORKER_THREADS=8


# we go to the folder where the folder script
cd ${SCRIPT_FOLDER}

# global online run with parallelization
python parallel_pcrglobwb_runner.py ${INI_FILE}
