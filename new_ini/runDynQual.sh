#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -p thin
#SBATCH -N 1
#SBATCH -n 128

cd /home/ejones/

#Calculate loads
python new_DynQualModel/deterministic_runner.py new_ini/setup_05min_GFDL_online_calcLoad.ini &
python new_DynQualModel/deterministic_runner_offline.py new_ini/setup_05min_GFDL_offline_calcLoad.ini &

#Force loads
python new_DynQualModel/deterministic_runner.py new_ini/setup_05min_GFDL_online_forceLoad.ini &
python new_DynQualModel/deterministic_runner_offline.py new_ini/setup_05min_GFDL_offline_forceLoad.ini &

wait