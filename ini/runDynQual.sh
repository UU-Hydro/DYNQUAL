#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -p thin
#SBATCH -N 1
#SBATCH -n 128

cd /home/ejones/

#Calculate loads
python new_DynQualModel/deterministic_runner.py new_ini_v2/setup_05min_online_calcLoad.ini &
python new_DynQualModel/deterministic_runner_offline.py new_ini/setup_05min_offline_calcLoad.ini &

#Force loads
python new_DynQualModel/deterministic_runner.py new_ini/setup_05min_online_forceLoad.ini &
python new_DynQualModel/deterministic_runner_offline.py new_ini/setup_05min_offline_forceLoad.ini &

wait