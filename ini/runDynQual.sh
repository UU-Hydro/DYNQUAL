#!/bin/bash
#SBATCH -t 10:00
#SBATCH -p thin
#SBATCH -N 1
#SBATCH -n 128

cd /home/ejones/

#Online DynQual runs with W5E5 forcing (calcLoads == True)
#python DynQualModel/deterministic_runner.py ini/DynQual_05min.ini & #single landmask
python DynQualModel/parallel_pcrglobwb_runner.py ini/DynQual_05min.ini & #parallel

#Offline DynQual runs with W5E5 forcing (calcLoads == False)
#python DynQualModel/deterministic_runner_offline.py ini/DynQual_05min_offline.ini & #single landmask
#python DynQualModel/parallel_pcrglobwb_runner_offline.py ini/DynQual_05min_offline.ini & #single landmask

wait