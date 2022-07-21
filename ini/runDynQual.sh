#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -p thin
#SBATCH -N 1
#SBATCH -n 32

cd /home/ejones/

###Online DynQual runs with W5E5 forcing (calcLoads == True or False)
#python DynQualModel/deterministic_runner.py ini/DynQual_05min.ini & #single landmask
#python DynQualModel/parallel_pcrglobwb_runner.py ini/DynQual_05min.ini & #parallel

###Offline DynQual runs with W5E5 forcing (calcLoads == False only.)
#python DynQualModel/deterministic_runner_offline.py ini/DynQual_05min_offline.ini & #single landmask
#python DynQualModel/parallel_pcrglobwb_runner_offline.py ini/DynQual_05min_offline.ini & #parallel

#wait