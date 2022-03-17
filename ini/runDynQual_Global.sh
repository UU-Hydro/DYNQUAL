#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -p thin
#SBATCH -N 1
#SBATCH -n 128

cd /home/ejones/

#Global DynQual runs with W5E5 forcing
python DynQualModel/parallel_pcrglobwb_runner.py ini/DynQual_05min_Global.ini &
#python DynQualModel/parallel_pcrglobwb_runner.py ini/DynQual_05min_offline_Global.ini &

wait