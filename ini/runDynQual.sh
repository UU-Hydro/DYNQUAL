#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH -p thin
#SBATCH -N 1
#SBATCH -n 64

cd /home/ejones/

#DynQual runs with W5E5 forcing
python DynQualModel/deterministic_runner.py ini/DynQual_05min.ini &
python DynQualModel/deterministic_runner_offline.py ini/DynQual_05min_offline.ini &

wait