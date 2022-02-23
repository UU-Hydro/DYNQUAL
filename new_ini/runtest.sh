#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -p thin
#SBATCH -N 1
#SBATCH -n 32

cd /home/ejones/

python new_DynQualModel/deterministic_runner_offline.py new_ini/test.ini