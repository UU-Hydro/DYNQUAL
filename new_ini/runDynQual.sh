#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -p normal
#SBATCH -N 1
cd /home/ejones/
python new_DynQualModel/deterministic_runner.py new_ini/setup_05min.ini