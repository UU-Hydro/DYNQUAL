#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -p normal
#SBATCH -N 1
cd /home/ejones/
python off_DynQualModel/deterministic_runner.py off_ini/setup_05min.ini