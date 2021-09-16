#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -p normal
#SBATCH -N 1
cd /home/ejones/
python DynQualModel/deterministic_runner.py historical_run/setup_test.ini