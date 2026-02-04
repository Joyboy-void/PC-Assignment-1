#!/bin/bash
#SBATCH --job-name=CILK_JOB    	      # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=email@iitdh.ac.in # Where to send mail	
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks-per-node=1           # Single task
#SBATCH --cpus-per-task=32            # 32 CPUs per task
#SBATCH --time=00:02:00               # Time limit hrs:min:sec
#SBATCH --output=cilk_test_%j.log   # Standard output and error log
#SBATCH --partition=cpu   # Queuename
pwd; hostname; date

set -ex
make clean
make all
make runSortCilk CILK_NWORKERS=32 SIZE=1000000000 CUTOFF=10
date
