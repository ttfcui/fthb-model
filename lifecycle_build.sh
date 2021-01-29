#!/bin/sh
#
#  Build lifecycle program using OpenMP with designated number of threads.
#  Compiler: Duke GRID interactive nodes now use gfortran 4.8.5 (2015XXXX).
#
#  Usage:
#  Run qhost to get threads on your node, and set as input argument to this script.
#  OPM_GET_NUM_THREADS() reports environment variable value, not actual number in use.
#
#  ./lifecycle_build.sh 1   :builds program to run with 32 threads

#SBATCH --job-name lifecycle     # create a short name for your job
#SBATCH --partition=scavenger
#SBATCH -c 8                     # number of threads requested for job
#SBATCH --mem-per-cpu=8G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=100:30:30
module load Python/2.7.11
module load Matlab/R2020a
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd model
make

echo "Made"

./lifecycle.out

# Discard executable.
#rm lifecycle.out
make clean
cd ..

echo "Program run complete. Output written to files."
echo " "
echo "Done!" > model/sbatch_out.txt
