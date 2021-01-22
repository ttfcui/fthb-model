#!/bin/bash
#
#  Build lifecycle program using OpenMP with designated number of threads.
#  Compiler: Booth GRID interactive nodes now use gfortran 4.8.5 (20150623).
#
#  Usage:
#  Run qhost to get threads on your node, and set as input argument to this script.
#  OPM_GET_NUM_THREADS() reports environment variable value, not actual number in use.
#
#  ./lifecycle_build.sh 32   :builds program to run with 32 threads

make

#  Run program with number of threads defined by user in shell argument.
#  Defaults to 40 threads (for Eric Zwick's server).
#
NUM_THREADS=${1:-40}
echo " "
echo "Run with $NUM_THREADS threads."
export OMP_NUM_THREADS=$NUM_THREADS
./lifecycle.out

# Discard executable.
#rm lifecycle.out

echo "Program run complete. Output written to files."
echo " "
