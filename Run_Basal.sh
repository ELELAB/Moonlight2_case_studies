#!/bin/bash
NTHREADS=4
export OMP_NUM_THREADS=$NTHREADS
tsp -N $NTHREADS Rscript scripts/run_DMA.R
#Rscript scripts/run_moonlight.R
#tsp -N $NTHREADS -L astrid Rscript scripts/run_moonlight.R
#tsp -N $NTHREADS Rscript scripts/00_init.R
