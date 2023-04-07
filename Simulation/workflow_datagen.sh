#!/bin/bash


# Simulation results
# Generate the synthetic data
for i in {1..9} 
do

# Experiment i
Rscript datagen.R $i 1 10
done
