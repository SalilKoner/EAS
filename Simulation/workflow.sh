#!/bin/bash

ndata=1000 # generate 1000 datasets for each of the experiments
echo $ndata

# Simulation results
# Generate the synthetic data and implement the estimation procedures
# This for-loop is most efficiently run in parallel on a computing cluster
# In parallel it will take approximately 1 day to finish the simulation 

# exper 1 to 6 is for the simulation replicated as in Bai and Ghosh (2018)
# For the design in Table 1 only.

#data generation
for i in {1..6} 
do
# Experiment i : This will created a subfolder named Experiment_$i and put all the 
Rscript datagen.R $i 1 $ndata
done

# Running all the 9 methods (EAS + MBSP + MLasso + MSGLasso + SPLS + SRRR + MSRL + MBGLSS)
for exper in {1..6}
do
for (( dataset=1; dataset<=$ndata; dataset++ ))
do
Rscript sim_EAS_BIC_and_CV.R $dataset 8 $exper 0.05 10 24 10 # EAS BIC and CV
Rscript sim_MBSP_MLasso_MSGLasso.R $dataset 1 $exper # MBSP + MLasso + MSGLasso
Rscript sim_SPLS_SRRR.R $dataset 1 $exper # SPLS + SRRR
Rscript sim_MSRL.R $dataset 1 $exper # MSRL
Rscript sim_MBGLSS.R $dataset 1 $exper # MBGLSS, only for Table 1
done
done

# Produce the summary statistics in Table 3 from the 1000 generated datasets
for exper in {1..6}
do
Rscript analysis_file_sim.R $exper $ndata # assuming that 1000 simulation datasets is used
done

# exper 7 to 9 is for the simulation replicated as requested by AE
# They corresponds to large q, i.e. Table 2 only.
# MBGLSS is not implemented for large q as it takes more than 4 days for a single data
# Caution: The entire simulation for Table 2 will take more than a month, if
# it is not run in parallel. 

#data generation
for i in {7..9} 
do
# Experiment i : This will created a subfolder named Experiment_$i and put all the 
Rscript datagen.R $i 1 $ndata
done

# Running all the 9 methods (EAS + MBSP + MLasso + MSGLasso + SPLS + SRRR + MSRL + MBGLSS)
for exper in {7..9}
do
for (( dataset=1; dataset<=$ndata; dataset++ ))
do
Rscript sim_EAS_BIC_and_CV.R $dataset 8 $exper 0.05 10 24 10 # EAS BIC and CV
Rscript sim_MBSP_MLasso_MSGLasso.R $dataset 1 $exper # MBSP + MLasso + MSGLasso
Rscript sim_SPLS_SRRR.R $dataset 1 $exper # SPLS + SRRR
Rscript sim_MSRL.R $dataset 1 $exper # MSRL
done
done

# Produce the summary statistics in Table 4 from the 1000 generated datasets
for exper in {7..9}
do
Rscript analysis_file_sim.R $exper $ndata # assuming that 1000 simulation datasets is used
done
