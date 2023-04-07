#!/bin/bash


# Simulation results
# Generate the synthetic data and implement the estimation procedures
# This for-loop is most efficiently run in parallel on a computing cluster
# In parallel it will take approximately 1 day to finish the simulation 
for i in {1..2} 
do
for j in {1..2}
do
Rscript sim_EAS_BIC_and_CV.R $i 8 $j 0.05 10 24 10
Rscript sim_MBSP_MLasso_MSGLasso.R $i 1 $j
Rscript sim_SPLS_SRRR.R $i 1 $j
Rscript sim_MSRL.R $i 1 $j
Rscript sim_MBGLSS.R $i 1 $j
done

#for j in {7..9}
#Rscript sim_EAS_BIC_and_CV.R $i 8 $j 0.05 10 24 10
#Rscript sim_MBSP_MLasso_MSGLasso.R $i 8 $j
#Rscript sim_SPLS_SRRR.R $i 8 $j
#Rscript sim_MSRL.R $i 8 $j
#done

done


# Produce the summary statistics in Table 1 from the 1000 generated datasets
#for i in {1..9}
#do
#Rscript analysis_file_sim.R $i
#done
