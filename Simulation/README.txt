This directory contains the code to reproduce the numerical results in the Section 5 of paper, The EAS
approach for grouped variable selection in multivariate linear regression model.

Please contact Salil Koner at salil.koner@duke.edu for any help or questions.

See workflow.sh for line-by-line Unix command line code for reproducing the 
results.  Running the command,

bash workflow.sh 

will reproduce all numerical results presented in the manuscript, but will take
many days.  If it is desired to reproduce the results, it is best to follow the
instructions in workflow.sh and parallelize accordingly on a computing cluster.
Setting fewer MCMC steps will also speed up computations, but will yield 
differing results.

Note: Following R packages must be installed in the system in order to successfully
run the code. 
"matrixsampling", "tidyverse", "Rfast", "MASS", "glmnet",
"MBSP", "rrpack", "spls", "MSGLasso", "MBSGS", Matrix

The directory Results contains six separate sub-directories named as 
Experiment_1-Experiment_6 where the output for each dataset will be stored.
Any such output will be named as output_{method}_dsn{dataset_number}.Rdata. 

When the output for all the datasets (we take 1000 datasets in our paper) is created,
the last part of workflow.sh will create the summary measures presented in Table 1 (or Table 2)
of the paper and save it at the same folder with the name Summary_all.Rdata.