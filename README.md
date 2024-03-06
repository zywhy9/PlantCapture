# Reproducibility Instructions

## Notes

### 1. Simulation study
* The simulation study is performed for three models in two sample size scenarios, small city (M,H) = (15, 150) for the first two models and (M,H) = (30, 300) for the last one, and large city (M, H) = (100, 1500).

### 2. Real data analysis
* We provide the analysis of the S-Night dataset which is reconstructed from literature.

### 3. Additional Notes
* Some codes need manually edit to fit small city scenario or large city scenario.

## Steps

### 1. Choose the model
The three models introduced in the study are located in separate folders, named "Basic Model", "Partial Identification" and "Heterogeneity Between Sites". Open the folder for the model of interest.

### 2. Simulate datasets
In the folder, there is a file called "Data_Simulation.R", which provides the codes to simulate dataset. There is a small city setting and a large city setting. You can choose the one you prefer or customize the values of parameters.

### 3. Analyze datasets
After simulating datasets, you can choose the methods to analyze datasets. "BNA.R" includes the code with Bayesian normal approximation method, "MLE.R" includes the code with MLE method, "UP_analysis.R" includes the code with uncertainty propagation method, and "MCMC_analysis.R" includes the code with MCMC method. "MCMC_evaluation.R" and "UP_evaluation.R" provide the summary statistics of the estimators.