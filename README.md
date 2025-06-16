The repository contains all the datafiles and scripts needed to replicate the simulation study. 

The scripts Simulation_sup1.R and Simulation_sup2.R contain the main code for simulation. They call the scripts GMSE.R, Vdesign.R and PEM_functions.R, which contain the functions needed to compute the accuracy estimators. Note that the script for GMSE is adapted from the code developed by the Italian Statistical agency (https://github.com/nina-DL/GMSE/tree/main). The scripts are specific to the simulation study and need further modification if one wishes to apply them to a different dataset. 

The datafile "samplonie.csv" contains the dataset for the simulation.

Scripts Simulation_analysis.R and Plots_sup1.R/Plots_sup2.R contain the code to compute the benchmark estimators and the plots in the results section respectively.
