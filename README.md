# Longitudinal-win-odds
R implementation for extension of probabilistic index model

lwo.R: main functions for the longitudinal win odds implementation, contain the modeling fitting function lwo() and summary(), predict() functions. 

functions.R: some user defined functions for simulation, including functions that simulates ordinal longitudinal data, calculates true win odds based on chosen simulation parameters, assesses simulation results.

LWO_analysis_example.qmd: R code for analyzing SID trial using the proposed longitudinal probabilistic index model. SID trial data is not publicly available therefore we provide a simulated dataset that mimics the SID trial in the script.

LWO_simulation.R: script for running the simulation.

LWO_results_assess.Rmd: The script sets up the simulation scenarios and evaluates the simulation results.

MC_simulation_shark.R: script for running the simulation in cloud computation.
