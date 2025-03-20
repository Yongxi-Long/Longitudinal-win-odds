# Longitudinal-win-odds
R implementation for extension of probabilistic index model

Analysis_example_LWO.qmd: R code for analyzing SID trial using the proposed longitudinal probabilistic index model. SID trial data is not publicly available therefore we provide a simulated dataset that mimics the SID trial in the script 

lwo.R: main functions for the longitudinal win odds implementation, contain the modeling fitting function lwo() and summary(), predict() functions 

functions.R: some user defined functions for simulation, including functions that simulates ordinal longitudinal data, calculates true win odds based on chosen simulation parameters, assesses simulation results.

MC_simulation.R: script for running the simulation

MC_simulation_shark.R: script for running the simulation in cloud computation
