# 3D Lake Nutrient Flux Network Modelling

This repository contains the code for the paper "Exploring Nonlinear Responses of Lake Nutrients to Restoration Measures: A Three-Dimensional Flux Network Modelling Approach".

## EFDC model
The parameter configurations of the EFDC model can be found in the SI.

## Meta file
All directories and coefficients should be assigned in the python script file ```meta.py```.

## Calculation of the TN, TP and Chla
Run the Python script ```water_quality_cal.py``` to calculate the TN, TP and Chla. 
Set the ```scenarios``` to compute values for certain scenarios. The default is to calculate all the scenarios.

## Calculation of the nutrient storages and fluxes
Run the Python script ```flux_cal.py``` to calculate the nutrient storages and fluxes. 
Set the ```scenarios``` to compute values for certain scenarios. The default is to calculate all the scenarios.

## Building the nutrient networks
Run the Python script ```flux_network_extract.py``` to build the nutrient networks based on nutrient storages and fluxes. 
Set the parameter ```frequency``` can control the network aggregation periods. Available values: "m" for monthly average and "d" for daily average.

## ENA analysis
Ecological network analysis (**ENA**) is a widely used method in the food web assessments. Here, it was used to provide systematic information on N and P cycling in lake ecosystems.
Set the correct directory of the network files and then run the R script ```flow_analysis.R``` to get the values of indicators and visualization results.

## Visualization
Files to create the figures in the main text:
- ```water_quality_calibration_maintext.py```: Fig. 2;
- ```flux_stack_maintext.py```: Fig. 3;
- ```water_quality_scenario_maintext_si.py```: Figs. 4a, 4b, 5a, and 5b;
- ```flux_network_maintext.py```: Figs. 4c, 4d, 5c, and 5d;
- ```nonlinearity_maintext.py```: Fig. 7.

Files to create the figures in the SI:
- ```water_quality_calibration_si.py```: Fig. S1-S7;
- ```water_quality_scenario_maintext_si.py```: Figs. S8, S9, and S10;
- ```flux_network_si.py```: Figs. S11-S17.


