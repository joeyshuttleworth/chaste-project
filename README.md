# Extrapolation To Steady State Project

## Purpose
This project investigates the detection and computation of limit cycles in single cell cardiomyocyte models in Chaste. TestErrorMeasures.hpp compares several candidate error measures for limit cycle detection and the class SmartSimulation contains a method to extrapolate to these limit cycles.

## Building and Running
It is recommended to use Docker to install Chaste and its dependencies. See https://github.com/Chaste/chaste-docker for more details.

Then, clone this repository into Chaste's projects folder (~/projects or ~/Chaste/projects/ if you're using the docker image).

Next, build Chaste by following these instructions https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/UserProjects .

You are now ready to run the project (and the Continuous test pack) with ```ctest -j4 Continuous -V```. This will produce output in /tmp/Chaste which can be copied and saved. plot_extrapolation.py assumes that /tmp/chaste/ has been copied to ~/Chaste/data/ and produces the plots in the plots folder.
