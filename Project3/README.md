# Project 3 FYS4150
## Building a model for the solar system using ordinary differential equations
A class-based code in C++ to simulate the solar system. Initial values are extracted from initial_raw.txt and placed in initial.txt with the python script make_initial.py. The program uses velocity verlet method and forward euler to evolve system.

To run simulation, run the makefile

> make all

then analyse data in results.py. In this script simply call functions of what you want to plot, then

> python3 results.py
