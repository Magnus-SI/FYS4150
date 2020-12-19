## Diffusion equation and modeling the litosphere

Here you can find the source codes for Project 5 in FYS4150.

If you want to run the programs simple type make all. Here, the functions are in diffusion.cpp, the 1d and 2d general case are run in main.cpp and the litosphere data in isomain.cpp

Then running plots.py will call functions for the plots in the report. Note that the litosphere plots are switched off by default
so that python does not generate too many plots at once. This can be changed at the bottom of plots.py. Also note that the litosphere plot
is set to read the 0.01 delta x used in the report, while in isomain.cpp delta x is set to 0.025 so the program can be run quickly.
Although 0.01 shouldn't take too much time either.

