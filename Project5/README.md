## Diffusion equation and modeling the litosphere

Here you can find the source codes for Project 5 in FYS4150.

If you want to run the programs simple type make all. Here, the functions are in diffusion.cpp, the 1d and 2d general case are run in main.cpp and the litosphere data in isomain.cpp

Then running plots.py will call functions for the plots in the report. Note that the litosphere plots are switched off by default
so that python does not generate too many plots at once. This can be changed at the bottom of plots.py.

Also note that the generated litosphere is set to run with 0.025 delta x instead of the 0.01 used in the report so the program can be run quickly for testing. The results are similar still.

If wanted, it shouldn't take too long to run with 0.01 instead. Then dx must also be changed to 0.01 in the plotlitosphere function in plots.py.
Although 0.01 shouldn't take too much time either.

