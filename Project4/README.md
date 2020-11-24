## Studies of phase transitions in magnetic systems

Here you can find the source codes for Project 4 in FYS4150.

If you want to run our programs, we recommend to study the 2x2 and 20x20 lattices, as this code should only take a few minutes to run. To do this, run

```bash
make all
```

which writes results to file and runs test functions. The test functions checks that the 2x2 lattice checks out with the analytic case, and tests the periodic boundary conditions.

To analyse the results from these lattice sizes, use the python script mean_vals by running the command 

```bash
python3 mean_vals.py
```

Now the terminal will ask you to type in various parameters which determine which file to open. First enter lattice size (either 2 or 20). Then temperature (1 for L = 2 and 1 or 2.4 for L = 20). The configuration to type in is either 0 (ordered initial configuration) or 0.5 (random intitial configuration). Only use 0.5 if you analyse the 20x20 lattice. When the script runs it will plot various quantities and write results to the terminal.


For the results from the paralellized programs, we have included the results in the data folder. To get the figures in the report, you only need to run plots.py.

To generate the results, openMP is required. The timer program can be run in less than a minute with
```bash
bash timerun.sh
```
but requires 8 cores by default to run the full comparison in the report. For the phase transition studies, this is the program mpimain.cpp which is compiled and run with 

```bash
mpic++ -O2 -o mmain.out mpimain.cpp
mpiexec mmain.out
```

It is currently set to the [2.2, 2.35] temperature interval with 4 million Monte Carlo cycles, which takes 2-3 hours to run on 8 cores. 
