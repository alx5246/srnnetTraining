
DESCRIPTION
One major problem using brian is that the wonderful standalone complilation option does not work for loops. That is if
one would like to simulate repeatedly, with different parameter values, you have to recompile every time. Also there
are some special things you must do even in this situation.

EXAMPLE 1: workAround_0_0.py
    Here we simply show a way to compile (deemed the 'standalone' in brian2 nomenclature) a simulation. The way it seems
    to work is that the code, written in python, is then translated to c++ and compiled with visual-studio, and then run,
    There are some caveats to using the method as indiciated in the brian2 manual online.

EXAMPLE 2: workAround_1_0.py & workAround_1_1.py
    Here there are two files. In one file we have a standalone brian2 simulation imbedded inside a function. In the
    second file we have a mechanim to call the other file and run the simulation. In these files we demonstarte two
    methods to run simulations iteratively (this is done by commenting out parts of the code in each file). The first
    way is to create a for-loop which calls the simualtion iteratively. The second method is to use parallel processing
    to run multiple iterations at once, which is the better option for us.