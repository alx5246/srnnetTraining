A. Lonsberry
May 2016

(b)rian (s)tand (a)lone (w)ork (a)rounds README (bsawaREAME)

DESCRIPTION

    One major problem using brian is that the wonderful standalone complilation option is somewhat fussy, and can be
    difficult to implement. One known problem is that we CANNOT at the point in time, compile code, and then send that
    code input to run. Rather, we will have to recomplie every time, which is still faster than using native python.

FILES AND FOLDERS

    (1) workAround_0_0.py
    Here we simply show a way to compile (deemed the 'standalone' in brian2 nomenclature) a simulation. The way it seems
    to work is that the code, written in python, is then translated to c++ and compiled with visual-studio, and then run,
    There are some caveats to using the method as indiciated in the brian2 manual online.

    (2) workAround_1_0.py and workAround_1_1.py
    Here there are two files. In one file we have a standalone brian2 simulation imbedded inside a function. In the
    second file we have a mechanim to call the other file and run the simulation. In these files we demonstarte two
    methods to run simulations iteratively (this is done by commenting out parts of the code in each file). The first
    way is to create a for-loop which calls the simualtion iteratively. The second method is to use parallel processing
    to run multiple iterations at once, which is the better option for us.