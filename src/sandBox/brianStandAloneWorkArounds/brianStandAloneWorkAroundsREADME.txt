
DESCRIPTION
One major problem using brian is that the wonderful standalone complilation option does not work for loops. That is if
I want to loop through variables or parameters and have the 'run()' inside said loop, the compiler will break. Thus
we need some reasonable work arounds so we do not always have to recomplie before running.