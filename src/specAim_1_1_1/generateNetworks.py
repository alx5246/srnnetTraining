# A.Lonsberry
# November 2016
#
# DESCRIPTION
#   Here I am going to have methods to create, test, and alter (evolutionary/genetic) networks that will be used in the
#   experiments. To be very clear, here I will essentially create and save the networks. These networks can then be
#   loaded into a script to do all the experiments.
#
#   We seperate network generation and experiments for a couple of reasons. Reason (1), the networks as randomly
#   instantiated may be unstable, and require iterations to change parameters to get within the acceptable range of
#   functionality. Reason (2), we may want to reuse networks for different experiments, so it makes sense to generically
#   make them and save them for reuse.

import matplotlib.pyplot as plt
import numpy

# I had a difficult time getting things to import, so now I have however, go about fixing this by going up a directoy
# which I think I can do in the following way.
import os
import sys
import pickle
sys.path.append(os.path.abspath("../"))