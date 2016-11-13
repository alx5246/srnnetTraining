A. Lonsberry
SEPTEMBER 2016

README for folder spikeDataGenerator

DESCRIPTION

    Here we have methods and classes necessary to take analog data and convert to spike trains. IMPORTANTLY, no source
    data (ie 1D moving mass) is created here. The methods within here are specifically for transferrence of analog data
    into different domains (ie smooth analog into multiple receptive fields into spike trains0

FILES AND FOLDERS

    (1) decompDimensions.py
        This python module/file has methods that help break down smooth input domain into a set of overlapping Gaussian
        or similar functions. This is necessary for example, in the instance we want to break down the input over a
        single DOF into multi overlapping fields where each field represents a neuron.

    (2) analogToSpikes.py
        This python module/file has methods that take in time-varying rates and produces spike-trains. Included are also
        methods to generate spike trains based on Poisson Processes.
