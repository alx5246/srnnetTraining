A. Lonsberry
June 2016

(p)y(s)pike (w)ork (a)round README (pswaREADME)

DESCRIPTION
    There is a nice python module for comparing spike trains. The primary problem with the original python module is the
    use of cython that would only work with very particular constraints different from brian2. The workaround is to use
    the python code (which they included in the source) and not install or compile the module using anacoda or pip.

    The python module is pyspike. The primary utitily of using the module is the calculation of the ISI-distnace
    (interspike interval), the SPIKE distance, and the Spike synchronization. These are three different measures
    recently developed and seem rather interesting.

    To see a generic overview of these approaches see the scholarpedia articles:
    "Measures of spike train synchrony" : http://www.scholarpedia.org/article/Measures_of_spike_train_synchrony
    "SPIKE - distance" : http://www.scholarpedia.org/article/SPIKE-distance

    The Original papers wehere these measures have been develeoped are enumerate as follows:
    (1) ISI-distance : Kreuz T, Haas JS, Morelli A, Abarbanel HDI, Politi A (2007a). "Measuring spike train synchrony". J Neurosci Methods 165:151–161.
    (2) SPIKE - distance : Kreuz T, Chicharro D, Houghton C, Andrzejak RG, Mormann F (2013). "Monitoring spike train synchrony". JNeurophysiol 109:1457-72.
    (3) Spike synchronization : ???

FILES AND FOLDERS

    (1) Folders "compiledCppSim/", "examples/", "pyspike/"
    These came with pyspike. The "compiledCppSim" is supposed to be where everything gets compiled if everything was not
    wonky and problematic for our conditions. The folder "pyspike/" has all the source code and functions that we will
    be using. The remianing folder "examples/" has examples that came with pyspike