# A. Lonsberry
# June 2017
#
# DESCRIPTION
#   This is a sandbox to help figure out how the cosine functions based conditional intensity functions work...

import numpy
import math
import matplotlib.pyplot as plt

#peakPositions = numpy.array([.01,.02,.03,.04,.05,.06,.07,.08])
#c = .5
peakPositions = numpy.array([0., 5., 10., 15., 20.])
c = 100.

print(peakPositions[-1])

a = math.pi * (peakPositions.size - 1.) * (2.*math.log((peakPositions[-1]+c)/(peakPositions[0]+c)))**-1.

timeArray = numpy.linspace(start=0.0, stop=25., num=1000)
cosineFuncValues = numpy.zeros((peakPositions.size, timeArray.size))

plt.figure(1)

for i in range(peakPositions.size):
    phi = a * math.log(peakPositions[0]+c) + (math.pi/2)*((i+1)-1) #MAKE SURE WE ADD 1 to i, because it starts at 0 and not 1
    for j, t in enumerate(timeArray):
        if math.fabs(a*math.log(t+c)-phi) <= math.pi:
            cosineFuncValues[i, j] = .5 * (1 + math.cos( a*math.log(t+c)-phi) )

    plt.plot(timeArray, cosineFuncValues[i, :])
plt.show()


# Make some overlapping Guassin Vectors

