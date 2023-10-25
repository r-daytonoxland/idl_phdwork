from matplotlib import pyplot as pyplot
import numpy as np 

def readall(fstart, data)

types = [intensity, blueshift, temperature]

fnames = []
[fnames.append('~/lib/' + fstart + types[i]) for i in len(2)]

data = np.array[3]
[data[i] = open(fnames[i], 'r') for i in len(2)]
