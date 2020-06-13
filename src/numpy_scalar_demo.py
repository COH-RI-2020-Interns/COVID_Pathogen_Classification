import numpy as np

x = np.arange(0,100)

def linear(slope, intercept):
    return slope*x + intercept

y = linear(2, 4)

y
