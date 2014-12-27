import os
import numpy as np

pref = 'N1000'

z = [e for e in os.listdir('../') if pref in e and 'list' in e and 'ens' in e]

ori = np.genfromtxt('../../tests/strengths.str',skiprows=1)
ll = np.genfromtxt('../'+z[0],skiprows=1)

import pylab as pl

