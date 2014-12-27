""" Checks cpc index """

import numpy as np
from scipy import stats

import sys

try:
    c1 = sys.argv[1]
    c2 = sys.argv[2]
except:
    print "Define some files to compare!"

f1 = np.genfromtxt(c1,skiprows=1)
f2 = np.genfromtxt(c2,skiprows=1)


N = np.max((np.unique(f1.T[0]),np.unique(f1.T[1])))+1
#N=94

m1 = np.zeros((N,N))
m2 = np.zeros((N,N))


for e in f1:
    m1[int(e[0]),int(e[1])] = int(e[-1])
for e in f2:
    m2[int(e[0]),int(e[1])] = int(e[-1])
    
fl1 = m1.flatten()
fl2 = m2.flatten()

rr = np.array([fl1,fl2]).T
r2 = stats.linregress(rr.T[0],rr.T[1])[2]

CPC = 2*rr.min(axis=1).sum()/rr.sum()

print " ### R^2:%f | CPC for %s and %s is : %r ###" % (r2*r2,c1,c2,CPC)

