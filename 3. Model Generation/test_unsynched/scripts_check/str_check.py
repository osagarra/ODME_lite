""" Checks strengths """

import numpy as np

import sys

try:
    c1 = sys.argv[1]
    c2 = sys.argv[2] 
    c3 = sys.argv[3] ## ensemble
except:
    print "Define some files to compare!"

f1 = np.genfromtxt(c1,skiprows=1)
f2 = np.genfromtxt(c2,skiprows=1)
f3 = np.genfromtxt(c3,skiprows=1)


from scipy import stats

r = stats.pearsonr(f1.T[1],f2.T[3])
rr = stats.pearsonr(f1.T[1],f3.T[5])

r2 = stats.pearsonr(f1.T[2],f2.T[10])
rr2 = stats.pearsonr(f1.T[2],f3.T[-10])

print "### Strength checker for %s %s ####" % (c1,c2)
print "R between str: out  %r | in %r for 1samp" % (r,r2)
print "R between str: out  %r | in %r for average" % (rr,rr2)