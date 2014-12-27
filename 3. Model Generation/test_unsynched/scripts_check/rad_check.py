""" Checks rad model """

import numpy as np
import numbapro as nb


import sys

fix_batty = True # Fix limited size?

try:
    c1 = sys.argv[1]
    c2 = sys.argv[2]
    c3 = sys.argv[3]
    c4 = sys.argv[4]
except:
    print "Define original .tr file, strength file, dist_file, rad_simulation!"


@nb.autojit
def compute_sij(d,s,N):
	ss = np.zeros((N,N))
	for i in range(N): # origin node
		for j in range(N): # dest node
			sij = 0
			for k in range(N): # dest node
				if d[i][k]<d[i][j] and k!=i and k!=j:
					sij+=s[k][1]
			ss[i][j]=sij
		if i%10 == 0:
			print "Done %d nodes %f frac" % (i+1,(i+1)*1./N)
	return ss


### imput ###
w = np.genfromtxt(c1,skiprows=1)
wm = np.genfromtxt(c4,skiprows=1)
s = np.genfromtxt(c2,skiprows=1)
s=s.T[1:].T
d = np.genfromtxt(c3,skiprows=1)


N = int(np.max((np.unique(w.T[0]),np.unique(w.T[1]))))+1

adj = np.zeros((N,N))
adj_m = np.zeros((N,N))
sij = np.zeros((N,N))
dist = np.zeros((N,N))

for e in w:
    adj[int(e[0]),int(e[1])] = int(e[-1])
for e in wm:
    adj_m[int(e[0]),int(e[1])] = int(e[-1])

for e in d:
    dist[int(e[0]),int(e[1])] = e[-1]

### compute radiation fake ###
sij = compute_sij(dist,s,N)

if fix_batty:
	fix = (1.-s.T[1]*1./s.T[1].sum())
else:
	fix=1.
adj_fake = np.einsum('i,j',s.T[0]*s.T[1]/fix,s.T[1])
den = (sij+(np.ones((N,N))*s.T[1]).T)*(sij+np.ones((N,N))*s.T[1]+(np.ones((N,N))*s.T[1]).T)

#adj_fake = np.einsum('i,j',s.T[0]*s.T[0]/fix,s.T[1])
#den = (sij+(np.ones((N,N))*s.T[0]).T)*(sij+np.ones((N,N))*s.T[1]+(np.ones((N,N))*s.T[0]).T)

adj_fake /=den
adj_fake = adj_fake - adj_fake*np.eye(N) # exclude self-loops





### compute sorensen ###
from scipy import stats
fl1 = adj.flatten()
fl2 = adj_fake.flatten()
fl3 = adj_m.flatten()

rr = np.array([fl1,fl2]).T

r2 = stats.linregress(rr.T[0],rr.T[1])[2]
r2 = r2*r2

CPC = 2*rr.min(axis=1).sum()/rr.sum()

print " ### R^2:%f CPC for %s and Radiation average is : %r ###" % (r2, c1,CPC)

rr = np.array([fl3,fl2]).T

r2 = stats.linregress(rr.T[0],rr.T[1])[2]
r2 = r2*r2

CPC = 2*rr.min(axis=1).sum()/rr.sum()

print " ### R^2:%f CPC for %s and Radiation average is : %r ###" % (r2,c4,CPC)


rr = np.array([fl1,fl3]).T

r2 = stats.linregress(rr.T[0],rr.T[1])[2]
r2 = r2*r2

CPC = 2*rr.min(axis=1).sum()/rr.sum()

print " ### R^2:%f CPC for %s and %s is : %r ###" % (r2,c1,c4,CPC)
