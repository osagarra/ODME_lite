"""
    Examples for optimization of gamma, x,y for the gravity case
"""
# imports #
import multi_edge_fitter as mw
import numpy as np
import sys
import os

pref = 'lagrange_mult'

try:
    os.mkdir(pref)
except OSError:
    pass

#########################
######## Inputs ########
#########################
## Import trip matrix ##
w = np.genfromtxt('../../tests/sample.tr',skiprows=1) # trip list (T length: origin dest t_ij)
s_list= np.genfromtxt('../../tests/strengths.str',skiprows=1) # Incoming and outgoing strength
sout = s_list.T[-2]
sin = s_list.T[-1]
N=len(sout)
## Import distance matrix ##
dd = np.genfromtxt('../../tests/cost_matrix.dists')  # cost matrix (NxN)
d = np.zeros((N,N))
for e in dd:
    d[e[0]][e[1]]=e[2]
    
C=np.sum([d[e[0]][e[1]]*e[-1] for e in w]) # total Cost


#######################################
######## Solving saddle points ########
#######################################
for selfs in [True,False]:

## Solve saddle points (check options) ##
    x,y,gamma= mw.fitter_grav.fit_gamma(sin,sout,C,d,tol_s=1e-9,tol_gamma=1e-14,maxreps=1000,
        max_rej=1000,fact=3.,verbose=True, print_tol=True, gamma_ini = sin.max()/C,selfs=selfs)


########################
##### Outputs ##########
########################

    if selfs:
        name = pref+'fixedgrav.xy'
    else:
        name = pref+'fixedgrav_noself.xy'


    lags = np.array([range(len(x)),x,y])
    np.savetxt(name,lags.T,header='#Gamma: %r Node_id x y' % gamma, fmt='%d %r %r')
    
    
    print "Stored result in %s" % name
    print "Accumulated error for sin:%r \t sout:%r" % mw.fitter_grav.dist_check_s(x,y,gamma,sout,sin,d,selfs=selfs)
    print "Accumulated error for cost : %r " % mw.fitter_grav.dist_check_C(x,y,gamma,C,d,selfs=selfs)
    
