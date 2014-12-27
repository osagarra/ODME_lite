"""
    Example for fixed binary degree sequence
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

#########################a
######## Inputs ########
#########################
## Import degree sequence ##
k_list = np.genfromtxt('../../tests/degrees.str',skiprows=1) # Incoming and outgoing degree
kout = k_list.T[-2]
kin = k_list.T[-1]

## preconditioning (initial conditions) ##



for selfs in [True,False]:

############################
######## Approach: Iterative fitting ########
###########################
    x=kout/np.sqrt(kout.sum())
    y=kin/np.sqrt(kin.sum())

    x,y = mw.fitter_k.balance_xy(kin,kout,tol=1e-9,verbose=True,selfs=selfs,x_ini=x,y_ini=y,print_tol=True)

#########################
######## Outputs ########
#########################
    lags = np.array([range(len(x)),x,y])


    if selfs:
        name = pref+'/fixedk.xy'
    else:
        name = pref+'/fixedk_noself.xy'
    
    np.savetxt(name,lags.T,header='#Node_id x y', fmt='%d %r %r')

    print "Stored result in %s" % name
    print "Accumulated error for kin:%r \t kout:%r" % mw.fitter_k.dist_check_k(x,y,kout,kin,selfs)