"""
    Example for fixed number of binary edges and fixed strength sequence optimization
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
## Import strength and degree sequence ##
s_list= np.genfromtxt('../../tests/strengths.str',skiprows=1) # Incoming and outgoing strength
sout = s_list.T[-2]
sin = s_list.T[-1]


## preconditioning (initial conditions) ##
x=sout/np.sqrt(sout.sum())
y=sin/np.sqrt(sin.sum())


for selfs in [True,False]:
    x,y = mw.fitter_s.balance_xy(sin,sout,tol=1e-9,verbose=True,x_ini=x,y_ini=y,selfs=selfs)
    lags = np.array([range(len(x)),x,y])


#########################
######## Outputs ########
#########################

    if selfs:
        name = pref+'/fixeds.xy'
    else:
        name = pref+'/fixeds_noself.xy'
    
    np.savetxt(name,lags.T,header='#Node_id x y', fmt='%d %r %r')

    print "Stored result in %s" % name
    print "Accumulated error for sin:%r \t sout:%r" % mw.fitter_s.dist_check_s(x,y,sout,sin,selfs)