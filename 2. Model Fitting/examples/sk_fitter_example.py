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
s_list= np.genfromtxt('../../tests/strengths2.str',skiprows=1) # Incoming and outgoing strength
sout = s_list.T[-2]
sin = s_list.T[-1]
k_list = np.genfromtxt('../../tests/degrees2.str',skiprows=1) # Incoming and outgoing degree
kout = k_list.T[-2]
kin = k_list.T[-1]
T = sout.sum()
N = len(kout)


x = sout/np.sqrt(T)
y = sin/np.sqrt(T)
#x = np.ones(N)
#y = np.ones(N)
z = np.ones(N)
w = np.ones(N)

for selfs in [False,True]:


############################
######## Approach 1: Iterative fitting ########
###########################
    """
    Preferred approach, convergence not guaranteed but good performance
    """
    tol_ini = 1e-12
    # Tolerance set very high (s-<s> < 10) due to convergence issues...
    x,y,z,w = mw.fitter_sk.balance_xyzw(sin,sout,kin,kout,tol=tol_ini,verbose=True,print_tol=True, 
        print_c=True, maxreps = 1e6, tol_c = 10, x_ini = x, y_ini = y, z_ini = z, w_ini = w)
    # fine tunning #
    z,w = mw.fitter_sk.balance_zw(x,y,kin,kout,tol_ini,print_tol=True,z_ini=z,w_ini=w) # fix degrees (can be done independently)
    lags = np.array([range(len(x)),x,y,z,w])


#########################
######## Outputs ########
#########################

    if selfs:
        name = pref+'/fixedsk.xyzw'
    else:
        name = pref+'/fixedsk_noself.xyzw'
    

    np.savetxt(name,lags.T,header='#Node_id x y z w', fmt='%d %r %r %r %r')

    print "Stored result in %s" % name
    print "Accumulated error for sin:%r \t sout:%r" % mw.fitter_sk.dist_check_s(x,y,z,w,sout,sin,selfs)
    print "Accumulated error for kin:%r \t out:%r" % mw.fitter_sk.dist_check_k(x,y,z,w,kout,kin,selfs)