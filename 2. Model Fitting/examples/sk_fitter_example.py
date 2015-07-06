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
s_list= np.genfromtxt('../../tests/strengths3.str',skiprows=1) # Incoming and outgoing strength
sout = s_list.T[-2]
sin = s_list.T[-1]
k_list = np.genfromtxt('../../tests/degrees3.str',skiprows=1) # Incoming and outgoing degree
kout = k_list.T[-2]
kin = k_list.T[-1]
T = sout.sum()
N = len(kout)

x = np.ones(N)*0.09
y = np.ones(N)*0.09
z = np.ones(N)*0.09
w = np.ones(N)*0.09

M = 37 # 

for agg in [True,False]: 
    for selfs in [True,False]:
        ############################
        ######## Approach 1: Iterative fitting ########
        ###########################
        """
        Preferred approach, convergence not guaranteed but good performance
        """
        print "####################################################"
        print "### Starting for agg:%r and selfs:%r (layers:%r) ###" % (agg,selfs,M)
        print "####################################################"
        tolc=1 # This refers to the accumulated absolute total error. # tolc = \sum_i | sout - <sout>| and \sum_i |kout - <kout>|
        tol_ini = 1e-13 # This refers to the convergence error # Tolerance set very high (s-<s> < 10) due to convergence issues...
        # preconditioning #
        #z,w = mw.fitter_sk.balance_zw(x,y,kin,kout,tol_ini,print_tol=True,selfs=selfs,M=M) # fix degrees (can be done independently)
        #x,y,z,w = mw.fitter_sk.balance_xyzw(sin,sout,kin,kout,tol=tol_ini,verbose=True,print_tol=False, 
        #    print_c=True, maxreps = 1e6, tol_c = tolc, selfs=selfs, agg=agg, M=M, act=100,x_ini=x,y_ini=y,z_ini=z,w_ini=w)
        x,y,z,w = mw.fitter_sk.balance_xyzw(sin,sout,kin,kout,tol=tol_ini,verbose=True,print_tol=False, 
            print_c=True, maxreps = 1e6, tol_c = tolc, selfs=selfs, agg=agg, M=M, act=100,alf=100)
        # fine tunning #
        lags = np.array([range(len(x)),x,y,z,w])
        #########################
        ######## Outputs ########
        #########################

        if selfs:
            name = pref+'/fixedsk'
        else:
            name = pref+'/fixedsk_noself'
        if agg:
            name = name + '_agg_layers%d' % M
        name = name+'.xyzw'

        np.savetxt(name,lags.T,header='#Node_id x y z w', fmt='%d %r %r %r %r')

        print "Stored result in %s" % name
        print "Accumulated error for sin:%r \t sout:%r" % mw.fitter_sk.dist_check_s(x,y,z,w,sout,sin,selfs,agg,M=M)
        print "Accumulated error for kin:%r \t out:%r" % mw.fitter_sk.dist_check_k(x,y,z,w,kout,kin,selfs,agg,M=M)
