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
k_list = np.genfromtxt('../../tests/degrees3.str',skiprows=1) # Incoming and outgoing degree
s_list = np.genfromtxt('../../tests/strengths3.str',skiprows=1) # Incoming and outgoing strneght
kout = k_list.T[-2]
kin = k_list.T[-1]

E=kout.sum()
T=s_list.T[-1].sum()

# average existing weights #
av_w =  T*1./E
M = 37 # number of layers

rhoME = mw.fitter_k.rho_calculator(av_w) # lagrange multiplier for rho ME case
rhoB = mw.fitter_k.rho_calculator(av_w,agg=True,M=1) # lagrange multiplier for rho Binary case
rhoW = mw.fitter_k.rho_calculator(av_w,indist=True,M=1) # lagrange multiplier for rho weighted case
rhoAW = mw.fitter_k.rho_calculator(av_w,indist=True,M=M) # lagrange multiplier for rho aggregated weighted case
rhoAB = mw.fitter_k.rho_calculator(av_w,agg=True,M=M) # lagrange multiplier for rho Aggregated Binary case


rhos = (rhoME,rhoAB,rhoB,rhoW,rhoAW)
## preconditioning (initial conditions) ##



for selfs in [True,False]:

############################
######## Approach: Iterative fitting ########
###########################
    #x=kout/np.sqrt(kout.sum())
    #y=kin/np.sqrt(kin.sum())
    x = np.ones(len(kin))*0.5
    y = np.ones(len(kin))*0.5
    #y=kin/np.sqrt(kin.sum())

    x,y = mw.fitter_k.balance_xy(kin,kout,tol=1e-14,verbose=True,selfs=selfs,x_ini=x,y_ini=y,print_tol=False)
    #x,y = mw.fitter_k.balance_xy(kin,kout,tol=1e-14,verbose=True,selfs=selfs,print_tol=True)

#########################
######## Outputs ########
#########################
    lags = np.array([range(len(x)),x,y])


    if selfs:
        name = pref+'/fixedk.xy'
    else:
        name = pref+'/fixedk_noself.xy'
    
    np.savetxt(name,lags.T,header='#Node_id x y  #rho: ME:%r AB:%r B:%r W:%r AW:%r' % rhos, fmt='%d %r %r')

    print "Stored result in %s" % name
    print "Accumulated error for kin:%r \t kout:%r" % mw.fitter_k.dist_check_k(x,y,kout,kin,selfs)
