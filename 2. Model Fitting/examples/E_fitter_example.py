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
k_list = np.genfromtxt('../../tests/degrees.str',skiprows=1) # Incoming and outgoing degree
kout = k_list.T[-2]
kin = k_list.T[-1]

E=kout.sum() # Number of edges
N=len(sout) # Number of nodes
T=sin.sum() # Number of events 
alf = np.sqrt(sout.max()*sin.max()/T)

M = 37 # 37 layers 

for agg in [True,False]: 
    for selfs in [True,False]:


        ## preconditioning (initial conditions) ##
        #lam = np.log(E)
        x=sout/np.sqrt(sout.sum())
        y=sin/np.sqrt(sin.sum())
        #x=np.ones(N)
        #y=np.ones(N)


        ############################
        ######## Approach 1: Iterative fitting ########
        ###########################
        """
        Preferred approach, convergence not guaranteed but good performance
        """
        # ideally do a first pass with low tolerance for s and lambda, and then do fine tunning #
        lam = E/T
        tol_ini = 1e-2
        x,y,lam = mw.fitter_E.fit_lambda(sin,sout,E,tol_gamma=tol_ini,tol_s=1e-8,verbose=False,
            lambda_ini=lam,x_ini=x,y_ini=y,print_tol=False,lambda_min=0.1,selfs=selfs,print_c=False,act=100,tol_c=1,agg=agg,M=M)
        delta0 = 10*tol_ini
        x,y,lam = mw.fitter_E.fit_lambda(sin,sout,E,tol_gamma=1e-9,tol_s=1e-13,verbose=False,
            lambda_ini=lam,x_ini=x,y_ini=y,print_tol=False,print_c=True,lambda_min=0.001,delta0=delta0,selfs=selfs,agg=agg,M=M) 
        # lambda_min is set for convergence, delta0 set to previous 10xtol_Gamma #

        lags = np.array([range(len(x)),x,y])

        if selfs and not agg:
            ############################
            ######## Approach2: Likelyhood (only implemented for the case with self-loops)  ########
            ###########################
            """
            This approach is much slower, use only if approach 1 does not converge
            """

            ## Compute log likelyhood constant term ## 
            #smalls = np.where(w.T[-1]<20)
            #bigs = np.where(w.T[-1]>=20)
            #Lcnt = np.einsum('i->',np.log(misc.factorial(w.T[-1][smalls]))) + np.einsum('i->',w.T[-1][bigs]*(np.log(w.T[-1][bigs])-1)) 
            Lcnt = 0;
            #constant term on loglikelyhood calculated using stirling#


            ## preconditioning with fixed lambda ###
            x2,y2 = mw.fitter_E.balance_xy(lam,sin,sout,verbose=True,x_ini=x,y_ini=y,tol=1e-9)
            #### minimization (using class) ####
            ## initialize class ##
            ME = mw.fitter_E.Multi_Edge_Maximizer_E(sin,sout,E,cnt=Lcnt,x=x2,y=y2,lam=lam)
            ## minimize ## (see options)
            #q=ME.minimize('L-BFGS-B',ftol=1e-14,gtol=1e-9) # use preferred method here...
            q=ME.minimize('TNC',ftol=1e-9,gtol=1e-9) # use preferred method here...
            lam2 =ME.lam
            x2= ME.x
            y2= ME.y
            lags = np.array([range(len(x2)),x2,y2])

        #########################
        ######## Outputs ########
        #########################

        if selfs:
            name = pref+'/fixedsE'
        else:
            name = pref+'/fixedsE_noself'
        if agg:
            name = name + '_agg_layers%d' % M
        name = name+'.xy'
        np.savetxt(name,lags.T,header='#Lamda: %r Node_id x y' % lam, fmt='%d %r %r')


        print "Stored result in %s" % name
        print "Accumulated error for sin:%r \t sout:%r" % mw.fitter_E.dist_check_s(x,y,lam,sout,sin,selfs=selfs,agg=agg,M=M)
        print "Accumulated error for E : %r " % mw.fitter_E.dist_check_E(x,y,lam,E,selfs=selfs,agg=agg,M=M)
        
