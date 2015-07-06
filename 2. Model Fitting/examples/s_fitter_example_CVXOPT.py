"""
    Example for fixed number of binary edges and fixed strength sequence optimization using CVXOPT
"""
# imports #

import multi_edge_fitter as mw
import numpy as np

import sys
import os

# options prefix #
pref = 'lagrange_mult'

try:
    os.mkdir(pref)
except OSError:
    pass
    

#########################
######## Inputs ########
#########################
## Import strength and degree sequence ##
s_list= np.genfromtxt('../../tests/strengths_hard.str',skiprows=1) # Incoming and outgoing strength
#s_list= np.genfromtxt('../../tests/strengths2.str',skiprows=1) # Incoming and outgoing strength
sout = s_list.T[-2]
sin = s_list.T[-1]
#layers = 37
N=len(sout)
inds_out = np.where(sout==0)
inds_in = np.where(sin==0)
layers = [1,10,37] 


abstol = 1e-10
reltol = 1e-10
feastol = 1e-10
case = 'W'

for l in layers:
    for selfs in [True,False]:
        print "######################################################"
        print "##### Indist case: Will do case:%r selfs:%r with layers:%d case ######" % (case,selfs,l)       
        ## generate object ##
        ME = mw.fitter_s_CVXOPT.Edge_Maximizer(sin,sout,selfs=selfs,M=l,case='W')
        ## precondition ##
        ME.precondition(tol=1e-5)
        ## Solve ##
        # Case 1: Unbounded, unconstrained
        print "=== Case 1: Unbounded, unconstrained ==="
        try:
            q = ME.solve_using_cvxopt(abstol = abstol,reltol = reltol,show_progress=True,feastol = feastol,refinement=2, maxiters=30, bound=True)
        except ValueError:
            print "!!!!! SOmething went south..."
        print "=== Case 2: Bounded, unconstrained ==="
        # Case 2: Unconstrained, bounded
        try:
            q2 = ME.solve_using_cvxopt(abstol = abstol,reltol = reltol,show_progress=True,feastol = feastol,refinement=2, bound=True, maxiters=30)
        except ValueError:
            print "!!!!! SOmething went south..."
        print "=== Case 3: UnBounded, constrained ==="
        # Case 3: Constrained, unbounded
        try:
            q3 = ME.solve_using_cvxopt(abstol = abstol,reltol = reltol,show_progress=True,feastol = feastol,refinement=2, cons=True, maxiters=100)        
        except ValueError:
            print "!!!!! SOmething went south..."

        #########
        # Note: 
        # If you get the error singular KKT matrix, the result still is updated if it is better than what was there before.
        # It is recommended to try the option ME.solve_as_system(tol=tolx) with tolx as desired tolerance, to tune the final result in cases with low accuracy 
        #########
        xy = ME.var_result()
        print "Xmax: %r Ymax: %r" % (xy.T[0].max(),xy.T[1].max())
        lags = np.array([range(ME.N),xy.T[0], xy.T[1]]).T
        if q['success'] or q2['success'] or q3['success'] or ME.check_s()[0]<1e-5:
            #########################
            ######## Outputs ########
            #########################

            if selfs:
                name = pref+'/fixeds_indist_CVXOPT_layers%d' % l
            else:
                name = pref+'/fixeds_noself_indist_CVXOPT_layers%d' % l
            
            name = name + '.xy'
            np.savetxt(name,lags,header='#Node_id x y', fmt='%d %r %r')

            print "Stored result in %s" % name
            print "===== Results ======"
            print "\t\t Accumulated error for sin:%r \t sout:%r" % ME.check_s()
            print "===== ======= ======"


