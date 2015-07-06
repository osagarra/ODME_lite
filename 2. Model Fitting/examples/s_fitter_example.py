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
#s_list= np.genfromtxt('../../tests/strengths_hard_ns.str',skiprows=1) # Incoming and outgoing strength
sout = s_list.T[-2]
sin = s_list.T[-1]
layers = [1,37]
N= len(sout)


"""
This problem can be done in several ways, one is by maximizing Loglikelyhood using any method, or by solving the system of equations 
defined by the gradient of the loglikeylhood. 
## case W ##
    SInce the loglikelyhood has a term  
    \sum_{ij} \ln (1-x_i y_j)
    Then a constraint of the type max(x_i y_j)<=1 must be imposed.
    When using the linear solver, these constraints are not enforced, but it provides a very good seed for the maximizer to start, in general.
    The maximizer uses an unconstrained, bounded method, so the only way to enforce the constraints is by setting,
    \max(x_i y_y) <= \max (x_i) \max(y_j) <= 1. This can be tunned using the parameter bmax on the call to the minimizer.
    Note that the above expression does not cover all the feasible phase space, hence the solutions found might be sub-optimal.

    Use each case according to your needs, the best solution is to apply both combined (as done in this example).

    Note as I said, that this method can be used to precondition the problem, in the event a suboptimal solution is found.

    If you want more precision, install CVXOPT (http://www.cvxopt.org) and use the alternative package provided.
## cases ME and B ##
    In these cases, the solution is not problematic and is easily achieved.
"""

print "\n\n"
print "###### Using Loglikelyhood & linear solver method #####"


cases = ['ME','W','B']
for case in cases:
    print "============================"
    print "#### Case: %r" % case
    print "============================"
    if case == 'B':
        layerss =  [37]
    else:
        layerss = layers
    for l in layerss:
        for selfs in [True,False]:
            print "\t##### Will do selfs:%r with layers:%d ######" % (selfs,l)       
            ME = mw.fitter_s.Edge_Maximizer(sin,sout,selfs=selfs,M=l,case=case)
            ### pre conditioning ###
            ME.precondition(verbose=False)
            tol_ini = 1e-12
            q1 = ME.solve_as_system(tol_c=tol_ini,verbose=True)
            q2 = {}
            if case != 'ME': # ME case is fully convex and balancing approach is enough
                ### One can do better playing around a bit with the parametters, depending on the case... ###
                # See: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.fmin_tnc.html
                eps0 = 1e-10
                epsi = eps0
                q1=ME.minimize('TNC',verbose=False,eta=0.5,xtol=tol_ini,tol=tol_ini,epsi=epsi,bmax=1.) # you may play with bmax
                #q = ME.solve_as_system('hybr') # apply linear solving, DOES NOT ENFORCE THE CONSTRAINTS
                q2=ME.minimize('TNC',verbose=False,eta=0.5,xtol=tol_ini,tol=tol_ini,epsi=epsi,bmax=1.5) # you may play with bmax
            if (q1['success'] or q2['success']) and not q1['violated']:
                xy = ME.var_result()
                lags = np.array([range(ME.N),xy.T[0],xy.T[1]])
                #########################
                ######## Outputs ########
                #########################

                if selfs:
                    name = pref+'/fixeds_indist_layers%d' % l
                else:
                    name = pref+'/fixeds_noself_indist_layers%d' % l
        
                name = name + '.xy'
                np.savetxt(name,lags.T,header='#Node_id x y', fmt='%d %r %r')

                print "\t=== Stored result in %s ===" % name

                print "\t===== Results ======"
                err = ME.check_s()
                print "\t\t Accumulated error for sin:%r \t sout:%r | Cons ok:%r | x_max:%.4f y_max:%.4f" % (err[0],err[1],
                    ME.check_domain(),ME.x.max(),ME.y.max())
                print "\t===== ======= ======"
                x = ME.x # for the next iteration, use as x0
                y = ME.y
            else:
                print "\t\t ==== Minimizaiton failed... ===="
        

