"""
Dummy test for CVXOPT and any maximizer
"""
## imports ##
import numpy as np
import multi_edge_fitter as mw

# opts #
selfs=True
M = 10 # layers

# DImension of network
N = 5
# Number of non-zero elements
nvarx = N-1
nvary = N-1
# Strength: Sequence of ones except the zero elements
sout = np.ones(N)
sin  = np.ones(N)
for i in range(N-nvarx):
    sout[i] = 0 # one empty
for i in range(N-1,N-nvary,-1):
    sin[i] = 0 # one empty
i = 0
if sin.sum()>sout.sum():
    while sin.sum()>sout.sum():
        sout[N-i-1] += 1
        i+=1
elif sout.sum()>sin.sum():
    while sin.sum()<sout.sum():
        sin[::-1][N-i-1] += 1
        i+=1
else: # all good
    pass


### Solve using one that works ###
W = mw.fitter_s.Weighted_Edge_Maximizer(sin,sout,selfs=selfs,M=M)
q=W.minimize(tol=1e-12)
sol_empty = W.var_adapt(False)
sol_full = W.var_adapt(True)

### Now load CVXOPT ###
W2 = mw.fitter_s_CVXOPT.Weighted_Edge_Maximizer_CVXOPT(sin,sout,selfs=selfs,M=M)
q=W2.solve_as_system()
xx_no = W2.var_adapt()
xx_cx = W2.to_CVXOPT_dense(W2.var_adapt())
z = [np.random.random()]
# check hessian #
extra_args_cx = W2.prepare_extra_args(True)
extra_args_no = W2.prepare_extra_args(False)
h_cx = mw.fitter_s_CVXOPT.loglikelyhood_hess_cvxopt(xx_cx,z,extra_args_cx,True)
h_no = mw.fitter_s_CVXOPT.loglikelyhood_hess_cvxopt(xx_no,z,extra_args_no)

h_cx_array = np.zeros(h_cx.size)
h_cx_array[(np.array(h_cx.I).flatten(),np.array(h_cx.J).flatten())] = np.array(h_cx.V).flatten()
h_no_array = np.zeros(h_no.size)
h_no_array[(np.array(h_no.I).flatten(),np.array(h_no.J).flatten())] = np.array(h_no.V).flatten()

from numpy import linalg
rh_h_cx = linalg.matrix_rank(h_cx_array)
rh_h_no = linalg.matrix_rank(h_no_array)
print "Rank hessian: %r" % rh_h_cx

###  ###