"""
    Optimize pij with strength sequence fixed and some trips also fixed
"""
# imports #
import multi_edge_fitter as mw
import numpy as np
import sys
import os
try:
    from numbapro import autojit
except ImportError:
    from numba import autojit


pref = 'lagrange_mult'

try:
    os.mkdir(pref)
except OSError:
    pass



@autojit
def fill_vect(pij,N):
    v = np.zeros((np.prod(pij.shape),3))
    q=0
    for i in range(N):
        for j in range(N):
            v[q][0]=i
            v[q][1]=j
            v[q][2]=pij[i][j]
            q+=1
    return v




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
## Import distance or cost matrix ##
dd = np.genfromtxt('../../tests/cost_matrix.dists',skiprows=1)
d = np.zeros((N,N))
for e in dd:
    d[e[0]][e[1]] = float(e[-1])
    d[e[1]][e[0]] = float(e[-1])

### Solve pij ### w_max is the minimal weight for trip to be fixed
w_max = 2
#pij,gamma = mw.fitter_pij.balance_qij_x_y(sin,sout,w,w_max,tol_s=1e-14,maxreps=1000,selfs=True) # if no distance is given, it just ignores the cost constraint

for selfs in [True,False]:
    pij,gamma = mw.fitter_pij.balance_qij_x_y(sin,sout,w,w_max,tol_s=1e-14,maxreps=1000,selfs=True, d=d)

#########################
######## Outputs ########
#########################
    vv = fill_vect(pij,N)
    
    if selfs:
        name = pref+'/pij_fixed_wmax%d.pij' % w_max
    else:
        name = pref+'/pij_fixed_wmax%d_noself.pij' % w_max

    
    np.savetxt(name,vv,header='#node_ori node_dest pij (gamma=%f)' % gamma, fmt='%d %d %r')
    print "Stored result in %s" % name
    