"""
Seq gravity model
"""

import numpy as np
from scipy import stats
import os
import sys

try:
    c1 = sys.argv[1] #str list
    c2 = sys.argv[2] #dist list
    gamma = np.float(sys.argv[3])
except:
    pass

###input ##    
s = np.genfromtxt(c1,skiprows=1)
N = len(s)
d = np.genfromtxt(c2,skiprows=1)

dist = np.zeros((N,N))

for l in d:
    dist[int(l[0])][int(l[1])] = l[-1]

### init conds ####
sin = [[e[0],e[-1]] for i,e in enumerate(s)]
sout = [[e[0],e[-2]] for i,e in enumerate(s)]
tij = np.zeros((N,N))
p = np.zeros(N)
## let's go ##
done = 0

import numbapro as nb
@nb.autojit
def compute_p(real_ori,p,sin,gamma,dist):    
    for i in range(len(sin)):
        if sin[i][0]!=real_ori:
            p[i] = sin[i][-1]*np.exp(-gamma*dist[real_ori][sin[i][0]])
        else:
            p[i] = 0
    return p
import time

while len(sout)>0:
    st = time.time()
    ori = np.random.randint(len(sout))
    p = compute_p(sout[ori][0],p,sin,gamma,dist)
    p/=np.sum(p)
    z = np.random.multinomial(1,p)
    dest = np.where(z==1)[0][0]
    ### update ###
    tij[sout[ori][0]][sin[dest][0]]+=1
    sout[ori][1]=sout[ori][1]-1
    sin[dest][1]=sin[dest][1]-1
    if sout[ori][1]<1:
        sout.pop(ori)
    if sin[dest][1]<1:
        p[len(sin)-1]=0
        sin.pop(dest)
    done+=1
    if done % 10000==0:
        print "DOne %d frac %f | took: %f, ETA: %f" % (done,done*1./s.T[-1].sum(),time.time()-st,(time.time()-st)*s.T[-1].sum()/10000)
        st = time.clock()

### output ###
nn = 'sim_ula_python_seq2.tr'
with open(nn,'w') as f:
    print >> f, "#NOde_ori Node_out tij"
    for i in range(N):
        for j in range(N):
            if tij[i][j]>0:
                print >> f, "%d %d %d" % (i,j,tij[i][j])
