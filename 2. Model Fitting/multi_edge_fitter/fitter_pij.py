""" Fitter for gravity law with given distance matrix and some trips fixed """
import numpy as np
import copy
import sys
try:
    from numbapro import autojit
except ImportError:
    from numba import autojit

### aux funcs ##
@autojit
def update_q(w,x,y):
    z = np.einsum('i,j',x,y)
    for ww in w:
        z[ww[0]][ww[1]] = ww[2]
    return z

@autojit
def update_q2(w,x,y,gamma,d):
    z = np.einsum('i,j,ij->ij',x,y,np.exp(-gamma*d))
    for ww in w:
        z[ww[0]][ww[1]] = ww[2]
    return z



@autojit
def fill_matrix(inds,n,selfs=True):
    qij =np.ones((n,n))
    for w in inds.T:
        qij[w[0]][w[1]] = 0
    if not selfs:
        for i in range(n):
            qij[i][i]=0
    return qij

def compute_excess_s(n,w,indd): # very slow!
    """ COmputes excess s
        input:
            n: number of nodes
            w: fixed trip seq (3 col i j t_ij)
            ind: indd to comptue (for out:0, for in:1)
        output:
            ss: strength correction sequence for each node
    """
    #delta = np.zeros(len(ss))
    #for i in xrange(n):
    #    indss=np.where(w.T[indd]==i)
    #    delta[i] = np.sum(w.T[-1][indd])
    #    #ww=np.delete(ww,ss)
    #return ss
    return np.array([w.T[-1][np.where(w.T[indd]==i)].sum() for i in range(n)])

#####################
### balance funcs ###
#####################
def balance_xy(ssout,ssin,inds,tol=1e-9,maxreps=10000,verbose=False,selfs=True,**kwargs):
    """ Balances the flow equations with given strength sequence in/out and given
        w values for some trips t_ij.
        
        Input:
            Sin_c: Corrected Incoming strength sequence (N)
                sin = sin - K_j^in
            Sout_c: Corrected Outgoing strength sequence (N)
                sout = sout - K_i^out
            inds: i,j pairs in set q
                Where qset is defined as the set of ij pairs that are fixed
            tol: Maximum tolerance [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
                q_ini: Initial guess for q_ij [matrix NxN]
                gamma: distance lagrange multiplier
                d: If include distance function
    """
    d_flag=False
    n = len(ssin)
    T = ssout.sum()
    a = kwargs.get('x_ini',ssout/np.sqrt(T))
    b = kwargs.get('y_ini',ssin/np.sqrt(T))
    gamma = kwargs.get('gamma',0)
    if gamma!=0:
        try:
            d = kwargs['dist']
            dd = np.exp(-gamma*d)
            d_flag=True
        except KeyError:
            print "Distance not defined, ignoring"
            d_flag=False
    ## inds ##
    #inds_in = np.where(ssin==0)
    #inds_out = np.where(ssout==0)
    afake = a
    bfake = b
    reps = 0
    qij = fill_matrix(inds,n,selfs)
    while True:
        if d_flag:
            tot=np.einsum('i,ij,ij->j',a,qij,dd)
        else:
            tot=np.einsum('i,ij',a,qij)
        b = ssin/tot
        #b[np.where(tot==0)]=0
        if d_flag:
            tot=np.einsum('j,ij,ij->i',b,qij,dd)
        else:
            tot=np.einsum('j,ij',b,qij)
        a = ssout/tot
        #a[np.where(tot==0)]=0
        tolb = np.max(np.abs(bfake-b))
        tola = np.max(np.abs(afake-a))
        #print "Tola %f Tolb %f" % (tola,tolb)
        #print "Dist s %f" % dist_check_s(a,b,ssout,ssin,inds_x,inds_y,selfs=selfs)[0].max()
        if  tolb< tol and tola<tol:
            break
            if verbose:
                print "took %d reps, tolb:%f, tola :%f tolq : %f" % (reps,tolb,tola,tolq)
        if reps>maxreps:
            print "Algorithm did not converge after %d reps" % maxreps
            break
        afake = a
        bfake = b
        #qij_fake = copy.deepcopy(qij)
        reps +=1
    if verbose:
        print "Algorithm converged after %d iters. Tola %r Tolb %r" % (reps,tola,tolb)
    return a,b



def balance_qij_x_y(sout,sin,w,w_max=10,tol_s=1e-9,maxreps=1000,selfs=True,verbose=True,**kwargs):
    """
         Balances x,y and qij with given strenght and some trips constraints
         input:
             sin,sout (strength sequences) [N arrays]
             w: Fixed trip sequence [3 col array, i j t_ij]
             w_max: W_max value to fix trips to (all trips below this value will not be fixed). [float]
             tol_s: Tolerance for strength balancing [float]
             tol_qij: Tolerance for qij [float]
             mareps: Maximum iterations [int]
             max_rej: maximum rejections
             selfs: Accept self loops? [Boolean]
             kwargs:
                 x_ini: initial guess for x
                 y_ini: initial guess for y
                 d :  if also consider distances and costs
        output:
            qij [NxN matrix with normalized pij probabilties]
            gamma value (lagrange multiplier, float)
    """
    if verbose:
        print "########################"
        print "##### Pij fitter #######"
        print "########################"
    ### initial opts ###
    n = len(sin)
    ## detect distances ##
    d = kwargs.get('d',int(0))
    if type(d)==int:
        d_flag = False
        d = None
    else:
        d_flag = True
        if verbose:
            print "... Distance detected, fixing also costs! ..."
    len_or = len(w)
    T = sin.sum()
    if d_flag:
        C = np.sum([e[-1]*d[e[0]][e[1]] for e in w]) # total cost
    w = w[np.where(w.T[-1]>w_max)]
    print "Fixing %d trips of out %d (%f frac) with w_max=%f. Fixing %d edges out of %d (%f frac)" % (w.T[-1].sum()
        ,T,w.T[-1].sum()*1./T,w_max,len(w),len_or,len(w)*1./len_or)
    x_ini = kwargs.get('x_ini',sout/np.sqrt(sin.sum()-w.T[-1].sum()))
    y_ini = kwargs.get('y_ini',sin/np.sqrt(sout.sum()-w.T[-1].sum()))
    inds =np.array(w.T[:-1],dtype=np.int)
    ####################
    if verbose:
        print ".... qij,x,y balancer ....."
        print ".. Computing excess strengths .."
    delta_sout = compute_excess_s(n,w,indd=0)
    delta_sin = compute_excess_s(n,w,indd=1)
    ssout = sout - delta_sout 
    ssin = sin- delta_sin
    if verbose:
        print ".. Balancing .."
    if not d_flag:
        x,y = balance_xy(ssout,ssin,inds,tol_s,maxreps,selfs=selfs,x_ini=x_ini,y_ini=y_ini,verbose=verbose)
        qij = update_q(w,x,y)
        gamma=0
    else:
        x,y,gamma = fit_gamma(w,ssout,ssin,C,d,inds,tol_s=tol_s,maxreps=maxreps,x_ini=x_ini,y_ini=y_ini,verbose=verbose)
        qij = update_q2(w,x,y,gamma,d)
    delta_w=np.sum(sout-np.einsum('ij->i',qij))
    if verbose:
        print "Total strength deviation: %lf" % delta_w
    return qij*1./qij.sum(),gamma


##############################     
# Fitter functions #
##############################     


def fit_gamma(w,ssout,ssin,C,d,inds,tol_s=1e-9,tol_gamma=1e-9,maxreps=1000,delta0 =-1,max_rej=100,fact=3.,selfs=True,verbose=False,**kwargs):
    """
         Minimizes gamma distance between C = t_ij d_ij and <C> = x_i y_j exp(-d_ij gamma) while keeping strength constraints
         and some edges fixed too.
         input:
             reduced sin,sout (strength sequences) [N arrays] --> minus fixed values
             d: cost matrix [NxN]
             Total C: total cost [float] --> not reduced!
             tol_s: Tolerance for strength balancing [float]
             tol_gamma: Tolerance for gamma [float]
             mareps: Maximum iterations [int]
             max_rej: maximum rejections
             fact: Factor to approach delta (the bigger, the more likely not to converge)
             selfs: Accept slef loops? [Boolean]
             kwargs:
                 gamma_ini: initial guess for distance [float]
                 x_ini: initial guess for x
                 y_ini: initial guess for y
                 delta0 :initial jump for gamma for the algorithm [float]
        output:
            qij (normalized) [NxN] matrixs
    """
    ### initial opts ###
    gamma_init = kwargs.get('gamma_ini',ssin.sum()/C)
    delta0 = kwargs.get('delta0',gamma_init/2.)
    x_ini = kwargs.get('x_ini',ssout/np.sqrt(ssout.sum()))
    y_ini = kwargs.get('y_ini',ssin/np.sqrt(ssin.sum()))
    T=ssin.sum()
    ####################
    print ".... Gamma fitter with fixed some trips ....."
    ####################
    gamma = gamma_init
    x,y = balance_xy(ssout, ssin, inds,gamma=gamma,dist=d,tol=tol_s,maxreps=maxreps,selfs=selfs,x_ini=x_ini,
        y_ini=y_ini,verbose=verbose)
    qij = update_q2(w,x,y,gamma,d)
    delta_C=dist_check_C(C,qij,d)
    print "C_grad/T:%f " % (delta_C/T)
    delta = delta0
    reps=0
    s=1
    while reps<maxreps and np.abs(delta_C)>tol_gamma:
        new_gamma=gamma+delta*s
        x,y = balance_xy(ssout, ssin, inds,gamma=new_gamma,dist=d,tol=tol_s,maxreps=maxreps,selfs=selfs,x_ini=x,y_ini=y,
            verbose=False)        
        qij  =update_q2(w,x,y,new_gamma,d)
        new_delta_C=dist_check_C(C,qij,d)/T        
        if np.abs(new_delta_C)<np.abs(delta_C): #accepted
            if np.sign(new_delta_C) != np.sign(delta_C): # if change of sign, change direction!
                s=-s
            gamma = new_gamma
            delta_C = new_delta_C
            #delta=np.abs(gamma-new_gamma)
            delta = delta*fact
            if verbose:
                print "Accepted! gamma:%r | Total strength deviation: %lf"  % (gamma,np.sum(ssout-np.einsum('ij->i',qij)) + w.T[-1].sum())
        else: #rejected, reduce step and change direction
            delta = delta/fact
            s=-s
            delta_gamma=new_gamma-gamma
            if verbose:
                print "Rejected! movement: %r delta_gamma:%r delta_C:%r" % (delta_C-new_delta_C,new_gamma-gamma,delta_C)
            if np.abs(delta_gamma)<1e-14:
                break
        reps+=1
    qij = update_q2(w,x,y,gamma,d)
    delta_C=dist_check_C(C,qij,d)/T        
    print "Delta Cost :%r for gamma:%r" % (delta_C,gamma)
    return x,y,gamma
    

#### distance checkers ####
def dist_check_C(C,qij,d):
    """
        Computes distance between prediction and reality for strengths
        qij matrix
        C: total cost
    """
    return C-np.einsum('ij,ij',qij,d)
    
