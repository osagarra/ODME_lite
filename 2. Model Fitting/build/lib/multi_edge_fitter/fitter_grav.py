""" Fitter for gravity law with given distance matrix """
import numpy as np
from scipy import optimize as opt


def balance_xy(gamma,sin,sout,d_mat,tol=1e-9,maxreps=10000,verbose=False,selfs=True,indist=False,**kwargs):
    """ Balances the flow equations
        d_mat is the cost matrix.
        Input:
            Gamma: Gamma value for distances [float]
            Sin: Incoming strength sequence (N)
            Sout: Outgoing strength sequence (N)
            d_mat: Distance matrix (NxN)
            tol: Maximum tolerance [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            Indist: True for indistinguishable weights
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
    """
    if gamma<0:
        raise ValueError("Gamma must be positive")
    n = len(sin)
    T = sout.sum()
    #sout = sout[sout>0]
    #sin = sin[sin>0]
    a=kwargs.get('x_ini',sout/np.sqrt(T))
    b=kwargs.get('y_ini',sin/np.sqrt(T))
    inds_in = np.where(sin==0)
    inds_out = np.where(sout==0)
    afake = a
    bfake = b
    reps = 0
    d_mat = np.exp(-gamma*d_mat)
    if not selfs:
        while True:
            b = sin/(np.einsum('i,ij',a,d_mat)-np.einsum('i,ii->i',a,d_mat))
            a = sout/(np.einsum('j,ij',b,d_mat)-np.einsum('j,jj->j',b,d_mat))
            a[inds_out]=0
            b[inds_in]=0
            #if reps>10:
            tolb = np.max(np.abs(bfake-b))
            tola = np.max(np.abs(afake-a))
            if  tolb< tol and tola<tol:
                if verbose:
                    print "took %d reps, tolb:%f, tola :%f" % (reps,tolb,tola)
                break
            if reps>maxreps:
                print "Algorithm did not converge after %d reps" % maxreps
                break
            afake = a
            bfake = b
            reps +=1
    else:
        while True:
            b = sin/np.einsum('i,ij',a,d_mat)
            a = sout/np.einsum('j,ij',b,d_mat)
            a[inds_out]=0
            b[inds_in]=0
            #if reps>10:
            tolb = np.max(np.abs(bfake-b))
            tola = np.max(np.abs(afake-a))
            #print "Deltax :%r Deltay: %r" % (tola,tolb)
            if  tolb< tol and tola<tol:
                if verbose:
                    print "took %d reps, tolb:%f, tola :%f" % (reps,tolb,tola)
                break
            if reps>maxreps:
                print "Algorithm did not converge after %d reps" % maxreps
                break
            afake = a
            bfake = b
            reps +=1
    return a,b
     
##############################     

################
### fine grain fitting
################
def fit_gamma(sin,sout,C,d,tol_s=1e-9,tol_gamma=1e-9,maxreps=1000,max_rej=100,fact=3.,selfs=True,verbose=False,**kwargs):
    """
         Minimizes gamma distance between C = t_ij d_ij and <C> = x_i y_j exp(-d_ij gamma) while keeping strength constraints
         input:
             sin,sout (strength sequences) [N arrays]
             d: cost matrix [NxN]
             C: total cost [float]
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
            x,y,gamma
    """
    ### initial opts ###
    gamma_init = kwargs.get('gamma_ini',sin.sum()*1./C)
    delta0 = kwargs.get('delta0',gamma_init/2.)
    x_ini = kwargs.get('x_ini',sout/np.sqrt(sout.sum()))
    y_ini = kwargs.get('y_ini',sin/np.sqrt(sin.sum()))
    T=sin.sum()
    ####################
    print ".... Gamma fitter ....."
    ####################
    gamma = gamma_init
    x,y = balance_xy(gamma,sin,sout,d,tol_s,maxreps,selfs=selfs,x_ini=x_ini,y_ini=y_ini)
    delta_C=dist_check_C(x,y,gamma,C,d,selfs=selfs)
    delta_s=dist_check_s(x,y,gamma,sin,sout,d,selfs=selfs)
    print "S_in grad:%f S_out grad:%f C_grad/T:%f C" % (delta_s[0],delta_s[1],delta_C)
    delta = delta0
    reps=0
    s=1
    while reps<maxreps and np.abs(delta_C)>tol_gamma:
        new_gamma=gamma+delta*s
        if new_gamma<0:
            i=1
            while new_gamma<0 and delta/i>tol_gamma:
                new_gamma=gamma+delta*s/i # lambda must be positive
                i+=1
            if new_gamma<0:
                print "Gamma <0,i:%d, new_gamma:%f. Aborting"% (i,new_gamma)
                break
        x,y = balance_xy(new_gamma,sin,sout,d,tol_s,maxreps,selfs=selfs,y_ini=y,x_ini=x)
        new_delta_C=dist_check_C(x,y,new_gamma,C,d,selfs=selfs)
        if np.abs(new_delta_C)<np.abs(delta_C): #accepted
            if np.sign(new_delta_C) != np.sign(delta_C): # if change of sign, change direction!
                s=-s
            gamma = new_gamma
            delta_C = new_delta_C
            #delta=np.abs(gamma-new_gamma)
            delta = delta*fact
            if verbose:
                print "Accepted! gamma:%r"  % gamma
        else: #rejected, reduce step and change direction
            delta = delta/fact
            s=-s
            delta_gamma=new_gamma-gamma
            if verbose:
                print "Rejected! movement: %r delta_gamma:%r" % (delta_C-new_delta_C,new_gamma-gamma)
            if np.abs(delta_gamma)<1e-14:
                break
        reps+=1
    x,y = balance_xy(gamma,sin,sout,d,tol_s,maxreps,selfs=selfs,y_ini=y,x_ini=x)
    delta_s = dist_check_s(x,y,gamma,sin,sout,d,selfs=selfs)
    print "All deltas: sin:%r sout:%r gamma:%r for gamma:%r" % (delta_s[0],delta_s[1],
        dist_check_C(x,y,gamma,C,d,selfs=selfs),gamma)
    return x,y,gamma
    
    
######### Loglikelyhood functions ######
# def loglikelyhood_grav(gamma,x,y,d,C,sin,sout):
#     """
#         Computes loglikelyhood function 
#         Inputs:
#             x,y [arrays] -> sequence of lagrange multipliers
#             d -> cost matrix
#             C -> total network cost
#     """
#     Lg = loglikelyhood_gamma(gamma,x,y,d,C)
#     indsx = np.where(x>0)
#     indsy = np.where(y>0)
#     #print "x len sout len: %d %d, ylen sin len: %d %d" % (len(indsx[0]),len(indsy[0]),len(np.where(sin>0)[0]),len(np.where(sout>0)[0]))
#     Lx = np.einsum('i,i',sout[indsx],np.log(x[indsx]))
#     Ly = np.einsum('i,i',sin[indsy],np.log(y[indsy]))
#     return -(Lg-Lx-Ly)
# 
# def loglikelyhood_gamma(gamma,x,y,d,C):
#     """
#         Computes loglikelyhood function for gamma only
#         (scalar function)
#         
#     """
#     Lc= C*gamma
#     #L = loglikelyhood_gamma_term(gamma,x,y,d)
#     L=np.einsum('i,j,ij',x,y,np.exp(-gamma*d))
#     return L+Lc
#     
    
#### distance checkers ####
def dist_check_s(x,y,gamma,sout,sin,d,selfs=True):
    """
        Computes distance between prediction and reality for distances
        x,y --> arrays of lagrange multpliers
        gamma: distance multiplier
        sout,sin -- > real strength sequences
        d: distance matrix
        self: True for self loops
    """
    if not selfs:
        extra1 = np.einsum('j,jj->j',x,np.exp(-d*gamma)) # diagonal
        extra2 = np.einsum('i,ii->i',y,np.exp(-d*gamma)) # diagonal
    else:
        extra1=0
        extra2=0
    gradin = np.sum(sin-y*(np.einsum('i,ij',x,np.exp(-d*gamma))-extra1))
    gradout = np.sum(sout-x*(np.einsum('j,ij',y,np.exp(-d*gamma))-extra2))
    return gradin,gradout
def dist_check_C(x,y,gamma,C,d,selfs=True):
    """
        Computes distance between prediction and reality for strengths
        x,y --> arrays of lagrange multpliers
        gamma: distance multiplier
        c: real cost [sum_{ij} t_ij * d_ij]
        d: distance matrix
        self: True for self loops
    """
    if not selfs:
        extra = np.einsum('j,j,jj',x,y,d*np.exp(-d*gamma))
    else:
        extra = 0
    return C-np.einsum('i,j,ij',x,y,d*np.exp(-d*gamma))+extra
