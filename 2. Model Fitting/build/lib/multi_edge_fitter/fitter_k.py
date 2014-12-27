""" Fitter for multi-edge network with given degree sequence """
#from numbapro import guvectorize,autojit
#from numba import autojit, double, int32, jit
import numpy as np
#import distances as dist
from scipy import optimize as opt


def balance_xy(kin,kout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,**kwargs):
    """ Balances the degree equations
        Input:
            kin: Incoming degree sequence (N)
            kout: Outgoing degree sequence (N)
            tol: Maximum tolerance [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
                print_tol: Print tolerance update
    """
    n = len(kin)
    E = kout.sum()
    #sout = sout[sout>0]
    #sin = sin[sin>0]
    a=kwargs.get('x_ini',kout/np.sqrt(E))
    b=kwargs.get('y_ini',kin/np.sqrt(E))
    print_tol = kwargs.get('print_tol',False)
    inds_in = np.where(kin==0)
    inds_out = np.where(kout==0)
    afake = a
    bfake = b
    reps = 0
    if not selfs:
        while True:
            aux = 1.+np.einsum('i,j',a,b)
            b = kin/(np.einsum('i,ij',a,1./aux)-a/(1.+a*b))
            aux = 1.+np.einsum('i,j',a,b)
            a = kout/(np.einsum('j,ij',b,1./aux)-b/(1.+a*b))
            a[inds_out]=0
            b[inds_in]=0
            #if reps>10:
            tolb = np.max(np.abs(bfake-b))
            tola = np.max(np.abs(afake-a))
            if print_tol:
				print "Delta_x:%r Delta_y:%r" % (tola,tolb)
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
            aux = 1.+np.einsum('i,j',a,b)
            b = kin/np.einsum('i,ij',a,1./aux)
            aux = 1.+np.einsum('i,j',a,b)
            a = kout/np.einsum('j,ij',b,1./aux)
            a[inds_out]=0
            b[inds_in]=0
            #if reps>10:
            tolb = np.max(np.abs(bfake-b))
            tola = np.max(np.abs(afake-a))
            if print_tol:
				print "Delta_x:%r Delta_y:%r" % (tola,tolb)
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



#### distance checkers ####
def dist_check_k(x,y,kout,kin,selfs=True):
    """
        Computes distance between prediction and reality for distances
        x,y --> arrays of lagrange multpliers
        kout,kin -- > real degree sequences
        self: True for self loops
    """
    if not selfs:
        extra1 = np.sum(np.einsum('j,j->j',x,y)/(1+np.einsum('j,j->j',x,y))) # diagonal
    else:
        extra1=0
    aux = np.einsum('i,j',x,y)/(1+np.einsum('i,j',x,y))
    gradin = np.sum(kin-aux.sum(axis=1))+extra1
    gradout = np.sum(kout-aux.sum(axis=0))+extra1
    return gradin,gradout
