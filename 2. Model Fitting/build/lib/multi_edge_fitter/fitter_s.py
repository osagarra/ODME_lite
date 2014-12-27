""" Fitter for multi-edge networks with given degree sequence """
import numpy as np
from scipy import optimize as opt


            


def balance_xy(sin,sout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,**kwargs):
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
    """
    n = len(sin)
    T = sout.sum()
    a=kwargs.get('x_ini',sout/np.sqrt(T))
    b=kwargs.get('y_ini',sin/np.sqrt(T))
    inds_in = np.where(sin==0)
    inds_out = np.where(sout==0)
    afake = a
    bfake = b
    reps = 0
    if not selfs:
        while True:
			b = sin/(np.einsum('i->',a)-a)
			a = sout/(np.einsum('i->',b)-b)
			a[inds_out]=0
			b[inds_in]=0
			#if reps>10:
			tolb = np.max(np.abs(bfake-b))
			tola = np.max(np.abs(afake-a))
			if  tolb< tol and tola<tol:
				break
				if verbose:
					print "took %d reps, tolb:%f, tola :%f" % (reps,tolb,tola)
			if reps>maxreps:
				print "Algorithm did not converge after %d reps" % maxreps
				break
			afake = a
			bfake = b
			reps +=1
	else: # it is analytical!
		a = sout/np.sqrt(T)
		b = sin/np.sqrt(T)
    return a,b
     
##############################     



#### distance checkers ####
def dist_check_s(x,y,sout,sin,selfs=True):
    """
        Computes distance between prediction and reality for distances
        x,y --> arrays of lagrange multpliers
        kout,kin -- > real degree sequences
        self: True for self loops
    """
    if not selfs:
		extra1 = (x*y).sum()
    else:
		extra1=0
    aux = np.einsum('i,j',x,y)
    gradin = np.sum(sin-aux.sum(axis=1))+extra1
    gradout = np.sum(sout-aux.sum(axis=0))+extra1        
    return gradin,gradout
    
######## deprecated for the moment ##########
#def dist_check_s(x,y,sout,sin,selfs=True,indist=False):
    #"""
        #Computes distance between prediction and reality for distances
        #x,y --> arrays of lagrange multpliers
        #kout,kin -- > real degree sequences
        #self: True for self loops
        #Indist: If indist true, then apply model of undistinguishable weights
    #"""
    #if not indist:
        #if not selfs:
            #extra1 = np.einsum('j,j',x,y) # diagonal
        #else:
            #extra1=0
        #aux = np.einsum('i,j',x,y)
    #else:
        #if not selfs:
            #extra1 = np.sum(np.einsum('j,j->j',x,y)/(1-np.einsum('j,j->j',x,y))) # diagonal
        #else:
            #extra1=0
        #aux = np.einsum('i,j',x,y)/(1-np.einsum('i,j',x,y))
    #gradin = np.sum(sin-aux.sum(axis=1))+extra1
    #gradout = np.sum(sout-aux.sum(axis=0))+extra1        
    #return gradin,gradout
    
    
