""" Fitter for multi-edge network with given degree sequence """
#from numbapro import guvectorize,autojit
#from numba import autojit, double, int32, jit
import numpy as np
#import distances as dist
from scipy import optimize as opt

def extra_selector(selfs):
    """
    Selects appropriate extra term for selfloops
    """
    if selfs:
        extra_k = lambda c,d,aux : 0
    else:
        extra_k = lambda c,d,aux: np.diag(c*d/aux)
    return extra_k


def balance_xy(kin,kout,tol=1e-9,tol_c=1e-9,maxreps=10000,verbose=False,selfs=True,print_c=True,**kwargs):
    """ Balances the degree equations
        Input:
            kin: Incoming degree sequence (N)
            kout: Outgoing degree sequence (N)
            tol: Maximum tolerance convergence [float]
            tol_c: Maximum tolerance on constraints [float]
            print_c : Print cost distance every "act" steps [Bool]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
                print_tol: Print tolerance update
                act: Compute cost gradient every act steps (default=100)
    """
    n = len(kin)
    E = kout.sum()
    act = kwargs.get('act',100)
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
    delta_k=dist_check_k(a,b,kout,kin,selfs)
    extra_k = extra_selector(selfs)
    print "### Fixing k ###"
    print "Initial errors | K_in:%f K_out:%f " % (delta_k[0],delta_k[1])
    if print_c:
        print "## Errors: \t\t\t\t || Convergence: ##"
        print "-----------------------------------------------------------"
	while True:
		aux = 1.+np.einsum('i,j',a,b)
		b = kin/(np.einsum('i,ij',a,1./aux)-extra_k(a,b,aux)/b)
		b[inds_in]=0
		aux = 1.+np.einsum('i,j',a,b)
		a = kout/(np.einsum('j,ij',b,1./aux)-extra_k(a,b,aux)/a)
		a[inds_out]=0
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
		if print_c and reps%act==0:
			da,db = dist_check_k(a,b,kout,kin,selfs)
			if np.abs(da) < tol_c and np.abs(db) < tol_c:
				if verbose:
					print "took %d reps, tola :%f, tolb:%f" % (reps,tola,tolb)
				break
			print "Delta k:%r || Delta_x:%r Delta_y:%r " % (da,tola,tolb)
		afake = a
		bfake = b
		reps +=1
    return a,b
     

def rho_calculator(av_w,agg=False,indist=False,M=1,**kwargs):
    """
        Computes the lagrange multiplier for the case of fixed k and a desired average existing weight av_w = T/<E> (where T is the total number of events).
        It does so by inverting the equation:
            <t | t>0 > (rho) -> rho(av_w) according to each case.
            
        Input:
            av_w: total number of events divided by total number of binary edges.
            agg: Set to true to analyze the case of aggregation of binary networks
            indist: Set to true to analyze the case of weighted networks 
            M: Number of aggregated layers (ignored for the ME case)
        Output:
            rho: Lagrange multiplier for the generation of random network ensembles.
    """
    av_w = 1.0*av_w
    if not indist and not agg: # case ME, fully analytical
        from scipy import special
        x = -av_w*np.exp(-av_w)
        rho = (special.lambertw(x,0)+av_w).real
    else:
        from scipy import optimize as opt
        if indist:
            if M==1: # fully analytical
                rho = 1. - 1/av_w
            else: # must solve equation numerically
                rho = opt.brentq(rho_ZINB,1e-16,1-1e-14,args=(av_w,M),**kwargs)
        else:
            if M==1:  # fully analytical
                rho = 1 #pretty absurd, it is a binary network!
            else:  # must solve equation numerically
                rho = opt.brentq(rho_ZIB,1e-12,2e20,args=(av_w,M),**kwargs)
    return rho

def rho_ZINB(x,av_w,M):
    return M*x/(1-x)*1/(1-(1-x)**M)-av_w
def rho_ZIB(x,av_w,M):
    return M*x/(1+x)*1/(1-(1+x)**(-M))-av_w
    
    

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
        #extra1 = np.sum(np.einsum('j,j->j',x,y)/(1+np.einsum('j,j->j',x,y))) # diagonal
        extra1 = x*y/(1.+x*y) # diagonal
    else:
        extra1=0
    aux = np.einsum('i,j',x,y)
    aux = aux/(1.+aux)
    gradin = np.abs(kin-aux.sum(axis=0)+extra1)
    gradout = np.abs(kout-aux.sum(axis=1)+extra1)
    return gradin.sum(),gradout.sum()


