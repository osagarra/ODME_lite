"""
 Fitter for fixed strength sequence and fixed number of edges.
"""

## Imports ##
from scipy import optimize as opt
import numpy as np




"""
#############################
#### Balancing approach: Preferred method #########
#############################
"""

import numpy as np
from scipy import optimize as opt


def balance_xy(lam,sin,sout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,print_tol=False,**kwargs):
    """ Balances the strength equations
        Input:
            Lam: Lambda value for lagrange multiplier on number of edges
            Sin: Incoming strength sequence (N)
            Sout: Outgoing strength sequence (N)
            tol: Maximum tolerance [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            print_tol [Bool]: Print tolerance convergence
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
    """
    # options
    if lam<=0:
        raise ValueError("Lambda must be positive")
    n = len(sin)
    T = sout.sum()
    a=kwargs.get('x_ini',sout/np.sqrt(T))
    b=kwargs.get('y_ini',sin/np.sqrt(T))
    alf = np.sqrt(sout.max()*sin.max()*1./sout.sum()) # constant to avoid overflow
    if 'x_ini' in kwargs.keys() or 'y_ini' in kwargs.keys():
        a = a*alf
        b = b*alf
    # pre - conditioning
    inds_in = np.where(sin==0)
    inds_out = np.where(sout==0)
    afake = a
    bfake = b
    reps = 0
    if not selfs:
        while True:
            aux = np.exp(np.einsum('i,j',a/alf,b/alf))
            b = sin/(np.einsum('i,ij',lam*a/(alf*alf),aux/(lam*(aux-1)+1))-lam*a/(alf*alf)*np.exp(a/alf*b/alf)/(lam*
                (np.exp(a/alf*b/alf)-1)+1))
            #np.einsum('i,ii->i',lam*a,aux/(lam*(aux-1)+1)))
            aux = np.exp(np.einsum('i,j',a/alf,b/alf))
            a = sout/(np.einsum('j,ij',lam*b/(alf*alf),aux/(lam*(aux-1)+1))-lam*b/(alf*alf)*np.exp(a/alf*b/alf)/(lam*
                (np.exp(a/alf*b/alf)-1)+1))
            #np.einsum('j,jj->j',lam*b,aux/(lam*(aux-1)+1)))
            a[inds_out]=0
            b[inds_in]=0
            #if reps>10:
            tolb = np.max(np.abs(bfake-b))
            tola = np.max(np.abs(afake-a))
            if print_tol: print "Delta_x:%r Delta_y:%r" % (tola,tolb)
            if  tolb< tol and tola<tol:
                break
                if verbose:
                    print "took %d reps, tolb:%f, tola :%f" % (reps,tolb,tola)
            if reps>maxreps:
                print "Algorithm did not converge after %d reps" % maxreps
                break
            if not all(np.isfinite(a)) or not all(np.isfinite(b)):
                raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
            afake = a
            bfake = b
            reps +=1
    else:
        while True:
            aux = np.exp(np.einsum('i,j',a/alf,b/alf))
            b = sin/(np.einsum('i,ij',lam*a/(alf*alf),aux/(lam*(aux-1)+1)))
            b[inds_in]=0
            aux = np.exp(np.einsum('i,j',a/alf,b/alf))
            a = sout/(np.einsum('j,ij',lam*b/(alf*alf),aux/(lam*(aux-1)+1)))
            a[inds_out]=0
            b[inds_in]=0
            #if reps>10:
            tolb = np.max(np.abs(bfake-b))
            tola = np.max(np.abs(afake-a))
            if print_tol: print "Delta_x:%r Delta_y:%r" % (tola,tolb)
            if  tolb< tol and tola<tol:
                break
                if verbose:
                    print "took %d reps, tolb:%f, tola :%f" % (reps,tolb,tola)
            if reps>maxreps:
                print "Algorithm did not converge after %d reps" % maxreps
                break
            if not all(np.isfinite(a)) or not all(np.isfinite(b)):
                raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
            afake = a
            bfake = b
            reps +=1
    return a/alf,b/alf
     
##############################     

################
### fine grain fitting
################
def fit_lambda(sin,sout,E,tol_s=1e-9,tol_gamma=1e-5,maxreps=1000,max_rej=100,fact=3.,selfs=True,verbose=False,**kwargs):
    """
         Minimizes lambda distance between E and <E> while keeping strength constraints
         input:
             sin,sout (strength sequences) [N arrays]
             E: total binary edges [int]
             tol_s: Tolerance for strength balancing [float]
             tol_gamma: Tolerance for gamma [float]
             mareps: Maximum iterations [int]
             max_rej: maximum rejections
             fact: Factor to approach delta (the bigger, the more likely not to converge)
             selfs: Accept slef loops? [Boolean]
             kwargs:
                 lambda_ini: initial guess for distance [float]
                 x_ini: initial guess for x
                 y_ini: initial guess for y
                 delta0 :initial jump for gamma for the algorithm [float]
                 lambda_min : minimum value of lambda [float] [default=0]
        output:
            x,y,gamma
    """
    ### initial opts ###
    gamma_init = kwargs.get('lambda_ini',sin.sum()*1./E)
    delta0 = kwargs.get('delta0',gamma_init/2.)
    gamma_min = kwargs.get('lambda_min',0)
    x_ini = kwargs.get('x_ini',sout/np.sqrt(sout.sum()))
    y_ini = kwargs.get('y_ini',sin/np.sqrt(sin.sum()))
    print_tol =kwargs.get('print_tol',False)
    T=sin.sum()
    ####################
    print ".... Lambda fitter ....."
    ####################
    gamma = gamma_init
    if verbose:
        print "Preconditioning x,y..."
    x,y = balance_xy(gamma,sin,sout,tol_s,maxreps,selfs=selfs,x_ini=x_ini,y_ini=y_ini,verbose=verbose,print_tol=print_tol)
    delta_E=dist_check_E(x,y,gamma,E,selfs=selfs)
    delta_s=dist_check_s(x,y,gamma,sout,sin,selfs=selfs)
    print "S_in grad:%f S_out grad:%f E_grad:%f " % (delta_s[0],delta_s[1],delta_E)
    delta = delta0
    reps=0
    s=1
    while reps<maxreps and np.abs(delta_E)>tol_gamma:
        new_gamma=gamma+delta*s # lambda must be positive
        if new_gamma<=gamma_min:
            i=1
            while new_gamma<=gamma_min and delta/i>tol_gamma:
                new_gamma=gamma+delta*s/i # lambda must be positive
                i+=1
            if new_gamma<=gamma_min:
                print "Lambda <=0,i:%d, new_lambda:%f. Aborting"% (i,new_gamma)
                break
        if verbose:
            print "Balancing x,y... for lambda : %r" % new_gamma
        x,y = balance_xy(new_gamma,sin,sout,tol_s,maxreps,selfs=selfs,y_ini=y,x_ini=x,verbose=verbose,print_tol=print_tol)
        new_delta_E=dist_check_E(x,y,new_gamma,E,selfs=selfs)
        delta_s = dist_check_s(x,y,gamma,sout,sin,selfs=selfs)
        if np.abs(new_delta_E)<np.abs(delta_E): #accepted
            if np.sign(new_delta_E) != np.sign(delta_E): # if change of sign, change direction!
                s=-s
            gamma = new_gamma                
            delta_E = new_delta_E
            #delta=np.abs(gamma-new_gamma)
            delta = delta*fact
            if verbose:
                print "Accepted! lambda:%r ||delta_E:%f delta_s:%f"  % (gamma,delta_E,delta_s[0])
        else: #rejected, reduce step and change direction
            delta = delta/fact
            s=-s
            delta_gamma=new_gamma-gamma
            if verbose:
                print "Rejected! movement: %r delta_lambda:%r" % (delta_E-new_delta_E,new_gamma-gamma)
            if np.abs(delta_gamma)<1e-14:
                break
        reps+=1
    x,y = balance_xy(gamma,sin,sout,tol_s,maxreps,selfs=selfs,y_ini=y,x_ini=x,verbose=verbose,print_tol=print_tol)
    delta_s = dist_check_s(x,y,gamma,sout,sin,selfs=selfs)
    print "### All deltas: sin:%r sout:%r E:%r for lambda:%r ###" % (delta_s[0].max(),delta_s[1].min(),
        dist_check_E(x,y,gamma,E,selfs=selfs),gamma)
    return x,y,gamma



#### Check functions (distance functions ) ####

def dist_check_s(x,y,lam,sout,sin,selfs=True):
    """
        Computes distance between prediction and reality for strengths
        x,y --> arrays of lagrange multpliers (kength N)
        lam --> Float, lagrange multiplier associated with number of binary edges
        sout,sin -- > real degree sequences (length N)
        self: True for self loops
    """
    xy_exp = np.exp(np.einsum('i,j',x,y))
    if not selfs:
        extra1 = lam*x*y*np.exp(x*y)/(lam*(np.exp(x*y)-1.)+1.) # diagonal
    else:
        extra1=0
    delta_x = sout - x*lam*np.einsum('j,ij',y,xy_exp/(lam*(xy_exp-1.)+1.)) + extra1
    delta_y = sin - y*lam*np.einsum('i,ij',x,xy_exp/(lam*(xy_exp-1.)+1.)) +  extra1
    return delta_x.sum(),delta_y.sum()

def dist_check_E(x,y,lam,E,selfs=True):
    """
        Computes distance between prediction and reality for binary edges
        x,y --> arrays of lagrange multpliers
        lam --> Float, lagrange multiplier associated with number of edges
        sout,sin -- > real degree sequences
        E --> total number of binary edges
        self: True for self loops
    """
    xy_exp = np.exp(np.einsum('i,j',x,y))
    if not selfs:
        extra1 = lam*np.einsum('ii->',(xy_exp-1.)/(lam*(xy_exp-1.)+1.))
    else:
        extra1=0
    delta_E = E- lam*np.einsum('ij->',(xy_exp-1)/(lam*(xy_exp-1)+1)) + extra1
    return delta_E



"""
#############################
#### Brute force maximization: Preferred as fine tunning method #########
#############################
"""
### Notes ###
"""
3 sets of values to be optimized:
    x [N array],y [N array] and lambda [single number]
    All floats.
    
Loglikelyhood to fix:
    L =  - Log denom_{ij} (x,y,lambda) + L(lambda) + L(sin,y) + L(sout,x)
"""


class Multi_Edge_Maximizer_E(object):
    """
        General maximizer class for fixed num of edges, sin,sout
            Allows to define global variables!
    """
    def __init__(self,sin,sout,E,selfs=True,**kwargs):
        """
            General initialization. Needs at least strength sequence
        """
        # Params #
        if not selfs:
            raise NotImplementedError("Not implemented without self-loops")
        if sin.sum()!=sout.sum() or len(sin)!=len(sout):
            raise TypeError("Incompatible strength lists provided")
        self.sin = sin
        self.sout = sout
        self.E = E
        self.T = sin.sum()
        self.N = len(sin)
        self.cnt = kwargs.get('cnt',0) #extra log term
        # Set init_conds #
        self.x = kwargs.get('x',self.init_conds()[0])
        self.y = kwargs.get('y',self.init_conds()[1])
        self.lam = kwargs.get('lam',self.init_conds()[2])
    def init_conds(self):
        """
            Sets init_conds if not given
        """
        x = self.sout*1./np.sqrt(self.T) + 1e-14
        y = self.sin*1./np.sqrt(self.T) + 1e-14
        #lam = 1./self.E
        lam=self.E/self.T
        return x,y,lam
    def var_adapt(self):
        """ Returns array to do computations on scipy all grouped """
        xx = np.zeros(2*self.N+1)
        xx[:self.N] = self.x
        xx[self.N:-1] = self.y
        xx[-1] = self.lam
        return xx
    def var_update(self,xx):
        """ Inverse of var_Adapt """
        self.x = xx[:self.N] 
        self.y = xx[self.N:-1]
        self.lam = xx[-1]

    def var_readapt(self,xx):
        """ Inverse of var_Adapt """
        return np.r_[self.x,self.y,self.lam]
        
    def check_s(self):
        """ Checks strength distance """
        #xx = np.zeros(2*self.N+1)
        #xx[:self.N] = self.x
        #xx[self.N:-1] = self.y
        #xx[-1] = self.lam
        #return check_s(xx,self.sin,self.sout,self.N)
        return dist_check_s(self.x,self.y,self.lam,self.sout,self.sin)

    def check_E(self):
        #xx = np.zeros(2*self.N+1)
        #xx[:self.N] = self.x
        #xx[self.N:-1] = self.y
        #xx[-1] = self.lam
        #return check_E(xx,self.sin,self.sout,self.E,self.N)
        return dist_check_E(x,y,lam,E,selfs=True)
        
    def minimize(self,meth='TNC',**kwargs): # meths: TNC, L-BFGS-B, COBYLA preferred
        ### better thing: use a cocktail
            # TNC
            # L-BFGS
            # COBYLA
        ### set options ###
        tol = kwargs.get('tol',1e-12)
        disp = kwargs.get('disp',True)
        bmin = kwargs.get('bmin',0)
        bmax = kwargs.get('bmax',None)
        ftol = kwargs.get('ftol',1e-12)
        gtol = kwargs.get('gtol',1e-9)
        bounds = np.reshape(([1e-14,None]*(2*self.N+1)),(2*self.N+1,2))
        bounds[-1][-1] = 1.
        opts = {'xtol':tol,'disp':disp,'ftol':ftol,'gtol':gtol}
        if meth=='COBYLA':
            opts['tol']=tol
        elif meth=='TNC':
            opts['xtol']=-1 # machine prec
        elif meth=='L-BFGS-B':
            opts['ftol']= ftol
            opts['gtol']= gtol
        else:
            print "Recommended only use: TNC, L-BFGS-B or COBYLA or SLSQP"
        #meths = {
        #    'BFGS':,
        #    'TNC':,
        #    'COBYLA':}
        #minimizer = opt.__getattribute__(meths[meth])
        minimizer=opt.minimize
        res = minimizer(loglikelyhood,self.var_adapt(),args=(self.sin,self.sout,self.E,self.N),method=meth,jac=loglikelyhood_grad,
            hess = loglikelyhood_hess, bounds=bounds, options = opts)
        dx,dy = check_s(res['x'],self.sin,self.sout,self.N)
        #dx/=self.sout
        #dy/=self.sin
        dx = dx[np.isfinite(dx)]
        dy = dy[np.isfinite(dy)]
        ch2 = check_E(res['x'],self.sin,self.sout,self.E,self.N)
        print "Result of minimization: sin:%f+-%f sout:%f+-%f E:%f" % (np.mean(dx),
            np.std(dx),np.mean(dy),np.std(dy),ch2)
        print "Final value of function (witch cosntant term): %f" % (loglikelyhood(res['x'],self.sin,self.sout,self.E,self.N)+self.cnt)
        if res['success']:
            self.var_update(res['x'])
        else:
            print "Minimization failed!"
        return res

#######  Log likelyhood functions ###


def loglikelyhood(xx,sin,sout,E,N):
    """ Computes minus loglikelyhood """
    x = xx[:N]
    y = xx[N:-1]
    lam = xx[-1]
    yy = np.log(y)
    inds1 = np.where(y==0)
    xx = np.log(x)
    inds2 = np.where(x==0)
    xx[inds2] = 0 # fix log(0)
    yy[inds1] = 0
    L1 = np.einsum('i,i',sin,yy)
    L2 = np.einsum('i,i',sout,xx)
    L3 = E*np.log(lam)
    L4 = np.einsum('ij->',np.log(1.+lam*np.einsum('i,j',x,y)))
    return L4-L1-L2-L3
    
def loglikelyhood_grad(xx,sin,sout,E,N):
    """ Computes gradient vector """
    x = xx[:N]
    y = xx[N:-1]
    lam = xx[-1]
    xy_exp = np.exp(np.einsum('i,j',x,y))
    denom =  1.+lam*(xy_exp-1)
    g1 = - sout/x + lam*np.einsum('j,ij',y,xy_exp/denom)
    g2 = - sin/y + lam*np.einsum('i,ij',x,xy_exp/denom)
    g3 = -E/lam + np.einsum('ij->', (xy_exp-1)/denom)
    return np.r_[g1,g2,g3]
    
def loglikelyhood_hess(xx,sin,sout,E,N):
    """ Computes Hessian matrix """
    x = xx[:N]
    y = xx[N:-1]
    lam = xx[-1]
    ## common terms ##
    xy = np.einsum('i,j',x,y)
    xy_exp = np.exp(xy)
    denom =  1.+lam*(xy_exp-1)
    denom = denom*denom
    ## diagonal ##
    Hxx = sout/(x*x) - (lam-1)*lam*np.einsum('j,ij',y*y,xy_exp/denom)
    Hxx = np.diag(Hxx)
    Hyy = sin/(y*y) - (lam-1)*lam*np.einsum('i,ij',x*x,xy_exp/denom)
    Hyy = np.diag(Hyy)
    Hll = E/lam/lam - np.einsum('ij->',(xy_exp-1)*(xy_exp-1)/denom)
    ## off diagonal ##
    Hxl = np.einsum('j,ij',y,xy_exp/denom)
    Hyl = np.einsum('i,ij',x,xy_exp/denom)
    Hxy = lam*xy_exp*((1-lam)*xy + lam*(xy_exp-1)+1)/denom
    Hll = np.r_[Hxl,Hyl,Hll]
    ## Stack together! ##
    return np.vstack((np.hstack((np.vstack((np.hstack((Hxx,Hxy)),np.hstack((Hxy,Hyy)))),np.reshape(Hll[:-1],(2*N,1)))),Hll))


#### Check functions (distance functions ) ####
def check_s(xx,sin,sout,N,selfs=True):
    x = xx[:N]
    y = xx[N:-1]
    lam = xx[-1]
    xy_exp = np.exp(np.einsum('i,j',x,y))
    delta_x = sout - x*lam*np.einsum('j,ij',y,xy_exp/(lam*(xy_exp-1)+1))
    delta_y = sin - y*lam*np.einsum('i,ij',x,xy_exp/(lam*(xy_exp-1)+1))
    return delta_x,delta_y

def check_E(xx,sin,sout,E,N,selfs=True):
    x = xx[:N]
    y = xx[N:-1]
    lam = xx[-1]
    xy_exp = np.exp(np.einsum('i,j',x,y))
    delta_E = E- lam*np.einsum('ij->',(xy_exp-1)/(lam*(xy_exp-1)+1))
    return delta_E


#### check functions ####
### Fitters to use ###
# BFGS : not using hessian
# Try SQLQST and finally try krkov scipy and Andersen (http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton_krylov.html#scipy.optimize.newton_krylov)
