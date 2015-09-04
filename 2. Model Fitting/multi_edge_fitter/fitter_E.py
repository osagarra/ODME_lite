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


def aux_selector(agg=False):
    """
    Selects appropriate aux function for the balancer
    """
    if agg:
        aux_f = lambda a,b,alf,M : np.power((np.einsum('i,j',a/alf,b/alf)+1),-M)
        corr = lambda a,b,alf : np.einsum('i,j',a/alf,b/alf)+1.
    else:
        aux_f = lambda a,b,alf,M : np.exp(np.einsum('i,j',-a/alf,b/alf)) # ignores M
        corr = lambda a,b,alf : 1.
    return aux_f,corr
        
def extra_selector(selfs):
    """
    Selects appropriate extra term for selfloops
    """
    if selfs:
        extra_s = lambda a,b,aux,lam,aux3 : 0
    else:
        extra_s = lambda a,b,aux,lam,aux3: np.diag(lam*a*b/(lam*(1.-aux)+aux)/aux3)
    return extra_s


def balance_xy(lam,sin,sout,tol=1e-9,tol_c=1e-9,maxreps=10000,verbose=False,selfs=True,print_tol=False,print_c=False,agg=False,M=1,**kwargs):
    """ Balances the strength equations
        Input:
            Lam: Lambda value for lagrange multiplier on number of edges
            Sin: Incoming strength sequence (N)
            Sout: Outgoing strength sequence (N)
            tol: Maximum tolerance on convergence [float]
            tol_c: Maximum tolerance permitted for the worse constraint-equation [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            print_tol [Bool]: Print tolerance convergence
            Agg: If the weighted network comes from aggreagation of binary networks set to True
            M: Number of layers (only in the case of aggregated network)
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
                act: Print error every act steps (if print_c is set to True)
                
    """
    ## add precision ##
    sin = np.array(sin,dtype=np.float128)
    sout = np.array(sout,dtype=np.float128)
    lam = np.float128(lam)
    ## print opts ##
    if verbose or print_tol:
        print "#### Total binary & Strength balancer| Agg: %r Layers:%d #####" % (agg,M)
    n = len(sin)
    if n>1000:
        act = 10
    elif n<1000 and n>500:
        act = 50
    else:
        act = 500
    ## options ##
    if lam<=0:
        raise ValueError("Lambda must be positive")
    if not agg:
        M=1
    act = kwargs.get('act',act)
    T = sout.sum()
    alf0 = np.sqrt(M*sout.max()*sin.max()*1./T)
    alf = kwargs.get('alf',alf0)
    #alf = kwargs.get('alf',np.sqrt(sout.max()*sin.max()*1./T))
    ## arranging variables & pre-conditioning ##
    #a=kwargs.get('x_ini',sout/np.sqrt(T)*0.1)
    #b=kwargs.get('y_ini',sin/np.sqrt(T)*0.1)
    a=kwargs.get('x_ini',np.ones(len(sout))*0.9)
    b=kwargs.get('y_ini',np.ones(len(sin))*0.9)
    inds_in = np.where(sin==0)
    inds_out = np.where(sout==0)
    a[inds_out]=0
    b[inds_in]=0
    a = np.array(a,dtype=np.float128)
    b = np.array(b,dtype=np.float128)
    afake = a
    bfake = b
    reps = 0
    if 'x_ini' in kwargs.keys() or 'y_ini' in kwargs.keys():
        a = a*alf
        b = b*alf
    delta_s = dist_check_s(a/alf,b/alf,lam,sout,sin,selfs=selfs,agg=agg,M=M)
    ## selfs or agg ##
    extra_s = extra_selector(selfs)
    aux_f,corr = aux_selector(agg)
    ## all ready, let's go! ##
    if verbose: 
        print "### Fixing s ###"
        print "Initial errors | S_in grad:%f S_out grad:%f " % (delta_s[1],delta_s[0])
    if np.max(np.abs(delta_s)) < tol_c:
        return a/alf,b/alf
    if print_c:
        print "## Errors: \t\t\t\t || Convergence: ##"
        print "-----------------------------------------------------------"
    while True:
        aux = aux_f(a,b,alf,M)
        b = sin*(alf*alf) /M/(np.einsum('i,ij',lam*a,1./corr(a,b,alf)/(lam*(1.-aux)+aux))-extra_s(a,b,aux,lam,corr(a,b,alf))/b)
        b[inds_in]=0
        aux = aux_f(a,b,alf,M)
        a = sout*(alf*alf)/M/(np.einsum('j,ij',lam*b,1./corr(a,b,alf)/(lam*(1.-aux)+aux))-extra_s(a,b,aux,lam,corr(a,b,alf))/a)
        a[inds_out]=0
        #if reps>10:
        tolb = np.max(np.abs(bfake-b))
        tola = np.max(np.abs(afake-a))
        if print_tol: print "Delta_x:%r Delta_y:%r" % (tola,tolb)
        if  tolb< tol and tola<tol:
            if verbose:
                print "took %d reps, tolb:%f, tola :%f" % (reps,tolb,tola)
            break
        if reps>maxreps:
            print "Algorithm did not converge after %d reps" % maxreps
            break
        if not all(np.isfinite(a)) or not all(np.isfinite(b)):
            raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
        if print_c and reps%act==0:
            da,db = dist_check_s(a/alf,b/alf,lam,sout,sin,selfs=selfs,agg=agg,M=M)
            if np.abs(da) < tol_c:
                if verbose:
                    print "took %d reps, tola :%f, tolb:%f" % (reps,tola,tolb)
                break
            print "Delta s:%r || Delta_x:%r Delta_y:%r " % (da,tola,tolb)
        afake = a
        bfake = b
        reps +=1

    return a/alf,b/alf
     
##############################     

################
### fine grain fitting
################
def fit_lambda(sin,sout,E,tol_s=1e-9,tol_gamma=1e-5,maxreps=1000,max_rej=100,fact=3.,selfs=True,verbose=False,agg=False,M=1,**kwargs):
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
             Agg: If the weighted network comes from aggreagation of binary networks set to True
             M: Number of layers (only in the case of aggregated network)
             kwargs:
                 lambda_ini: initial guess for distance [float]
                 x_ini: initial guess for x
                 y_ini: initial guess for y
                 delta0 :initial jump for gamma for the algorithm [float]
                 lambda_min : minimum value of lambda [float] [default=0]
                 tol_c : Tolerance with respect to constraints
        output:
            x,y,gamma
    """
    ### initial opts ###
    gamma_init = kwargs.get('lambda_ini',sin.sum()*1./E)
    delta0 = kwargs.get('delta0',gamma_init/2.)
    gamma_min = kwargs.get('lambda_min',1e-15)
    x_ini = kwargs.get('x_ini',sout/np.sqrt(sout.sum()))
    y_ini = kwargs.get('y_ini',sin/np.sqrt(sin.sum()))
    print_tol =kwargs.get('print_tol',False)
    print_c =kwargs.get('print_c',False)
    tol_c =kwargs.get('tol_c',1e-9)
    act = kwargs.get('act',100)
    T=sin.sum()
    ####################
    print ".... Lambda fitter ....."
    ####################
    gamma = gamma_init
    if verbose:
        print "Preconditioning x,y..."
    x,y = balance_xy(gamma,sin,sout,tol_s,maxreps=maxreps,selfs=selfs,x_ini=x_ini,y_ini=y_ini,verbose=verbose,print_c=print_c,tol_c=tol_c,print_tol=print_tol,act=act,agg=agg,M=M)
    if verbose:
        print "Optimizing lambda..."
    delta_E=dist_check_E(x,y,gamma,E,selfs=selfs,agg=agg,M=M)
    delta_s=dist_check_s(x,y,gamma,sout,sin,selfs=selfs,agg=agg,M=M)
    print "S_in grad:%r S_out grad:%r E_grad:%r " % (delta_s[0],delta_s[1],delta_E)
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
        x,y = balance_xy(new_gamma,sin,sout,tol=tol_s,tol_c=tol_c,maxreps=maxreps,verbose=verbose,selfs=selfs,print_tol=print_tol,print_c=print_c,agg=agg,M=M,x_ini=x,y_ini=y,act=act)
        new_delta_E=dist_check_E(x,y,new_gamma,E,selfs=selfs,agg=agg,M=M)
        if np.abs(new_delta_E)<np.abs(delta_E): #accepted
            if np.sign(new_delta_E) != np.sign(delta_E): # if change of sign, change direction!
                s=-s
            gamma = new_gamma                
            delta_E = new_delta_E
            #delta=np.abs(gamma-new_gamma)
            delta = delta*fact
            delta_s = dist_check_s(x,y,gamma,sout,sin,selfs=selfs,agg=agg,M=M)
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
    
    #x,y = balance_xy(gamma,sin,sout,tol=tol_s,tol_c=tol_c,maxreps=maxreps,verbose=verbose,selfs=selfs,print_tol=print_tol,print_c=print_c,agg=agg,M=M,x_ini=x,y_ini=y,act=act)
    x,y = balance_xy(gamma,sin,sout,tol=tol_s,tol_c=tol_c,maxreps=maxreps,verbose=verbose,selfs=selfs,print_tol=print_tol,print_c=print_c,agg=agg,M=M,x_ini=x,y_ini=y,act=act)
    if verbose:
        delta_s=dist_check_s(x,y,gamma,sout,sin,selfs=selfs,agg=agg,M=M)
        print "### All deltas: sin:%r sout:%r E:%r for lambda:%r ###" % (delta_s[0],delta_s[1],
            dist_check_E(x,y,gamma,E,selfs=selfs,agg=agg,M=M),gamma)
    return x,y,gamma



#### Check functions (distance functions ) ####

def dist_check_s(x,y,lam,sout,sin,selfs=True,agg=False,M=1,absolute=True):
    """
        Computes absolute distance between prediction and reality for strengths
        x,y --> arrays of lagrange multpliers (kength N)
        lam --> Float, lagrange multiplier associated with number of binary edges
        sout,sin -- > real degree sequences (length N)
        self: True for self loops
        absolute: True for absolute error
        Agg: If the weighted network comes from aggreagation of binary networks set to True
        M: Number of layers (only in the case of aggregated network)

    """
    if not agg:
        xy_exp = np.exp(-np.einsum('i,j',x,y))
        if not selfs:
            extra1 = lam*x*y/(lam*(1.-np.exp(-x*y))+np.exp(-x*y)) # diagonal
        else:
            extra1=0
        delta_x = sout - x*lam*np.einsum('j,ij',y,1./(lam*(1.-xy_exp)+xy_exp)) + extra1
        delta_y = sin  - y*lam*np.einsum('i,ij',x,1./(lam*(1.-xy_exp)+xy_exp)) + extra1
    else:
        xy = 1+np.einsum('i,j',x,y)
        if not selfs:
            extra1 = lam*M*x*y/(1+x*y)/((1+x*y)**(-M)+lam*(1.-(1+x*y)**(-M)))
        else:
            extra1 = 0
        delta_x = sout - M*x*lam*np.einsum('j,ij',y,1./xy/(lam*(1.-xy**(-M))+xy**(-M))) + extra1
        delta_y = sin  - M*y*lam*np.einsum('i,ij',x,1./xy/(lam*(1.-xy**(-M))+xy**(-M))) + extra1        
    if absolute:
        delta_x = np.abs(delta_x)
        delta_y = np.abs(delta_y)
    return delta_x.sum(),delta_y.sum()

def dist_check_E(x,y,lam,E,selfs=True,agg=False,M=1):
    """
        Computes distance between prediction and reality for binary edges
        x,y --> arrays of lagrange multpliers
        lam --> Float, lagrange multiplier associated with number of edges
        sout,sin -- > real degree sequences
        E --> total number of binary edges
        self: True for self loops
        Agg: If the weighted network comes from aggreagation of binary networks set to True
        M: Number of layers (only in the case of aggregated network)
    """
    if not agg:
        xy_exp = np.exp(-np.einsum('i,j',x,y))
        if not selfs:
            extra1 = lam*np.einsum('ii->',(1.-xy_exp)/(lam*(1.-xy_exp)+xy_exp))
        else:
            extra1=0
        delta_E = E- lam*np.einsum('ij->',(1.-xy_exp)/(lam*(1.-xy_exp)+xy_exp)) + extra1
    else:
        xy = np.power(1+np.einsum('i,j',x,y),-M)
        if not selfs:
            extra1 = lam*np.einsum('ii->',(1.-xy)/(lam*(1.-xy)+xy))
        else:
            extra1=0
        delta_E = E- lam*np.einsum('ij->',(1.-xy)/(lam*(1.-xy)+xy)) + extra1
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
