""" 
    Fitter for multi-edge networks with given degree sequence and aggregated binary networks
"""
import numpy as np
from scipy import optimize as opt
import copy

#epsilon = 1e-15 # minimal value
epsilon = 0

def aux_selector(case='ME'):
    """
    Selects appropriate aux function for the balancer
    """
    if case=='B':
        corr = lambda a,b,alf : 1./(np.einsum('i,j',a/alf,b/alf)+1.)
    elif case=='ME':
        corr = lambda a,b,alf : np.ones((len(a),len(a)))
    else:
        raise ValueError("Not defined for caes other than ME or B")
    return corr
                    


def balance_xy(sin,sout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,case='ME',tol_c=1e-9,print_tol=False,print_c=False,**kwargs):
    """ Balances the degree equations
        Input:
            sin: Incoming strength sequence (N)
            sout: Outgoing strength sequence (N)
            tol: Maximum tolerance on convergence[float]
            tol_c: Maximum tolerance on constraints [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            Case: Considered case
            kwargs:
                M : Number of aggregated layers from which the weighted network is obtained [only if agg=True]
                x_ini: Initial guess for x
                y_ini: Initial guess for y
                fact: Factor to ensure convergence (by default 1)
                act: Compute cost gradient every act steps (default=100)
                inds_selfs: If selfs is true, index of the self loops
    """
    ## print opts ##
    if case == 'B':
        M = kwargs.get('M',0)
        if M<=0:
            raise ValueError('M needs to be positive integer if case is B!')
    elif case == 'ME':
        M=1
        if selfs:
			# Analytic solution
			return sout/np.sqrt(sout.sum()),sin/np.sqrt(sin.sum())
    else:
        raise ValueError("Not defined for caes other than ME or B")
    if not selfs:
        inds_selfs = kwargs.get('inds_selfs',False)
        if not inds_selfs:
            raise ValueError("If no selfloops accepted, need to pass also their indices!")
    if verbose or print_tol:
        print "#### Strength balancer| Case: %r Layers:%d #####" % (case,M)
    indist = kwargs.get('indist',False)
    act = kwargs.get('act',100)
    n = len(sin)
    if n>1000:
        act = 10
    elif n<1000 and n>500:
        act = 50
    else:
        act = 500
    T = sout.sum()
    a=kwargs.get('x_ini',sout/np.sqrt(T))
    b=kwargs.get('y_ini',sin/np.sqrt(T))
    iin = np.where(sin==0)
    iout = np.where(sout==0)
    afake = a
    bfake = b
    reps = 0
    alf = 1.
    if 'x_ini' in kwargs.keys() or 'y_ini' in kwargs.keys():
        a = a*alf
        b = b*alf
    delta_s=dist_check_s(a/alf,b/alf,sout,sin,selfs=selfs,case=case,M=M)    
    if delta_s[0] < tol_c:
        return a/alf,b/alf
    print "### Fixing s ###"
    print "Initial errors | S_in grad:%f S_out grad:%f" % (delta_s[0],delta_s[1])
    if print_c:
        print "## Errors: \t\t\t\t || Convergence: ##"
        print "-----------------------------------------------------------"
    kin = sin*1./M
    kout = sout*1./M
    ## selfs or B ##
    corr = aux_selector(case)
    while True:
        aux = corr(a,b,alf)
        if not selfs:
            aa = aux.flatten()
            aa[inds_selfs] = 0
            aux = aa.reshape(len(a),len(b))
        b = kin*(alf*alf) /np.einsum('i,ij',a,aux)
        b[iin]=0
        aux = corr(a,b,alf)
        if not selfs:
            aa = aux.flatten()
            aa[inds_selfs] = 0
            aux = aa.reshape(len(a),len(b))
        a = kout*(alf*alf)/np.einsum('j,ij',b,aux)
        a[iout]=0
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
            da,db = dist_check_s(a/alf,b/alf,sout,sin,selfs,case,M)
            if np.abs(da) < tol_c and np.abs(dc) < tol_c:
                if verbose:
                    print "took %d reps, tola :%f, tolb:%f" % (reps,tola,tolb)
                break
            print "Delta s:%r || Delta_x:%r Delta_y:%r " % (da,tola,tolb)
        afake = a
        bfake = b
        reps +=1
    return a/alf,b/alf
     
##############################     



#### distance checkers ####
def dist_check_s(x,y,sout,sin,selfs=True,case='ME',M=1):
    """
        Computes absolute distance between prediction and reality for distances
        x,y --> arrays of lagrange multpliers
        kout,kin -- > real degree sequences
        self: True for self loops
        M: Number of layers [set only if case=B]
        Case: Particular case considered
    """
    aux = np.einsum('i,j',x,y)
    if case == 'B':
        aux = M*aux/(1.+aux)
    elif case == 'W':
        aux = M*aux/(1.-aux)        
    elif case == 'ME':
        pass
    else:
        raise ValueError("Not defined for caes other than ME or B or W")
    if not selfs:
        aux[np.diag_indices(len(x))] = 0
    gradin = np.abs(sin-aux.sum(axis=0))
    gradout = np.abs(sout-aux.sum(axis=1))      
    return gradin.sum(),gradout.sum()
    

"""
#############################
#### All embeded in a class #########
#############################
"""
### Notes ###
"""
2 sets of values to be optimized:
    x [N array],y [N array]
    All floats.
    
"""


class Edge_Maximizer(object):
    """
        General maximizer class for fixed sin,sout
    """
    def __init__(self,sin,sout,selfs=False,case='W', **kwargs):
        """
            General initialization. Needs at least strength sequence
        """
        if case != 'B' and case != 'W' and case != 'ME':
            raise ValueError("Not defined for caes other than ME or B or W")
        # main variables
        self.case = case
        self.N = len(sin)
        self.selfs=selfs
        # check on layers
        self.M = kwargs.get('M',1)
        if case == 'ME':
            self.M = 1 # ignoring layers
        if case == 'B':
            max_l = np.int(np.max((sout.max(),sin.max()))*1./len(sout))
            if max_l > self.M:
                raise ValueError("For binary case, NM_min = max(sout,sin)")
        if self.M<1:
            raise ValueError("Number of layers must be at least 1")
        # Checks on s #
        if sin.sum()!=sout.sum() or len(sin)!=len(sout):
            raise ValueError("Incompatible strength lists provided")
        # check on directed #
        if np.abs(sout-sin).sum()<epsilon:
            self.directed=False
        else:
            self.directed=True
		# sort strengths #
        if np.max(sout)>np.max(sin):
            sort_id = np.argsort(sout)[::-1]
        else:
            sort_id = np.argsort(sin)[::-1]
        self.sort_id = sort_id
        sin  = sin[sort_id]
        sout = sout[sort_id]
        # Variables, eliminating 0 entries... #
        self.iout = np.where(sout==0)
        self.iin = np.where(sin==0)
        if set(self.iout[0]).intersection(self.iin[0]):
            print "\t...Some (0,0) strength pairs detected. Solving issue, but check your strength list!..."
        self.sin = copy.deepcopy(sin[sin>0])
        self.sout = copy.deepcopy(sout[sout>0])
        self.sout_full = copy.deepcopy(sout)
        self.sin_full  = copy.deepcopy(sin)
        # Number of reduced variables #
        self.nvarx = self.N-self.iout[0].shape[0]
        self.nvary = self.N-self.iin[0].shape[0]
        self.nvar = self.nvarx + self.nvary
        # generate inds for 0 strengths #
        self.generate_inds()
        # Set init_conds #
        x_ini = kwargs.get('x',None)
        y_ini = kwargs.get('y',None)
        if x_ini is None or y_ini is None:
            x_ini,y_ini = self.init_conds()
        else:
            if x_ini.shape != self.sout_full.shape or y_ini.shape != self.sin_full.shape:
                x_ini = self.init_conds()[0]
                y_ini = self.init_conds()[1]
                print "Non Valid initial conditions! Setting default"
            else:
                # first sort
                x_ini = x_ini[sort_id]
                y_ini = y_ini[sort_id]
                # then prune 0 values
                x_ini = x_ini[np.where(self.sout_full>0)]
                y_ini = y_ini[np.where(self.sin_full>0)]
            if not self.check_domain(np.r_[x_ini,y_ini]):
                x_ini,y_ini = self.init_conds()
                print "Initial conditions outside domain! Setting default"

        self.x = x_ini
        self.y = y_ini
    ## utils ##
    def generate_inds(self):
        """ Generates indics for non-zero entries """
        maskx = np.ones(self.N)
        maskx[self.iout] = 0
        self.indsx = np.where(maskx)
        masky = np.ones(self.N)
        masky[self.iin] = 0
        self.indsy = np.where(masky)
        self.indsxy = (np.r_[self.indsx[0],self.N+self.indsy[0]],)
        # also generates indices for self-looped entries if needed
        if not self.selfs:
            self.inds_selfs = generate_inds_selfs(self.nvarx,self.nvary)
        else:
            self.inds_selfs = False # dummy value
    ## variable adaptation ##
    def var_transform(self,xx):
        """ Adapts xx alien types to array """
        if type(xx) != np.array:
            try:
                xx = np.array(xx).flatten()
            except TypeError:
                raise TypeError("Are you sure xx is correct? Does not seem to be convertible to array...")
        return xx
    def var_update(self,xx):
        """ Inverse of var_Adapt """
        xx = self.var_transform(xx)
        if np.sum(self.check_s(xx)) < np.sum(self.check_s()): # norms are defined positive!
            self.x = copy.deepcopy(xx[:self.nvarx])
            self.y = copy.deepcopy(xx[self.nvarx:])
    def var_adapt(self,fill=False):
        """ Return variables in an array """
        if fill:
            xf = np.zeros(self.N)
            xf[self.indsx] = np.array(self.x).flatten()
            yf = np.zeros(self.N)
            yf[self.indsy] = np.array(self.y).flatten()
        else:
            xf = np.array(self.x).flatten()
            yf = np.array(self.y).flatten()
        return np.r_[xf,yf]
    def var_result(self):
        """ Returns the result with the positions as the original given list """
        trans 	= self.sort_id
        xx		= self.var_adapt(fill=True)
        resx = np.zeros(len(trans))
        resy = np.zeros(len(trans))
        for i,e in enumerate(xx[:self.N]):
			resx[trans[i]]	= e
			resy[trans[i]]	= xx[i+self.N]
        return np.array([resx,resy]).T
    ## extra args preparation ##
    def prepare_extra_args(self):
        extra = (self.case,self.sin,self.sout,self.selfs,self.M,self.nvarx,self.nvary,
            self.inds_selfs,generate_inds_full_ij(self.nvarx,self.nvary))
        return extra
  
    ## checks ##
    def check_s(self,xx=None):
        """ Checks strength distance """
        if xx is None:
            x,y = self.x,self.y
        else:
            xx = self.var_transform(xx)
            x = xx[:self.nvarx]
            y = xx[self.nvarx:]
        return check_strengths(x,y,self.sout,self.sin,self.M, self.selfs, self.inds_selfs, case = self.case)
        
    def check_domain(self,xx=None):
        """ Check constraint, if all_good returns True """
        if xx is None:
            x,y = self.x,self.y
        else:
            xx = self.var_transform(xx)
            x = xx[:self.nvarx]
            y = xx[self.nvarx:]
        return check_domain(x,y,self.selfs, self.inds_selfs, case=self.case)
        
    def check_L(self):
        return loglikelyhood(self.var_adapt(),*self.prepare_extra_args())
    
    def check_L_grad(self):
        return loglikelyhood_grad(self.var_adapt(),*self.prepare_extra_args())
    
    def check_L_hess(self):
        return loglikelyhood_hess(self.var_adapt(),*self.prepare_extra_args())

    ## preconditioning ##
    def init_conds(self,rand=False):
        """
            Sets init_conds if not given
        """
        case = self.case
        if case == 'ME' or case =='B':
            x = self.sout/np.sqrt(self.sout.sum())
            y = self.sin/np.sqrt(self.sout.sum())
        else: # weighted
            x = np.ones(self.nvarx)*0.9
            y = np.ones(self.nvary)*0.9
        if rand:
            #x = np.random.random(self.nvarx)*0.9 + 1e-8
            #y = np.random.random(self.nvary)*0.9 + 1e-8
            x = x*(1+0.001*np.abs(np.random.random(self.nvarx)))
            y = y*(1+0.001*np.abs(np.random.random(self.nvary)))
            if case == 'W':
                x/=x.max()*0.95
                y/=y.max()*0.95
        return x,y
    def precondition(self,order=True,solve_system=False,**kwargs):
        """
            Calls TNC method to precondition x,y variables + linear solver (optionally)
            or the other way around if ord set to False.
        """
        print "\t# - # Preconditioning # - #"
        if self.check_s()[0]<0.1:
            return
        tol_ini             = 1e-5
        kwargs['tol']       = tol_ini
        if solve_system:
            if order:
                q = self.solve_as_system(tol= tol_ini,**kwargs) # you may play with existing scipy methods
                q = self.minimize(**kwargs)
            else:
                q = self.minimize(**kwargs)
                q = self.solve_as_system(tol= tol_ini,**kwargs) # you may play with existing scipy methods
        else:
            q = self.minimize(**kwargs)
        if not q['success']:
             print "\t\t# - # Preconditniong failed... # - #"
    ########################################################
    #############  Solvers    ############################## 
    ########################################################
    def solve_as_system(self,meth='hybr', force_update=True, **kwargs): # you may play with existing scipy methods
        # see http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html#scipy.optimize.root
        """
            Cases ME and B: Uses custom solver
            Case W:
                Tries to solve using scipy methods for large problems
                You can set the extra options of the solvers by passing arguments to kwargs...
            force_update forces to update result even if no convergence is attained (just in case it is better than the actual result)
        """
        case = self.case
        print "\t# ================= System solver==================== #"
        if case == 'ME' or case =='B':
            if not self.selfs:
                kwargs['inds_selfs'] = self.inds_selfs
            x,y = balance_xy(self.sin,self.sout,x_ini=self.x,y_ini=self.y,case=self.case,M=self.M,selfs=self.selfs,**kwargs)
            xx = np.r_[x,y]
            res = {}
            res['x'] = xx
            if np.sum(np.isfinite(xx)) == len(xx):
                res['success'] = True
            else: res['success'] = False
        else:
            meth = kwargs.get('meth',meth)
            try:
                solver = opt.root
            except:
                print "\t\tMethod not available"
                return
            extra_args = self.prepare_extra_args()
            res = solver(loglikelyhood_grad,self.var_adapt(),method=meth,jac=loglikelyhood_hess_array,
                args=extra_args,options=kwargs)
        res['violated'] = not self.check_domain(res['x'])        
        if not res['violated']: # if constraints are not violated
            if res['success'] or force_update:
                self.var_update(res['x'])
        if res['success']:
            if not res['violated']:
                self.var_update(res['x'])
                print "\t\tFinal value of function: %f" % (self.check_L())
                print "\t\tFinal value of delta constraints: %f" % (self.check_s()[0])
            else:
                print "\t\tSystem solved, but the constraints are violated"
        else:
            #if kwargs.get('verbose',False):
            print "\t\tMinimization failed!"
        return res
        

    def minimize(self,meth='TNC', force_update=True,**kwargs): # meths: TNC prefered
        """
            Unconstrained minimization problem on a squared N-dimensional domain.
            You can set the extra options of the solvers by passing arguments to kwargs...
            force_update forces to update result even if no convergence is attained (just in case it is better than the actual result)
        """
        case = self.case
        epsi = kwargs.get('epsi',1e-12)
        print "\t# ================= Likelyhood Maximizer ============= #"
        meth = kwargs.get('meth',meth)
        ### set options ###
        # prepare kwargs #
        kwargs['tol'] = kwargs.get('tol',1e-5) # default tol
        kwargs['ftol']=0. # its a loglikelyhood, it can give 0 at max.
        opts = kwargs
        opts['disp']=kwargs.get('verbose',False)
        # check if method is allowed #
        allowed_methods = ['L-BFGS-B','TNC','SLSQP']
        if meth not in allowed_methods:
            print "\t\tRecommended use: TNC"
        minimizer=opt.minimize
        extra_args = self.prepare_extra_args()
        ## set bounds ##
        bmin = kwargs.get('bmin',epsi)
        bounds = np.ones((self.nvar,2))
        bounds.T[0]=bmin
        if case == 'ME' or case == 'B':
            bounds.T[-1] = 2**100
        else:
            bmax = kwargs.get('bmax',1.)
            if bmax == 1:
                bmax = bmax-epsi
                bmax_x = bmax
                bmax_y = bmax
            elif bmax>1:                
                if self.sout.max()>self.sin.max():
                    bmax_x = bmax - epsi
                    if self.selfs:
                        bmax_y = 1./bmax
                    else:
                        bmax_y = 1./(bmax-0.1*bmax)
                else:
                    bmax_y = bmax - epsi
                    if self.selfs:
                        bmax_x = 1./bmax
                    else:
                        bmax_x = 1./(bmax-0.1*bmax)                 
            else:
                bmax_x=bmax
                bmax_y=bmax
            self.x[self.x>bmax_x] = bmax_x
            self.y[self.y>bmax_y] = bmax_y
            self.y[self.y<=0]=0
            self.x[self.x<=0]=0            
            bounds.T[1][:self.nvarx] = bmax_x
            bounds.T[1][self.nvarx:] = bmax_y            
        if self.check_s()[0]<kwargs['tol']:
            res = {}
            res['x'] = self.var_adapt()
            res['success'] = True
        else:
            res = minimizer(loglikelyhood,self.var_adapt(),args=extra_args,method=meth,jac=loglikelyhood_grad,
                hess = loglikelyhood_hess, options = opts, bounds = bounds)
        res['violated'] = not self.check_domain(res['x'])
        if not res['violated']: # if constraints are not violated
            if res['success'] or force_update:
                self.var_update(res['x'])
        if res['success']:
            if not res['violated']:
                print "\t\tFinal value of function: %f" % (self.check_L())
            else:
                print "\t\tSystem solved, but the constraints are violated"
        else:
            if kwargs.get('verbose',False):
                print "\t\tMinimization failed!"        
        return res

#######  Log likelyhood functions ###
import itertools
from scipy import sparse

def loglikelyhood(xxf,case,sin,sout,selfs,M,nvarx,nvary,inds_selfs,*args,**kwargs):
    """ Computes minus loglikelyhood """
    # vars #
    x = xxf[:nvarx]
    y = xxf[nvarx:]
    # logs #
    xl = np.log(x)
    yl = np.log(y)
    if case == 'W':
        aux = 1.-np.einsum('i,j',x,y)
    elif case == 'B':
        aux = 1./(1.+np.einsum('i,j',x,y))
    else:
        aux = np.einsum('i,j',x,y)
    if not selfs:
        aa = aux.flatten()
        if case == 'ME':
            aa[inds_selfs] = 0.
        else:
            aa[inds_selfs] = 1.
        aux = aa.reshape(nvarx,nvary)
    L1 = np.einsum('i,i',sin,yl)
    L2 = np.einsum('i,i',sout,xl)
    if case == 'ME':
        L4 = np.einsum('ij->',-aux)
    else:
        L4 = np.einsum('ij->',np.log(aux))
    L = -(M*L4+L1+L2)
    return L
def loglikelyhood_grad(xxf,case,sin,sout,selfs,M,nvarx,nvary,inds_selfs,*args,**kwargs):
    """ Computes gradient vector """
    # sin,sout,selfs,M
    # nvarx, nvary, inds_selfs
    x = xxf[:nvarx]
    y = xxf[nvarx:]
    ## common terms ##
    if case == 'W':
        aux = 1.-np.einsum('i,j',x,y)
    elif case == 'B':
        aux = 1.+np.einsum('i,j',x,y)
    else:
        aux = np.ones((nvarx,nvary))
    aux = 1./aux
    if not selfs:
        aa = aux.flatten()
        aa[inds_selfs[0]] = 0
        aux = aa.reshape(nvarx,nvary)
    g1 =  sout/(x+epsilon) - M*np.einsum('j,ij',y,aux)
    g2 =  sin/ (y+epsilon) - M*np.einsum('i,ij',x,aux)
    grad = -np.r_[g1,g2]
    return grad
    
def loglikelyhood_hess(xxf,case,sin,sout,selfs,M,nvarx,nvary,inds_selfs,inds_full_ij,*args,**kwargs):
    """ Computes Hessian matrix in sparse ccoo_scipy format """
    x = xxf[:nvarx]
    y = xxf[nvarx:]
    ## common terms ##
    if case == 'W':
        aux = 1.-np.einsum('i,j',x,y)
    elif case == 'B':
        aux = 1.+np.einsum('i,j',x,y)
    else:
        aux = np.ones((nvarx,nvary))
    aux = 1./(aux*aux)
    if not selfs:
        aa = aux.flatten()
        aa[inds_selfs] = 0
        aux = aa.reshape(nvarx,nvary)
    if case == 'B':
        ## diagonal ##
        Hxx = sout/(x*x+epsilon) - M*np.einsum('j,ij',y*y,aux)
        Hyy = sin /(y*y+epsilon) - M*np.einsum('i,ij',x*x,aux)        
        ## off diagonal, lower and upper corner ##
    elif case =='W':
        ## diagonal ##
        Hxx = sout/(x*x+epsilon) + M*np.einsum('j,ij',y*y,aux)
        Hyy = sin /(y*y+epsilon) + M*np.einsum('i,ij',x*x,aux)
    else:
        ## diagonal ##
        Hxx = sout/(x*x+epsilon)
        Hyy = sin /(y*y+epsilon)
    ## off diagonal, lower and upper corner ##
    Hxy = M*aux
    # to sparse matrix #
    # rows and columns of diagonal
    rows_d = np.arange(nvarx+nvary)
    cols_d = rows_d 
    # rows and columns of the upper and lower block
    rows_ndl = inds_full_ij[1]+nvarx #lower (invert indices)
    cols_ndl = inds_full_ij[0] 
    rows_ndu = inds_full_ij[0] 
    cols_ndu = inds_full_ij[1]+nvarx #upper 
    ## Stack together! ##
    vals = np.r_[Hxx,Hyy,Hxy.flatten(),Hxy.flatten()] # diag_x, diag_y, upper, lower
    inds = (np.r_[rows_d,rows_ndl,rows_ndu],np.r_[cols_d,cols_ndl,cols_ndu])
    hess = sparse.coo_matrix((vals, inds), shape=(nvarx+nvary, nvarx+nvary))
    return hess

def loglikelyhood_hess_array(xxf,sin,sout,selfs,M,nvarx,nvary,inds_selfs,inds_full_ij,*args,**kwargs):
    return loglikelyhood_hess(xxf,sin,sout,selfs,M,nvarx,nvary,inds_selfs,inds_full_ij,*args,**kwargs).toarray()
    # very crappy patch to return an array, should be changed for a decorator #


## Domain checker & cons checker ##    
def check_strengths(x,y,sout,sin,M,selfs,inds_selfs,case='W'):
    """
        Returns norm between constraints and actual values
    """
    nvarx = len(sout)
    nvary = len(sin)
    aux = np.einsum('i,j',x,y)
    if not selfs:
        aa = aux.flatten()
        aa[inds_selfs] = 0
        aux = aa.reshape(nvarx,nvary)
    if case == 'B':
        aux = M*aux/(1.+aux)
    elif case == 'W':
        aux = M*aux/(1.-aux)        
    else:
        pass
    g1 = aux.sum(axis=1) - sout
    g2 = aux.sum(axis=0) - sin
    return np.linalg.norm(g1),np.linalg.norm(g2)


def check_domain(x,y,selfs,inds_selfs,case = 'W'):
    """
        Checks wheter {x,.y} is in the domain, returns true or false according to case
    """
    nvarx = len(x)
    nvary = len(y)
    aux = np.einsum('i,j',x,y)
    if not selfs:
        aa = aux.flatten()
        aa[inds_selfs] = 0
        aux = aa.reshape(nvarx,nvary)
    if case == 'B' or case == 'ME':
        return (aux.min()>=0)
    else:
        return (aux.max()<1-epsilon and aux.min()>=0)

### handy stuff ###
def generate_inds_selfs(nvarx,nvary):
    """
        Generate indices where self-loop elements must be zero in  a matrix (nvarx,nvary)
        On a flattened array ordered as x_1 y_1, x_1 y_2..x_1y_nvary... x2_y_1....x_nvarx y_nvarj
        The elements where i == j.
    """
    iself = np.fromiter(((nvary+1)*i for i in xrange(nvarx)),int)
    return (iself[iself<nvarx*nvary],) # numpy index format

def generate_inds_full_ij(nvarx,nvary):
    """
        Generates all indices of a matrix (nvarx x nvary)
    """
    rows = np.array([np.ones(nvary,dtype=int)*i for i in range(nvarx)]).flatten()
    cols = np.array([np.arange(nvary)]*nvarx).flatten()
    return (rows,cols)
    
