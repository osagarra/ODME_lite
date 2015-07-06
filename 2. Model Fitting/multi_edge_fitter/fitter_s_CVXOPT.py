"""
    CVXOPT implementation
"""
import multi_edge_fitter as mw
from scipy import sparse
import numpy as np
import cvxopt
##  from : http://cvxopt.org/userguide/solvers.html#s-cp

#epsilon = 1e-15
epsilon = 0

class Edge_Maximizer(mw.fitter_s.Edge_Maximizer):
    """
        General maximizer class for fixed sin,sout
    """
    def to_CVXOPT_dense(self,x,from_arr = True):
        if not from_arr:
            return cvxopt.matrix(x,x.shape)
        else:
            return cvxopt.matrix(x,(x.shape[0],1))
    #### Re-definition of some methods ###       
    def prepare_extra_args(self,to_CVXOPT = False):
        extra = super(Edge_Maximizer, self).prepare_extra_args()
        if not to_CVXOPT:
            return extra
        else:
            sin_dense  = self.to_CVXOPT_dense(self.sin)
            sout_dense = self.to_CVXOPT_dense(self.sout)
            extra = list(extra)
            extra[1] = sin_dense
            extra[2] = sout_dense
            extra.append(generate_inds_jac_cons(self.nvarx,self.nvary))
            return tuple(extra)
    def init_conds(self,rand=False):
        """
            Sets init_conds if not given
        """
        #return super(Edge_Maximizer, self).init_conds(rand=True)
        return super(Edge_Maximizer, self).init_conds(rand=rand)


    ##### solver #######
    def solve_using_cvxopt(self, cons = False, bound = False, force_update=True, **kwargs):
        """
        Solve the convex optimization problem
            min f_0(xx) = min Loglikelyhood(xx)
            subject to:
                Domain constraint x_max y_max <= 1
        Implementing only assuming x need to be inside the domain
        Comments from: http://cvxopt.org/userguide/solvers.html#s-cp
        
        3 alternative implementations:
            1. Unconstrained, unbounded problem (only requires that x,y >0) [default]
            2. Unconstrained, bounded problem (x,y < 1) [bound=True]
            3. Constrained, unbounded problem (x,y>0, xy<1 enforced as constraints) [cons=True]
        
        force_update forces to update result even if no convergence is attained (just in case it is better than the actual result)
        """
        print "\t# ================= CVXOPT solver==================== #"
        # number of variables #
        nvar = self.nvar
        if self.case == 'ME' or self.case =='B':
            cons  = False
            bound = False
            to_CVXOPT = False
        else:
            if cons:
                ncons = self.nvarx*self.nvary # constraints
                to_CVXOPT = True
                bound = False
                #ncons = 1 # constraints
            else:
                to_CVXOPT = False
                ncons = 0
        m,n = ncons,nvar
        extra_args = self.prepare_extra_args(to_CVXOPT=to_CVXOPT)
        def F(xx=None, z=None,extra_args=extra_args, to_CVXOPT = to_CVXOPT, cons=cons, bound = bound):
            """
                to_CVXOPT choses to use native LAPACK stuff
                cons choses to manually enforce the nonlinear constraint x_iy_j<1
                bound choses to keep all variables x,y in [0,1]
            """
            """ 
                if x is none, return returns a tuple (m, x0), where m is the number of nonlinear constraints 
                and x_0 is a point in the domain of f. x0 is a dense real matrix of size (n, 1).
            """
            if xx is None:
                if self.check_domain(): 
                    x0 = self.to_CVXOPT_dense(self.var_adapt())
                else:
                    #print "x Not in domain..."
                    x0 = cvxopt.matrix(np.ones(nvar)*0.9,(nvar,1),'d')
                return ncons, x0
            """
                If x is not in the domain of f, F(x) returns None or a tuple (None, None).
            """
            if not self.check_domain(np.array(xx).flatten()): return None
            """ 
                F(x), with x a dense real matrix of size (n, 1), returns a tuple (f, Df). 
                f is a dense real matrix of size (m+1, 1), with f[k] equal to f_k(x). 
                (If m is zero, f can also be returned as a number.) 
            """
            if to_CVXOPT:
                if cons:
                    f =  loglikelyhood_and_cons_cvxopt(xx,extra_args)
                else:
                    f =  loglikelyhood_cvxopt(xx,extra_args,to_CVXOPT=True)
            else:
                f =  loglikelyhood_cvxopt(xx,extra_args,to_CVXOPT=False)
            """
                Df is a dense or sparse real matrix of size (m + 1, n) with Df[k,:] equal to the transpose of the gradient \nabla f_k(x).
            """
            if to_CVXOPT:
                if cons:
                    Df = loglikelyhood_and_cons_grad_cvxopt(xx,extra_args).T
                    #print Df.size
                else:
                    Df = loglikelyhood_grad_cvxopt(xx,extra_args,to_CVXOPT=True).T             
            else:
                Df = loglikelyhood_grad_cvxopt(xx,extra_args,to_CVXOPT=False).T             
            if z is None: return f, Df
            """
                F(x,z), with x a dense real matrix of size (n, 1) and z a positive dense real matrix of size (m + 1, 1) returns a tuple (f, Df, H). 
                f and Df are defined as above. H is a square dense or sparse real matrix of size (n, n), 
                whose lower triangular part contains the lower triangular part of
                    z_0 \nabla^2f_0(x) + z_1 \nabla^2f_1(x) + \cdots + z_m \nabla^2f_m(x).

                If F is called with two arguments, it can be assumed that x is in the domain of f.
            """
            if to_CVXOPT:
                if cons:
                    H = loglikelyhood_and_cons_hess_cvxopt(xx,z,extra_args)
                    ###
                    #from numpy import linalg
                    #h_cx_array = np.zeros(H.size)
                    #h_cx_array[(np.array(H.I).flatten(),np.array(H.J).flatten())] = np.array(H.V).flatten()
                    #g_cx_array = np.zeros(Df.size)
                    #g_cx_array[(np.array(Df.I).flatten(),np.array(Df.J).flatten())] = np.array(Df.V).flatten()
                    #print linalg.matrix_rank(g_cx_array),linalg.matrix_rank(h_cx_array)
                    ###                 
                else:
                    H = loglikelyhood_hess_cvxopt(xx,z,extra_args,to_CVXOPT=True)
            else:
                H = loglikelyhood_hess_cvxopt(xx,z,extra_args,to_CVXOPT=False)
            return f, Df, H
        # Cone programming #
        onnes = np.ones(nvar)
        inds  = np.arange(nvar)
        G     = cvxopt.spmatrix(onnes,inds,inds)
        h     = cvxopt.matrix(onnes,(nvar,1),'d')
        #dims = {'l':,'q':0,'s':0}
        # Set opts #
        key_kwargs = ('show_progress','maxiters','abstol','reltol','feastol','refinement','maxiters')
        defaults = {'show_progress':False,'maxiters':100,'abstol':1e-7,'reltol':1e-6,'feastol':1e-7,'refinement':1}
        for k,v in defaults.iteritems(): # set defaults (if there are multiple calls)
            cvxopt.solvers.options[k] = v
        for k in kwargs.iterkeys():
            opt = kwargs.get(k,None)
            if opt:
                if k in key_kwargs:
                    cvxopt.solvers.options[k]=opt
        # Solve #
        if bound:
            sol = cvxopt.solvers.cp(F,G,h)          
        else:
            sol = cvxopt.solvers.cp(F)
        # Wrap solution #
        sol['violated'] = not self.check_domain()
        if not sol['violated']:
            if force_update or sol['status']=='optimal':
                self.var_update(np.array(sol['x']).flatten())
        if sol['status']=='optimal' or  self.check_s(sol['x'])<self.check_s(): #if it is better
            sol['success'] = True
        else:
            sol['success'] = False
        return sol


#### handies ####
import itertools

def generate_inds_jac_cons(nvarx,nvary):
    """
        Generates indices for the jacobian sparse structure of the constraints
        cols: number of constraints
        rows: variables
    """  
    cols = np.r_[np.arange(nvarx*nvary),np.arange(nvarx*nvary)]
    rows_y = np.array([i*np.ones(nvary,dtype=int) for i in xrange(nvarx)]).flatten()
    rows_x = np.array([np.arange(nvary,dtype=int)]*nvarx).flatten() + nvarx
    return (np.r_[rows_y,rows_x],cols) # numpy index format




#######################################################
#######  Log likelyhood functions adapted to CVXOPT ###
#######################################################

def loglikelyhood_cvxopt(xxf,extra_args,to_CVXOPT=False):
    if not to_CVXOPT:    
        return mw.fitter_s.loglikelyhood(np.array(xxf).flatten(),*extra_args)
    case,sin,sout,selfs,M,nvarx,nvary,inds_selfs = extra_args[:8]
    x = xxf[:nvarx]
    y = xxf[nvarx:]
    xl = cvxopt.log(x)
    yl = cvxopt.log(y)
    if case == 'W':
        aux = 1.-np.einsum('ik,jk',x,y)
    elif case == 'B':
        aux = 1./(1.+np.einsum('ik,jk',x,y))
    else:
        aux = np.einsum('ik,jk',x,y)
    if not selfs:
        aa = aux.flatten()
        if case == 'ME':
            aa[inds_selfs] = 0.
        else:
            aa[inds_selfs] = 1.
        aux = aa.reshape(nvarx,nvary)
    L1 = np.einsum('ik,ik',sin,yl)
    L2 = np.einsum('ik,ik',sout,xl)
    if case == 'ME':
        L4 = np.einsum('ij->',-aux)
    else:
        L4 = np.einsum('ij->',np.log(aux))
    L = -(M*L4+L1+L2)
    return -(M*L4+L1+L2)
    
def loglikelyhood_grad_cvxopt(xxf,extra_args,to_CVXOPT=False):
    # implement BLAS and LAPACK stuff!
    """ Computes gradient vector """
    if not to_CVXOPT:    
        grad= mw.fitter_s.loglikelyhood_grad(np.array(xxf).flatten(),*extra_args)
    else:
        case,sin,sout,selfs,M,nvarx,nvary,inds_selfs = extra_args[:8]
        x = xxf[:nvarx]
        y = xxf[nvarx:]
        if case == 'W':
            aux = 1.-np.einsum('ik,jk',x,y) # lapack! (to do)
        elif case == 'B':
            aux = 1.+np.einsum('ik,jk',x,y) # lapack! (to do)        
        else:
            aux = np.ones((nvarx,nvary,1))
        #assert(np.all(aux>epsilon))
        aux =  1./aux
        if not selfs:
            aa = aux.flatten()
            aa[inds_selfs] = 0
            aux = aa.reshape(nvarx,nvary)
        g1 =  cvxopt.div(sout,(x+epsilon)) - M*np.einsum('jk,ij',y,aux)
        g2 =  cvxopt.div(sin ,(y+epsilon)) - M*np.einsum('ik,ij',x,aux)
        grad = -np.r_[g1,g2]
    return cvxopt.matrix(grad,(grad.shape[0],1),'d')

    
def loglikelyhood_hess_cvxopt(xxf,z,extra_args,to_CVXOPT=False):
    """ Computes Hessian matrix, using sparse matrix format CVXOPT """
    # implement BLAS and LAPACK stuff!
    case,sin,sout,selfs,M,nvarx,nvary,inds_selfs,inds_full_ij = extra_args[:9]
    if not to_CVXOPT:    
        hess= mw.fitter_s.loglikelyhood_hess(np.array(xxf).flatten(),*extra_args)
        inds = np.tril_indices(nvarx+nvary)
        hess = hess.toarray()[inds] #2Nx2N, keep lower diagonal
    else:
        x = xxf[:nvarx]
        y = xxf[nvarx:]
        x2 = cvxopt.mul(x,x)
        y2 = cvxopt.mul(y,y)
        ## common terms ##
        if case == 'W':
            aux = 1.-np.einsum('ik,jk',x,y)
        elif case == 'B':
            aux = 1.+np.einsum('ik,jk',x,y)
        else:
            aux = np.ones((nvarx,nvary))
        #assert(np.all(aux>epsilon))
        aux = 1./(aux*aux)
        if not selfs:
            aa = aux.flatten()
            aa[inds_selfs] = 0
            aux = aa.reshape(nvarx,nvary)
        ## diagonal ##
        if case == 'B':
            ## diagonal ##
            Hxx = cvxopt.div(sout,(x2+epsilon)) -  M*np.einsum('jk,ij',y2,aux)
            Hyy = cvxopt.div(sin ,(y2+epsilon)) -  M*np.einsum('ik,ij',x2,aux)
            ## off diagonal, lower and upper corner ##
        elif case =='W':
            ## diagonal ##
            Hxx = cvxopt.div(sout,(x2+epsilon)) +  M*np.einsum('jk,ij',y2,aux)
            Hyy = cvxopt.div(sin ,(y2+epsilon)) +  M*np.einsum('ik,ij',x2,aux)
        else:
            ## diagonal ##
            Hxx = cvxopt.div(sout,(x2+epsilon))
            Hyy = cvxopt.div(sin ,(y2+epsilon))
        ## off diagonal, lower and upper corner ##
        Hxy = M*aux
        ## Return as CVXOPT matrix ## -> http://cvxopt.org/userguide/matrices.html#sparse-matrices
        # triplet description (value,col,row)
        # rows and columns of diagonal
        rows_d = np.arange(nvarx+nvary)
        cols_d = rows_d 
        # rows and columns of the lower block
        rows_ndl = inds_full_ij[1]+nvarx #lower (invert indices)
        cols_ndl = inds_full_ij[0] 
        inds = (np.r_[rows_d,rows_ndl],np.r_[cols_d,cols_ndl])
        hess = np.r_[Hxx.reshape((nvarx,)),Hyy.reshape((nvary,)),Hxy.reshape((nvarx*nvary,))]
    return cvxopt.spmatrix(z[0]*hess,inds[0],inds[1])


### Adapted functions to constraints ###

def loglikelyhood_and_cons_cvxopt(xxf,extra_args):
    f1 = loglikelyhood_cvxopt(xxf,extra_args,to_CVXOPT=True)
    extra_args2 = tuple(extra_args[i] for i in [3,5,6,7,8,9])
    xx = np.array(xxf).flatten()    
    fk = constraints(xx,*extra_args2)
    #fk = constraints_max(xx,*extra_args2)
    nvarx = extra_args[1]
    nvary = extra_args[2]
    return cvxopt.matrix(np.r_[f1,fk])
    
def loglikelyhood_and_cons_grad_cvxopt(xxf,extra_args):
    g1 = loglikelyhood_grad_cvxopt(xxf,extra_args,to_CVXOPT=True)
    extra_args2 = tuple(extra_args[i] for i in [3,5,6,7,8,9])
    xx = np.array(xxf).flatten()    
    gk = constraints_jac(xx,*extra_args2)
    #gk = constraints_jac_max(xx,*extra_args2)
    vals = np.r_[np.array(g1).flatten(),gk.data]
    rows = np.r_[np.arange(xxf.size[0],dtype=int),gk.row] 
    cols = np.r_[np.zeros(xxf.size[0],dtype=int),1+gk.col,]
    ########
    # very dirty hack to avoid CVXOPT problem!!!! <-- ignore
    ########
    #vals = np.r_[vals,0.9] 
    #rows = np.r_[rows,rows[np.random.randint(1,high=len(rows))]] 
    #cols = np.r_[cols,cols[[np.random.randint(1,high=len(cols))]]]
    ###
    return cvxopt.spmatrix(vals,rows,cols) # needs to be transposed
    
def loglikelyhood_and_cons_hess_cvxopt(xxf,z,extra_args):
    h1 = loglikelyhood_hess_cvxopt(xxf,z,extra_args,to_CVXOPT=True)
    extra_args2 = tuple(extra_args[i] for i in [3,5,6,7,8,9])
    xx = np.array(xxf).flatten()    
    hk = constraints_hess(xx,z,*extra_args2)
    #hk = constraints_hess_max(xx,z,*extra_args2)
    hk = cvxopt.spmatrix(hk.data,hk.row,hk.col,h1.size)
    return h1+hk


### constraints ####

def constraints(xx,selfs,nvarx,nvary,inds_selfs,*args):
    """
        Returns constraints vector Nvarx Nvary in length, using numpy arrays as input
    """
    x = xx[:nvarx]
    y = xx[nvarx:]
    ## common terms ##
    aux =  np.einsum('i,j',x,y).flatten()
    if not selfs:
        aux[inds_selfs[0]] = 0
    return aux-1+np.sqrt(epsilon) # we want the constraints smaller or equal 0


def constraints_max(xx,selfs,nvarx,nvary,inds_selfs,*args):
    """
        Returns constraints vector Nvarx Nvary in length, using numpy arrays as input
    """
    x = xx[:nvarx]
    y = xx[nvarx:]
    if selfs:
        v = x.max()*y.max()
    else:
        xm = np.sort(x)[::-1][:2] # 2 largest values
        ym = np.sort(y)[::-1][:2] # 2 largest values
        v = np.max(xm*ym[::-1])
    return v - 1 + np.sqrt(epsilon)
    
def constraints_jac(xx,selfs,nvarx,nvary,inds_selfs,inds_full_ij,inds_jac_cons,*args): # inds_full_ij added so the same arguments are called on hessian and jac
    """
        Returns jacobian of constraints
    """
    x = xx[:nvarx]
    y = xx[nvarx:]
    # block of horizontal y's on top side # block of diagonal x's on the bottom side
    vals = np.r_[np.array([y]*nvarx).flatten(),np.array([np.ones(nvary)*x[i] for i in xrange(nvarx)]).flatten()]
    # take care of selfs#
    if not selfs:
        vals[inds_selfs[0]] = 0
        vals[inds_selfs[0]+nvarx*nvary] = 0
    # its a very sparse matrix #
    jac_c = sparse.coo_matrix((vals, inds_jac_cons), shape=(nvarx+nvary,nvarx*nvary))
    return jac_c

def constraints_jac_max(xx,selfs,nvarx,nvary,inds_selfs,inds_full_ij,inds_jac_cons,*args): # inds_full_ij added so the same arguments are called on hessian and jac
    """
        Returns jacobian of constraints
    """
    x = xx[:nvarx]
    y = xx[nvarx:]
    vals = np.zeros(2)
    if selfs:
        xm = x.max()
        ym = y.max()
    else:
        xm = np.sort(x)[::-1][:2] # 2 largest values
        ym = np.sort(y)[::-1][:2] # 2 largest values
        v = xm*ym[::-1]
        if v[0]>v[1]: # means combination is x_m and second largest y
            xm = xm[0]
            ym = ym[-1]
        else: # means comb is y_max and second largest x
            xm = xm[1]
            ym = ym[0]
    xm_ind = np.where(x==xm)
    ym_ind = np.where(y==ym)           
    vals[0] = ym
    vals[1] = xm
    rows = np.r_[xm_ind[0],ym_ind[0]+nvarx]
    # its a very sparse matrix #
    jac_c = sparse.coo_matrix((vals, (rows,np.zeros(2))), shape=(nvarx+nvary,1))
    return jac_c


def constraints_hess(xx,z,selfs,nvarx,nvary,inds_selfs,inds_full_ij,inds_jac_cons,*args): # z are the lagrange multipliers
    """
        Returns hessian constraint function with lagrange multipliers z
    """
    # return sparse matrix #
    rows = inds_full_ij[1]+nvarx #we transpose the full_ij indices, because we want the lower triangle
    cols = inds_full_ij[0] 
    inds = (rows,cols)
    # multipliers #
    #print "%r" % z
    vals =   -np.array(z).flatten()[1:] # z[0] is the main function (we can impose that it is negative)
    if not selfs:
        vals[inds_selfs[0]] = 0
    hess_c = sparse.coo_matrix((vals, inds), shape=((nvarx+nvary), (nvarx+nvary)))
    return hess_c


def constraints_hess_max(xx,z,selfs,nvarx,nvary,inds_selfs,inds_full_ij,inds_jac_cons,*args): # z are the lagrange multipliers
    """
        Returns hessian constraint function with lagrange multipliers z
    """
    x = xx[:nvarx]
    y = xx[nvarx:]
    vals = np.zeros(2)
    if selfs:
        xm = x.max()
        ym = y.max()
    else:
        xm = np.sort(x)[::-1][:2] # 2 largest values
        ym = np.sort(y)[::-1][:2] # 2 largest values
        v = xm*ym[::-1]
        if v[0]>v[1]: # means combination is x_m and second largest y
            xm = xm[0]
            ym = ym[-1]
        else: # means comb is y_max and second largest x
            xm = xm[1]
            ym = ym[0]
    xm_ind = np.where(x==xm)
    ym_ind = np.where(y==ym)           
    vals[0] = z[1]
    vals[1] = z[1]
    rows = np.r_[xm_ind[0],ym_ind[0]+nvarx]
    cols = np.r_[ym_ind[0]+nvarx,xm_ind[0]]
    inds = (rows,cols)
    hess_c = sparse.coo_matrix((vals, inds), shape=((nvarx+nvary), (nvarx+nvary)))
    return hess_c


