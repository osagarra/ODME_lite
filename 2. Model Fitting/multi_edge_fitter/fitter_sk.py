"""
 Fitter for fixed strength sequence and fixed degree sequence
"""

## Imports ##
import numpy as np



"""
#############################
#### Sequential balancing: all at once (much better, preferred approach) #########
#############################
"""

def aux_selector(agg=False):
    """
    Selects appropriate aux function for the balancer
    """
#    if agg:
#        aux_f = lambda a,b,alf,M : (1+np.einsum('i,j',a/alf,b/alf))**(M)
#        corr = lambda a,b,alf : np.einsum('i,j',a/alf,b/alf)+1.
#    else:
#        aux_f = lambda a,b,alf,M : np.exp(np.einsum('i,j',a/alf,b/alf)) # ignores M
#        corr = lambda a,b,alf : 1.
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
#    if selfs:
#        extra_k = lambda c,d,aux,aux2 : 0
#        extra_s = lambda a,b,c,d,aux,aux2,aux3 : 0
#    else:
#        extra_k = lambda c,d,aux,aux2: np.diag(c*d*(aux-1.)/(aux2*(aux-1.)+1.))
#        extra_s = lambda a,b,c,d,aux,aux2,aux3: np.diag(c*d*a*b*aux/(aux2*(aux-1.)+1.)/aux3)
    if selfs:
        extra_k = lambda c,d,aux,aux2 : 0
        extra_s = lambda a,b,c,d,aux,aux2,aux3 : 0
    else:
        extra_k = lambda c,d,aux,aux2: np.diag(c*d*(1.-aux)/(aux2*(1.-aux)+aux))
        extra_s = lambda a,b,c,d,aux,aux2,aux3: np.diag(c*d*a*b/(aux2*(1.-aux)+aux)/aux3)
    return extra_k,extra_s



def balance_xyzw(sin,sout,kin,kout,tol=1e-9,tol_c=10,maxreps=10000,verbose=False,selfs=True,print_tol=False,print_c=False,agg=False,M=1,**kwargs):
    """ Balances the strength and degree  equations
        Input:
            Sin: Incoming strength sequence (N)
            Sout: Outgoing strength sequence (N)
            Kin: Incoming degree sequence (N)
            Kout: Outgoing degree sequence (N)
            tol: Maximum iteration tolerance [float]
            tol_c: Maximum tolerance permitted for the worse constraint-equation [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            print_tol:  Print convergence of algorithm
            Agg: If the weighted network comes from aggreagation of binary networks set to True
            M: Number of layers (only in the case of aggregated network)
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
                z_ini: Initial guess for z
                w_ini: Initial guess for w
                act: Print error every act steps (if print_c is set to True)
                alf: Reduction factor on variables (by default max(xy))
    """
    ## add precision ##
    sout = np.array(sout,dtype=np.float128)
    kout = np.array(kout,dtype=np.float128)
    sin = np.array(sin,dtype=np.float128)
    kin = np.array(kin,dtype=np.float128)
    ## print opts ##
    if verbose or print_tol:
        print "#### Degree & Strength balancer #####"
    n = len(sin)
    if n>1000:
        act = 10
    elif n<1000 and n>500:
        act = 50
    else:
        act = 500
    ## options ##
    if not agg:
        M=1
    act = kwargs.get('act',act)
    T = sout.sum()
    alf0 = np.sqrt(M*sout.max()*sin.max()*1./T)
    alf = kwargs.get('alf',alf0)
    #alf = kwargs.get('alf',np.sqrt(sout.max()*sin.max()*1./T))
   
    ## arranging variables & pre-conditioning ##
    a=kwargs.get('x_ini',sout/np.sqrt(T))
    b=kwargs.get('y_ini',sin/np.sqrt(T))
    c=kwargs.get('z_ini',np.ones(n)*0.9)
    d=kwargs.get('w_ini',np.ones(n)*0.9)
    #a=kwargs.get('x_ini',0.09*np.ones(n))
    #b=kwargs.get('y_ini',0.09*np.ones(n))
    #c=kwargs.get('z_ini',np.ones(n)*0.9)
    #d=kwargs.get('w_ini',np.ones(n)*0.09)
    if 'z_ini' not in kwargs.keys() or 'w_ini' not in kwargs.keys():
        if verbose: print "Pre conditioning..."
        c,d = balance_zw(a/alf,b/alf,kin,kout,tol=tol,tol_c=tol_c,maxreps=maxreps,verbose=verbose,selfs=selfs,print_tol=print_tol,print_c=print_c,agg=agg,M=M,act=act)
    inds_in = np.where(sin==0)
    inds_out = np.where(sout==0)
    a[inds_out]=0
    b[inds_in]=0
    c[inds_out]=0
    d[inds_in]=0
    a = np.array(a,dtype=np.float128)
    b = np.array(b,dtype=np.float128)
    c = np.array(c,dtype=np.float128)
    d = np.array(d,dtype=np.float128)
    afake = a
    bfake = b
    cfake = c
    dfake = d
    reps = 0
    if 'x_ini' in kwargs.keys() or 'y_ini' in kwargs.keys():
        a = a*alf
        b = b*alf
    delta_k=dist_check_k(a/alf,b/alf,c,d,kout,kin,selfs,agg,M)
    delta_s=dist_check_s(a/alf,b/alf,c,d,sout,sin,selfs,agg,M)    
    ## selfs or agg ##
    extra_k,extra_s = extra_selector(selfs)
    aux_f,corr = aux_selector(agg)
    ## all ready, let's go! ##
    print "### Fixing s,k ###"
    print "Initial errors | S_in grad:%f S_out grad:%f K_in:%f K_out:%f " % (delta_s[1],delta_s[0],delta_k[1],delta_k[0])
<<<<<<< HEAD
    if np.max(np.abs(delta_s)) < tol_c and np.max(np.abs(delta_k))<tol_c:
        return a/alf,b/alf,c,d
    if print_c:
        print "## Errors: \t\t\t\t || Convergence: ##"
        print "-----------------------------------------------------------"
    while True:
        ### balancing ####
        aux  = aux_f(a,b,alf,M)
        aux2 = np.einsum('i,j',c,d)
        #strengths
        #b = sin*(alf*alf)/M/(d*np.einsum('i,ij',a*c,aux/corr(a,b,alf)/(aux2*(aux-1.)+1.)) - extra_s(a,b,c,d,aux,aux2,corr(a,b,alf))/b)
        b = sin*(alf*alf)/M/(d*np.einsum('i,ij',a*c,1./corr(a,b,alf)/(aux2*(1.-aux)+aux)) - extra_s(a,b,c,d,aux,aux2,corr(a,b,alf))/b)
        b[inds_in]=0
        aux  = aux_f(a,b,alf,M)
        #a = sout*(alf*alf)/M/(c*np.einsum('j,ij',b*d,aux/corr(a,b,alf)/(aux2*(aux-1.)+1.)) - extra_s(a,b,c,d,aux,aux2,corr(a,b,alf))/a)
        a = sout*(alf*alf)/M/(c*np.einsum('j,ij',b*d,1./corr(a,b,alf)/(aux2*(1.-aux)+aux)) - extra_s(a,b,c,d,aux,aux2,corr(a,b,alf))/a)
        a[inds_out]=0
        aux  = aux_f(a,b,alf,M)
        #degrees
        #d = kin/( np.einsum('i,ij',c,(aux-1.)/(aux2*(aux-1.)+1.)) - extra_k(c,d,aux,aux2)/d)
        d = kin/( np.einsum('i,ij',c,(1.-aux)/(aux2*(1.-aux)+aux)) - extra_k(c,d,aux,aux2)/d)
        d[inds_in]=0                                              
        aux2 = np.einsum('i,j',c,d)                               
        #c = kout/( np.einsum('j,ij',d,(aux-1.)/(aux2*(aux-1.)+1.)) - extra_k(c,d,aux,aux2)/c)
        c = kout/( np.einsum('j,ij',d,(1.-aux)/(aux2*(1.-aux)+aux)) - extra_k(c,d,aux,aux2)/c)
        c[inds_out]=0
        aux2 = np.einsum('i,j',c,d)
        ### checking convergence ####               
        tola = np.max(np.abs(afake-a))
        tolb = np.max(np.abs(bfake-b))
        tolc = np.max(np.abs(cfake-c))
        told = np.max(np.abs(dfake-d))
        if print_tol:
            print "Delta_x:%r Delta_y:%r Delta_z:%r Delta_w:%r" % (tola,tolb,tolc,told)
        if print_c and reps%act==0:
            da,db = dist_check_s(a/alf,b/alf,c,d,sout,sin,selfs,agg,M)
            dc,dd = dist_check_k(a/alf,b/alf,c,d,kout,kin,selfs,agg,M)
            if np.abs(da) < tol_c and np.abs(dc) < tol_c:
=======
    if np.max(np.abs(delta_k)) < tol_c or np.max(np.abs(delta_k))<tol_c:
        return a,b,c,d
    if not selfs:
        while True:
            aux = np.exp(np.einsum('i,j',-a/alf,b/alf))
            aux2 = np.einsum('i,j',c,d)
            #strengths
            b = sin*(alf*alf)/(d*np.einsum('i,ij',a*c,1./(aux2*(1.-aux)+aux)) -  d*a*c/(c*d*(1.-np.exp(-a/alf*b/alf))+np.exp(-a/alf*b/alf)))
            b[inds_in]=0            
            aux = np.exp(-np.einsum('i,j',a/alf,b/alf))
            a = sout*(alf*alf)/(c*np.einsum('j,ij',b*d,1./(aux2*(1.-aux)+aux)) - c*b*d/(c*d*(1.-np.exp(-a/alf*b/alf))+np.exp(-a/alf*b/alf)))
            a[inds_out]=0
            aux = np.exp(-np.einsum('i,j',a/alf,b/alf))
            #degrees
            d = kin/(  np.einsum('i,ij',c,(1.-aux)/(aux2*(1.-aux)+aux)) - np.einsum('i,i->i',c,(1.-np.exp(-a/alf*b/alf))/(c*d*(1.-np.exp(-a/alf*b/alf))+np.exp(-a/alf*b/alf))))
            d[inds_in]=0
            aux2 = np.einsum('i,j',c,d)
            c = kout/( np.einsum('j,ij',d,(1.-aux)/(aux2*(1.-aux)+aux)) - np.einsum('j,j->j',d,(1.-np.exp(-a/alf*b/alf))/(c*d*(1.-np.exp(-a/alf*b/alf))+np.exp(-a/alf*b/alf))))
            c[inds_out]=0
            #if reps>10:
            tola = np.max(np.abs(afake-a))
            tolb = np.max(np.abs(bfake-b))
            tolc = np.max(np.abs(cfake-c))
            told = np.max(np.abs(dfake-d))
            if print_tol: 
				print "Delta_x:%r Delta_y:%r Delta_z:%r Delta_w:%r" % (tola,tolb,tolc,told)
            if print_c and reps%act==0:
				da,db = dist_check_s(a/alf,b/alf,c,d,sout,sin,selfs)
				dc,dd = dist_check_k(a/alf,b/alf,c,d,kout,kin,selfs)
				if np.abs(da) < tol_c and np.abs(dc) < tol_c:
					if verbose:
						print "took %d reps, tola :%f, tolb:%f, tolc :%f, told:%f" % (reps,tola,tolb,tolc,told)
					break
				print "## Ds:%r Dk:%r" % (da,dd)
            if  tolb< tol and tola<tol and tolc<tol and told<tol:
                if verbose:
                    print "took %d reps, tola :%f, tolb:%f, tolc :%f, told:%f" % (reps,tola,tolb,tolc,told)
                break                
            if reps>maxreps:
                print "Algorithm did not converge after %d reps" % maxreps
                break
            if not all(np.isfinite(a)) or not all(np.isfinite(b)) or not all(np.isfinite(c)) or not all(np.isfinite(d)):
                raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
            afake = a
            bfake = b
            cfake = c
            dfake = d
            reps +=1
		#raise NotImplementedError("Feature not implemented")
    else:
        while True:
            aux = np.exp(np.einsum('i,j',-a/alf,b/alf))
            aux2 = np.einsum('i,j',c,d)
            #strengths
            b = sin* (alf*alf)/(d*np.einsum('i,ij',a*c,1./(aux2*(1.-aux)+aux)))
            b[inds_in]=0            
            aux = np.exp(-np.einsum('i,j',a/alf,b/alf))
            a = sout*(alf*alf)/(c*np.einsum('j,ij',b*d,1./(aux2*(1.-aux)+aux)))
            a[inds_out]=0
            aux = np.exp(-np.einsum('i,j',a/alf,b/alf))
            #degrees
            d = kin/(  np.einsum('i,ij',c,(1.-aux)/(aux2*(1.-aux)+aux)))
            d[inds_in]=0
            aux2 = np.einsum('i,j',c,d)
            c = kout/( np.einsum('j,ij',d,(1.-aux)/(aux2*(1.-aux)+aux)))            
            c[inds_out]=0
            #if reps>10:
            tola = np.max(np.abs(afake-a))
            tolb = np.max(np.abs(bfake-b))
            tolc = np.max(np.abs(cfake-c))
            told = np.max(np.abs(dfake-d))
            if print_tol: 
				print "Delta_x:%r Delta_y:%r Delta_z:%r Delta_w:%r" % (tola,tolb,tolc,told)
            if print_c and reps%act==0:
				da,db = dist_check_s(a/alf,b/alf,c,d,sout,sin)
				dc,dd = dist_check_k(a/alf,b/alf,c,d,kout,kin)
				if np.abs(da) < tol_c and np.abs(dc) < tol_c:
					if verbose:
						print "took %d reps, tola :%f, tolb:%f, tolc :%f, told:%f" % (reps,tola,tolb,tolc,told)
					break
				print "## Ds:%r Dk:%r" % (da,dd)
            if  tolb< tol and tola<tol and tolc<tol and told<tol:
>>>>>>> e4d232dadd44cc2975bd6a834f9b613f25a7c94b
                if verbose:
                    print "took %d reps, tola :%f, tolb:%f, tolc :%f, told:%f" % (reps,tola,tolb,tolc,told)
                break
            print "s:%r k:%r || Delta_x:%r Delta_y:%r Delta_z:%r Delta_w:%r" % (da,dd,tola,tolb,tolc,told)
        if  (tolb< tol and tola<tol) and (tolc<tol and told<tol):
            if verbose:
                print "took %d reps, tola :%f, tolb:%f, tolc :%f, told:%f" % (reps,tola,tolb,tolc,told)
            break                
        if reps>maxreps:
            print "Algorithm did not converge after %d reps" % maxreps
            break
        if not all(np.isfinite(a)) or not all(np.isfinite(b)) or not all(np.isfinite(c)) or not all(np.isfinite(d)) or np.max([a,b,c,d])>1e100:
            raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
        afake = a
        bfake = b
        cfake = c
        dfake = d
        reps +=1
    return a/alf,b/alf,c,d


## pre-conditioning function ##
def balance_zw(x,y,kin,kout,tol=1e-9,tol_c=1e-7,maxreps=10000,verbose=False,selfs=True,print_tol=False,print_c=False,agg=False,M=1,**kwargs):
    """ Balances the degree  equations
        Input:
            x,y: Lagrange multipliers from stength (N)
            Kin: Incoming degree sequence (N)
            Kout: Outgoing degree sequence (N)
            tol: Maximum iteration tolerance [float]
            tol_c: Maximum tolerance permitted for the worse constraint-equation [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            print_tol:  Print convergence of algorithm
            Agg: If the weighted network comes from aggreagation of binary networks set to True
            M: Number of layers (only in the case of aggregated network)
            kwargs:
                z_ini: Initial guess for z
                w_ini: Initial guess for w
                act: Print error every act steps (if print_c is set to True)
    """
    x = np.array(x,dtype=np.float128)
    y = np.array(y,dtype=np.float128)
    kin = np.array(kin,dtype=np.float128)
    kout = np.array(kout,dtype=np.float128)
    if print_tol or verbose:
        print "#### Degree balancer #####"
    n = len(kin)
    c=kwargs.get('z_ini',np.ones(n)*0.9)
    d=kwargs.get('w_ini',np.ones(n)*0.9)
    act = kwargs.get('act',100)
    alf = 1. # imposed
    inds_in = np.where(kin==0)
    inds_out = np.where(kout==0)
    c[inds_out]=0
    d[inds_in]=0
    c = np.array(c,dtype=np.float128)
    d = np.array(d,dtype=np.float128)
    cfake = c
    dfake = d
    reps = 0
    ## selfs or agg ##
    extra_k,extra_s = extra_selector(selfs)
    aux_f = aux_selector(agg)[0]
    aux = aux_f(x,y,alf,M)
    while True:
        ### balancing ####
        aux2 = np.einsum('i,j',c,d)
        #degrees
        #d = kin/( np.einsum('i,ij',c,(aux-1.)/(aux2*(aux-1.)+1.)) - extra_k(c,d,aux,aux2)/d)
        d = kin/( np.einsum('i,ij',c,(1.-aux)/(aux2*(1.-aux)+aux)) - extra_k(c,d,aux,aux2)/d)
        d[inds_in]=0                                              
        aux2 = np.einsum('i,j',c,d)                               
        #c = kout/( np.einsum('j,ij',d,(aux-1)/(aux2*(aux-1)+1)) - extra_k(c,d,aux,aux2)/c)
        c = kout/( np.einsum('j,ij',d,(1.-aux)/(aux2*(1.-aux)+aux)) - extra_k(c,d,aux,aux2)/c)
        c[inds_out]=0       
        ### checking convergence ####               
        tolc = np.max(np.abs(cfake-c))
        told = np.max(np.abs(dfake-d))        
        if print_tol:
            print "Delta_z:%r Delta_w:%r" % (tolc,told)
        if print_c and reps%act==0:
            dc,dd = dist_check_k(x,y,c,d,kout,kin,selfs,agg,M)
            if np.abs(dc) < tol_c:
                if verbose:
                    print "took %d reps, tolc :%f, told:%f" % (reps,tolc,told)
                break
            print "k:%r || Delta_z:%r Delta_w:%r" % (dd,tolc,told)
        if (tolc<tol and told<tol):
            if verbose:
                print "took %d reps, tolc :%f, told:%f" % (reps,tolc,told)
            break                
        if reps>maxreps:
            print "Algorithm did not converge after %d reps" % maxreps
            break
        if not all(np.isfinite(c)) or not all(np.isfinite(d)):
            raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
        cfake = c
        dfake = d
        reps +=1
    return c,d


###############################################
#### Check functions (distance functions ) ####
###############################################

def dist_check_s(x,y,z,w,sout,sin,selfs=True,agg=False,M=1):
    """
        Computes distance between prediction and reality for strengths
        x,y,z,w --> arrays of lagrange multpliers (kength N)
        sout,sin -- > real strength sequences (length N)
        self: True for self loops
        Agg: If the weighted network comes from aggreagation of binary networks set to True
        M: Number of layers (only in the case of aggregated network)
    """
<<<<<<< HEAD
    if not agg:
        #xy_exp = np.exp(-np.einsum('i,j',x,y))
        xy_exp = np.exp(np.einsum('i,j',x,y))
        zw = np.einsum('i,j',z,w)
        if not selfs:
            #extra1 = np.einsum('i,i,ii->i',x*z,y*w,1./(zw*(1.-xy_exp)+xy_exp)) # diagonal
            extra1 = np.einsum('i,i,ii->i',x*z,y*w,xy_exp/(zw*(xy_exp-1)+1)) # diagonal
        else:
            extra1=0
        #delta_x = np.abs(sout - x*z*np.einsum('j,ij',y*w,1./(zw*(1.-xy_exp)+xy_exp)) + extra1).sum()
        #delta_y = np.abs(sin - y*w*np.einsum('i,ij',x*z,1./(zw*(1.-xy_exp)+xy_exp)) + extra1).sum()
        delta_x = np.abs(sout - x*z*np.einsum('j,ij',y*w,xy_exp/(zw*(xy_exp-1)+1)) + extra1).sum()
        delta_y = np.abs(sin - y*w*np.einsum('i,ij',x*z,xy_exp/(zw*(xy_exp-1)+1)) + extra1).sum()
    else:
        #xy = (np.einsum('i,j',x,y)+1)**(-M)
        xy = (np.einsum('i,j',x,y)+1)**(M)
        zw = np.einsum('i,j',z,w)
        if not selfs:
            #extra1 = M*x*y/(1+x*y) * z*w/((1+x*y)**(-M) + z*w*(1-(1+x*y)**(-M)))
            extra1 = M*x*y/(1+x*y) * z*w*(1+x*y)**(M)/(1+ z*w*((1+x*y)**(M)-1))
        else:
            extra1 = 0
        #delta_x = np.abs(sout - M*x*z*np.einsum('j,ij',y*w,1./(1+np.einsum('i,j',x,y))/(xy + zw *(1-xy))) + extra1).sum()
        #delta_y = np.abs(sin - M*y*w*np.einsum('i,ij',x*z,1./(1+np.einsum('i,j',x,y))/(xy + zw *(1-xy))) + extra1).sum()
        delta_x = np.abs(sout - M*x*z*np.einsum('j,ij',y*w,xy/(1+np.einsum('i,j',x,y))/(1 + zw *(xy-1))) + extra1).sum()
        delta_y = np.abs(sin - M*y*w*np.einsum('i,ij',x*z,xy /(1+np.einsum('i,j',x,y))/(1 + zw *(xy-1))) + extra1).sum()
    return delta_y,delta_x 
=======
    xy_exp = np.exp(-np.einsum('i,j',x,y))
    zw = np.einsum('i,j',z,w)
    if not selfs:
        extra1 = np.einsum('i,i,ii->i',x*z,y*w,1./(zw*(1.-xy_exp)+xy_exp)) # diagonal
    else:
        extra1=0
    #delta_x = sout - x*z*np.einsum('j,ij',y*w,xy_exp/(zw*(xy_exp-1.)+1.)) + extra1
    #delta_y = sin - w*y*np.einsum('i,ij',x*z,xy_exp/(zw*(xy_exp-1.)+1.)) +  extra1
    #return delta_x.sum(),delta_y.sum()
    delta_x = np.abs(sout - x*z*np.einsum('j,ij',y*w,1./(zw*(1.-xy_exp)+xy_exp)) + extra1).sum()
    return delta_x,delta_x
>>>>>>> e4d232dadd44cc2975bd6a834f9b613f25a7c94b

def dist_check_k(x,y,z,w,kout,kin,selfs=True,agg=False,M=1):
    """
        Computes distance between prediction and reality for degrees in absolute value
        x,y,z,w --> arrays of lagrange multpliers
        kout,kin -- > real degree sequences
        self: True for self loops
        Agg: If the weighted network comes from aggreagation of binary networks set to True
        M: Number of layers (only in the case of aggregated network)
    """
<<<<<<< HEAD
    if not agg:
        #xy_exp = np.exp(-np.einsum('i,j',x,y))
        xy_exp = np.exp(np.einsum('i,j',x,y))
        zw = np.einsum('i,j',z,w)
        if not selfs:
            #extra1 = z*w*(1-np.exp(-x*y))/(z*w*(1.-np.exp(-x*y)) + np.exp(-x*y))
            extra1 = z*w*(np.exp(x*y)-1)/(z*w*(np.exp(x*y)-1) + 1)
            #np.einsum('i,i->i',z*w,(1.-np.exp(-x*y))/(z*w*(1.-np.exp(-x*y))+np.exp(-x*y))) # diagonal
        else:
            extra1=0
        #delta_z = np.abs(kout - z*np.einsum('j,ij',w,(1.-xy_exp)/(zw*(1.-xy_exp)+xy_exp)) + extra1).sum()
        #delta_w = np.abs(kin - w*np.einsum('i,ij',z,(1.-xy_exp)/(zw*(1.-xy_exp)+xy_exp)) + extra1).sum()
        delta_z = np.abs(kout - z*np.einsum('j,ij',w,(xy_exp-1)/(zw*(xy_exp-1)+1)) + extra1).sum()
        delta_w = np.abs(kin - w*np.einsum('i,ij',z,(xy_exp-1)/(zw*(xy_exp-1)+1)) + extra1).sum()

    else:
        #xy = (1+np.einsum('i,j',x,y))**(-M)
        xy = (1+np.einsum('i,j',x,y))**(M)
        zw = np.einsum('i,j',z,w)       
        if not selfs:
            #extra1 = z*w*(1-(1+x*y)**(-M))/((1+x*y)**(-M)+z*w*(1-(1+x*y)**(-M)))
            extra1 = z*w*((1+x*y)**(M)-1)/(1+z*w*((1+x*y)**(M)-1))
        else:
            extra1 = 0
        #delta_z = np.abs(kout - z*np.einsum('j,ij',w,(1-xy)/(xy+zw*(1-xy))) + extra1).sum()     
        #delta_w = np.abs(kin - w*np.einsum('i,ij',z,(1-xy)/(xy+zw*(1-xy))) + extra1).sum()     
        delta_z = np.abs(kout - z*np.einsum('j,ij',w,(xy-1)/(1+zw*(xy-1))) + extra1).sum()     
        delta_w = np.abs(kin - w*np.einsum('i,ij',z,(xy-1)/(1+zw*(xy-1))) + extra1).sum()     
    return delta_w,delta_z




"""
##############################
##### Part balancing: degrees, then strengths, not good convergence #########
#WARNING: NOT GOOD CONVERGENCE, Deprecated.
##############################
"""



#def balance_zw(x,y,kin,kout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,print_tol=False,**kwargs):
    #""" Balances the degree equations
        #Input:
            #x,y: Lagrange multipliers from stength (N)
            #Kin: Incoming degree sequence (N)
            #Kout: Outgoing degree sequence (N)
            #tol: Maximum tolerance [float]
            #Maxreps: Maximum number of iterations [float]
            #Verbose [Bool]
            #Selfs: Allow self loops? [True for yes]
            #kwargs:
                #z_ini: Initial guess for z
                #w_ini: Initial guess for w
    #"""
    #if print_tol or verbose:
        #print "#### Degree balancer #####"
    #n = len(kin)
    #c=kwargs.get('z_ini',np.ones(n)*0.9)
    #d=kwargs.get('w_ini',np.ones(n)*0.9)
    #inds_in = np.where(kin==0)
    #inds_out = np.where(kout==0)
    #c[inds_out]=0
    #d[inds_in]=0
    #cfake = c
    #dfake = d
    #reps = 0
    #aux = np.exp(np.einsum('i,j',x,y))
    #if not selfs:
        #raise NotImplementedError("Feature not implemented")
    #else:
        #while True:
            #aux2 = np.einsum('i,j',c,d)
            ##degrees
            #d = kin/np.einsum('i,ij',c,(aux-1.)/(aux2*(aux-1.)+1.))
            #d[inds_in]=0
            #aux2 = np.einsum('i,j',c,d)
            #c = kout/np.einsum('j,ij',d,(aux-1.)/(aux2*(aux-1.)+1.))
            #c[inds_out]=0
            #tolc = np.max(np.abs(cfake-c))
            #told = np.max(np.abs(dfake-d))
            #if print_tol: print "Delta_z:%r Delta_w:%r" % (tolc,told)
            #if  tolc<tol and told<tol:
                #break
                #if verbose:
                    #print "took %d reps, tolc :%f, told:%f" % (reps,tolc,told)
            #if reps>maxreps:
                #print "Algorithm did not converge after %d reps" % maxreps
                #break
            #if not all(np.isfinite(c)) or not all(np.isfinite(d)):
                #raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
            #cfake = c
            #dfake = d
            #reps +=1
    #return c,d



#def balance_xy(z,w,sin,sout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,print_tol=False,**kwargs):
    ###### does not work nor converge properly...  #####
    #""" Balances the strength equations
        #Input:
            #z,w: Lagrange multipliers from degree (N)
            #Sin: Incoming strength sequence (N)
            #Sout: Outgoing strength sequence (N)
            #tol: Maximum tolerance [float]
            #Maxreps: Maximum number of iterations [float]
            #Verbose [Bool]
            #Selfs: Allow self loops? [True for yes]
            #kwargs:
                #x_ini: Initial guess for x
                #y_ini: Initial guess for y
    #"""
    #if verbose or print_tol:
        #print "#### Strength balancer #####"
    #n = len(sin)
    #alf = np.sqrt(sout.max()*sin.max()/sout.sum())
    #a=kwargs.get('x_ini',sout/np.sqrt(T))
    #b=kwargs.get('y_ini',sin/np.sqrt(T))
    #if 'x_ini' in kwargs.keys() or 'y_ini' in kwargs.keys():
        #a = alf*a
        #b = alf*b
    #inds_in = np.where(sin==0)
    #inds_out = np.where(sout==0)
    #a[inds_out]=0
    #b[inds_in]=0
    #afake = a
    #bfake = b
    #reps = 0
    #aux2 = np.einsum('i,j',z,w)
    #if not selfs:
        #raise NotImplementedError("Feature not implemented")
    #else:
        #while True:
            #aux = np.exp(np.einsum('i,j',a/alf,b/alf))
            #b = sin*(alf*alf)/(w*np.einsum('i,ij',a*z,aux/(aux2*(aux-1.)+1.)))
            #b[inds_in]=0
            #aux = np.exp(np.einsum('i,j',a/alf,b/alf))
            #a = sout*(alf*alf)/(z*np.einsum('j,ij',b*w,aux/(aux2*(aux-1.)+1.)))
            #a[inds_out]=0
            #tola = np.max(np.abs(afake-a))
            #tolb = np.max(np.abs(bfake-b))
            #if print_tol: print "Delta_x:%r Delta_y:%r" % (tola,tolb)
            #if  tolb< tol and tola<tol:
                #break
                #if verbose:
                    #print "took %d reps, tola :%f, tolb:%f" % (reps,tola,tolb)
            #if reps>maxreps:
                #print "Algorithm did not converge after %d reps" % maxreps
                #break
            #if not all(np.isfinite(a)) or not all(np.isfinite(b)):
                #raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
            #afake = a
            #bfake = b
            #reps +=1
    #return a/alf,b/alf


#def balance_xy_zw(sin,sout,kin,kout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,print_tol=False,**kwargs):
    #""" Balances the strenght and degree equations in parts 
        #Input:
            #Sin: Incoming strength sequence (N)
            #Sout: Outgoing strength sequence (N)
            #Kin: Incoming degree sequence (N)
            #Kout: Outgoing degree sequence (N)
            #tol: Maximum tolerance [float]
            #Maxreps: Maximum number of iterations [float]
            #Verbose [Bool]
            #Selfs: Allow self loops? [True for yes]
            #kwargs:
                #x_ini: Initial guess for x
                #y_ini: Initial guess for y
                #z_ini: Initial guess for z
                #w_ini: Initial guess for w 
    #"""
    #if verbose or print_tol:
        #print "#### Strength & degree part balancer #####"
    #n = len(sin)
    #T = sout.sum()
    #a=kwargs.get('x_ini',np.ones(n)*0.9)
    #b=kwargs.get('y_ini',np.ones(n)*0.9)
    #c=kwargs.get('z_ini',np.ones(n)*0.9)
    #d=kwargs.get('w_ini',np.ones(n)*0.9)
    #inds_in = np.where(sin==0)
    #inds_out = np.where(sout==0)
    #b[inds_in]=0
    #d[inds_in]=0
    #c[inds_out]=0
    #a[inds_out]=0
    #afake = a
    #bfake = b
    #cfake = c
    #dfake = d
    #reps = 0
    #if not selfs:
        #raise NotImplementedError("Feature not implemented")
    #else:
        #while True:
            #c,d = balance_zw(a,b,kin,kout,tol,maxreps,verbose,selfs,print_tol,z_ini=c,w_ini=d) # fix degrees
            #a,b = balance_xy(c,d,sin,sout,tol,maxreps,verbose,selfs,print_tol,x_ini=a,y_ini=b) # fix strengths
            #tola,tolb = dist_check_s(a,b,c,d,sout,sin)
            #tolc,told = dist_check_k(a,b,c,d,kout,kin)
            #if print_tol: print "Delta_x:%r Delta_y:%r Delta_z:%r Delta_w:%r" % (tola,tolb,tolc,told)
            #if  tolb< tol and tola<tol and tolc<tol and told<tol:
                #break
                #if verbose:
                    #print "took %d reps, tol_sout :%f, tol_sin:%f, tol_kout :%f, tol_kin:%f" % (reps,tola,tolb,tolc,told)
            #if reps>maxreps:
                #print "Algorithm did not converge after %d reps" % maxreps
                #break
            #if not all(np.isfinite(a)) or not all(np.isfinite(b)) or not all(np.isfinite(c)) or not all(np.isfinite(d)):
                #raise ValueError("Something is wrong, algorithm did not converge. Check constraints or change method")
            #afake = a
            #bfake = b
            #cfake = c
            #dfake = d
            #reps +=1
    #return a,b,c,d
     
###############################     
=======
    xy_exp = np.exp(-np.einsum('i,j',x,y))
    zw = np.einsum('i,j',z,w)
    if not selfs:
        extra1 = np.einsum('i,i->i',z*w,(1.-np.exp(-x*y))/(z*w*(1.-np.exp(-x*y))+np.exp(-x*y))) # diagonal
    else:
        extra1=0
    #delta_z = (kout - z*np.einsum('j,ij',w,(xy_exp-1.)/(zw*(xy_exp-1.)+1.)) + extra1)
    #delta_w = kin - w*np.einsum('i,ij',z,(xy_exp-1.)/(zw*(xy_exp-1.)+1.)) +  extra1
    #return delta_z.sum(),delta_w.sum()
    delta_z = np.abs(kout - z*np.einsum('j,ij',w,(1.-xy_exp)/(zw*(1.-xy_exp)+xy_exp)) + extra1).sum() # equal for both
    return delta_z,delta_z
>>>>>>> e4d232dadd44cc2975bd6a834f9b613f25a7c94b
