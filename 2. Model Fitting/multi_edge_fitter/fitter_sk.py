"""
 Fitter for fixed strength sequence and fixed degree sequence
"""

## Imports ##
from scipy import optimize as opt
import numpy as np



"""
#############################
#### Balancing approach:  #########
#############################
"""


"""
#############################
#### Sequential balancing: all at once (much better, preferred approach) #########
#############################
"""


def balance_xyzw(sin,sout,kin,kout,tol=1e-9,tol_c=10,maxreps=10000,verbose=False,selfs=True,print_tol=False,print_c=False,**kwargs):
    """ Balances the strength and degree  equations
        Input:
            Sin: Incoming strength sequence (N)
            Sout: Outgoing strength sequence (N)
            Kin: Incoming degree sequence (N)
            Kout: Outgoing degree sequence (N)
            tol: Maximum iteration tolerance [float]
            tol_c: Maximum tolerance equations [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            print_tol:  Print convergence of algorithm
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
                z_ini: Initial guess for z
                w_ini: Initial guess for w
    """
    if verbose or print_tol:
		print "#### Degree & Strength balancer #####"
    n = len(sin)
    T = sout.sum()
    #a=kwargs.get('x_ini',sout/np.sqrt(T))
    #b=kwargs.get('y_ini',sin/np.sqrt(T))
    a=kwargs.get('x_ini',np.ones(n)*0.9)
    b=kwargs.get('y_ini',np.ones(n)*0.9)
    c=kwargs.get('x_ini',np.ones(n)*0.9)
    d=kwargs.get('y_ini',np.ones(n)*0.9)
    #c=kwargs.get('w_ini',sout)
    #d=kwargs.get('z_ini',sin)
    #c=kwargs.get('z_ini',np.ones(n)*0.9)
    #d=kwargs.get('w_ini',np.ones(n)*0.9)
    #if 'z_ini' not in kwargs.keys() or 'w_ini' not in kwargs.keys():
	#	c,d = balance_zw(a,b,kin,kout,tol,maxreps,verbose,selfs,print_tol)
    #else:
	#	c=kwargs['z_ini']
	#	d=kwargs['w_ini']
    inds_in = np.where(sin==0)
    inds_out = np.where(sout==0)
    b[inds_in]=0
    d[inds_in]=0
    c[inds_out]=0
    a[inds_out]=0
    afake = a
    bfake = b
    cfake = c
    dfake = d
    reps = 0

    delta_k=dist_check_k(a,b,c,d,kout,kin,selfs)
    delta_s=dist_check_k(a,b,c,d,sout,sin,selfs)
    print "### Fixing s,k ###"
    print "Initial errors | S_in grad:%f S_out grad:%f K_in:%f K_out:%f " % (delta_s[1],delta_s[0],delta_k[1],delta_k[0])
    if np.max(np.abs(delta_k)) < tol_c or np.max(np.abs(delta_k))<tol_c:
        return a,b,c,d
    if not selfs:
        while True:
            aux = np.exp(np.einsum('i,j',a,b))
            aux2 = np.einsum('i,j',c,d)
            #strengths
            b = sin/(d*np.einsum('i,ij',a*c,aux/(aux2*(aux-1.)+1.)) - d*np.einsum('i,i->i',a*c,np.exp(a*b)/(c*d*(np.exp(a*b)-1.)+1.)))
            b[inds_in]=0
            aux = np.exp(np.einsum('i,j',a,b))
            a = sout/(c*np.einsum('j,ij',b*d,aux/(aux2*(aux-1.)+1.)) - c*np.einsum('j,j->j',b*d,np.exp(a*b)/(c*d*(np.exp(a*b)-1.)+1.)))
            a[inds_out]=0
            aux = np.exp(np.einsum('i,j',a,b))
            #degrees
            d = kin/( np.einsum('i,ij',c,(aux-1.)/(aux2*(aux-1.)+1.)) - np.einsum('i,i->i',c,(np.exp(a*b)-1.)/(c*d*(np.exp(a*b)-1.)+1.)))
            d[inds_in]=0
            aux2 = np.einsum('i,j',c,d)
            c = kout/( np.einsum('j,ij',d,(aux-1.)/(aux2*(aux-1.)+1.)) - np.einsum('j,j->j',d,(np.exp(a*b)-1.)/(c*d*(np.exp(a*b)-1.)+1.)))
            c[inds_out]=0
            #if reps>10:
            tola = np.max(np.abs(afake-a))
            tolb = np.max(np.abs(bfake-b))
            tolc = np.max(np.abs(cfake-c))
            told = np.max(np.abs(dfake-d))
            if print_tol: 
				print "Delta_x:%r Delta_y:%r Delta_z:%r Delta_w:%r" % (tola,tolb,tolc,told)
            if print_c and reps%10==0:
				da,db = dist_check_s(a,b,c,d,sout,sin,selfs)
				dc,dd = dist_check_k(a,b,c,d,kout,kin,selfs)
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
            afake = a
            bfake = b
            cfake = c
            dfake = d
            reps +=1
		#raise NotImplementedError("Feature not implemented")
    else:
        while True:
            aux = np.exp(np.einsum('i,j',a,b))
            aux2 = np.einsum('i,j',c,d)
            #strengths
            b = sin/(d*np.einsum('i,ij',a*c,aux/(aux2*(aux-1.)+1.)))
            b[inds_in]=0
            aux = np.exp(np.einsum('i,j',a,b))
            a = sout/(c*np.einsum('j,ij',b*d,aux/(aux2*(aux-1.)+1.)))
            a[inds_out]=0
            aux = np.exp(np.einsum('i,j',a,b))
            #degrees
            d = kin/np.einsum('i,ij',c,(aux-1.)/(aux2*(aux-1.)+1.))
            d[inds_in]=0
            aux2 = np.einsum('i,j',c,d)
            c = kout/np.einsum('j,ij',d,(aux-1.)/(aux2*(aux-1.)+1.))
            c[inds_out]=0
            #if reps>10:
            tola = np.max(np.abs(afake-a))
            tolb = np.max(np.abs(bfake-b))
            tolc = np.max(np.abs(cfake-c))
            told = np.max(np.abs(dfake-d))
            if print_tol: 
				print "Delta_x:%r Delta_y:%r Delta_z:%r Delta_w:%r" % (tola,tolb,tolc,told)
            if print_c and reps%10==0:
				da,db = dist_check_s(a,b,c,d,sout,sin)
				dc,dd = dist_check_k(a,b,c,d,kout,kin)
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
            afake = a
            bfake = b
            cfake = c
            dfake = d
            reps +=1
    return a,b,c,d



"""
#############################
#### Part balancing: degrees, then strengths, not good convergence #########
WARNING: NOT GOOD CONVERGENCE, USE ABOVE METHOD IN PREFERENCE!!!!
#############################
"""
import numpy as np
from scipy import optimize as opt


def balance_zw(x,y,kin,kout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,print_tol=False,**kwargs):
    """ Balances the degree equations
        Input:
			x,y: Lagrange multipliers from stength (N)
            Kin: Incoming degree sequence (N)
            Kout: Outgoing degree sequence (N)
            tol: Maximum tolerance [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            kwargs:
                z_ini: Initial guess for z
                w_ini: Initial guess for w
    """
    if print_tol or verbose:
		print "#### Degree balancer #####"
    n = len(kin)
    c=kwargs.get('z_ini',np.ones(n)*0.9)
    d=kwargs.get('w_ini',np.ones(n)*0.9)
    inds_in = np.where(kin==0)
    inds_out = np.where(kout==0)
    c[inds_out]=0
    d[inds_in]=0
    cfake = c
    dfake = d
    reps = 0
    aux = np.exp(np.einsum('i,j',x,y))
    if not selfs:
		raise NotImplementedError("Feature not implemented")
    else:
        while True:
            aux2 = np.einsum('i,j',c,d)
            #degrees
            d = kin/np.einsum('i,ij',c,(aux-1.)/(aux2*(aux-1.)+1.))
            d[inds_in]=0
            aux2 = np.einsum('i,j',c,d)
            c = kout/np.einsum('j,ij',d,(aux-1.)/(aux2*(aux-1.)+1.))
            c[inds_out]=0
            tolc = np.max(np.abs(cfake-c))
            told = np.max(np.abs(dfake-d))
            if print_tol: print "Delta_z:%r Delta_w:%r" % (tolc,told)
            if  tolc<tol and told<tol:
                break
                if verbose:
                    print "took %d reps, tolc :%f, told:%f" % (reps,tolc,told)
            if reps>maxreps:
                print "Algorithm did not converge after %d reps" % maxreps
                break
            cfake = c
            dfake = d
            reps +=1
    return c,d



def balance_xy(z,w,sin,sout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,print_tol=False,**kwargs):
	##### does not work nor converge properly...  #####
    """ Balances the strength equations
        Input:
			z,w: Lagrange multipliers from degree (N)
            Sin: Incoming strength sequence (N)
            Sout: Outgoing strength sequence (N)
            tol: Maximum tolerance [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
    """
    if verbose or print_tol:
		print "#### Strength balancer #####"
    n = len(sin)
    a=kwargs.get('x_ini',sout/np.sqrt(T))
    b=kwargs.get('y_ini',sin/np.sqrt(T))
    inds_in = np.where(sin==0)
    inds_out = np.where(sout==0)
    a[inds_out]=0
    b[inds_in]=0
    afake = a
    bfake = b
    reps = 0
    aux2 = np.einsum('i,j',z,w)
    if not selfs:
		raise NotImplementedError("Feature not implemented")
    else:
        while True:
            aux = np.exp(np.einsum('i,j',a,b))
            b = sin/(w*np.einsum('i,ij',a*z,aux/(aux2*(aux-1.)+1.)))
            b[inds_in]=0
            aux = np.exp(np.einsum('i,j',a,b))
            a = sout/(z*np.einsum('j,ij',b*w,aux/(aux2*(aux-1.)+1.)))
            a[inds_out]=0
            tola = np.max(np.abs(afake-a))
            tolb = np.max(np.abs(bfake-b))
            if print_tol: print "Delta_x:%r Delta_y:%r" % (tola,tolb)
            if  tolb< tol and tola<tol:
                break
                if verbose:
                    print "took %d reps, tola :%f, tolb:%f" % (reps,tola,tolb)
            if reps>maxreps:
                print "Algorithm did not converge after %d reps" % maxreps
                break
            afake = a
            bfake = b
            reps +=1
    return a,b


def balance_xy_zw(sin,sout,kin,kout,tol=1e-9,maxreps=10000,verbose=False,selfs=True,print_tol=False,**kwargs):
    """ Balances the strenght and degree equations in parts 
        Input:
            Sin: Incoming strength sequence (N)
            Sout: Outgoing strength sequence (N)
            Kin: Incoming degree sequence (N)
            Kout: Outgoing degree sequence (N)
            tol: Maximum tolerance [float]
            Maxreps: Maximum number of iterations [float]
            Verbose [Bool]
            Selfs: Allow self loops? [True for yes]
            kwargs:
                x_ini: Initial guess for x
                y_ini: Initial guess for y
                z_ini: Initial guess for z
                w_ini: Initial guess for w	
    """
    if verbose or print_tol:
        print "#### Strength & degree part balancer #####"
    n = len(sin)
    T = sout.sum()
    a=kwargs.get('x_ini',np.ones(n)*0.9)
    b=kwargs.get('y_ini',np.ones(n)*0.9)
    c=kwargs.get('z_ini',np.ones(n)*0.9)
    d=kwargs.get('w_ini',np.ones(n)*0.9)
    inds_in = np.where(sin==0)
    inds_out = np.where(sout==0)
    b[inds_in]=0
    d[inds_in]=0
    c[inds_out]=0
    a[inds_out]=0
    afake = a
    bfake = b
    cfake = c
    dfake = d
    reps = 0
    if not selfs:
		raise NotImplementedError("Feature not implemented")
    else:
        while True:
            c,d = balance_zw(a,b,kin,kout,tol,maxreps,verbose,selfs,print_tol,z_ini=c,w_ini=d) # fix degrees
            a,b = balance_xy(c,d,sin,sout,tol,maxreps,verbose,selfs,print_tol,x_ini=a,y_ini=b) # fix strengths
            tola,tolb = dist_check_s(a,b,c,d,sout,sin)
            tolc,told = dist_check_k(a,b,c,d,kout,kin)
            if print_tol: print "Delta_x:%r Delta_y:%r Delta_z:%r Delta_w:%r" % (tola,tolb,tolc,told)
            if  tolb< tol and tola<tol and tolc<tol and told<tol:
                break
                if verbose:
                    print "took %d reps, tol_sout :%f, tol_sin:%f, tol_kout :%f, tol_kin:%f" % (reps,tola,tolb,tolc,told)
            if reps>maxreps:
                print "Algorithm did not converge after %d reps" % maxreps
                break
            afake = a
            bfake = b
            cfake = c
            dfake = d
            reps +=1
    return a,b,c,d
     
##############################     

#### Check functions (distance functions ) ####

def dist_check_s(x,y,z,w,sout,sin,selfs=True):
    """
        Computes distance between prediction and reality for strengths
        x,y,z,w --> arrays of lagrange multpliers (kength N)
        sout,sin -- > real strength sequences (length N)
        self: True for self loops
    """
    xy_exp = np.exp(np.einsum('i,j',x,y))
    zw = np.einsum('i,j',z,w)
    if not selfs:
        extra1 = np.einsum('i,i,ii->i',x*z,y*w,xy_exp/(zw*(xy_exp-1.)+1.)) # diagonal
    else:
        extra1=0
    #delta_x = sout - x*z*np.einsum('j,ij',y*w,xy_exp/(zw*(xy_exp-1.)+1.)) + extra1
    #delta_y = sin - w*y*np.einsum('i,ij',x*z,xy_exp/(zw*(xy_exp-1.)+1.)) +  extra1
    #return delta_x.sum(),delta_y.sum()
    delta_x = (sout - x*z*np.einsum('j,ij',y*w,xy_exp/(zw*(xy_exp-1.)+1.)) + extra1).sum()
    return delta_x,delta_x


def dist_check_k(x,y,z,w,kout,kin,selfs=True):
    """
        Computes distance between prediction and reality for degrees
        x,y,z,w --> arrays of lagrange multpliers
        kout,kin -- > real degree sequences
        self: True for self loops
    """
    xy_exp = np.exp(np.einsum('i,j',x,y))
    zw = np.einsum('i,j',z,w)
    if not selfs:
        extra1 = np.einsum('i,i->i',z*w,(np.exp(x*y)-1.)/(z*w*(np.exp(x*y)-1.)+1.)) # diagonal
    else:
        extra1=0
    #delta_z = (kout - z*np.einsum('j,ij',w,(xy_exp-1.)/(zw*(xy_exp-1.)+1.)) + extra1)
    #delta_w = kin - w*np.einsum('i,ij',z,(xy_exp-1.)/(zw*(xy_exp-1.)+1.)) +  extra1
    #return delta_z.sum(),delta_w.sum()
    delta_z = (kout - z*np.einsum('j,ij',w,(xy_exp-1.)/(zw*(xy_exp-1.)+1.)) + extra1).sum() # equal for both
    return delta_z,delta_z
