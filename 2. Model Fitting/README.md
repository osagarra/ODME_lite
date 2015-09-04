Origin-Destination Multi-Edge Package. Module 2: Max Entropy Fitter
========================================================================

 Copyright 2014 Oleguer Sagarra. All rights reserved. Code under License GPLv3.
______________________________________________________________________________________

The present program may be used to fit the lagragian multipliers that are needed to generate maximum entropy ensembles for a variety of constraints and types of weighted networks. 
This lagrangian multipliers allow to later generate networks belonging to ensembles of many different types of networks (weighted, multi-edge, accumulated weighted and accumulated binary networks) with a set of prescribed properties (such as degree sequence, strength sequence, total cost...). This procedure can be either used to generate networks to model several phenomena or to assess relevance of features detected in real data. It is very usefull for hypotehsis testing.

If you use our software, please do cite us.

## References 

[1] Statistical mechanics of multiedge networks
	Sagarra O., Pérez-Vicente C. and Díaz-Guilera, A.  Phys. Rev. E 88, 062806 (2013)
    [Physical Review E](http://pre.aps.org/abstract/PRE/v88/i6/e062806)

<<<<<<< HEAD
[2] Supersampling and network reconstruction of urban mobility.
	Sagarra, O., Szell, M., Santi, P., and Ratti, C. Arxiv Preprint (2015)
	[arXiv:1504.01939v1](http://arxiv.org/abs/1504.01939)
=======
[2] The configuration multi-edge model: Assessing the effect of fixing node strengths on weighted network magnitudes
	Sagarra O., Font-Clos F., Pérez-Vicente C. and Díaz-Guilera, A.
	[Europhysics Letters](http://iopscience.iop.org/0295-5075/107/3/38002/article;jsessionid=D22CAAF312F43653DA0C1279853CBF0C.c3)

[3] Supersampling and network reconstruction of urban mobility.
	Sagarra, O., Szell, M., Santi, P., and Ratti, C. Arxiv Preprint (2015)
	[arXiv:1504.01939v1](http://arxiv.org/abs/1504.01939)

>>>>>>> e4d232dadd44cc2975bd6a834f9b613f25a7c94b

## Requirements and Installation


To install, simply use the usual python command,
```
	$ python setup.py install
```		

Make sure the python dependencies are installed (subdependencies of these packages not listed):
```
	- Numpy
	- Scipy
	- Numba (preferably numbapro)
```

If one wants to use the refined solver CVXOPT (see examples), then one can also OPTIONALLY install CVXOPT suite from [http://cvxopt.org](http://cvxopt.org).

The installation instructions can be checked [here](http://cvxopt.org/install/index.html).
 
## Brief description of the Max Entropy Fitter

This python package may be used to obtain the lagrange multipliers which then can be used by module 3 of the ODME package to generate maximum entropy network models according to some fixed constraints.
To do so, a set of saddle point equations need to be solved (for theoretical details, see reference [1]).

The present software considers up to 4 different types of weighted networks:
	
	- Multi-Edge Networks (ME): Single networks composed of distinguishable events.
	- Weighted Networks (W): Single networks composed of indistinguishable events.
	- Accumulated Binary Networks (AB): Networks composed by the aggregation of M layers of binary networks.
	- Accumulated Weighted Networks (AW): Networks composed by the aggregation of M layers of weighted networks.

The present software is able to solve the saddle point equations for 5 different cases (in parethesis, the types of networks available):
	
	0. Fixing the strength sequence. (ME,W,AB,AW)
	1. FIxing the strength sequence (incoming and outgoing - ME,AB,W,AW) and an average trip cost constraint. (only ME)
	2. Fixing the strength sequence and some individual trip entries t_{ij}. (ME)
	3. Fixing the strength sequence and the total number of binary edges. (ME,AB)
	4. Fixing the degree and strength sequence. (ME,AB) [EXPERIMENTAL]
	5. Fixing the degree sequence. (ME,W,AB,AW, note that in this case the equations for all cases look the same, hence the lagrangian multipliers are also the same).

	
It allows for considering self loops or not.
	
Each of the cases uses a different module:

	0. Uses fitter_s.py (and fitter_s_CVXOPT.py). Accepts also the case of indistiguishable events and aggregated binary events.
		Inputs:
			Incoming strength sequence [array of length N number of nodes] 
			Outgoing strength sequence [array of length N number of nodes] 
		Relevant function to use from the package:
			balance_xy --> Returns x,y. Lagrange multipliers such that:
				- if indist (W,AW cases) is set to true: <t_{ij}> = M x y /(1-xy) where M is the number of considered layers. 
				(Note: This particular case and implementation is done by brute force maximization).
				- if agg (AB case) is set to true: <t_{ij}> = M x y /(1+xy) where M is the number of considered layers
				- any other case (ME), normal Multi Edge network: <t_{ij}> = x y
				and keeps the strength sequence provided).
	1. Uses fitter_grav.py
		Inputs:
			Incoming strength sequence [array of length N number of nodes] 
			Outgoing strength sequence [array of length N number of nodes]
			Distance matrix [matrix NxN]
			Total cost C [float number]
		Relevant function to use from the package:
			fit_gamma --> Returns x,y,gamma (lagrange multipliers such that <t_{ij}> = x y exp(-gamma d_{ij}))
		
	2. Uses fitter_pij.py
		Inputs:
			Incoming strength sequence [array of length N number of nodes] 
			Outgoing strength sequence [array of length N number of nodes]
			Trips one wants to fix. One row per trip: [origin_node_id destination_node_id t_ij]
		Relevant function to use from the package:
			balance_qij_x_y --> Returns p_ij matrix (already normalized) (NxN matrix such that <t_{ij}> = T p_{ij}  where T is the total number of desired events).
	
	3. Uses fitter_E.py (**Experimental, convergence not assured!**)
		Two implementations are provided both with convergence non guaranteed: a balancing approach and a likelyhood maximization approach.
		We recommend to use first the balancing approach and then the likelyhood maximization approach (See the examples).		
		Input:
			Incoming strength sequence [array of length N number of nodes] 
			Outgoing strength sequence [array of length N number of nodes]
			Total number of binary edges (E) [float].
		In this case it uses a custom class or the function fit_lambda. Returns sets of lagrange multipliers x,y,lambda (Lagrange multipliers such that:
				- IF ME case: <t_{ij}> = lam* x y exp(xy) / (lam*(exp(xy)-1)+1) )
				- IF AB case: <t_{ij}> = lam* M* x y /(1 + x*y) *(1+x*y)**M / (lam*((1+x*y)**M-1)+1) 
	4. Uses fitter_sk.py (**Experimental, convergence not assured!**)
		Inputs:
			Incoming strength sequence [array of length N number of nodes] 
			Outgoing strength sequence [array of length N number of nodes]
		Relevant function to use from the package:
			balance_xyzw --> Returns x,y,z,w (lagrange multipliers such that:
				- If AB case: <t_{ij}> = M xy /(1+xy) zw (1+xy)^M / (1+zw((1+xy)^M -1 ))
				- If ME case: <t_{ij}> = x y lam* zw exp(zw) / (x y lam*(exp(zw)-1)+1) )
        (Note: This particular algorithm is stable, but its convergence is lengthly)
	5. Uses fitter_k.py
		Inputs:
			Incoming and outgoing degree [array of lenght N number of nodes]
		Relevant function to use from the package:
			balance_xy --> Returns x,y (lagrange multipliers such that <t_{ij}> = xy/(1+xy))

All the functions are fully documented. Some examples are provided in folder "examples".


## Examples

In the folder *examples* some example code for the different cases (0,1,2,3,4,5) listed above can be found.


## Performance

Since the code uses squared non-sparse numpy arrays (NxN), some memory problems might be encountered when considering large networks (number of nodes N>10000) 
depending on the memory of the available computer.

## License

(C) Copyright 2014 Oleguer Sagarra.
All rights reserved. 
Code under License GPLv3.

This module is a part of the ODME package, and thus it is distributed under its same license (see main README).

This module depends on scipy, CVXOPT and numba which are not distributed with this package. Please check the kind of license that applies in each case.

## Roadmap

The software is for the moment ready and works smoothely under apropiate conditions.
If you wish to contribute please contact osagarra@ub.edu.

It may be expanded to cover additional cases.

Its bottelnecks may be also translated in fortran / C for faster convergence.

Finally, the cases of W and AW rely on a brute-force, massive optimization software for non-linear problems. 
For the moment it uses CVXOPT and Scipy but could be extended to other cases.

## Known issues 

In all cases, the convergence of the algorithm is not assured (we are dealing with a large optimization problem). 
It might or might not converge. In any case, all results need to be checked with simulations to see if the level of accuracy is according to the preference of the user.

Anyhow, if the algorithm does not converge, make sure the constraints you are imposing are feasible.

This package has been tested with networks up to N=10000 nodes.

## Additional notes

Check the DOCS/requirements to see the dependencies

This soft has been tested on Linux Ubuntu 12.10 and MacOS X 6.x
