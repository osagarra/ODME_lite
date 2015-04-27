Origin-Destination Multi-Edge Package. Module 3: Multi Edge Network generator
========================================================================

 Copyright 2014 Oleguer Sagarra. All rights reserved. Code under License GPLv3.
______________________________________________________________________________________

The present program may be used to generate maximally random multi-edge networks with prescribed properties using different ensembles: Canonical (if available) or Grand-Canonical. It can also be used to generate ME networks according to the Radiation model and Seq. Gravity model proposed by other authors.
Multi-edge networks are weighted networks composed by distinguishable weights (weights which may be univocally identified such as time-stamped interactions or users riding on transportation systems).

Generating randomized instances of a network keeping some properties of its original network can be useful for studying the effect of a certain feature on any topological or dynamical property of such network.
Using different ensembles may be useful for performance and coding issues: Rewiring networks preserving some types of constraints may be a difficult endeavour both coputationally and algoritmically.
Comparing real networks with their randomized counterparts is very useful for data analysis and detection of statistically significant patterns. With this software you migh do just that, but if you use our software please do cite us.

In case you are interested in randomizing binary networks in a variety of situations take a look at [RandNetGen](https://github.com/polcolomer/RandNetGen) by Pol Colomer-de-Simón.

## References 

[1] Statistical mechanics of multiedge networks
	Sagarra O., Pérez-Vicente C. and Díaz-Guilera, A.  Phys. Rev. E 88, 062806 (2013)
    [Physical Review E](http://pre.aps.org/abstract/PRE/v88/i6/e062806)
    
[2] Supersampling and network reconstruction of urban mobility.
	Sagarra, O., Szell, M., Santi, P., and Ratti, C. Arxiv Preprint (2015)
	[arXiv:1504.01939v1](http://arxiv.org/abs/1504.01939)

External references:

[3] A universal model of commuting networks
	Lenormand, M., Huet, S., Gargiulo, F., & Deffuant, G. (2012).  PloS One, 7(10), e45985.

[4] A universal model for mobility and migration patterns
	Simini, F., González, M. C., Maritan, A., & Barabási, A.-L. (2012).  Nature, 484(7392), 96–100

[5] The configuration multi-edge model: Assessing the effect of fixing node strengths on weighted network magnitudes
    Sagarra O., Font-Clos F., Pérez-Vicente C. J. and Díaz-Guilera, A. EPL (Europhysics Lett. 107, 38002 (2014).

## Requirements and Installation

  You only need to install the GSL library: http://www.gnu.org/software/gsl/
  
  How to tutorial: http://www.brianomeara.info/tutorials/brownie/gsl

## Compilation

  Simply run:

    $ make 

  A single executable is created *MultiEdgeGen*.

## Execution

The executables is called **MultiEdgeGen** and has a set of options. To introduce on option you must introduce first its option key (preceded by a dash), one white space and then the argument value.
Arguments can appear in any order. If an argument does not appear the program gets the default value except for the compulsory items.

Examples
______________

 Other models: Simulate Multinomial Radiation model computing average sorensen index, providing cost_matrix, original network and original strength sequence.
 
    $ ./MultiEdgeGen -N 1000 -d 1 -f ../tests/strengths.str -C 0 -e 0 -g 0 -z 1 -a ../tests/cost_matrix.dists -m 25000 -v 1 -S 1 -Q ../tests/sample.tr
 
 Custom pij matrix: Simulte average values for a network belonging to the ensemble with provided probability matrix pij and a sampling T = 10000

    $ ./MutliEdgeGen -N 1000 -d 1 -f ../tests/sample.tr -C 1 -e 2 -T 10000

 Fixed strength: Simulte average values for a network belonging to the ensemble with fixed strength nodes ignoring the first line of the files for the header with provided lagrange multiplier list xxx.xy.

    $ ./MultiEdgeGen -N 1000 -d 1 -f ../tests/lagrange_mult/fixeds.xy -C 2 -e 2 -g 0 -h 1

 Fixed strength and cost: Simulte average values for a network belonging to the ensemble with fixed strength and fixed average cost, *gamma value = 0.00023* and providing a cost_matrix.

    $ ./MultiEdgeGen -N 1000 -d 1 -f ../tests/lagrange_mult/fixeds.xy -C 2 -e 2 -g 0.000023 -z 1 -a ../tests/cost_matrix.dists

 
 Fixed strength and total number of binary edges: Simulte average values for a network belonging to the ensemble with fixed strength and total number of binary edges and a provided cost matrix *cost_matrix.dists*, providing also the original network (*sample.tr*) to compare sorensen values at each run, not accepting self loops and *gamma value =0.77032662076117631*.

    $ ./MultiEdgeGen -N 1000 -d 1 -f ../tests/lagrange_mult/fixedE_noself.xy -C 3 -e 2 -g 0.77032662076117631 -z 1 -a ../tests/cost_matrix.dists -l 0 -S 1 -Q ../tests/sample.tr

 Fixed strength and degree sequence: Simulte average values for a network belonging to the ensemble with fixed strength and degree sequence, with verbose and 10 reps for averaging.

    $ ./MultiEdgeGen -N 4085 -d 1 -f ../tests/lagrange_mult/fixedsk.xyzw -C 4 -e 2 -v 1 -r 10

 Fixed degree: Simulte average values for a network belonging to the ensemble of fixed degrees and multi-edge network with *T=200000* events.
 	
    $ ./MultiEdgeGen -N 1000 -d 1 -f ../tests/lagrange_mult/fixedk.xy -C 5 -e 2 -T 200000


## Options and Arguments

All the programs share a general set of arguments (some might not include all of them).

#### Options

```
    Compulsory:    
		-N      	Number of nodes [int]
		-d			Directed (1) or undirected (0) option [int]
		-f			Path to file with appropiate Lagrange Multipliers in each case or srtength sequence (see README) [array of floats]
        -C          Case to apply to the model, see README for more info and description of each case [int]
        -g          Gamma value. Deferrence value for distances or number of edges (according to each case) 	
    Optional (depending on case)
        -e          Ensemble method to use in each case [int]
        -T          Number of trips (if needed) [int]
		-z			Distance Option: Include distances in analysis? [int] (default=0 no, 1 for yes)
  	        -a      Path to cost matrix list in format node_i node_j d_ij [must be triangular matrix symmetric]
		    -L		Log-distance option (to compute the logarithm of the cost matrix) [int] (default=0 np, >0 for yes) 
            -m		Maximum distance for binning [float] (default= 20000) [in meters]"
        -S          Sorensen option: Compute sorensen with given original network as benchmark [int] (default=0)
            -Q      Path to original network to compare to, if wanted in edge list format: node_origin node_dest t_ij (int int int)
 		-l 			Self-loop option (>0 for accepting them) (default =1)
 		-c 			Clustering option (1 for yes) (warning: Depending on av_s makes simulations orders of magnitude slower) [int] (default=0)
 		-v 			Verbose (1 for on, 0 for off) [int] (default 0)
        -r          Number of reps for averaging (default=100)
        -p          Printing option of sample adjacency list (default=0, >0 for yes)\n"
		-s 			Initial seed for random generator [int] (default=1)
 		-x 			Exponent for log-binning on network stats (-1 for no log binning) [int] (default=-1)
        -h 			Number of header lines on file_s [int] (default=1)
        -R          Number of maximum trials reps for binary linking (default=100) [only for binary cases, see README]\n"
```

## Brief description of the Multi-Edge Network Generator

This program generates ensemble expectations of various properties for networks generated according to a variety of models. These models include the well known Radiation Model, Sequential Gravity, doubly constrained gravity model or weighted configuration model. In a nuthsell, it implements many of the cases presented in reference [1]. The program is able to provide single graph realizations as well as averages over the ensemble of considered graphs for various network magnitudes.


## Description of different cases
This part of the documentation describes the particularities that each case has (the numeration corresponds to the option -C in the code). 

#### Case 0: Other models
This case refers to the implementation of both the Radiation model (ref [3]) and the Sequential gravity model (ref [4]). In this case, 2 versions are implemented for each of the models (which are asymptotically equivalent), which may be specified with the flag *-e*. The only compulsory inputs for these models are the incoming and ougtoing strength sequences as well as the cost matrix considered. For the sequential gravity model, if the parametter *gamma* (*\beta* in the original paper) is not given, then it is approximated from the data (assuming the cost_matrix is provided in kilometers). Available methods (*-e* flag):

*Radiation Model*


- 0: Multinomial Radiation model. Fast version of the algorithm, which uses a multinomial assignation of events from source.
- 1: Stochastic Radiation model. Slow version of the algorithm, implementing the micro-simulation of job assignation per each traveller or trip (event).

The expectation in this case for the average number of trips between two nodes is,

	<t_{ij}> = s^{out}_i  s^{in}_i/(1 - s^{in}_i/T) s^{in}_j / [ (s_ij + s^{in}_i)(s_ij + s^{in}_i + s^{in}_j)  

Note that this implementation is based on the network formula adapted of the model (see equation 5 in [Masucci et atr. Arxiv](http://arxiv.org/abs/1206.5735)), which differs from the one in reference [3]. You may adapt the code changing the code (see apropriate lines in file *other_null_models.c**).

*Seq. Gravity Model*

- 2: Multinomial Sequential Gravity model. Fast multinomial implementation of the algorithm by Lenormand et altr.
- 3: Bernouilli Sequential Gravity model. Rejection based implementation of the model.

This model has no analytical expectation for the total number of trips.

#### Case 1: Custom matrix *p_ij* given
This case simply implements the situation described in ref [2]. The flag -T fixing the total number of events needs to be given as well as the matrix of probabilities for the allocation of events (does not need to be normalized). In this case, the expectation for the events between nodes ij is simply,

    <t_{ij}> = T p_{ij}

and the events can be distributed according either to a set of independent Poisson processes with mean <t_{ij}> (grand canonical ensemble, -e flag set to 2) or a multinomial approach with fixed probabilities (canonical ensemble) with a fixed total sampling T (flag -e set to 0) or varying sampling according to a Poisson process of mean <T> (alternative grand-canonical case, flag -e set to 1).

#### Case 2: Fixed cost and strength sequence
This case implements the situation where the strength as well as the total cost is fixed. In this case, it needs as input a lagrange multiplier sequence (2 lagrange multipliers per node in the directed case, 1 for the undirected case). The expectation for the average events between nodes is,

    <t_{ij}> = x_i y_j exp(-gamma c_{ij}).

Note that the case *gamma=0* is the weigthed configuration model (ref [5]). Also note that depending on the definition of costs, one recovers the doubly constrained gravity model in its exponential form or power-law deference form for distances. The events can be distributed according either to a set of independent Poisson processes with mean <t_{ij}> (grand canonical ensemble, -e flag set to 2) or a multinomial approach with fixed probabilities (canonical ensemble) with a fixed total sampling T (flag -e set to 0) or varying sampling according to a Poisson process of mean <T> (alternative grand-canonical case, flag -e set to 1).

#### Case 3: Fixed strength sequence and total number of binary edges
This case fixes the strength sequence of each node and also the total number of binary edges. The inputs are the *gamma* parameter, and a lagrange multiplier sequence for each node. The expectation for the average events between nodes is,

    <t_{ij}> = gamma x_i y_j exp(x_i y_j) / (gamma(exp(x_i y_j)-1)+1).

The only method implemented is the grand-canonical case (-e:2), which generates a set of independent ZIP processes.

#### Case 4: Fixed strength and degree sequence
This case fixes the strength sequence of each node and its degree also. The inputs is a lagrange multiplier sequence for each node (4 values in the directed case, 2 values for the undirected one). The expectation for the average events between nodes is,

    <t_{ij}> = z_i w_j  x_i y_j exp(x_i y_j) / (z_i w_j (exp(x_i y_j)-1)+1).

The only method implemented is the grand-canonical case (-e:2), which generates a set of independent ZIP processes.

#### Case 5: Fixed degree sequence
This case fixes only the degree sequence of each node. In this case the input is the total number of events T as well as a lagrange multiplier sequence for each node (2 values in the directed case, 1 values for the undirected one). The expectation for the average events between nodes is,

    <t_{ij}> = T/<E> z_i w_j / (z_i w_j +1)
    <E>     = \sum_{ij} z_i w_j / (z_i w_j +1)

Note that T must be bigger than <E>. The only method implemented is the grand-canonical case (-e:2), which generates a set of independent ZIP processes.

## Outputs
For both directed and undirected cases:

The generated outputs are files with extension *.hist* if they correspond to network statistics (binned) and *.list* for the node list with different features.


The files with the suffix *ens* reffer to averages over ensembles. The ones without the *ens* suffix refer to a single realization.


**.list files**
   The file *node_list.list* contains different node features.
   For the undirected case: (Note that the clustering only appears if -c flag is set)
   

```
Node_num k k_analitic s Y2 k_nn k^w_nn s^w_nn (optionally  clust c^w)
```
Explicit formulas: (in latex)
```
	k_i = \sum \Theta(t_ij)
	k^{anal}_i = \sum (1-e^{-s_is_j/T})
	s_i = \sum t_{ij}
	Y2 = \sum t_ij^2 / s_i
	k_nn_i = \sum \Theta(t_{ij}) k_j / k_i
	k^w_nn_i = \sum t_{ij} k_j /s_i
	s^w_nn = \sum t_{ij} s_j / s_i
	c = \sum_{jk} \Theta(t_ij)\Thetat(t_jk)\Theta(t_ki) / k_i(k_i-1)
	c^w = \sum_{jk} (t_ij + t_ik) \Theta(t_ij) \Thetat(t_jk)\Theta(t_ki) / [2s_i(k_i-1)]
```

For the directed case: (In this case the clustering is not defined)

```
Node_num k k_analitic s Y2 k_nn k^w_nn s^w_nn 
```
First for the out case, then for the in case. Note that the averages over Y2,s^w_nn,k^w_nn and k_nn are only defined if s>0 (otherwise they are 0).

Explicit formulas (in latex):
```
Out:
	k_i^{out} = \sum_j \Theta(t_ij)
	k^{anal}_i^{out} = \sum_j (1-e^{-s_i^out}s_j^{in}/T})
	s_i^{out} = \sum_j t_{ij}, 
	Y2^{out} = \sum_j t_ij^2 / s_i^{out}
	k_nn_i^{out} = \sum_j \Theta(t_{ij}) k_j^{in} / k_i^{out}
	k^w_nn_i = \sum_j t_{ij} k_j^{in} /s_i^{out}, 
	s^w_nn = \sum_j t_{ij} s_j^{in} / s_i^{out}
In:
	k_j^{in} = \sum_j \Theta(t_ij)
	k^{anal}_j = \sum_j (1-e^{-s_i^out s_j^in/T})
	s_j^in = \sum_i t_{ij}
	Y2^{in} = \sum_i t_ij^2 / s_j^{in}
	k_nn^{in} = \sum_i \Theta(t_{ij}) k_i^{out} / k_j^{in}
	k^w_nn = \sum_i t_{ij} k_i^{out} /s_j^{in}, 
	s^w_nn^{in} = \sum_i t_{ij} s_i^{out} / s_j^{in}
```

**ens_XX.list files**
   The file *ens_XXnode_list.list* contains different node features averaged over the ensemble for XX reps.
   For the undirected case: (Note that the clustering only appears if -c flag is set)
   

```
Node_num <k> std_k <s> std_s <Y2> std_Y2 <k_nn> std_k_nn <k^w_nn> std_k^w_nn <s^w_nn> std_s^w_nn (optionally  <clust> std_clust <c^w> std_c^w)
```
Explicit formulas: (in latex)
```
	k_i = \sum \Theta(t_ij)
	k^{anal}_i = \sum (1-e^{-s_is_j/T})
	s_i = \sum t_{ij}
	Y2 = \sum t_ij^2 / s_i
	k_nn_i = \sum \Theta(t_{ij}) k_j / k_i
	k^w_nn_i = \sum t_{ij} k_j /s_i
	s^w_nn = \sum t_{ij} s_j / s_i
	c = \sum_{jk} \Theta(t_ij)\Thetat(t_jk)\Theta(t_ki) / k_i(k_i-1)
	c^w = \sum_{jk} (t_ij + t_ik) \Theta(t_ij) \Thetat(t_jk)\Theta(t_ki) / [2s_i(k_i-1)]
```

For the directed case: (In this case the clustering is not defined)

```
Node_num <k> std_k <s> std_s <Y2> std_Y2 <k_nn> std_k_nn <k^w_nn> std_k^w_nn <s^w_nn> std_s^w_nn
```
First for the *out* case, then for the *in* case. Note that the averages over Y2,s^w_nn,k^w_nn and k_nn are only taken over realizations where the nodes have strength >0 (due to the not definition of the case 0/0).


**w_s_io.hist** files

This files contain the graph-average existing occupation number as function of the product of strenghts in the format:
```
x^{out}_i*x^{in}_j \bar{t_ij} \std{\bar{t_ij}}
```
Where x stands for degree and strengths respectively and the averages $\bar{x}$ are taken over a single graph realization (and the standard deviations also). The cases *w_s_i.hist* and *w_s_o.hist* only present the average weight of links according to the strength of node origin or destination respectively.

**_w.hist files**
This files contain the distribution of existing occupation numbers in the format:
```
bin_id bin_min bin_max Bin_count Bin_std CCDF
```
The histogram is not normalized.


**entropies.hist files**
If more than 10 repetitions of the ensemble are asked, the entropy histogram for the different realizations is also computed. Note that the entropy definition for each ensemble is different, check reference [1] in case of doubt.

#### Additional files
If the flag *-z* is set to *>0*, then also some distance related quantities are computed.

**d_edges.hist and d_trips.hist and p_ij.hist files**

Distribution of binary edge costs costs and trip costs. Probability of binary connection as a function of cost between nodes.

**w_dij.hist**

Graph average value of existing weights as a function of cost.

**Example**:

The following command

```
time ./net_analysis -N 94 -d 1 -f tests/sample.tr -z 0 -a tests/cost_matrix.dists  -h 1
```

takes on a intel® Core™ i5-2500 CPU @ 3.30GHz × 4 with 12 GM RAM the following time

```
real	0m0.025s
user	0m0.021s
sys	0m0.003s
```

on a and generates the following files:

```
N94avs37746.57447_d_edges.hist		- Binary edge cost histogram
N94avs37746.57447_d_trips.hist		- Event cost histogram
N94avs37746.57447_pij.hist			- Binary connection probability
N94avs37746.57447_w_dij.hist		- Average weight as function of cost
N94avs37746.57447_w_s_i.hist		- Average Weight as a function of incoming strength
N94avs37746.57447_w_s_o.hist		- Average Weight as a function of outgoing strength
N94avs37746.57447_w_s_oi.hist		- Average Weight as a function of product of strengths
N94avs37746.57447node_list.list		- List of node attributes
N94avs37746.57447_w.hist			- Weight distribution P(w)
```

The undirected version of the algorithm produces equivalent files for the *in* and *out* values.


## Acknowledgements

I would like to thank Pol Colomer and Sergio Oller for their very useful comments and suggestions. I would also like to thank Maxime Lenormand for making his original code accessible for comparison purposes.



## License

Copyright 2014 Oleguer Sagarra.
All rights reserved. 
Code under License GPLv3.

## Roadmap

The software is for the moment ready and works smoothely under apropiate conditions.

It could be extended to include indistinguishable cases.

Also the Radiation model with selection could be included.

If you wish to contribute please contact osagarra@ub.edu.



## Known issues 

When the Log option is activated, special care in the maximum distance parametter need to be taken, or the bins will not be correctly selected.
The enumeration of nodes between distance and adjacency list needs to be one to one. Additionally, the counting of nodes needs to start at 0 (the id of the first node to appear in the adjacency list and distance list needs to be 0).

Some memory issues may arise if used on very big networks as the scaling with the total number of nodes N goes as N^2.

## Additional notes


Check the DOCS/requirements to see the dependencies

This soft has been tested on Linux Ubuntu 12.10 and MacOS X 6.x

Recommended compiler for mac is Clang while icc for Linux.
