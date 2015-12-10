Origin-Destination Multi-Edge Package.
========================================================================

 Copyright 2014 Oleguer Sagarra. All rights reserved. Code under License GPLv3.
______________________________________________________________________________________


## Introduction 

The Origin-Destination Multi-Edge network package (ODME) is a set of 
codes that I have written during the development of my PhD thesis. It 
is a package intended to be used for the analysis of many kinds of 
networks, focusing specially on spatial, multi-edge networks.

A multi-edge network is a non binary network composed of 
distinguishable, integer weights. A paradigmatic example of it are 
Origin-Destination matrices used to represent flows between nodes or 
locations. The present package is able to analyze many types of any 
general form of non binary networks (being binary networks a particular case of them), 
be them spatial or not.

The package is structured in three separated modules which perform 
different operations. The first one is an analysis tool, while the 
second is used to solve the parametters needed for the generation of 
null models. The third one is used for the generation of such models.

The package may be used to study spatial networks or non-spatial ones, 
and the analysis tool can be used in any weighted network, as long as 
the weights are integer numbers.

## References 


[1] Statistical mechanics of multiedge networks
	Sagarra O., Pérez-Vicente C. and Díaz-Guilera, A.  Phys. Rev. E 88, 062806 (2013)
    [Physical Review E](http://pre.aps.org/abstract/PRE/v88/i6/e062806)

[2] The configuration multi-edge model: Assessing the effect of fixing node strengths on weighted network magnitudes
	Sagarra O., Font-Clos F., Pérez-Vicente C. and Díaz-Guilera, A.
	[Europhysics Letters](http://iopscience.iop.org/0295-5075/107/3/38002/article;jsessionid=D22CAAF312F43653DA0C1279853CBF0C.c3)

[3] Supersampling and network reconstruction of urban mobility.
	Sagarra, O., Szell, M., Santi, P., and Ratti, C. PLoS One 10, e0134508 (2015)
	[arXiv:1504.01939v1](http://arxiv.org/abs/1504.01939)
    [PLoS One 10, e0134508 (2015)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0134508)

[4] Role of adjacency-matrix degeneracy in maximum-entropy-weighted network models.
    Sagarra, O., Pérez-Vicente C. and Díaz-Guilera, A. Phys. Rev. E 92, 052816
    [arxiv:1509.01383](http://arxiv.org/abs/1509.01383)
    [Physical Review E 92, 052816 (2015)](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.92.052816)




## Contents of the package

This package is structured in three separated modules. Each module has 
its own README and documentation and the organization goes as follows:
- Module 1: Network analysis. Codes are written in pure C and the 
module constitutes a tool for the analysis of general purpose weighted networks (be them 
spatial or not), as long as their weights are integer numbers.
- Module 2: Model Fitting. Codes written in python (numpy and other 
packages needed). This module serves to fit the parametters used for 
the generation of maximum entropy null models.
- Module 3: Model generation. Codes written in C. This module 
generates expectations from different maximum entropy models, as well 
as two well-known models: Radiation and sequential gravity model.




## Acknowledgements

I would like to thank Pol Colomer and Sergio Oller for their very useful comments and suggestions.



## License

(C) Copyright 2014 Oleguer Sagarra.
All rights reserved. 
Code under License GPLv3.

Each file in this folder is part of the ODME package. This code has no warranty whatsoever nor any kind of support is provided. 
You are free to do what you like with this code as long as you leave this copyright in place. Please cite us if you use our software.

ODME is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

ODME is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the package. If not, see http://www.gnu.org/licenses/.


## Additional notes


Check the DOCS/ folder in each module as well as their separate README files for details.

