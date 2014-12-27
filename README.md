Origin-Destination Multi-Edge Package.
========================================================================

 Copyright 2014 Oleguer Sagarra. All rights reserved. Code under License GPLv3.
______________________________________________________________________________________


## Introduction 

The Origin-Destination Multi-Edge network package (ODME) is the set of 
codes that I have written during the development of my PhD thesis. It 
is a package intended to be used for the analysis of an kind of 
multi-edge network, focusing specially on spatial networks.

A multi-edge network is a weighted network composed of 
distinguishable, integer weights. A paradigmatic example of it are 
Origin-Destination matrices used to represent flows between nodes or 
locations.

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

[3] To be added




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

We would like to thank Pol Colomer and Sergio Oller for their very useful comments and suggestions.



## License

Copyright 2014 Oleguer Sagarra.
All rights reserved. 
Code under License GPLv3.


## Additional notes


Check the DOCS/ folder in each module as well as their separate README files for details.
