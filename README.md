# The 'onlineChange' R package
The onlineChange R package is designed to quickest detect any change in distributions for online streaming data, supporting any user-specified distributions. It supports sequential Monte Carlo to perform Bayesian change point analysis, where the parameters before or after the change can be unknown. It also supports likelihood ratio based test to determine the stopping time (i.e. when change occurs).
This is still an ongoing project by a group of researchers from the University of Minnesota. The published software here is an initial version that is ready to use.
For questions and references, please contact Jie Ding at dingj@umn.edu

## Getting Started

First install the devtools package

install.packages("devtools")

library("devtools")

Then install this package

install_github('JieGroup/onlineChange')

## Using This Package

To see the available function to use, type 

ls("package:onlineChange")

A quick guide of package can be found [here](https://github.com/JieGroup/onlineChange/blob/master/vignettes/user-guide.pdf) 

## Acknowledgment

This research is funded by the Defense Advanced Research Projects Agency (DARPA) under grant number HR00111890040.

Part of this package is based on the Restricted Boltzmann Machine project at https://github.com/TimoMatzen/RBM
