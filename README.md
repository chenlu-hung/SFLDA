# SFLDA
The files demonstrate the algorithm in our paper [Sensible functional linear discrminant analysis](https://arxiv.org/abs/1606.03844)

## Installation

1. Download and install [PACE](http://www.stat.ucdavis.edu/PACE/)

2. Edit *example_xxx.m*: uncomment the 3rd line and replace *path/to/PACE* with your PACE directory.

3. Run *example_xxx.m*

## PACE 
1. PACE 2.16 (http://www.stat.ucdavis.edu/PACE/) included in this computer code is recommended, where the file cvfda_lwls.m has been slightly modified for functional data (‘regular=2’). 

2. If you use PACE 2.16 from its website, you might encounter some trouble with the bandwidth selection for mean function if the leave-one-subject-out CV is selected ('bwmu_gcv=0'). 

3. If you use PACE 2.17, in addition to the problem of cvfda_lwls.m, you might encounter another issue: mysample.m will reset the random seed to rng(123,’twister’) after estimating the covariance function if 'bwxcov_gcv=0' is used. Thus, it might cause problem when you perform simulation studies.