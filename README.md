# Split Knockoffs

This is a Matlab package to reproduce the experiments in our paper.



## Dependency

- [Glmnet in Matlab](https://web.stanford.edu/~hastie/glmnet_matlab/)



## Installation Details

This package is tested in Matlab R2020b on Windows 10 with [the updated version of Glmnet](https://web.stanford.edu/~hastie/glmnet_matlab/glmnet_matlab_new.zip). A fortran compiler is required for using Glmnet in Matlab. We tested the [Intel Fortran Compiler](https://software.intel.com/content/www/us/en/develop/articles/oneapi-standalone-components.html#fortran) (Version 2021.2.0) with [Visual Studio 2019](https://visualstudio.microsoft.com/). 



## Usage

To finish the installation, please run 'startup.m'. For usage of this package, type 'help split_knockoffs.filter' and 'help split_knockoffs.cv_filter' in the command line of Matlab.



## Acknowledgement

This package adapts the '+knockoffs' folder from [Knockoffs for matlab](https://web.stanford.edu/group/candes/knockoffs/software/knockoffs/) for common functions and comparisons.



## Reproducible Experiments

- To reproduce the results of simulation experiments in our paper, please check the respective scripts in the 'simu_exp' folder for details.
- To reproduce the results of experiments on Alzheimer‘s Disease in our paper, please check the respective scripts in the 'AD_exp' folder for details.
