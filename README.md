This is the package for split Knockoffs. The main content of this package is the following two functions:

1. split_knockoffs.filter(X, D, y, nu_s, option); the split Knockoff filter, details will be explained later
2. split_knockoffs.cv_filter(X, D, y, nu_s, option); the split Knockoff filter with cross validation, details will be explained later.

To reproduce the figures and tables in the paper, please go to "simulation" folder and "AD_experiments" folder and run the respective file. The recommended setting for the range of lambda is 10<sup>0</sup> - 10<sup>-6</sup> with step size 0.01. Before using this package, please add the folder to path by changing the "addpath('C:\Users\Yang\desktop\Split_Knockoffs')" in each document. This package is built based on [Knockoffs for matlab](https://web.stanford.edu/group/candes/knockoffs/software/knockoffs/).



Usage of functions:

[results, Z, t_Z] = split_knockoffs.filter(X, D, y, nu_s, option): split Knockoff filter for structural sparsity problem.

input argument
X : the design matrix
y : the response vector
D : the linear transform
nu_s: a set of nu, appointed by the user
nu: the parameter for variable splitting
option: options for creating the Knockoff statistics
	option.eta : the choice of eta for creating the knockoff copy
	option.q: the desired FDR control bound
	option.method: "knockoff" or "knockoff+"
	option.stage0: choose the method to conduct split knockoff
		"fixed": fixed intercept assignment for PATH ORDER method
			option.beta : the choice of fixed beta for step 0: 
				"mle": maximum likelihood estimator; 
				"ridge": ridge regression choice beta with lambda = 1/nu
				"cv_split": cross validation choice of split LASSO over nu and lambda
				"cv_ridge": cross validation choice of ridge regression over lambda
		"path": take the regularization path of split LASSO as the intercept assignment for PATH ORDER method
		"magnitude": using MAGNITUDE method
	option.lambda_s: the choice of lambda for path
	option.normalize: whether to normalize the data

output argument
results: a cell of length(nu_s) * 1 with the selected variable set in each cell w.r.t. nu.
Z: a cell of length(nu_s) * 1 with the feature significance Z in each cell w.r.t. nu.
t_Z: a cell of length(nu_s) * 1 with the knockoff significance tilde_Z in each cell w.r.t. nu.

[result, CV_loss, nu_optimal] = split_knockoffs.cv_filter(X, D, y, nu_s, option): split Knockoff filter for structural sparsity problem, using cross validation.

input argument
X : the design matrix
y : the response vector
D : the linear transform
nu_s: a set of nu, appointed by the user
nu: the parameter for variable splitting
option: options for creating the Knockoff statistics
	option.eta : the choice of eta for creating the knockoff copy
	option.q: the desired FDR control bound
	option.method: "knockoff" or "knockoff+"
	option.stage0: choose the method to conduct split knockoff
		"fixed": fixed intercept assignment for PATH ORDER method
			option.beta : the choice of fixed beta for step 0: 
				"mle": maximum likelihood estimator; 
				"ridge": ridge regression choice beta with lambda = 1/nu
				"cv_split": cross validation choice of split LASSO over nu and lambda
				"cv_ridge": cross validation choice of ridge regression over lambda
		"path": take the regularization path of split LASSO as the intercept assignment for PATH ORDER method
		"magnitude": using MAGNITUDE method
	option.lambda_s: the choice of lambda for path
	option.k_fold: the fold used in cross validation
	option.cv_rule: the rule used in CV
		"min": choose nu with minimal CV loss
		"complexity": choose nu with minimal model complexity in the range of 0.99 * CV_loss <= min(CV_loss)

output argument
result: CV optimal selection set
CV_loss: the CV loss w.r.t. nu.
nu_optimal: CV optimal nu.
