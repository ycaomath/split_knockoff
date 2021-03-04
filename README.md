# split_knockoff

For simulations, please go to simulation folder and run simulation.m

The main functions for split knockoff is:

[result, Z, t_Z] = split_knockoff_original(X, D, y, nu, q, s_size, method)
This function runs the original definition of step 2. Here s_size = 2 - eta controls the size of S, method = "knockoff" or "knockoff +". Returns the selected sets, Z and tilde_Z.

[result, Z, t_Z] = split_knockoff_alter(X, D, y, nu, q, s_size, method)
This function runs the alternative definition of step 2.

[CV_error, nu] = CV(X, D, y, nu_s, q, s_size, k_fold, method)
This function runs the CV for split knockoff. Return chosen nu and respective CV loss for each nu. k_fold is to choose how many folds that CV would use.

[result] = gen_knockoff(X, D, y, q, method)
This function returns the selected set from knockoff package. This function uses U X pinv(D) as the design matrix, where U is the orthogonal complement.
