function simu_data = simu_unit_cv(n, p, D, A, c, k, option)
% the simulation unit for simulation experiments.
%
% input arguments
% n: the sample size
% p: the dimension of variables
% D: the linear transform
% A: the SNR for gaussian noise
% c: the feature correlation
% k: the number of nonnulls in beta
% option: options for split knockoffs
%
% output arguments
% simu_data: a structure contains simulation results

sigma = 1; % noise level
tests = option.tests; % number of experiments

% generate X
Sigma = zeros(p, p);
for i = 1: p
    for j = 1: p
        Sigma(i, j) = c^(abs(i - j));
    end
end

rng(100);
X = mvnrnd(zeros(p, 1), Sigma, n); % generate X
m = size(D, 1);


% generate beta and gamma
beta_true = zeros(p, 1);
for i = 1: k
    beta_true(i, 1) = A;
    if rem(i, 3) == 1
        beta_true(i, 1) = 0;
    end
end
gamma_true = D * beta_true;

% create matrices to store results
fdr_split = zeros(tests, 2);
power_split = zeros(tests, 2);
raw = cell(tests, 1);
fdr_k = zeros(tests, 2);
power_k = zeros(tests, 2);

%%%%%%%%%%%%%%% begin simulation %%%%%%%%%%%%%

for test = 1: tests

    % generate varepsilon
    rng(test);
    varepsilon = randn(n, 1) * sqrt(sigma);

    y = X * beta_true + varepsilon;
    
    [results, raw{test}] = split_knockoffs.filter(X, D, y, option);

    result = results.sk;
    [fdr_split(test, 1), power_split(test, 1)] ...
        = split_knockoffs.private.simu_eval(gamma_true, result);
    result = results.sk_plus;
    [fdr_split(test, 2), power_split(test, 2)] ...
        = split_knockoffs.private.simu_eval(gamma_true, result);
    if m <= p && n >= p
        result = split_knockoffs.private.convert_knockoff(X, D, y, option);
        
        [fdr_k(test, 1), power_k(test, 1)] ...
            = split_knockoffs.private.simu_eval(gamma_true, result.k);
        [fdr_k(test, 2), power_k(test, 2)] ...
            = split_knockoffs.private.simu_eval(gamma_true, result.k_plus);
    end
    fprintf('%.1f%%...',test/tests*100);
    if test == tests
        fprintf('done.\n');
    end
    if mod(test, 10) == 0 && test ~= tests
        fprintf('\n');
    end
end

% output the results

simu_data = struct;


simu_data.fdr_split = fdr_split;
simu_data.power_split = power_split;
simu_data.fdr_k = fdr_k;
simu_data.power_k = power_k;
simu_data.raw = raw;

end