function simu_data = simu_unit(n, p, D, A, c, k, option)
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
% simu_data: a structure contains simulation results, especially Z, t_Z

sigma = 1; % the noise level
tests = option.tests; % the number of experiments
num_nu = length(option.nu);

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
fdr_split = zeros(tests, num_nu, 2);
power_split = zeros(tests, num_nu, 2);
raw = cell(tests, 1);

%%%%%%%%%%%%%%% begin simulation %%%%%%%%%%%%%

for test = 1: tests

    % generate varepsilon
    rng(test);
    
    % generate noise and y
    varepsilon = randn(n, 1) * sqrt(sigma);
    y = X * beta_true + varepsilon;
    
    [results, raw{test}] = split_knockoffs.filter(X, D, y, option);

    for i = 1: num_nu
        result = results{i}.sk;
        [fdr_split(test, i, 1), power_split(test, i, 1)] ...
            = split_knockoffs.private.simu_eval(gamma_true, result);
        result = results{i}.sk_plus;
        [fdr_split(test, i, 2), power_split(test, i, 2)] ...
            = split_knockoffs.private.simu_eval(gamma_true, result);
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
simu_data.raw = raw;

end