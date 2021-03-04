% add paths, can be changed to other directions
addpath(genpath('C:\Users\Yang\desktop\knockoffs_matlab'));
addpath(genpath("C:\Users\Yang\Documents\MATLAB\SLEP_package_4.1"));
% this simulation uses parallel computing
% parameter settings
n = 150;
p = 50;
k = 15; % sparsity level
A = 0.5; % magnitude
c = 0.5; % feature correlation
cv_on = 0; % whether to use CV for alternative step 2 knockoff 
% calculation of CV uses parallel computing
k_fold = 5; % cv folds
s_size = 1.5; % adjust s_size = 2 - eta
method = "knockoff"; % can chosen from knockoff and knockoff+
q = 0.2; % target FDR
sigma = 1; % noise level
tests = 10; % number of experiments

% settings for nu
nu_s = [0.1, 1, 10];
num_nu = length(nu_s);


% generate D
D_G = zeros(p-1, p);

for i = 1:(p-1)
    D_G(i, i) = 1;
    D_G(i, i+1) = -1;
end

% D = eye(p);
% D = [eye(p); D_G];
D = D_G; % various choice of D
m = size(D, 1);

% generate X
Sigma = zeros(p, p);
for i = 1: p
    for j = 1: p
        Sigma(i, j) = c^(abs(i - j));
    end
end

rng(100);
X = mvnrnd(zeros(p, 1), Sigma, n); % generate X


% generate beta and gamma
beta_true = zeros(p, 1);
for i = 1: k
    beta_true(i, 1) = A;
    if rem(i, 3) == 1
        beta_true(i, 1) = -A;
    end
end
gamma_true = D * beta_true;

% create matrices to store results
FDR_s_alter = zeros(tests, num_nu);
POWER_s_alter = zeros(tests, num_nu);
FDR_s_original = zeros(tests, num_nu);
POWER_s_original = zeros(tests, num_nu);
if cv_on == 1
    FDR_CV = zeros(tests, 1);
    POWER_CV = zeros(tests, 1);
    CV_list = zeros(num_nu, tests);
end
FDR_knockoff = zeros(tests, 1);
POWER_knockoff = zeros(tests, 1);

%%%%%%%%%%%%%%% begin simulation %%%%%%%%%%%%%

parfor test = 1: tests

    % generate varepsilon
    rng(test);
    
    % generate noise and y
    varepsilon = randn(n, 1) * sqrt(sigma);
    y = X * beta_true + varepsilon;
    
    % running CV for alternative step2 split knockoff
    if cv_on == 1
        [CV_list(:, test) , chosen_nu] = CV(X, D, y, nu_s, q, s_size, k_fold, method);
        [results, ~, ~] = split_knockoff_alter(X, D, y, chosen_nu, q, s_size, method);
        [POWER_CV(test, 1), FDR_CV(test, 1)] = Power_FDR_for_simulation(gamma_true, results);
    end
    
    % running knockoff as a comparison
    results = gen_knockoff(X, D, y, q, method);
    [POWER_knockoff(test, 1), FDR_knockoff(test, 1)] = Power_FDR_for_simulation(gamma_true, results);
    
    % runing original and alternative step 2 knockoff
    for i = 1: num_nu
        nu = nu_s(1, i);
        [results, ~, ~] = split_knockoff_alter(X, D, y, nu, q, s_size, method);
        [POWER_s_alter(test, i), FDR_s_alter(test, i)] = Power_FDR_for_simulation(gamma_true, results);
        
        [results, ~, ~] = split_knockoff_original(X, D, y, nu, q, s_size, method);
        [POWER_s_original(test, i), FDR_s_original(test, i)] = Power_FDR_for_simulation(gamma_true, results);
    end
end

% compute the means

mean_FDR_alter = mean(FDR_s_alter);
mean_POWER_alter = mean(POWER_s_alter);

mean_FDR_original = mean(FDR_s_original);
mean_POWER_original = mean(POWER_s_original);

mean_FDR_knockoff = mean(FDR_knockoff);
mean_POWER_knockoff = mean(POWER_knockoff);

if cv_on == 1
    FDR = mean(FDR_CV);
    POWER = mean(POWER_CV);
end