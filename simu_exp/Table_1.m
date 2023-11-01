% This script reproduces Table 1 of our paper. The performance of Split
% Knockoffs can be found in '*_sk' variables, while that of Knockoffs be in
% '*_knockoff' variables.

%% choose parameters

k = 20; % the sparsity level
A = 0.25; % the signal noise ratio
n = 500;% the sample size
p = 100;% the number of features
c = 0.5; % the feature correlation

option = struct;
option.q = 0.2;
option.beta = 'cv_all';
option.frac = 2/5;
option.tests = 200; % the number of simulation instances
option.normalize = true;
option.lambda = 10.^[0: -0.01: -6];

% set random seed
option.rng = 1;

%% generate 3 types of transformation

% generate D1, D2, and D3
D_G = zeros(p-1, p);

for i = 1:(p-1)
    D_G(i, i) = 1;
    D_G(i, i+1) = -1;
end

D_1 = eye(p);
D_2 = D_G;
D_3 = [eye(p); D_G];
D_s = {D_1, D_2, D_3};

fdr_knock = zeros(2, option.tests, 2);
power_knock = zeros(2, option.tests, 2);
rawvalue = cell(3, 1);


%% start the simulation experiments
for i = 1: 3
    % choose the respective D for each example
    fprintf('Running simulation experiments for D_%d.\n', i);
    D = D_s{i};
    simu_data = split_knockoffs.private.simu_unit_cv(n, p, D, A, c, k, option);
    rawvalue{i, 1} = simu_data.raw;
    if i <= 2
        fdr_knock(i,  :, :) = simu_data.fdr_k;
        power_knock(i,  :, :) = simu_data.power_k;
    end
    save(sprintf('%s/result/temp',pwd));
end


%% start calculating results on 3 types of W statistics

beta_true = zeros(p, 1);
for i = 1: k
    beta_true(i, 1) = A;
    if rem(i, 3) == 1
        beta_true(i, 1) = 0;
    end
end

fdr_sk = zeros(3, option.tests, 3, 2);
power_sk = zeros(3, option.tests, 3, 2);
for di = 1: 3
    D = D_s{di};
    gamma_true = D * beta_true;
    dataframe = rawvalue{di};
    for i = 1: option.tests
        stats = dataframe{i};
        Z = stats.Z;
        t_Z = stats.t_Z;
        r = stats.r;
        t_r = stats.t_r;
        Z_tilde = t_Z .* (real(r) == real(t_r));
        
        % Ws
        W = Z .* sign(Z - t_Z);
        S = knockoffs.select(W, option.q, 'knockoff');
        [fdr_sk(di, i, 1, 1), power_sk(di, i, 1, 1)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
        S = knockoffs.select(W, option.q, 'knockoff+');
        [fdr_sk(di, i, 1, 2), power_sk(di, i, 1, 2)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
        
        % Wst
        W = Z .* sign(Z - Z_tilde);
        S = knockoffs.select(W, option.q, 'knockoff');
        [fdr_sk(di, i, 2, 1), power_sk(di, i, 2, 1)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
        S = knockoffs.select(W, option.q, 'knockoff+');
        [fdr_sk(di, i, 2, 2), power_sk(di, i, 2, 2)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
        
        % Wbc
        W = max(Z, t_Z) .* sign(Z - t_Z);
        S = knockoffs.select(W, option.q, 'knockoff');
        [fdr_sk(di, i, 3, 1), power_sk(di, i, 3, 1)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
        S = knockoffs.select(W, option.q, 'knockoff+');
        [fdr_sk(di, i, 3, 2), power_sk(di, i, 3, 2)] ...
            = split_knockoffs.private.simu_eval(gamma_true, S);
    end
end

mean_fdr_knockoff = mean(fdr_knock, 2);
mean_power_knockoff = mean(power_knock, 2);
sd_fdr_knockoff = std(fdr_knock, 0,  2);
sd_power_knockoff = std(power_knock, 0,  2);

mean_fdr_sk = reshape(mean(fdr_sk, 2), [3, 3, 2]);
mean_power_sk = reshape(mean(power_sk, 2), [3, 3, 2]);
sd_fdr_sk = reshape(std(fdr_sk, 0,  2), [3, 3, 2]);
sd_power_sk = reshape(std(power_sk, 0,  2), [3, 3, 2]);

save(sprintf('%s/result/Table_1', pwd));

clearvars -except mean_fdr* mean_power* sd_fdr* sd_power*
