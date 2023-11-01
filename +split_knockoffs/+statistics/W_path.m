function[W, stats] = W_path(X, D, y, nu, option)
% split_knockoffs.statistics.pathorder.W_path generates the split knockoff
% statistics W based on the beta(lambda) from a split LASSO path in the
% intercepetion assignment step.
%
% input arguments:
% X : the design matrix
% y : the response vector
% D : the linear transformation
% nu: the parameter for variable splitting
% option: options for creating the Knockoff statistics
%
% output arguments
% W: the split knockoff statistics
% stats: the intermedia statistics

[n, ~] = size(X);

% split the dataset
rand_rank = option.rand_rank;

ind1 = rand_rank(1: floor(n * option.frac));
ind2 = rand_rank(floor(n * option.frac+1): end);
X_1 = X(ind1, :);
y_1 = y(ind1, :);
X_2 = X(ind2, :);
y_2 = y(ind2, :);


% %%%%%%%%%%%%% step 0 %%%%%%%%%%%%%%

[m, p] = size(D);

% generate the design matrix

opts.copy = false;
[A_beta,A_gamma,tilde_y,~] = split_knockoffs.create(X_1, y_1, D, nu, opts);

% set lambda
if isfield(option, 'lambda') == true
    lambda_vec = option.lambda;
else
    lambda_vec = 10.^[0:-0.01:-6];
end
nlambda = length(lambda_vec);

% set the penalty
penalty = ones(m+p,1);
for i = 1: p
    penalty(i, 1) = 0;
end

% lasso path settings for glmnet
opts = struct; 
opts.lambda = lambda_vec;
opts = glmnetSet(opts);
opts.penalty_factor = penalty;
fit_step0 = glmnet([A_beta, A_gamma], tilde_y, [], opts);


coefs = fit_step0.beta;

% store beta(lambda)
betas = coefs(1: p, :);

%%%%%%%%%%%%% step 1 %%%%%%%%%%%%%%

% generate the design matrix

opts.copy = true;
[A_beta,A_gamma,tilde_y,tilde_A_gamma] = split_knockoffs.create(X_2, y_2, D, nu, opts);

coef1 = zeros(m, nlambda);
parfor i = 1: nlambda
    % take beta_lambda, gamma_lambda as calculated in step 1
    y_new = tilde_y - A_beta * betas(:, i);
    % calculate LASSO
    opts = struct; 
    opts.lambda = lambda_vec(i);
    opts = glmnetSet(opts);
    fit_step1 = glmnet(A_gamma, y_new, [], opts);
    coef1(:, i) = fit_step1.beta;
end  

% calculate r and Z
r = zeros(m, 1);
Z = zeros(m, 1);
for i = 1: m
    [Z(i), r(i)] = split_knockoffs.private.hittingpoint(coef1(i, :), lambda_vec);
end

% deal with the case where some features have alrealdy been screened off
if option.m > m
    temp = zeros(option.m, 1);
    temp(option.gamma_supp) = Z;
    Z = temp;

    temp = zeros(option.m, 1);
    temp(option.gamma_supp) = r;
    r = temp;
end

%%%%%%%%%%%%% step 2 %%%%%%%%%%%%%% 
coef2 = zeros(m, nlambda);
parfor i = 1: nlambda
    % take beta_lambda, gamma_lambda as calculated in step 1
    y_new = tilde_y - A_beta * betas(:, i);
    % calculate LASSO
    opts = struct; 
    opts.lambda = lambda_vec(i);
    opts = glmnetSet(opts);
    fit_step2 = glmnet(tilde_A_gamma, y_new, [], opts);
    coef2(:, i) = fit_step2.beta;
end  
% calculate tilde_Z tilde_r and W
t_Z = zeros(m, 1);
t_r = zeros(m, 1);

for i = 1: m
    [t_Z(i), t_r(i)] = split_knockoffs.private.hittingpoint(coef2(i, :), lambda_vec); 
end

% deal with the case where some features have alrealdy been screened off
if option.m > m
    temp = zeros(option.m, 1);
    temp(option.gamma_supp) = t_Z;
    t_Z = temp;

    temp = zeros(option.m, 1);
    temp(option.gamma_supp) = t_r;
    t_r = temp;
end

stats.Z = real(Z);
stats.t_Z = real(t_Z);
stats.r = real(r);
stats.t_r = real(t_r);

if isfield(option, 'gamma_supp') == true
    stats.gamma_supp = option.gamma_supp;
end

%%%%%%%%%%%%% W %%%%%%%%%%%%%% 
if isfield(option, 'W') == false
    option.W = 'st';
end
Z_tilde = t_Z .* (r == t_r);
switch option.W
    case 's'
        W = Z .* sign(Z - t_Z);
    case 'st'
        W = Z .* sign(Z - Z_tilde);
    case 'bc'
        W = max(Z, t_Z) .* sign(Z - t_Z);
    case 'bct'
        W = max(Z, Z_tilde) .* sign(Z - Z_tilde);
end
W = real(W);
stats.Ws = Z .* sign(Z - t_Z);
stats.Wst = Z .* sign(Z - Z_tilde);
stats.Wbc = max(Z, t_Z) .* sign(Z - t_Z);
stats.Wbct = max(Z, Z_tilde) .* sign(Z - Z_tilde);
end