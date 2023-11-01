function[W, stats] = W_fixed(X, D, y, nu, option)
% split_knockoffs.statistics.pathorder.W_path generates the split knockoff
% statistics W based on a fixed \beta(\lambda) = \hat\beta in the
% intercepetion assignment step.
%
% input arguments:
% X : the design matrix
% y : the response vector
% D : the linear transformation
% nu: the parameter for variable splitting
% option: options for creating theSplit Knockoff statistics
%
% output arguments
% W: the split knockoff statistics
% stats: the intermedia statistics

[m, ~] = size(D);

if m == 0
    W = zeros(option.m, 1);
    stats.Z = zeros(option.m, 1);
    stats.t_Z = zeros(option.m, 1);
    stats.r = zeros(option.m, 1);
    stats.t_r = zeros(option.m, 1);
    return
end


beta_hat = option.beta_hat;

%%%%%%%%%%%%% step 1 %%%%%%%%%%%%%%

% generate the design matrix

opts.copy = true;
[A_beta,A_gamma,tilde_y,tilde_A_gamma] = split_knockoffs.create(X, y, D, nu, opts);

% set lambda
if isfield(option, 'lambda') == true
    lambda_vec = option.lambda;
else
    lambda_vec = 10.^[0:-0.01:-6];
end

% lasso path settings for glmnet
opts = struct; 
opts.lambda = lambda_vec;
opts = glmnetSet(opts);

y_new = tilde_y - A_beta * beta_hat;

fit_step1 = glmnet(A_gamma, y_new, [], opts);
coef1 = fit_step1.beta;
% calculate r and Z
r = zeros(m, 1);
Z = zeros(m, 1);
for i = 1: m
    [Z(i), r(i)] = split_knockoffs.private.hittingpoint(coef1(i, :), lambda_vec);
end

if isfield(option, 'gamma_supp') == true
    temp = zeros(option.m, 1);
    temp(option.gamma_supp) = Z;
    Z = temp;

    temp = zeros(option.m, 1);
    temp(option.gamma_supp) = r;
    r = temp;
end

%%%%%%%%%%%%% step 2 %%%%%%%%%%%%%% 
fit_step2 = glmnet(tilde_A_gamma, y_new, [], opts);
coef2 = fit_step2.beta;
t_Z = zeros(m, 1);
t_r = zeros(m, 1);

for i = 1: m
    [t_Z(i), t_r(i)] = split_knockoffs.private.hittingpoint(coef2(i, :), lambda_vec); 
end

if isfield(option, 'gamma_supp') == true
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