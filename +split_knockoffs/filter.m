function [results, stats] = filter(X, D, y, option)
% Split Knockoff filter for structural sparsity problem.
% 
% Input Arguments
% X : the design matrix.
% y : the response vector.
% D : the linear transformation.
% option: options for creating the Split Knockoff statistics.
%	option.q: the desired FDR control target.
%   option.beta: choices on \beta(\lambda), can be: 'path', \beta(\lambda)
%       is taken from a regularization path; 'cv_beta', \beta(\lambda) is taken
%       as the cross validation optimal estimator \hat\beta; or 'cv_all',
%       \beta(\lambda) as well as \nu are taken from the cross validation
%       optimal estimators \hat\beta and \hat\nu.The default setting is
%       'cv_all'.
%	option.lambda_cv: a set of lambda appointed for cross validation in
%       estimating \hat\beta, default 10.*[0:-0.4:-8].
%	option.nu_cv: a set of nu appointed for cross validation in
%       estimating \hat\beta and \hat\nu, default 10.*[0:0.4:2].
%   option.nu: a set of nu used in option.beta = 'path' or 'cv_beta' for
%       Split Knockoffs, default 10.*[0:0.2:2].
%	option.lambda: a set of lambda appointed for Split LASSO path
%       calculation, default 10.*[0:-0.01:-6].
%	option.normalize: whether to normalize the data, default true.
%   option.W: the W statistics used for Split Knockoffs, can be 's', 'st',
%       'bc', 'bct', default 'st'.
% 
% Output Arguments
% results: cells with the selected variable set in each cell w.r.t. nu
%   when appropriate.
% stats: various statistics used by Split Knockoffs, e.g. Z, t_Z

if option.normalize == true || isfield(option, 'normalize') == false
    X = split_knockoffs.private.normc(X); % normalize(X)
    y = split_knockoffs.private.normc(y); % normalize(y)
end

q = option.q;
[n, p] = size(X);
m = size(D, 1);


% set default choice of beta and nu as cv optimal 
if isfield(option, 'beta') == false
    option.beta = 'cv_all';
end

if isfield(option, 'frac') == true
    option.n_2 = n - floor(n * option.frac);
end

switch option.beta
    % choose cv optimal \hat\beta and \hat\nu
    case 'cv_all'
        if isfield(option, 'rng') == true
            rng(option.rng);
        end
        rand_rank = randperm(n);
        
        % record random order for function W_fixed
        option.rand_rank = rand_rank;
        ind1 = rand_rank(1: floor(n * option.frac));
        ind2 = rand_rank(floor(n * option.frac+1): end);
        X_1 = X(ind1, :);
        y_1 = y(ind1, :);
        X_2 = X(ind2, :);
        y_2 = y(ind2, :);
        
        n_1 = size(X_1, 1);
        if n_1<p || (n-n_1)<p+m
            % give beta_hat nu_hat and support set estimation with dataset
            % 1 in high dimensional cases
            option.n2 = n-n_1;
            [beta_hat, stat_cv] = split_knockoffs.cv.cv_screen(X_1, y_1, D, option);
            X_new = X_2(:, stat_cv.beta_supp);
            D_new = D(stat_cv.gamma_supp, stat_cv.beta_supp);
            option.gamma_supp = stat_cv.gamma_supp;
        else
            [beta_hat, stat_cv] = split_knockoffs.cv.cv_all(X_1, y_1, D, option);
            X_new = X_2;
            D_new = D;
        end
        option.beta_hat = beta_hat;
        option.m = m;
        nu = stat_cv.nu;
        [W, stats] = split_knockoffs.statistics.W_fixed(X_new, D_new, y_2, nu, option);
        method = 'knockoff';
        results.sk = knockoffs.select(W, q, method);   
        method = 'knockoff+';
        results.sk_plus = knockoffs.select(W, q, method);
        stats.nu = stat_cv.nu;
        if n_1<p || (n-n_1)<p+m
            stats.gamma_supp = stat_cv.gamma_supp;
        end
    % choose \beta(\lambda) from split lasso path with a input sequence of nu
    case 'path'
        if isfield(option, 'rng') == true
            rng(option.rng);
        end
        rand_rank = randperm(n);
        % record random order for function W_fixed
        option.rand_rank = rand_rank;
        ind1 = rand_rank(1: floor(n * option.frac));
        X_1 = X(ind1, :);
        y_1 = y(ind1, :);
        n_1 = size(X_1, 1);
        option.m = m;
        if n_1<p || (n-n_1)<p+m
            % give support set estimation with dataset 1 in high
            % dimensional cases
            option.n2 = n-n_1;
            [~, stat_cv] = split_knockoffs.cv.cv_screen(X_1, y_1, D, option);
            X_new = X(:, stat_cv.beta_supp);
            D_new = D(stat_cv.gamma_supp, stat_cv.beta_supp);
            option.gamma_supp = stat_cv.gamma_supp;
        else
            X_new = X;
            D_new = D;
        end
        if isfield(option, 'nu') == true
            nu_s = option.nu;
        else
            nu_s = 10.^[0:0.2:2];
        end
        num_nu = length(nu_s);
        results = cell(num_nu, 1);
        stats = cell(num_nu, 1);
        for i = 1: num_nu
            nu = nu_s(i);
            [W, stats{i}] = split_knockoffs.statistics.W_path(X_new, D_new, y, nu, option);
            method = 'knockoff';
            results{i}.sk = knockoffs.select(W, q, method);
            method = 'knockoff+';
            results{i}.sk_plus = knockoffs.select(W, q, method);
        end
    % choose cv optimal \hat\beta with a input sequence of nu
    case 'cv_beta'
        if isfield(option, 'rng') == true
            rng(option.rng);
        end
        rand_rank = randperm(n);
        option.rand_rank = rand_rank;
        ind1 = rand_rank(1: floor(n * option.frac));
        ind2 = rand_rank(floor(n * option.frac+1): end);
        X_1 = X(ind1, :);
        y_1 = y(ind1, :);
        X_2 = X(ind2, :);
        y_2 = y(ind2, :);
        n_1 = size(X_1, 1);
        if n_1<p || (n-n_1)<p+m
            % give beta_hat and support set estimation with dataset
            % 1 in high dimensional cases
            option.n2 = n-n_1;
            [beta_hat, stat_cv] = split_knockoffs.cv.cv_screen(X_1, y_1, D, option);
            X_new = X_2(:, stat_cv.beta_supp);
            D_new = D(stat_cv.gamma_supp, stat_cv.beta_supp);
            option.gamma_supp = stat_cv.gamma_supp;
        else
            [beta_hat, stat_cv] = split_knockoffs.cv.cv_all(X_1, y_1, D, option);
            X_new = X_2;
            D_new = D;
        end
        option.beta_hat = beta_hat;
        option.m = m;
        if isfield(option, 'nu') == true
            nu_s = option.nu;
        else
            nu_s = 10.^[0:0.2:2];
        end
        num_nu = length(nu_s);
        results = cell(num_nu, 1);
        stats = cell(num_nu, 1);
        for i = 1: num_nu
            nu = nu_s(i);
            [W, stats{i}] = split_knockoffs.statistics.W_fixed(X_new, D_new, y_2, nu, option);
            method = 'knockoff';
            results{i}.sk = knockoffs.select(W, q, method);
            method = 'knockoff+';
            results{i}.sk_plus = knockoffs.select(W, q, method);
            stats{i}.nu = stat_cv.nu;
            if n_1<p || (n-n_1)<p+m
                stats{i}.gamma_supp = stat_cv.gamma_supp;
            end
        end
    % fill in random noise when p<n_2<p+m
    case 'fill'
        if isfield(option, 'rng') == true
            rng(option.rng);
        end
        rand_rank = randperm(n);
        option.rand_rank = rand_rank;
        ind1 = rand_rank(1: floor(n * option.frac));
        ind2 = rand_rank(floor(n * option.frac+1): end);
        X_1 = X(ind1, :);
        y_1 = y(ind1, :);
        X_2 = X(ind2, :);
        y_2 = y(ind2, :);
        n_2 = size(X_2, 1);
        if p<=n_2 && n_2<p+m
            beta_hat = (X' * X)^(-1) * X' * y;
            sigma_hat = sqrt(sum((y - X * beta_hat).^2) / (n - p));
            add_number = p+m-n_2;
            add_noise = randn(add_number, 1) * sigma_hat;
            y_2 = [y_2; add_noise];
            X_2 = [X_2; zeros(add_number, p)];
        end
        [beta_hat, stat_cv] = split_knockoffs.cv.cv_all(X_1, y_1, D, option);
        X_new = X_2;
        D_new = D;
        option.beta_hat = beta_hat;
        option.m = m;
        if isfield(option, 'nu') == true
            nu_s = option.nu;
        else
            nu_s = 10.^[0:0.2:2];
        end
        num_nu = length(nu_s);
        results = cell(num_nu, 1);
        stats = cell(num_nu, 1);
        for i = 1: num_nu
            nu = nu_s(i);
            [W, stats{i}] = split_knockoffs.statistics.W_fixed(X_new, D_new, y_2, nu, option);
            method = 'knockoff';
            results{i}.sk = knockoffs.select(W, q, method);
            method = 'knockoff+';
            results{i}.sk_plus = knockoffs.select(W, q, method);
            stats{i}.nu = stat_cv.nu;
        end
    % directly appoint beta(lambda)
    case 'appoint'
        option.beta_hat = option.beta_appointed;
        if isfield(option, 'nu') == true
            nu_s = option.nu;
        else
            nu_s = 10.^[0:0.2:2];
        end
        num_nu = length(nu_s);
        results = cell(num_nu, 1);
        stats = cell(num_nu, 1);
        for i = 1: num_nu
            nu = nu_s(i);
            [W, stats{i}] = split_knockoffs.statistics.W_fixed(X, D, y, nu, option);
            method = 'knockoff';
            results{i}.sk = knockoffs.select(W, q, method);
            method = 'knockoff+';
            results{i}.sk_plus = knockoffs.select(W, q, method);
        end
end
end