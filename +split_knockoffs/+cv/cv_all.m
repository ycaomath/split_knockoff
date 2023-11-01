function [beta, stat]= cv_all(X, y, D, option)
% split_knockoffs.cv.cv_all calculate the CV optimal beta
% in the problem 1/n |y - X beta|^2 + 1/nu |D beta - gamma|^2 + lambda |gamma|_1.
% 
% input arguments
% X : the design matrix
% y : the response vector
% D : the linear transform
% 
% output arguments
% beta: CV optimal beta
% stat: various intermedia statistics


k_fold = 5;
[n, ~] = size(X);
[m, p] = size(D);

% appoint the set of \nu
if isfield(option, 'nu_cv') == false
    nu_s = 10.^[0:0.4:2];
else
    nu_s = option.nu_cv;
end

% appoint a set of lambda
if isfield(option, 'lambda_cv') == false
    lambda_s = 10.^[0:-0.4:-8];
%     lambda_s = 1;
else
    lambda_s = option.lambda_cv;
end
nlambda = length(lambda_s);

test_size = floor(n / k_fold);

% randomly split data
rng(1);
rand_rank = randperm(n);

% create matrix to store result for split lasso test loss
loss_sl = zeros(k_fold, length(nu_s), nlambda);


for k = 1: k_fold
    
    % generate test set
    test_index = [test_size * (k-1) + 1: test_size * k];
    test = rand_rank(test_index);
    
    % training set
    X_train = X;
    X_train(test, :) = [];
    y_train = y;
    y_train(test, :) = [];
    
    % test set
    X_test = X(test, :);
    y_test = y(test, :);
    
    % normalization
    X_train = bsxfun(@minus,X_train,mean(X_train,1));
    y_train = bsxfun(@minus,y_train,mean(y_train,1));
    X_test = bsxfun(@minus,X_test,mean(X_test,1));
    y_test = bsxfun(@minus,y_test,mean(y_test,1));
    
    % test loss for split lasso
    for i = 1: length(nu_s)
        nu = nu_s(i);
        
        % generate split LASSO design matrix
        [n_train, ~] = size(X_train);
        A_beta = [X_train/sqrt(n_train);D/sqrt(nu)];
        A_gamma = [zeros(n_train,m);-eye(m)/sqrt(nu)];
        tilde_y = [y_train/sqrt(n_train);zeros(m,1)];
        
        % set the penalty for glmnet
        penalty = ones(m+p,1);
        for tmp = 1: p
            penalty(tmp) = 0;
        end

        opts.standardize = false; % cancel standardize for glmnet
        opts.lambda = lambda_s; % choose lambda
        opts.intr = false;
        opts.penalty_factor = penalty; % choose the penalty
        options = glmnetSet(opts);

        % fit glmnet
        fit = glmnet([A_beta, A_gamma], tilde_y, [], options);
        % extract beta
        coefs = fit.beta;
        
        for j =  1: length(lambda_s)
            coef = coefs(:, j);
            beta = coef(1:p);
            % calculate loss
            y_sl = X_test * beta;
            loss_sl(k, i, j) = norm(y_sl-y_test)^2/test_size;
        end
    end
    
end

mean_loss_sl = mean(loss_sl, 1);
mean_loss_sl = reshape(mean_loss_sl, [length(nu_s), nlambda]);

% find minimal
[nu_number, lambda_number] = find(mean_loss_sl == min(min(mean_loss_sl)), 1);
nu_sl = nu_s(nu_number);
lambda_sl = lambda_s(lambda_number);

stat.nu = nu_sl;
stat.lambda = lambda_sl;
% calculate beta
A_beta = [X/sqrt(n);D/sqrt(nu_sl)];
A_gamma = [zeros(n,m);-eye(m)/sqrt(nu_sl)];
tilde_y = [y/sqrt(n);zeros(m,1)];

opts.standardize = false; % cancel standardize for glmnet
opts.lambda = lambda_sl; % choose lambda
opts.intr = false;
opts.penalty_factor = penalty; % choose the penalty
options = glmnetSet(opts);
            
% fit glmnet
fit = glmnet([A_beta, A_gamma], tilde_y, [], options);
% extract beta
coef = fit.beta;
beta = coef(1:p);
end