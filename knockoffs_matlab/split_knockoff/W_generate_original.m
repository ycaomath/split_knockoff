function[Z, t_Z, W] = W_generate_original(X, D, y, nu, s_size)
    % this is the main function for original step 2 knockoff, returns Z,
    % tilde_Z and W
    [~, p] = size(X);
    m = size(D, 1);
    
    % generate the design matrix
    [A_beta,A_gamma,tilde_A_gamma,tilde_y] = knockoff_construct(nu,X,D,y, s_size);
    
    %%%%%%%%% step 1 %%%%%%%%
    
    % set penalty, since SLEP does not allow the penalty to be 0, we set
    % the penalty that we want to be 0 as 10^-10
    penalty = ones(m+p,1);
    for i = 1: p
        penalty(i, 1) = 10^-10;
    end
    
    % lasso path settings for SLEP
    opts = struct; 
    opts.q = 1; % choose L_1 penalty
    opts.ind = (0: (m+p)); % choose m+p groups
    opts.gWeight = penalty;
    
    % set lambda from 0 to maxium
    lambda_max = max(abs([A_beta, A_gamma]'*tilde_y));
    lambda_min = 0;
    nlambda = 500;
    lambda_vec = [lambda_max:-(lambda_max-lambda_min)/nlambda:lambda_min];
    lambda_vec(nlambda+1) = [];
    coefs = zeros((m+p), length(lambda_vec));
    
    % lasso for step 1
    for i = 1: length(lambda_vec)
        [coefs(:, i), ~, ~] = glLeastR([A_beta, A_gamma], tilde_y, lambda_vec(i), opts);
    end
    % store beta and gamma(coef1) for step2
    betas = coefs(1: p, :);
    coef1 = coefs((p+1): (p+m), :);
    
    % calculate r and Z
    r = zeros(m, 1);
    Z = zeros(m, 1);
    for i = 1: m
        [Z(i, 1), r(i, 1)] = calculate_Z_r(coef1(i, :), lambda_vec);
    end
    % first let W=Z, when tilde_Z>Z, set W = -tilde_Z
    W = Z;
    
    %%%%%%%%%%%%% step 2 %%%%%%%%  
    coef2 = zeros(m, nlambda);
    for i = 1: nlambda
        % take beta_lambda, gamma_lambda as calculated in step 1
        y_new = tilde_y - A_beta * betas(:, i) - A_gamma * coef1(:, i);
        % calculate LASSO
        opts = struct; 
        opts.q = 1;
        opts.ind = (0: m);
        [coef2(:, i), ~, ~] = glLeastR(tilde_A_gamma, y_new, lambda_vec(i), opts);
    end  
    % calculate tilde_Z tilde_r and W
    t_Z = zeros(m, 1);
    t_r = zeros(m, 1);
    
    for i = 1: m
        [tilde_Z, tilde_r] = calculate_Z_r(coef2(i, :), lambda_vec); 
        
        % store tilde_r
        t_r(i, 1) = tilde_r;
        
        if tilde_r == r(i, 1)
            % store tilde_Z when it is positive
            t_Z(i) = tilde_Z;
            % change W to -tilde_Z when tilde_Z is larger
            if tilde_Z >= Z(i, 1)
                W(i, 1) = - tilde_Z;
            end
        end
    end
end