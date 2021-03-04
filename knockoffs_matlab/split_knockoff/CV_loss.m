function[CV_error] = CV_loss(X, D, y, nu, q, s_size, k_fold, method)
    % this function calculate CV loss for any given nu
    [n, ~] = size(X);
    test_size = n / k_fold;
    CV_error = 0;
    
    parfor k = 1: k_fold
        % generate random order of sample here for parallel computation
        rng(1);
        rand_rank = randperm(n);
        
        % determine the index for test set
        test_index = [floor(test_size * (k-1) + 1): floor(test_size * k)];
        test = rand_rank(test_index);
        
        % calculate training set and test set
        X_train = X;
        X_train(test, :) = [];
        y_train = y;
        y_train(test, :) = [];
        
        X_test = X(test, :);
        y_test = y(test, :);
        
        % find the estimated non-zero set index
        [index, ~, ~] = split_knockoff_original(X, D, y, nu, q, s_size, method);
        
        % fit beta gamma_{hat{S}} to regression
        [A_beta,A_gamma,~,tilde_y] = knockoff_construct(nu,X_train,D,y_train, s_size);
        A_gamma = A_gamma(:, index);
        lm = fitlm([A_beta,A_gamma], tilde_y, "intercept", false);
        coef = lm.Coefficients.Estimate;
        
        % calculate design matrix for beta gamma in test set
        [n_test, ~] = size(X_test);
        m = size(D, 1);
        A_beta = [X_test/sqrt(n_test);D/sqrt(nu)];
        A_gamma = [zeros(n_test,m);-eye(m)/sqrt(nu)];
        A_gamma = A_gamma(:, index);

        % calculate y for test set
        tilde_y = [y_test/sqrt(n_test);zeros(m,1)];
        
        % make prediction and store the CV loss
        y_predict = [A_beta,A_gamma] * coef;
        CV_error = CV_error + (norm(y_predict - tilde_y))^2 / k_fold / 2 / n;
    end
end