function[CV_error, nu] = CV(X, D, y, nu_s, q, s_size, k_fold, method)
    % this function out put cross validation error and best nu
    
    num_nu = length(nu_s);
    CV_error = zeros(num_nu, 1);
    
    % calculate CV_error
    parfor i = 1: num_nu
        nu = nu_s(i);
        CV_error(i, 1) = CV_loss(X, D, y, nu, q, s_size, k_fold, method);
    end
    
    % find best nu
    index = find(CV_error == min(CV_error));
    nu = nu_s(index);
end