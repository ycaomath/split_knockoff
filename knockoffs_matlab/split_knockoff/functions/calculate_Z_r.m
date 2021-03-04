function[Z, r] = calculate_Z_r(coef, lambdas)
    % this function calculate Z and r respectively for any given lasso path
    n_lambda = size(coef, 2);
    
    Z = 0;
    r = 0;
    
    % calculate Z and r
    for j = 1: n_lambda
        if abs(coef(1, j)) >= 10^-6
            Z = lambdas(1, j);
            r = sign(coef(1, j));
            break
        end
    end
end