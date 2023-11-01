function[Z, r] = hittingpoint(coef, lambdas)
% private.hittingpoint calculates the hitting time and the sign of
% respective variable in a path. 
%
% input arguments
% coef: the regularization path for one variable
% lambdas: the respective values of lambda in the path
%
% output arguments
% Z: the hitting time
% r: the sign of respective variable at the hitting time


n_lambda = length(lambdas);

Z = 0;
r = 0;

% calculate Z and r
for j = 1: n_lambda
    if abs(coef(j)) ~= 0
        Z = lambdas(j);
        r = sign(coef(j));
        break
    end
end
end