function result = convert_knockoff(X, D, y, option)
% split_knockoffs.private.convert_knockoff converts transformational
% sparsity to direct sparsity when the transformation is of full row rank,
% and present the selection results of Knockoffs in the trasferred
% problems.
% 
%
% input arguments
% X : the design matrix
% y : the response vector
% D : the linear transformation
% option.q : the desired FDR control target
%
% output arguments
% result: the estimated support set

[m, p] = size(D);

% calculate X_new
X_new = X * pinv(D);
y_new = y;

% multiply the  orthogonal complement
if m < p
    X_0 = X * null(D);
    Ortho = null(X_0');
    X_new = Ortho' * X_new;
    y_new = Ortho' * y_new;
end

% output selection sets
result = struct;
method = 'knockoff';
result.k = knockoffs.filter(X_new, y_new, option.q, {'fixed'}, 'Method', 'equi', 'Threshold', method);
method = 'knockoff+';
result.k_plus = knockoffs.filter(X_new, y_new, option.q, {'fixed'}, 'Method', 'equi', 'Threshold', method);
end