function[result] = gen_knockoff(X, D, y, q, method)
    % method = "knockoff" or "knockoff+"
    % convert the problem to regression
    [X_new, y_new] = Transfer_design(X, D, y);
    [result, ~] = knockoffs.filter(X_new, y_new, q, {'fixed'}, 'Method', 'equi', 'Threshold', method);
end