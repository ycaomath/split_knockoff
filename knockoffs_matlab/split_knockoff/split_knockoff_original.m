function[result, Z, t_Z] = split_knockoff_original(X, D, y, nu, q, s_size, method)
    % record Z, tilde_Z as well as W for further calculation
    [Z, t_Z, W] = W_generate_original(X, D, y, nu, s_size);
    % selection on W
    result = knockoffs.select(W, q, method);
end