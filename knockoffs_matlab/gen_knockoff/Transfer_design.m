function [X_new, y_new] = Transfer_design(X, D, y)
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
end