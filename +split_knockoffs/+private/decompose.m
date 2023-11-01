function [U,S,V,U_perp] = decompose(X, D)
% split_knockoffs.private.decompose  Decompose the design matrix X.

[n, p] = size(X);
[m, ~] = size(D);

% Factorize X as X = USV' (reduced SVD).
[U,S,V] = knockoffs.private.canonicalSVD(X);

% Construct an orthogonal matrix U_perp such that U_perp'*X = 0.
[Q,~] = qr([U zeros(n,m)], 0); % Skinny QR.
U_perp = Q(:,p+1:end);

end