function deltaX = delta(X)
% Compute the squared Euclidean distance matrix from a Gram (inner product) matrix.
%
% Input:
%   X - Symmetric square matrix of inner products (i.e., a Gram matrix).
%
% Output:
%   deltaX - Matrix of squared pairwise Euclidean distances.
%
    [n, ~] = size(X);
    deltaX = diag(X)*ones(1, n) + ones(n, 1)*diag(X)' - X - X';
end

