function [deltasX] = deltastar(X)
% Adjoint of the Euclidean distance matrix map (delta).
%
% This function computes the adjoint (with respect to the Frobenius inner product)
% of the linear map that takes a Gram matrix to a squared Euclidean distance matrix.
%
% Input:
%   X - An n-by-n matrix (typically a matrix of squared distances).
%
% Output:
%   deltasX - The adjoint image, an n-by-n symmetric matrix.
%
    [n, ~] = size(X);
    deltasX = diag(X*ones(n, 1)) + diag(X'*ones(n, 1)) - X - X';
end

