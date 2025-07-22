function A = erdosrenyimask(n, p)
% Generate a random symmetric binary mask (Erdős–Rényi model) with i.i.d. edges.
%
% Inputs:
%   n - Number of nodes (size of the matrix is n x n).
%   p - Probability that an edge (i, j) is present (i.e., A(i,j) = 1).
%
% Output:
%   A - n x n symmetric binary adjacency mask.
%
    A = triu(rand(n, n) <= p, 1);
    A = A + A';
end
