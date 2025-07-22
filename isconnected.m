function b = isconnected(mask)
% Check whether an undirected graph is connected.
%
% Input:
%   mask - n x n binary adjacency matrix of the graph (assumed undirected).
%
% Output:
%   b - Boolean indicating whether the graph is connected.
%
    bins = conncomp(graph(mask));
    b = all(bins == 1);
end
