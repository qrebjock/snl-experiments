function [Y, cost, info, options, problem] = findminimum(n, k, D, mask, Y0, options)
% Attempt to find a minimizer of a SNL problem instance.
%
% Inputs:
%   n       - Number of points (nodes) in the problem.
%   k       - Embedding dimension for optimization (e.g., 2 or 3).
%   D       - n x n matrix of squared distances (may contain noise) that defines the SNL problem.
%   mask    - n x n binary matrix indicating which pairwise distances are known.
%   Y0      - Initial configuration (n x k matrix) for the optimization algorithm.
%   options - Struct of optimization options (e.g., tolerance, max iterations).
%
% Outputs:
%   Y       - Final estimated configuration (n × k).
%   cost    - Final value of the cost function.
%   info    - Struct containing optimization metadata.
%   options - The options struct, possibly with additional diagnostic information.
%   problem - Struct describing the optimization problem (e.g., objective function, gradient).
%
    if nargin < 4
        mask = [];
    end

    if nargin < 5
        Y0 = center(randn(n, k));
    end

    if nargin < 6
        options = struct();
    end

    defaultoptions.verbosity = 0;
    defaultoptions.maxiter = 5000;
    defaultoptions.tolgradnorm = 1e-13;
    defaultoptions.Delta_bar = 30;

    options = fillstructwithdefault(options, defaultoptions);

    problem.M = euclideanfactory(n, k);

    problem.cost  = @(Y) f(Y, D, mask);
    problem.grad = @(Y) gradf(Y, D, mask);
    problem.hess = @(Y, U) hessf(Y, D, U, mask);

    [Y, cost, info] = trustregions(problem, Y0, options);
    Y = center(Y);
end
