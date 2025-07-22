function fY = f(Y, D, mask)
% Evaluate the SNL cost function.
%
% Inputs:
%   Y    - Configuration at which to evaluate the cost.
%   D    - Matrix of squared distances (that defines the SNL problem).
%   mask - Binary matrix indicating which distances are known (same size as D).
%
% Output:
%   fY   - Value of the cost function at configuration Y.
%
    X = Y*Y';
    M = delta(X) - D;

    if (nargin == 3) & ~isempty(mask)
        M = mask.*M;
    end

    fY = 1/4*inner(M, M);

end

