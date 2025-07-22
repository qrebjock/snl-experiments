function gradY = gradf(Y, D, mask)
% Compute the gradient of the SNL cost function.
%
% Inputs:
%   Y     - Configuration matrix (n × k), where each row is a point in k-dimensional space.
%   D     - n x n matrix of squared distances defining the SNL problem.
%   mask  - n x n binary matrix indicating which pairwise distances are known.
%
% Output:
%   gradY - Gradient of the SNL cost function at Y (same size as Y).
%
    if (nargin < 3) | isempty(mask)
        gradY = deltastar(delta(Y*Y') - D)*Y;
    else
        gradY = deltastar(mask.*(delta(Y*Y') - D))*Y;
    end

end

