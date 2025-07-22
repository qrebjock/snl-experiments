function hessY = hessf(Y, D, V, mask)
% Compute an Hessian-vector product of the SNL cost function.
%
% Inputs:
%   Y     - Configuration matrix (n × k), the point at which the Hessian is evaluated.
%   D     - n x n matrix of squared distances defining the SNL problem.
%   V     - n x k matrix, direction in which the Hessian is evaluated.
%   mask  - n x n binary matrix indicating which pairwise distances are known.
%
% Output:
%   hessY - The product of the Hessian of the SNL cost function at Y with the direction V.
%
    if (nargin < 4) | isempty(mask)
        hessY = deltastar(delta(Y*V' + V*Y'))*Y + deltastar(delta(Y*Y') - D)*V;
    else
        hessY = deltastar(delta(Y*V' + V*Y').*mask)*Y + deltastar((delta(Y*Y') - D).*mask)*V;
    end

end

