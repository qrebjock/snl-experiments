function s = inner(A, B)
% Compute the Frobenius inner product between two matrices.
%
% Inputs:
%   A - First matrix.
%   B - Second matrix, same size as A.
%
% Output:
%   s - Frobenius inner product of A and B.
%
    s = sum(A.*B, "all");
end
