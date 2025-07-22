function centeredY = center(Y, where)
% Center a vector or matrix using the centering matrix I - ones(n, n)/n.
%
% Inputs:
%   Y      - The vector or matrix to be centered.
%   where  - Specifies how to apply the centering matrix. 
%            Options are:
%              "left"  — apply on the left: (I - 11^\top/n) * Y
%              "right" — apply on the right: Y * (I - 11^\top/n)
%              "both"  — apply on both sides.
%
% Output:
%   centeredY - Centered version of Y.
%
  if nargin == 1
      where = "left";
  end

  [m, n] = size(Y);

  if where == "left"
      centeredY = Y - ones(m, 1) * (ones(1, m) * Y)/m;
  elseif where == "right"
      centeredY = Y - (Y * ones(n, 1)) * ones(1, n)/n;
  elseif where == "both"
      centeredY = Y - ones(m, 1) * (ones(1, m) * Y)/m;
      centeredY = centeredY - (centeredY * ones(n, 1)) * ones(1, n)/n;
  else
      warning("center: wrong argument for variable where.")
  end

end
