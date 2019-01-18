function a = rbfvector(x, n, N, c, rad, normalize)
% Computes normalized or unnormalized radial basis functions values for a continuous x
%   A = RBFVECTOR(X, n, N, C, RAD, [NORMALIZE])
%
% Parameters:
%   X       - vector of continuous variables. Matrix of dimension p x n
%       where each column is a continuous (p-dimensional) variables
%   n       - number of points for which to compute the RBF values (number of columns in X)
%   N       - number of RBFs
%   C       - centers of RBFs, one on each column, size p x N
%   RAD     - radii of basis functions, one on each column, size p x N
%   NORMALIZE   - whether to normalize RBF values. Optional, default 1 (normalize).
% Returns:
%   A       - normalized or unnormalized (depending on NORMALIZE) 
%           RBF values for given X. A matrix of size n x N, each row the activation
%           vector of one point.

% Compute matrix of distances normalized by radii for all points,
% by flattening x and repeating the centers and radii accordingly
a = ((repmat(x(:), 1, N) - repmat(c, n, 1)) ./ repmat(rad, n, 1)) .^ 2;
% Reshape such that variable dimension is on the first dimension of the matrix
a = reshape(a, [], n, N);
% Sum across variable dimesion, apply exponential
a = exp(-squeeze(sum(a)));
% Normalize by default or if NORMALIZE = 1
if nargin < 6 || normalize,
    a = a ./ repmat(sum(a, 2), 1, N);
end;

% END rbf returning a ------------------------------------------------------

