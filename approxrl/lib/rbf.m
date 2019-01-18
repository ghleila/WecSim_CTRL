function a = rbf(x, N, c, rad)
% Computes unnormalized radial basis functions values for a continuous x
%   A = RBF(X, N, C, RAD)
%
% Parameters:
%   X       - continuous variable (column vector of dimension p)
%   N       - number of basis functions
%   C       - centers of basis functions, one on each column, size p x N
%   RAD     - radii of basis functions, one on each column, size p x N
% Returns:
%   A       - unnormalized RBF values for given X

a = exp(-sum( ((repmat(x, 1, N)  - c) ./ rad).^2, 1 ));

% END rbf returning a ------------------------------------------------------
