function [ind, a] = nrbftrunc(x, N, c, rad, thresh)
% Computes normalized radial basis functions values for a continuous x
%   [IND, A] = NRBF(X, N, C, RAD, THRESH)
%
% Parameters:
%   X       - continuous variable (column vector of dimension p)
%   N       - number of basis functions
%   C       - centers of basis functions, one on each column, size p x N
%   RAD     - radii of basis functions, one on each column, size p x N
%   THRESH  - truncation threshold
% Returns:
%   IND     - indices of RBF values > 0, length(IND) <= N
%   A       - normalized RBF values > 0 for given X, 1xlength(IND)

a = exp(-sum( ((repmat(x, 1, N)  - c) ./ rad).^2 ));
% first truncate
ind = find(a > thresh);
a = a(ind);
% then normalize
a = a/sum(a);

end
% END nrbf returning a ------------------------------------------------------
