function [x, d] = unflat(X)
% Rebuilds the set of grids used to generate a flat enumeration of points with flat()
%   x = UNFLAT(X)
% Parameters:
%   X       - flat array of grid points
% Returns:
%   x       - cell array of one-dimensional grids
%   d       - an array containing the grid sizes

% quick return when there is nothing to flatten
n = size(X, 1);
if n == 1, 
    x = {X}; % quick return
    d = length(X);
else
    x = cell(1, n);
    d = zeros(1, n);
    for i = 1:n,
        x{i} = sortx(X(i, :));
        d(i) = length(x{i});
    end;
end;

