function [xu, yp] = removeduplicates(x, y, tol, aggregatex, aggregatey)
% works on ROW points
% Notes finds sets of points {X | for any xi in X, there is xj in X with |xi-xj|_inf <= tol}
% and NOT sets {X | |xi - xj| <= tol for any xi, xj in X}
% It can happen that the distance between points in a set is larger than the tolerance

% default aggregation methods
if nargin < 4, aggregatex = @(x) x(1, :); end;     % pick the first value
if nargin < 5, aggregatey = @nanmean; end;         % compute the mean, ignoring missing (nan) values

[x2, i2] = sortrows(x);
y2 = y(i2, :);
% find indices where points are NOT close -- these are the starting points for the index ranges
ind = find([1; any(abs(diff(x2, 1)) > tol, 2)]);

N = length(ind);
xu = zeros(N, size(x, 2));
yp = zeros(N, size(y, 2));
ind = [ind; size(x, 1)+1];  % for the last interval
for i = 1:N,
    xu(i, :) = feval(aggregatex, x2(ind(i):ind(i+1)-1, :));
    yp(i, :) = feval(aggregatey, y2(ind(i):ind(i+1)-1, :));
end;
