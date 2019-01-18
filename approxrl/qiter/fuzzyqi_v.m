function varargout = fuzzyqi_v(x, H)
% Compute value function for a given x, or return entire grid of values
%   [X, V] = FUZZYQI_V(H)
%   [V] = FUZZYQI_V(X, H)
% Parameters:
%   X   - state vector
%   H   - fuzzy data structure. Use FUZZYQI_PREPAREH to obtain it
% Returns:
%   X   - state grids, cell array of grids, when called only with H
%   V   - matrix of state values, when called only with H. 
%       size(V) = [length(X{1}), ... length(X{end})]
%       Scalar value of state X, when called with X and H
%
% See also FUZZYQI_PREPAREH, FUZZYQI_H

% Example:
% [X, V] = fuzzyqi_v(H); mesh(X{1}, X{2}, V'); title('V(x)');

if nargin == 1,     % just H
    H = x;
    % return the value function table
    varargout = {H.X, reshape(max(H.theta, [], 2), H.DIMS.dimx)};
    return;
end;    % otherwise, go on with x and H

% compute mdegs of current state
[ind, mu] = mdegs_p(x, H.X, H.roll, H.DIMS.dimx, H.DIMS.p, H.tab);

Qa = mu' * H.theta(ind, :);
varargout = {max(Qa)};
