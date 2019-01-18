function varargout = fuzzyqi_q(x, u, H)
% Compute Q-value for a given (x, u) pair, or return entire Q-function
%   [X, U, V] = FUZZYQI_V(H)
%   [V] = FUZZYQI_V(X, H)
% Parameters:
%   X   - state vector
%   U   - action vector (must be one of the discrete actions)
%   H   - fuzzy data structure. Use FUZZYQI_PREPAREH to obtain it
% Returns:
%   X   - state grids, cell array of grids, when called only with H
%   U   - action grids, cell array of grids, when called only with H
%   V   - matrix of state values, when called only with H. 
%       size(V) = [length(X{1}), ... length(X{end}), length(U{1}), ... length(U{end})]
%       Scalar value of state X, when called with X and H
%
% See also FUZZYQI_PREPAREH, FUZZYQI_H, FUZZYQI_V

% Example:
% [X, U, Q] = fuzzyqi_q(H); mesh(X{1}, X{2}, Q(:, :, 1)'); title('Q(x, u_0)');

if nargin == 1,     % just H
    H = x;
    % return the value function table
    varargout = {H.X, H.U, reshape(H.theta, [H.DIMS.dimx H.DIMS.dimu])};
    return;
end;    % otherwise, go on with x and H

% compute mdegs of current state
[ind, mu] = mdegs_p(x, H.X, H.roll, H.DIMS.dimx, H.DIMS.p, H.tab);

Qa = mu' * H.theta(ind, :); % Q-values for all the actions
varargout = {Qa(findflat(flat(H.U), u))};
