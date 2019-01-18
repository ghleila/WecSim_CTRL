function u = fuzzyqi_h(x, t, H)
% Compute policy for a given x
%   U = FUZZYQI_H(X, T, H)
% Parameters:
%   X   - state vector
%   T   - time argument. Ignored, added for compatibility with policy function format
%   H   - fuzzy data structure. Use FUZZYQI_PREPAREH to obtain it
% Returns:
%   U   - action dictated by policy at x
%
% See also FUZZYQI_PREPAREH, FUZZYQI_V

% t argument is ignored

% compute mdegs of current state
[ind, mu] = mdegs_p(x, H.X, H.roll, H.DIMS.dimx, H.DIMS.p, H.tab);

% compute optimal action
if H.interph,     % interpolated action
    u = (mu' * H.hstar(ind, :))';
else                % or crisp (ties broken randomly)
    Qa = mu' * H.theta(ind, :);
    ui = find(Qa == max(Qa)); ui = ui(ceil(rand * length(ui)));
    ui = lin2ndi(ui, H.DIMS.dimu);
    % compose u, variable-wise
    for q = H.DIMS.q:-1:1,
        u(q) = H.U{q}(ui(:, q));
    end;
end;
