function [eJ, J, K] = mc_fuzzyq(XMFS, theta, model, X0, cfg)
% Monte-Carlo evaluation of fuzzy Q-iteration parameter matrix
% This is used in two ways:
%   - as a part of cross-entropy optimization of Q-iteration.
%   - called by mc_rbfqi
% Assumes X0 is a uniform distribution over a discrete set of states;
% each such state is a column of X0
% Assumes 'flat' action space in cfg.U
% In its current form evaluates only deterministic systems.

n0 = size(X0, 2);
J = zeros(n0, 1); K = zeros(n0, 1);
% precompute discounting vector, pre-initialize reward vector
disc = cumprod([1 cfg.gamma+zeros(1, cfg.mc_maxsteps)]);
r = zeros(cfg.mc_maxsteps, 1);
p = model.p;
X = cell(p, 1);
Np = zeros(1, p);       % this is a row vector so that it works with ndi2lin
for i = 1:length(XMFS),
    X{i} = XMFS{i}.c;
    Np(i) = length(X{i});
end;

% shorthand vars to avoid structure accesses in the loop
interph = cfg.interph;
mc_maxsteps = cfg.mc_maxsteps;
U = cfg.U;

if cfg.interph,         % compute best local actions
    [Qmax imax] = max(theta, [], 2);
    h = U(:, imax);
end;

% auxiliary vars
roll = 0 * Np;
tab = dec2base(0:2^model.p-1, 2) - 47;       % auxiliary indices table

% initial states loop
for i = 1:n0,
    x = X0(:, i);
    for k = 1:mc_maxsteps,
        % compute crisp (approximate) optimal action and apply to system
        [ind, mu] = mdegs_p(x, X, roll, Np, p, tab);
        if interph,
            [x r(k) terminal] = feval(model.fun, model, x, h(:, ind) * mu);
        else
            Qa = mu' * theta(ind, :);
            [x r(k) terminal] = feval(model.fun, model, x, U(:, find(Qa == max(Qa), 1)));
        end;
        if terminal, break; end;    % no other rewards contribute (equivalently, they are all 0)
    end;
    % return of x0 is inner product of discounting and rewards
    J(i) = disc(1:k) * r(1:k);
    K(i) = k;   % length of trajectory until terminal state was reached (or mc_maxsteps was exhausted)
end;

% expected value is the mean
eJ = mean(J);

% END rbfmceval() RETURNING expected cost J ===========================

