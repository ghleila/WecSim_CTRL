function [Jmean, J, K, trajhist] = mc_approx(model, approx, theta, X0, cfg)
% Evaluates control performance for an approximator
% Performance is evaluated by computing (approximate) returns from the set X0 of initial
% states.

% default config
CFG.gamma = [];         % discount factor
CFG.mc_eps = [];        % precision in return evaluation
CFG.mc_N = [];
CFG.Kmax = [];          % takes precedence over mc_eps
CFG.quantize = 0;       % whether to quantize the control actions
CFG.U = {};             % discrete actions (required when quantize = 1)
cfg = parseconfig(cfg, CFG);

if isempty(cfg.Kmax),
    Kmax = ceil(log(cfg.mc_eps * (1-cfg.gamma) / abs(model.maxr)) / log(cfg.gamma));
else
    Kmax = cfg.Kmax;
end;

if ~isempty(cfg.mc_N) && cfg.mc_N > 1, error('mc_approx_vectorized: stochastic systems not supported yet'); end;

p = model.p; q = model.q;

% compute X0 size
if iscell(X0),
    for i = length(X0):-1:1, N0(i) = length(X0{i}); end;
    X0 = flat(X0);
else
    N0 = [];
end;
n0 = size(X0, 2);

if cfg.quantize, 
    disp('mc_approx using quantized u');
    U = cfg.U;  % shorthand variable
    for i = length(cfg.U):-1:1, Un(i) = length(cfg.U{i}); end;  % find length of discretizations
end;

Xh = nan(p, Kmax+1, n0);
Uh = nan(q, Kmax, n0);
Rh = nan(Kmax+1, n0);
T = zeros(n0, 1);
J = nan(n0, 1); K = nan(n0, 1);

% precompute discounting vector
disc = cumprod([1 cfg.gamma+zeros(1, Kmax)]);
for i0 = 1:n0,
    Xh(:, 1, i0) = X0(:, i0);
    % simulate the system for at most Kmax, while possibly quantizing control actions
    for k = 1:Kmax,
        if cfg.quantize,
            Uh(:, k, i0) = quantize(approx.h(approx, theta, Xh(:, k, i0)), U, Un);
        else
            Uh(:, k, i0) = approx.h(approx, theta, Xh(:, k, i0));
        end;
        [Xh(:, k+1, i0), Rh(k+1, i0), T(i0)] = model.fun(model, Xh(:, k, i0), Uh(:, k, i0));
        if T(i0), break; end;    % no other rewards contribute (equivalently, they are all 0)
    end;
    
    % return of x0 is inner product of discounting and rewards
    J(i0) = disc(1:k) * Rh(2:k+1, i0);
    K(i0) = k;   % length of trajectory until terminal state was reached (or Kmax was exhausted)
end;

Jmean = mean(J);

% reshape in the shape of X0 if X0 was given in the form of a cell array of grids
if ~isempty(N0),
    J = reshape(J, N0);
    K = reshape(K, N0);
end;

if nargout >= 4,
    trajhist = varstostruct('Xh', 'Uh', 'Rh', 'T');
end;
