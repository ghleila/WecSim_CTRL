function [eJ, J] = mc_fuzzyqiter(c, theta, model, X0, cfg)
% Monte-Carlo evaluation of fuzzy Q-iteration parameter matrix
% Intended for use with cefuzzyqi, use mc_fuzzyq to evaluate solutions
% computed with fuzzyqi
% Assumes X0 is a uniform distribution over a discrete set of states;
% each such state is a column of X0
% Assumes 'flat' action space in cfg.U
% In its current form evaluates only deterministic systems.

Np = cfg.Np; roll = cfg.roll; p = model.p;
% helper index table for mdegs computation
tab = cfg.tab;
n0 = size(X0, 2);

J = zeros(n0, 1);
% precompute discounting vector, pre-initialize reward vector
disc = cumprod([1 cfg.gamma+zeros(1, cfg.mc_maxsteps)]);
r = zeros(cfg.mc_maxsteps, 1);
for i = 1:n0,
    x = X0(:, i);
    for k = 1:cfg.mc_maxsteps,
        % compute crisp (approximate) optimal action and apply to system
        [ind, mu] = mdegs_p(x, c, roll, Np, p, tab);
        Qa = mu' * theta(ind, :);
        [x r(k) terminal] = feval(model.fun, model, x, cfg.U(:, find(Qa == max(Qa), 1)));
        if terminal, break; end;    % no other rewards contribute (equivalently, they are all 0)
    end;    
    % return of x0 is inner product of discounting and rewards
    J(i) = disc(1:k) * r(1:k);
end;

% expected value is the mean
eJ = mean(J);

% END rbfmceval() RETURNING expected cost J ===========================