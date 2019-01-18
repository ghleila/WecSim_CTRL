function [eJ, J] = hiv_mc_rbfvotingpolicy_varode(c, rad, uind, m, X0, cfg, W)
% Optimized Monte Carlo evaluation of voting-RBF policy score, for HIV
% system and variable-step integration

n0 =size(X0, 2);    % number of states to start from
N = size(c, 2);     % number of RBFs

% precompute discounting vector, pre-initialize reward vector
disc = cumprod([1 cfg.gamma+zeros(1, cfg.mc_maxsteps)]);
r = zeros(cfg.mc_maxsteps, 1);
[junique, selectors] = ind2selectors(uind);
% Make phi a column vector such that the sum across colums works properly even when all the RBFs
% have the same assigned discrete action
phi = zeros(N+1, 1); % pad with a zero "dummy" for the sums
onetoN = 1:N;

Q = m.Q; R1 = m.R1; R2 = m.R2; S = m.S; % to avoid repeated structure accesses
J = zeros(n0, 1);
for i = 1:n0,
    x = X0(:, i);
    for k = 1:cfg.mc_maxsteps,
        % compute action
        if cfg.xfiltering,
            phi(onetoN) = exp(-sum( ((repmat(m.xfilter(x, m), 1, N)  - c) ./ rad).^2 ));  % direct RBF computation
        else
            phi(onetoN) = exp(-sum( ((repmat(x, 1, N)  - c) ./ rad).^2 ));  % direct RBF computation
        end;
        [actmax imax] = max(sum(phi(selectors), 1));
        u = cfg.U(:, junique(imax));

        % transition
    	[odet odey] = m.odesolver(@hiv_trans, m.odet, x, m.odeopt, m, u);
        x = max(m.minx, odey(end, :)'); % negative values are nonsensical
        % reward
        r(k) = - Q * x(5) - R1 * u(1) * u(1) - R2 * u(2) *u(2) + S * x(6);
    end;
    % return of x0 is inner product of discounting and rewards
    J(i) = disc(1:k) * r(1:k);
end;
% expected value is the mean
eJ = mean(J);

% END RETURNING expected cost eJ, cost vector J ===========================
