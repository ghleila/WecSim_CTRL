function [eJ, J] = hiv_mc_rbfnearestpolicy_varode(c, rad, uind, m, X0, cfg, W)
% Optimized Monte Carlo evaluation of nearest-RBF policy score, for HIV
% system and variable-step integration

n0 =size(X0, 2);    % number of states to start from
N = size(c, 2);     % number of RBFs

% precompute discounting vector, pre-initialize reward vector
disc = cumprod([1 cfg.gamma+zeros(1, cfg.mc_maxsteps)]);
r = zeros(cfg.mc_maxsteps, 1);
% form explicit policy for ease of use
h = cfg.U(:, uind);

Q = m.Q; R1 = m.R1; R2 = m.R2; S = m.S; % to avoid repeated structure accesses
J = zeros(n0, 1); 
for i = 1:n0,
    x = X0(:, i);
    for k = 1:cfg.mc_maxsteps,
        if cfg.xfiltering,
            [actmax imax] = max( exp(-sum( ((repmat(m.xfilter(x, m), 1, N)  - c) ./ rad).^2 )) );
        else
            [actmax imax] = max( exp(-sum( ((repmat(x, 1, N)  - c) ./ rad).^2 )) );
        end;
        u = h(:, imax);
        % Equivalent to: (implemented here directly to save a function call)
        % [actmax imax] = max(rbf(m.xfilter(x, m), N, c, rad));
        
        % transition
    	[odet odey] = m.odesolver(@hiv_trans, m.odet, x, m.odeopt, m, u);
        x = max(m.minx, odey(end, :)'); % negative values are nonsensical
        % reward
        r(k) = - Q * x(5) - R1 * u(1) * u(1) - R2 * u(2) *u(2) + S * x(6);
    end;
    % return of x0 is inner product of discounting and rewards
    J(i) = disc(1:k) * r(1:k);
end;

if isfield(cfg, 'P0') && ~isempty(cfg.P0) && ~any(isnan(cfg.P0)),
    eJ = sum(cfg.P0(:) .* J(:));
else
    % score is the average of the returns over the initial states
    eJ = mean(J);
end

% END rbfmceval() RETURNING expected cost eJ, cost vector J
% ===========================