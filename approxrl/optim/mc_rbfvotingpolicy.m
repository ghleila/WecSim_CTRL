function [eJ, J] = mc_rbfvotingpolicy(c, rad, uind, model, X0, cfg, W)
% Monte-Carlo evaluation of parameter vector for cross-entropy voting policy search with RBFs
%   [EJ, J] = MC_RBFVOTINGPOLICY(C, RAD, UIND, MODEL, X0, CFG, [W])
% Monte-Carlo evaluation of RBF parameters and action assignment for one parameter sample
% of cross-entropy policy search. The policy is computing by summing up the activations of
% all the RBFs assigned to the same action, and then choosing the action with the largest vote.
%
% Parameters:
%   C       - RBF centers, p x N where p state variable size, N number of basis functions
%   RAD     - RBF radii, p x N where p state variable size, N number of basis functions
%   UIND    - assignment of discrete actions to centers, in 1-based integer indices
%   MODEL   - problem model
%   X0      - distribution of initial states. For now, assumed to be uniform over a 
%       set of discrete values, each given in a column of X0
%   CFG     - configuration parameters. Requires fields (subject to change):
%       U, mc_maxsteps, mc_nsim - see cerbfqi for description
%   W       - (optional) pre-generated noise sequences, if they are to be used
%
% Returns:
%   EJ      - mean cost
%   E       - cost of every initial state sample
% 2009-02-03: Commented out fixed noise code

n0 =size(X0, 2);    % number of states to start from
N = size(c, 2);     % number of RBFs

% precompute discounting vector, pre-initialize reward vector
disc = cumprod([1 cfg.gamma+zeros(1, cfg.mc_maxsteps)]);
r = zeros(cfg.mc_maxsteps, 1);
% Compute the selectors corresponding to the indices in uind, see
% ind2selectors comments
[junique, selectors] = ind2selectors(uind);
% Make phi a column vector such that the sum across colums works properly even when all the RBFs
% have the same assigned discrete action
phi = zeros(N+1, 1); % add a zero "dummy" which is referred to by the padded elements in selectors
onetoN = 1:N;

if isfield(cfg, 'xfiltering'),
    xfiltering = cfg.xfiltering;
else xfiltering = 0;
end;
if model.det,         % deterministic system
    J = zeros(n0, 1);
    for i = 1:n0,
        x = X0(:, i);
        for k = 1:cfg.mc_maxsteps,
            if xfiltering,
                phi(onetoN) = exp(-sum( ((repmat(model.xfilter(x, model), 1, N)  - c) ./ rad).^2 ));  % direct RBF computation
            else 
                phi(onetoN) = exp(-sum( ((repmat(x, 1, N)  - c) ./ rad).^2 ));  % direct RBF computation
            end;
            % 1. sum phi(selectors) on rows, obtaining the cumulated
            % activation for every unique index in uind (junique)
            % since phi has one extra zero, the padded indices equal to
            % card(U)+1 in selectors have no effect
            % 2. Then, find the maximum sum and its corresponding index
            [actmax imax] = max(sum(phi(selectors), 1));
            % use junique to recover the actual action index that has the
            % maximum summed activation
            [x r(k) terminal] = feval(model.fun, model, x, cfg.U(:, junique(imax)));
            if terminal, break; end;
        end;
        % return of x0 is inner product of discounting and rewards
        J(i) = disc(1:k) * r(1:k);
    end;
    % expected value is the mean
    eJ = mean(J);
else                                            % stochastic system, run mc_nsim simulations for each point
%     useW = nargin >= 7;
    J = zeros(n0, cfg.mc_nsim);
    for i = 1:n0,
        for s = 1:cfg.mc_nsim,
            x = X0(:, i);
%             % verify if noise sequences supplied
%             if useW, w = reshape(W(:, i, s, :), model.nw, cfg.mc_maxsteps); end;
            for k = 1:cfg.mc_maxsteps,
                if xfiltering,
                    phi(onetoN) = exp(-sum( ((repmat(model.xfilter(x, model), 1, N)  - c) ./ rad).^2 ));  % direct RBF computation
                else 
                    phi(onetoN) = exp(-sum( ((repmat(x, 1, N)  - c) ./ rad).^2 ));  % direct RBF computation
                end;
                [actmax imax] = max(sum(phi(selectors), 1));
%                 if useW,
%                     [x r(k) terminal] = feval(model.fun, model, x, cfg.U(:, junique(imax)), w(:, k));
%                 else
                    [x r(k) terminal] = feval(model.fun, model, x, cfg.U(:, junique(imax)));
%                 end;
                if terminal, break; end;
            end;
            % return of x0 is inner product of discounting and rewards
            J(i, s) = disc(1:k) * r(1:k);
        end;
    end;
    % first compute expectation over stochasticity
    % then over initial state to get final value
    eJ = mean(mean(J, 2));
end;

% END rbfmceval() RETURNING expected cost eJ, cost vector J ===========================

