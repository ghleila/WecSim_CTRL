function [eJ, J] = mc_rbfnearestpolicy(c, rad, uind, model, X0, cfg) %, W)
% Monte-Carlo evaluation of parameter vector for cross-entropy policy search with RBFs
%   [EJ, J] = MC_RBFNEARESTPOLICY(C, RAD, UIND, MODEL, X0, CFG, [W])
% Monte-Carlo evaluation of RBF parameters and action assignment for one parameter sample
% of cross-entropy policy search. The policy is computing by choosing the action assigned to
% the nearest RBF.
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
%
% Returns:
%   EJ      - mean cost across initial states (weighted cost not supported!)
%   J       - cost of every initial state sample
%
% 2009-02-03: commented out fixed-noise code

n0 =size(X0, 2);    % number of states to start from
N = size(c, 2);     % number of RBFs

% precompute discounting vector, pre-initialize reward vector
disc = cumprod([1 cfg.gamma+zeros(1, cfg.mc_maxsteps)]);
r = zeros(cfg.mc_maxsteps, 1);
% form explicit policy for ease of use
h = cfg.U(:, uind);

% Use this with old data
% if ~isfield(model, 'det') || model.det,         % deterministic system
%     % this is a quick fix for older models that do not initialize the 'det' field
%     model.det = 1;
% end;

if isfield(cfg, 'xfiltering'),
    xfiltering = cfg.xfiltering;
else xfiltering = 0;
end;

if model.det,         
    % deterministic system
    J = zeros(n0, 1);
    for i = 1:n0,
        x = X0(:, i);
        for k = 1:cfg.mc_maxsteps,
            % Equivalent to: (implemented here directly to save a function call)
            % [actmax imax] = max(rbf(x, N, c, rad));
            if xfiltering,
                [actmax imax] = max( exp(-sum( ((repmat(model.xfilter(x, model), 1, N)  - c) ./ rad).^2 )) );
            else
                [actmax imax] = max( exp(-sum( ((repmat(x, 1, N)  - c) ./ rad).^2 )) );
            end;
            [x r(k) terminal] = feval(model.fun, model, x, h(:, imax));
            if terminal, break; end;
        end;
        % return of x0 is inner product of discounting and rewards
        J(i) = disc(1:k) * r(1:k);
    end;
    % mean score
    eJ = mean(J);

else
    % stochastic system, run mc_nsim simulations for each point
    J = zeros(n0, cfg.mc_nsim);
%     useW = nargin >= 7;     % use fixed noise sequences
    for i = 1:n0,
        for s = 1:cfg.mc_nsim,
            x = X0(:, i);
%             % W supplied, extract noise vector for current simulation
%             if useW, w = reshape(W(:, i, s, :), model.nw, cfg.mc_maxsteps); end;
            for k = 1:cfg.mc_maxsteps,
                if xfiltering,
                    [actmax imax] = max( exp(-sum( ((repmat(model.xfilter(x, model), 1, N)  - c) ./ rad).^2 )) );
                else 
                    [actmax imax] = max( exp(-sum( ((repmat(x, 1, N)  - c) ./ rad).^2 )) );
                end;
%                 if useW,    
%                     [x r(k) terminal] = feval(model.fun, model, x, h(:, imax), w(:, k));
%                 else
                    [x r(k) terminal] = feval(model.fun, model, x, h(:, imax));
%                 end;
                if terminal, break; end;
            end;
            % return of x0 is inner product of discounting and rewards
            J(i, s) = disc(1:k) * r(1:k);
        end;
    end;
    % first compute expectation over stochasticity
    % then average over initial states to get final value
    eJ = mean(mean(J, 2));
    
end;


% END rbfmceval() RETURNING expected cost eJ, cost vector J ===========================